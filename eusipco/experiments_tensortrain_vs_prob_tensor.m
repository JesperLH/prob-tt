function experiments_tensortrain_vs_prob_tensor()
%% Tensor Train Experiments
addpath('./requiredfunctions/')
addpath(genpath('./'))
addpath(genpath('../ncptensor'))
addpath(genpath('../tools'))
addpath(strjoin(strcat('../thirdparty-matlab/oseledets-TT-Toolbox-8332a6c/',{'','core', 'exp', 'cross', 'fmex', 'misc','solve'}),';'))

rng(131462234)
%% Setup problem and generate data
num_repeats = 10;
maxiter = 50;
tt_threshold = logspace(-16,-1,16); % Vary threshold...)
snr_list = [-20:2.5:20, 30:10:100];

N = 20:-1:16;
D = [1,length(N)+1:-1:3, 1];

[X_clean, G_true, ES, EV] = generateTensorTrain(N, D);
X_clean = X_clean/sqrt(var(X_clean(:)));


%% Run experiments


n_methods = length(tt_threshold)... % Oseledets TT
    +1 +1; % One fixed-order TT and prob_tt with random init
name_methods = strcat('TT (\epsilon=',strsplit(strtrim(sprintf('%2.0e ',tt_threshold)),' '),')');
name_methods = {name_methods{:}, 'TT (fixed)', 'Prob. TT (rand)'};

%%

all_rmse_noise = nan(length(snr_list), n_methods); % This could be varied over datasets..
all_rmse_true = nan(length(snr_list), n_methods);
all_elbo = nan(length(snr_list), n_methods);
all_tt_comp = nan(length(snr_list), n_methods);
all_numel = nan(length(snr_list), n_methods);


for snr_i = 1:length(snr_list)
    fprintf('Running SNR = %3.2f. (level %i of %i)\n', snr_list(snr_i), snr_i, length(snr_list))
    
    X=addTensorNoise(X_clean, snr_list(snr_i));
    method_i = 1;
    %% Calculate Oseledets TT
    fprintf('\tRunning Oseledets TT-toolbox\n')
    for t = 1:length(tt_threshold)
        fprintf('\t\t epsilon = %3.2e .... ', tt_threshold(t)); t00 = tic;
        tt = tt_tensor(X, tt_threshold(t));
        cr= tt.core ; ps= tt.ps ;
        core=cell(ndims(X),1);
        for k=1:ndims(X)
            core{k} =cr(ps(k): ps(k +1) -1); %#ok<SAGROW>
            core{k} = reshape(core{k}, tt.r(k), tt.n(k), tt.r(k+1)); %#ok<SAGROW>
        end
        toc(t00);
        
        [e_rmse, e_rmse_clean, e_numel, ~] = calculateErrorStuff(X, X_clean, core);% G_true)
        all_rmse_noise(snr_i, method_i) = e_rmse;
        all_rmse_true(snr_i, method_i) = e_rmse_clean;
        all_numel(snr_i, method_i) = e_numel;
        method_i = method_i+1;
    end
    
    %% Calculate TT fixed rank
    fprintf('\tRunning Oseledets TT-toolbox with fixed rank. ...'); t00 = tic;
    G_est = tt_tensor_fixed_rank(X,D);
    [e_rmse, e_rmse_clean, e_numel, e_tt_comp] = calculateErrorStuff(X, X_clean, G_est, G_true);
    all_rmse_noise(snr_i, method_i) = e_rmse;
    all_rmse_true(snr_i, method_i) = e_rmse_clean;
    all_tt_comp(snr_i, method_i) = e_tt_comp;
    all_numel(snr_i, method_i) = e_numel;
    method_i = method_i +1;
    toc(t00);
    
    
    %% Calculate tt_prob
    G_best = [];
    elbo_best = -inf;
    found_a_solution = false;
    fprintf('\tRunning Probabilistic TT-toolbox\n')
    for j = 1:num_repeats
        fprintf('\t\t Repeat %i of %i ...', j, num_repeats); t00 = tic;
        try
            elbo=[];
            [G_est, S_est, V_est, tau_est, elbo] = tt_prob_tensor(X, [], D,...
                'maxiter',maxiter,'verbose','no');

            if ~isempty(elbo) && elbo(end) > elbo_best
                elbo_best = elbo(end);
                G_best = G_est;
                found_a_solution = true;
            end
        catch e
            warning(sprintf('Something went wrong... Error message was:\n%s\n',e.message))
        end
        toc(t00);
    end
    if found_a_solution
        [e_rmse, e_rmse_clean, e_numel, e_tt_comp] = calculateErrorStuff(X, X_clean, G_best, G_true);
        all_rmse_noise(snr_i, method_i) = e_rmse;
        all_rmse_true(snr_i, method_i) = e_rmse_clean;
        all_elbo(snr_i, method_i) = elbo_best;
        all_tt_comp(snr_i, method_i) = e_tt_comp;
        all_numel(snr_i, method_i) = e_numel;
    end
    method_i = method_i +1;
    
    
end

save('./eusipco/synthetic_tensortrain_vs_probtt','all_rmse_noise', 'all_rmse_true',...
    'all_elbo', 'all_tt_comp', 'all_numel','snr_list','N','D','tt_threshold','name_methods')
% save('./eusipco/synthetic_tensortrain_cleandata', 'final_elbo', 'final_rrmse', 'final_tau',...
%     'final_tt_recon', 'final_tt_rank','final_numel_G')
end

function [e_rmse, e_rmse_clean, e_numel, e_tt_comp] = calculateErrorStuff(X_noise, X_clean, G_est, G_true)
% Calculates relevant error measures and solution statistics.
    if nargin < 4
        G_true = [];
    end
    
    X_recon = constructTensorTrain(G_est);
    
    e_rmse = norm(X_noise(:)-X_recon(:),'fro')^2/norm(X_noise(:),'fro')^2;
    e_rmse_clean = norm(X_clean(:)-X_recon(:),'fro')^2/norm(X_clean(:),'fro')^2;
    e_numel = sum(cellfun(@numel, G_est));
    
    if ~isempty(G_true)
        e_tt_comp = tt_compare(G_true, G_est);
    else
        e_tt_comp = nan;
    end
end