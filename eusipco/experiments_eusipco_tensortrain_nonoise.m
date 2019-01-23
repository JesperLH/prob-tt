%% Tensor Train Experiments
addpath('./requiredfunctions/')
addpath(strjoin(strcat('../thirdparty-matlab/oseledets-TT-Toolbox-8332a6c/',{'','core', 'exp', 'cross', 'fmex', 'misc','solve'}),';'))

rng(131462234)
%% Setup problem and generate data
num_repeats = 10;
maxiter = 50;

N = 20:-1:16;
D = [1,length(N)+1:-1:3, 1];
D_est = D;

[X_clean, G, ES, EV] = generateTensorTrain(N, D);
X_clean = X_clean/sqrt(var(X_clean(:)));

% All permutations
all_permuations = perms(1:ndims(X_clean));
all_permuations = all_permuations(end:-1:1,:); %Make sure the first one is 1:5
num_permutes = size(all_permuations,1);

%%
num_noise_level = length(list_snrdb);
%% Run experiments
final_elbo = nan(num_permutes,num_repeats);
final_rrmse = nan(num_permutes,num_repeats);
final_tau = nan(num_permutes,num_repeats);
final_tt_recon = nan(num_permutes,1);
final_tt_rank = nan(num_permutes,length(D));
final_numel_G = nan(num_permutes,1);

X=X_clean;%
%X=addTensorNoise(X_clean, 100); % Added 5 dB noise

for i = 1:num_permutes
    t0 = tic;
    fprintf('Permutation %i of %i...', i, num_permutes)
    perm_idx = all_permuations(i,:);
    
    org_idx = 1:ndims(X);
    perm2org_idx = sum(bsxfun(@times,(org_idx == perm_idx')', org_idx),2);
    
    D_est_perm = [];
    
    %% Calculate TT
    tt = tt_tensor(permute(X,perm_idx), 1e-6);
    cr= tt.core ; ps= tt.ps ;
    for k=1:5
        core{k} =cr(ps(k): ps(k +1) -1); %#ok<SAGROW>
        core{k} = reshape(core{k}, tt.r(k), tt.n(k), tt.r(k+1)); %#ok<SAGROW>
    end
    G_init = core;
    
    X_recon = constructTensorTrain(G_init);
    X_recon = permute(X_recon, perm2org_idx);
    final_tt_recon(i) = norm(X(:)-X_recon(:),'fro')^2/norm(X(:),'fro')^2;
    final_tt_rank(i,:) = tt.r;
    final_numel_G(i) = sum(cellfun(@numel, G_init));
    %% Calculate prob from TT init
    for j = 1:num_repeats
        
        try
            [G_est, S_est, V_est, tau_est, elbo] = tt_prob_tensor(permute(X,perm_idx), G_init, [],...
                'maxiter',maxiter,'verbose','no', 'fixed_tau',0);
            
            X_recon = constructTensorTrain(G_est);
            X_recon = permute(X_recon, perm2org_idx);
            
            final_rrmse(i,j) = norm(X(:)-X_recon(:),'fro')^2/norm(X(:),'fro')^2;
            final_elbo(i,j) = elbo(end);
            final_tau(i,j) = tau_est;
        catch e
            warning(sprintf('Something went wrong... Error message was:\n%s\n',e.message))
        end
    end
    
    %%
    %%
    toc(t0)
end
save('./eusipco/synthetic_tensortrain_cleandata', 'final_elbo', 'final_rrmse', 'final_tau',...
    'final_tt_recon', 'final_tt_rank','final_numel_G')