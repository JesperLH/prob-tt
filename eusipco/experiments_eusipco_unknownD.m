%% Tensor Train Experiments
% Known Dmax and (un)known model order

addpath('./requiredfunctions/')
addpath('../thirdparty-matlab/oseledets-TT-Toolbox-8332a6c/')

rng(131462234)
%% Setup problem and generate data
%snr_db = -5;
list_snrdb = [-10:2.5:10];
num_repeats = 10;
maxiter = 50;

N = 20:-1:16;
D = [1,length(N)+1:-1:3, 1];
D_est = [1, repmat(max(D),1,length(N)-1),1];

[X_clean, G, ES, EV] = generateTensorTrain(N, D);
X_clean = X_clean/sqrt(var(X_clean(:)));

% All permutations
all_permuations = perms(1:ndims(X_clean));
all_permuations = all_permuations(end:-1:1,:); %Make sure the first one is 1:5
num_permutes = size(all_permuations,1);


model_G = true;
model_S = true;
model_tau = true;
model_V = true;

%%
num_noise_level = length(list_snrdb);
%% Run experiments.
final_elbo = nan(num_permutes,num_repeats, num_noise_level);
final_rrmse = nan(num_permutes,num_repeats, num_noise_level);
final_tau = nan(num_permutes,num_repeats, num_noise_level);

for snr = 1:num_noise_level
    
    X=addTensorNoise(X_clean, list_snrdb(snr)); % Added 5 dB noise
    
    for i = 1:num_permutes
        t0 = tic;
        fprintf('SNR=%4.2f\t Permutation %i of %i...', list_snrdb(snr), i, num_permutes)
        perm_idx = all_permuations(i,:);
        
        org_idx = 1:ndims(X);
        perm2org_idx = sum(bsxfun(@times,(org_idx == perm_idx')', org_idx),2);

        D_est_perm = D_est;%[1, D_est(perm_idx(1:end-1)+1), 1];
        for j = 1:num_repeats

            try
                [G_est, S_est, V_est, tau_est, elbo] = tt_prob_tensor(permute(X,perm_idx), [], D_est_perm,...
                    'model_tau', model_tau,...
                    'model_S', model_S,...
                    'model_V', model_V,...
                    'model_G', model_G, 'maxiter',maxiter,'verbose','no');


                %G_est = G_est(perm2org_idx);
                X_recon = constructTensorTrain(G_est);
                X_recon = permute(X_recon, perm2org_idx);
%                 fprintf('Reconstruction error is %6.4f\n',norm(X(:)-X_recon(:),'fro')^2/norm(X(:),'fro')^2)

                final_rrmse(i,j,snr) = norm(X(:)-X_recon(:),'fro')^2/norm(X(:),'fro')^2;
                final_elbo(i,j,snr) = elbo(end);
                final_tau(i,j,snr) = tau_est;
            catch e
               warning(sprintf('Something went wrong... Error message was:\n%s\n',e.message)) 
            end
        end
        toc(t0)
    end
end
save('./eusipco/synthetic_unknownD', 'final_elbo', 'final_rrmse', 'final_tau',...
    'all_permuations', 'list_snrdb')
