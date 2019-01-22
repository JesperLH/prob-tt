%% Tensor Train Experiments
% Unknown D (but known max(D) and (un)known model order

addpath('./requiredfunctions/')
addpath('../thirdparty-matlab/oseledets-TT-Toolbox-8332a6c/')

rng(131462234)
%% Setup problem and generate data
snr_db = 5;
num_repeats = 4;
num_permutes = 5;
maxiter = 5;
N = 15:-1:11;
D = [1,length(N)+1:-1:1+1];
D_est = [1,repmat(max(D),1,ndims(X))];
% D_est(2:end) = D(2:end)*2;
% D_est(2) = D_est(2)-1;

[X, G, ES, EV] = generateTensorTrain(N, D);
X=addTensorNoise(X/sqrt(var(X(:))), snr_db); % Added 5 dB noise

model_G = true;
model_S = true;
model_tau = true;
model_V = true;

%%


%% Run experiments
% - TODO loop over different SNR ..
final_elbo = zeros(num_permutes,num_repeats);
final_rrmse = zeros(num_permutes,num_repeats);
final_tau = zeros(num_permutes,num_repeats);
final_tt_comp = nan(num_permutes,num_repeats);
final_permutation = zeros(num_permutes, ndims(X));

for i = 1:num_permutes
    if i == 1
        perm_idx = 1:ndims(X);
%         G_init = G;
%     elseif i==2
%         perm_idx = 1:ndims(X);
%         G_init = [];
%     else
        perm_idx = randperm(ndims(X));
        G_init = [];
    end
    final_permutation(i,:) = perm_idx;
    org_idx = 1:ndims(X);
    perm2org_idx = sum(bsxfun(@times,(org_idx == perm_idx')', 1:ndims(X)),2);
    
    D_est_perm = [1, D_est((perm_idx)+1)];
    for j = 1:num_repeats
    
    
        [G_est, S_est, V_est, tau_est, elbo] = tt_prob_tensor(permute(X,perm_idx), G_init, D_est_perm,...
            'model_tau', model_tau,...
            'model_S', model_S,...
            'model_V', model_V,...
            'model_G', model_G, 'maxiter',maxiter,'verbose','no');


        %G_est = G_est(perm2org_idx);

        if all(D==D_est_perm)
             final_tt_comp(i,j) = tt_compare(G, G_est);
        end
        X_recon = constructTensorTrain(G_est);
        X_recon = permute(X_recon, perm2org_idx);
        fprintf('Reconstruction error is %6.4f\n',norm(X(:)-X_recon(:),'fro')^2/norm(X(:),'fro')^2)

        final_rrmse(i,j) = norm(X(:)-X_recon(:),'fro')^2/norm(X(:),'fro')^2;
        final_elbo(i,j) = elbo(end);
        final_tau(i,j) = tau_est;
    end
end

save('./eusipco/synthetic_knownD', 'final_elbo', 'final_rrmse', 'final_tau', 'final_tt_comp', 'final_permutation')

%% Analyse results
[min_elbo, min_elbo_idx] = min(final_elbo, [], 1)
[min_rmse, min_rmse_idx] = min(final_rrmse, [], 1)
[min_tau, min_tau_idx] = min(final_tau, [], 1)
