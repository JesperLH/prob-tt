function experiments_aminoacid_modelorder(permute_order)

if nargin < 1
    permute_order = 1:3;
end

addpath('../thirdparty-matlab/nway331/')
load('../Data/Amino-Acid/claus.mat')
scale_x = sqrt(var((X(:))));
X=permute(X,permute_order);

X = X/scale_x;

D_1mode = 2:max(size(X,1));
D_2mode = 2:2:min(size(X,2), size(X,3));
n_repeats = 5;
constr = [0,0,0]; %Unconstrained

final_elbo = zeros(length(D_1mode), length(D_2mode), n_repeats)*nan;
final_rmse = zeros(length(D_1mode), length(D_2mode), n_repeats)*nan;
final_tau = zeros(length(D_1mode), length(D_2mode), n_repeats)*nan;
final_cp = zeros(length(D_1mode), length(D_2mode), n_repeats)*nan; 

for d1 = 1:length(D_1mode)
    for d2 = 1:length(D_2mode)
        D_est = [1, D_1mode(d1), D_2mode(d2) 1];
        t0 = tic;
        for rep = 1:n_repeats
            try 
                [G,S,V, tau, elbo] = tt_prob_tensor(X, [], D_est,'conv_crit',1e-6);
                X_tt = constructTensorTrain(G);
            
                final_elbo(d1,d2,rep) = elbo(end);
                final_rmse(d1,d2,rep) = sqrt(mean((X(:)-X_tt(:)).^2));
                final_tau(d1,d2,rep) = tau;

                [~,A_clean] = evalc('parafac(constructTensorTrain(G), 3, [], constr)');
                [avg_est_clean, idx_clean] = optimal_component_match(abs(corr(A_clean{find(size(X)==5)},y)));
                final_cp(d1,d2,rep) = avg_est_clean;
            catch
               warning('Something did not go as expected..')
            end
            
        end
        toc(t0)
        
    end
end

save(['./eusipco/amino_acid_modelorder', strrep(num2str(permute_order),'  ','_') '.mat'],...
    'D_1mode', 'D_2mode', 'n_repeats', 'final_elbo', 'final_rmse', 'final_tau', 'final_cp')
%%
% figure; 
% imagesc(max(final_cp,[],3)');  colorbar
% title('ELBO for different model orders (higher is better')
% set(gca,'XTick',1:length(D_1mode), 'XTickLabel', D_1mode);
% xlabel('#Components in mode 1')
% set(gca,'YTick',1:length(D_2mode), 'YTickLabel', D_2mode);
% ylabel('#Components in mode 2')
%ind2sub(size(final_elbo), idx_lin)
%[x,y,z] = ind2sub(size(final_elbo), idx_lin)