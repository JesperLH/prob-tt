%% Run Experiments for EUSIPCO-2019
% Runs the experiments described in [1].
%
% [1] Hinrich, J. L. and Mørup, M., Probabilistic Tensor Train Decomposition.
% https://github.com/JesperLH/prob-tt
%

addpath(genpath('./requiredfunctions'))
%addpath(genpath('../ncptensor'))
%addpath(genpath('../tools'))

elapsed_time = zeros(4,1)*nan;
for j = 1:4
    t0_ = tic;
    switch j
        case 1
            experiments_eusipco_knownD
        case 2
            experiments_eusipco_unknownD
        case 3
            experiments_tensortrain_vs_prob_tensor
        case 4
            perm_idx = perms(1:3);
            for p = 1:size(perm_idx,1)
                experiments_aminoacid_modelorder(perm_idx(p,:))
            end    
    end
    elapsed_time(j) = toc(t0_);
end

%%
return
%% Visualization (requires display).

% Figure 1
visualize_exp_synthetictensor

% Figure 2
visualize_exp_probTTvsOseTT

% Table 1
visualize_exp_aminoorder

% Figure 3
experiments_aminoacid 


