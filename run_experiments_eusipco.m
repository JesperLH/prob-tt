%% Run Experiments for EUSIPCO
addpath(genpath('./'))
addpath(genpath('../ncptensor'))
addpath(genpath('../tools'))

elapsed_time = zeros(6,1)*nan;
for j = 2:2
    t0 = tic;
    switch j
        case 1
            experiments_aminoacid
        case 2
            perm_idx = perms(1:3);
            for p = 1:size(perm_idx,1)
                experiments_aminoacid_modelorder(perm_idx(p,:))
            end
        case 3
            experiments_eusipco_tensortrain_nonoise
        case 4
            experiments_eusipco_tensortrain_amino
        case 5
            experiments_eusipco_knownD
        case 6
            experiments_eusipco_unknownD
    end
    elapsed_time(j) = toc(t0);
end
