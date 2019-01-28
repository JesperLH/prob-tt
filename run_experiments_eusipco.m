%% Run Experiments for EUSIPCO
addpath(genpath('./'))
addpath(genpath('../ncptensor'))
addpath(genpath('../tools'))

elapsed_time = zeros(6,1)*nan;
for i = 5:6
    t0 = tic;
    switch i
        case 1
            experiments_aminoacid
        case 2
            experiments_aminoacid_modelorder
        case 3
            experiments_eusipco_tensortrain_nonoise
        case 4
            experiments_eusipco_tensortrain_amino
        case 5
            experiments_eusipco_knownD
        case 6
            experiments_eusipco_unknownD
    end
    elapsed_time(i) = toc(t0);
end
