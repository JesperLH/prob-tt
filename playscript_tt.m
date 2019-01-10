
addpath('./requiredfunctions/')

%%
N = 15:-1:11;
D = [1,length(N)+1:-1:1+1];
D_est = D;
%D(2:end-1) = D(2:end-1) +1;
% D_est(4)=3;
%D_est(2:3) = D_est(2:3)-1;

[X, G, ES, EV] = generateTensorTrain(N, D);

X=addTensorNoise(X/sqrt(var(X(:))), 10); % Added 5 dB noise

%%
clc
for model_G = [false, true]
    
    for model_S = [false, true]
        for model_V = [false, true]
            %model_V = true;%model_S;
            for model_tau = [false, true]
                model_G = true;
                model_S = true;
                model_tau = true;
                model_V = true;
                
                s_ = {'tau', 'S', 'V', 'G'};
                idx = [model_tau, model_S, model_V, model_G];
                fprintf('\n\nModeling.. \n%s\n\n', strjoin(s_(idx),' and '))
                
                [G_est, S_est, V_est, tau_est] = tt_prob_tensor(X, [], D_est,...
                    'model_tau', model_tau,...
                    'model_S', model_S,...
                    'model_V', model_V,...
                    'model_G', model_G);
                
                
                if all(D==D_est)
                    tt_compare(G, G_est)
                end
            end
        end
    end
end

