function [X, G, S, V] = generateTensorTrain(N, D)
%% GENERATETENSORTRAIN generates a tensor train with 
    G=cell(length(N),1);
    for d = 1:length(N)-1
        G{d} = randn(D(d),N(d),D(d+1));
        G{d} = reshape(orth(reshape(G{d}, D(d)*N(d), D(d+1))), ...
            D(d),N(d),D(d+1));
    end
    
    % Last mode
    S = diag(rand(D(end-1),1)*2+0.5);
    V = orth(randn(D(end-1), N(end))');
    G{end} = S*V';
    
    % Construct data
    X = constructTensorTrain(G);
    
end