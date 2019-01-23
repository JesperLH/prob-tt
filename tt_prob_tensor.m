function [G, ES, EV, Etau, elbo] = tt_prob_tensor(X,G, D, varargin)
%%TT_PROB_TENSOR Finds a probabilistic Tensor Train for Data X
%% Read input and Setup Parameters
paramNames = {'conv_crit', 'maxiter', ...
    'model_tau', 'fixed_tau', 'tau_a0', 'tau_b0',...
    'model_S', 'fixed_S',...
    'model_V', 'fixed_V',...
    'model_G', 'fixed_G','verbose'};
defaults = {1e-8, 100, ...
    true, 5, 1e-3, 1e-3,...
    true, 0,...
    true, 0,...
    true, 0,...
	'yes'};
% Initialize variables to default value or value given in varargin
[conv_crit, max_iter, ...
    model_tau, fixed_tau, tau_alpha0, tau_beta0,...
    model_S, fixed_S, model_V, fixed_V, model_G, fixed_G, verbose]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});

%% Check D is valid
N = size(X);
if length(D) == ndims(X)
    if D(1) == 1
        D = [D,1];
    elseif D(end) == 1
        D = [1,D];
    else
        error('Input error: specified model order not understood')
    end
    
elseif length(D) == ndims(X)-1
    D = [1,D,1];
elseif length(D) == ndims(X)+1
    assert(D(1)==1 && D(end)==1, 'Input error, D(1) and D(end) must be 1')
end

assert(all(D(2:end-1)<=N(1:end-1)), 'Model order cannot be higher than mode observations.')
assert(D(end-1) <= N(end), 'Model order in last mode cannot be higher than mode end or end-1')
%% Initialize
SST = sum(sum(X(:).^2));
prior_s = ones(D(end-1),1)*sqrt(SST/numel(X));
s_prior_log_value = ones(D(end-1),1)*0;
s_lb = zeros(D(end-1),1);
s_ub = ones(D(end-1),1)*inf;


if ~isempty(G)
    G_sizes = cellfun(@size, G,  'UniformOutput', false);
    D = ones(ndims(X)+1,1);
    for d = 1:ndims(X)-1
        D(d+1) = G_sizes{d}(end);
    end
    
    EV = G{end};
    ES = eye(size(G{end},1));
else
    
    G=cell(length(N),1);
    for d = 1:length(N)-1
        G{d} = randn(D(d),N(d),D(d+1));
        G{d} = reshape(orth(reshape(G{d}, D(d)*N(d), D(d+1))), ...
            D(d),N(d),D(d+1));
    end
    
    % Last mode
    ES = diag(prior_s);
    ES = sqrt(norm(X(:),'fro')^2/length(prior_s))*eye(length(prior_s));
    %ES = diag(rand(D(end-1),1)*2+0.5);
    EV = orth(randn(N(end), D(end-1)));
    G{end} = ES*EV';
    
%     for d = 1:length(N)
%         G{d} = G{d} + randn(size(G{d}))*1e-2*sqrt(var(G{d}(:)));
%     end
    
end

Etau=tau_alpha0/tau_beta0;%*1/mean(X(:).^2);

if sum(cellfun(@numel, G)) >= numel(X)
    warning('The number of elements in G is larger than X. This does not provide compression of X and may lead to unexpected performance')
end

%%
%% Setup how information should be displayed.
if ~strcmpi(verbose,'no')
    disp([' '])
    disp(['Probabilistic Tensor Train'])        %TODO: Better information about what is running..
    disp(['A D=' strcat('(',strjoin(strsplit(num2str(D(:)')),', '), ')') ' component model will be fitted']);
    disp([' ']);
    dheader = sprintf(' %16s | %16s | %16s | %16s | %16s | %16s | %16s |','Iteration', ...
        'Cost', 'rel. \Delta cost', 'Noise (s.d.)', 'Var. Expl.', 'Time (sec.)', 'Time CPU (sec.)');
    dline = repmat('------------------+',1,7);
    fprintf('%s\n%s\n%s\n',dline, dheader, dline);
end

%
iter = 0;
elbo = zeros(max_iter,1,'like', SST);
rmse = zeros(max_iter,1,'like', SST);
total_realtime = tic;
total_cputime = cputime;

%% Run algo
iter = 0;
dELBO = inf;


while iter < max_iter && dELBO > conv_crit || iter <= fixed_tau
    iter = iter+1;
    time_tic_toc = tic;
    time_cpu = cputime;
    %% Update each factor...
    if iter == 1 || model_G
        for i = 1:ndims(X)-1
%         for i = randperm(ndims(X)-1)
            
            % Contract modes
            sG = size(G{i});
            x_con = contract_fixed_modes(X,G,i);
            assert(sum(abs(x_con(:)))>0, 'No residual..')
            
            % Update von-Mises-Fisher Distribution
            F=reshape(x_con, prod(sG(1:end-1)), sG(end))*Etau;
            
            % F is obs x comp, where it is desired that comp x comp is
            % orthogonal
            [UU,SS,VV]=svd(F, 'econ');
            [f,V,lF]=hyperg(size(F,1),diag(SS),3);
%             assert(max(f)>eps, 'Subspace is pruned...')
            
            % Get expected value of G{i} and the entropy of vMF
            G{i}=reshape(UU*diag(f)*VV', sG);
            H_G(i)= lF-sum(sum(F(:).*G{i}(:)));  %#ok<AGROW>
            
            assert(sum(abs(G{i}(:)))>0, 'An entire cart was turned off')
        end
    end
    %% Update the last factor
    x_con = contract_fixed_modes(X,G,ndims(X));
    
    % Update V
    if iter == 1 || model_V
        % Update von-Mises-Fisher Distribution
        F=(ES*x_con)'*Etau;
        [UU,SS,VV]=svd(F, 'econ');
        [f,V,lF]=hyperg(N(end),diag(SS),3);
        assert(~isinf(lF))
        
        % Get expected value of G{i} and the entropy of vMF
        EV=(UU*diag(f)*VV');
        H_V= lF-sum(sum(F.*EV));
        P_V = 0;
    end
    
    % Update S
    if iter == 1 || model_S
        
        ss_sig2 = 1./(1*Etau + prior_s);
        ss_mu = ss_sig2.*sum(EV'.*x_con, 2)*Etau;
                
        [logZhatOut, ~ , s_mu, s_var ] = ...
            truncNormMoments_matrix_fun(s_lb, s_ub, ss_mu , ss_sig2);
        
        N_s = N(end);
        
        [P_S, H_S] = getPriorAndEntropy_tnorm(s_mu, s_var, logZhatOut, N_s, ...
            s_prior_log_value, prior_s);
        
        ES = diag(s_mu);
        ESS = sum(s_mu.^2 + s_var);
    end
    G{end} = ES*EV';
    
    %% Update the noise
    x_con = contract_fixed_modes(X,G,ndims(X));
    
    SSE = (SST+ESS-2*sum(sum(x_con .* G{end})));
    assert(SSE>=0, 'SSE was negative!')
    
    if model_tau && iter > fixed_tau
        est_alpha = tau_alpha0 + numel(X)/2;
        est_beta = tau_beta0 + SSE/2;
    else
        est_alpha = tau_alpha0;
        est_beta = tau_beta0;
    end
    Etau = est_alpha/est_beta;
    Elog_tau = psi(est_alpha)-log(est_beta);
    H_tau = gamma_entropy(est_alpha, est_beta);
    P_tau = -gammaln(tau_alpha0)+tau_alpha0*log(tau_beta0)...
        +(tau_alpha0-1)*Elog_tau - tau_beta0*Etau;
    
    assert(Etau>0, 'Noise precision was negative!')
    assert(isreal(Elog_tau), 'Not real!')
    
    %% Calculate Elbo
    elbo(iter) = 0.5*numel(X)*(-log(2*pi)+Elog_tau)-0.5*SSE*Etau...
        +P_tau+H_tau...
        +sum(H_G)... %Note, P_G is 0, when P(G) is a uniform prior on the sphere.
        +P_V+H_V...
        +P_S+sum(H_S);
  
    rmse(iter) = 0.5*SSE/numel(X);
    %% Display
    if mod(iter,50)==0 && ~strcmpi(verbose,'no')
        fprintf('%s\n%s\n%s\n',dline, dheader, dline);
    end
    
    time_tic_toc = toc(time_tic_toc);
    time_cpu = cputime-time_cpu;
    if iter > 1
        dCost = (rmse(iter)-rmse(iter-1))/abs(rmse(iter-1)); % Should decrease
        dELBO = (elbo(iter)-elbo(iter-1))/abs(elbo(iter-1)); % Shoud increase

        if ~strcmpi(verbose,'no')
            fprintf(' %16i | %16.4e | %16.4e | %16.4e | %16.4f | %16.4f | %16.4f |\n',...
                iter, elbo(iter), dELBO, sqrt(1./Etau), 1-SSE/SST, time_tic_toc, time_cpu);
        end
        
        if abs(dELBO) > conv_crit
           assert(dELBO>=0, 'ELBO DECREASED!')
        elseif dELBO < 0 
            warning('ELBO decreased, but was below convergence threshold')
        end
            
        
    elseif ~strcmpi(verbose,'no')
        fprintf(' %16i | %16.4e | %16.4e | %16.4e | %16.4f | %16.4f | %16.4f |\n',...
            iter, elbo(iter), nan, sqrt(1./Etau), 1-SSE/SST, time_tic_toc, time_cpu);
    end
   
    
    
end

if abs((2*rmse(iter)/SST)-1) < eps
    warning('Algorithm converged, but the Tensor Train does not describe the data!!!!!')
end

elbo = elbo(1:iter);
rmse = rmse(1:iter);

end

function X_contr = contract_fixed_modes(X_contr,G, idx)

n_modes = ndims(X_contr);

for i = 1:n_modes
    if i == 1 && i < idx
        X_contr = tprod(X_contr, [-1, i+1:ndims(X_contr)], squeeze(G{i}), [-1, 1]);
    elseif i < idx
        X_contr = tprod(X_contr, [-1,-2, 2:ndims(X_contr)-1], G{i}, [-1, -2, 1]);
    elseif i > idx
        if i == idx+1
            if idx == 1
                X_contr = tprod(X_contr, [1,-1, 4:ndims(X_contr)+1], G{i}, [2, -1, 3]);
            else
                X_contr = tprod(X_contr, [1,2,-1, 5:ndims(X_contr)+1], G{i}, [3, -1, 4]);
            end
        else
            if idx == 1
                X_contr = tprod(X_contr, [1,2,-1,-2, 4:ndims(X_contr)-1], G{i}, [-1, -2, 3]);
            else
                X_contr = tprod(X_contr, [1,2,3,-1,-2, 5:ndims(X_contr)-1], G{i}, [-1, -2, 4]);
            end
        end
    else
        
    end
%     if idx ~=1 || i > 1
%         X_contr=X_contr;
%     end
end

end


function entropy = gamma_entropy(a,b)
if size(b,2) > 1
    entropy = bsxfun(@minus,(gammaln(a)-(a-1).*psi(a)+a)',log(b));
else
    entropy = gammaln(a)-(a-1).*psi(a)-log(b)+a;
end
%entropy = gammaln(a)-a*psi(a)-log(b)+a;
end

function [prio, entr] = getPriorAndEntropy_tnorm(mu, sig2, logZhatOut, N, ...
    prior_log_value, prior_value)

lowerbound = 0;
upperbound = inf;

log_psi_func = @(t) -0.5*t.^2-0.5*log(2*pi);

sig = sqrt(sig2);
alpha=(lowerbound-mu)./sig;
beta=(upperbound-mu)./sig;

if upperbound==inf
    %                 if isa(alpha,'gpuArray')
    %                     r=exp(log(abs(alpha))+log_psi_func(alpha)-logZhatOUT);
    %                 else
    r=real(exp(log(alpha)+log_psi_func(alpha)-logZhatOut));
    %                 end
end
entr = log(sqrt(2*pi*exp(1)))+log(sig)+logZhatOut+0.5*r;

assert(~any(isnan(entr(:))),'Entropy was NaN')
assert(~any(isinf(entr(:))),'Entropy was Inf')

% Gets the prior contribution to the cost function.
% Note. The hyperparameter needs to be updated before calculating
% this.
prio = N*(-log(1/2)-1/2*log(2*pi))...
    +N/2*sum(prior_log_value(:))...
    -1/2*sum(sum(bsxfun(@times, prior_value , mu.^2 + sig2)));


end
