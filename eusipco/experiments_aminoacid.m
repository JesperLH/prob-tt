function experiments_aminoacid()
addpath('../thirdparty-matlab/nway331/')
load('../Data/Amino-Acid/claus.mat')
scale_x = sqrt(var((X(:))));
X = X/scale_x;

%Best solution
D_est = [1, 5, 5, 1];
perm_idx = [2,3,1];
reverse_perm_idx = [3,1,2];

%Worst solution
% D_est = [1, 3, 5, 1];
% perm_idx = [2,3,1];
% reverse_perm_idx = [3,1,2];


n_repeats = 50;
best_model = [];
final_elbo = zeros(n_repeats,1)*-inf;

for rep = 1:n_repeats
    try
        [G,S,V, tau, elbo] = tt_prob_tensor(permute(X,perm_idx), [], D_est,'conv_crit',1e-6);
    catch e
        elbo = -inf;
        G=[];
        warning(e.message)
%         [G,S,V, tau, elbo] = tt_prob_tensor(X, [], D_est,'conv_crit',1e-6);
    end
    final_elbo(rep) = elbo(end);
    
    [~, idx] = max(final_elbo);
    if idx == rep
        best_model = G;
    end
end
G = best_model;
figure; bar(final_elbo)

X_recon = permute(constructTensorTrain(G), reverse_perm_idx);

plotDataCleanDataAndDiff(X*scale_x, X_recon*scale_x)


constr = [0,0,0]; %Unconstrained

% A_raw = parafac(X, 3, [], constr);
% A_nonneg = parafac(X, 3, [], constr+2);
% A_clean = parafac(constructTensorTrain(G), 3, [], constr);
[~,A_raw] = evalc('parafac(X, 3, [], constr)');
[~,A_nonneg] = evalc('parafac(X, 3, [], [2,2,2])');
[~,A_clean] = evalc('parafac(X_recon, 3, [], constr)');

[avg_est_raw, idx_raw] = optimal_component_match(abs(corr(A_raw{1},y)));
[avg_est_nonneg, idx_nonneg] = optimal_component_match(abs(corr(A_nonneg{1},y)));
[avg_est_clean, idx_clean] = optimal_component_match(abs(corr(A_clean{1},y)));

avg_est_raw, avg_est_nonneg, avg_est_clean

% Get worst solution
%X_recon_worst=getWorstSolution(X, n_repeats);
%plotDataCleanDataAndDiffBestWorst(X*scale_x, X_recon*scale_x, X_recon_worst*scale_x)
plotDataCleanDataAndDiffBestWorst(X*scale_x, X_recon*scale_x, nmodel(A_nonneg)*scale_x)



end

function plotAminoAcidData(X)%, Xoverlay)
figure('units', 'normalized', 'position', [0.1,0.3,0.7,0.3])
N=size(X,1);
for i = 1:N
    %subplot(2,ceil(N/2),i)
    subplot(1,N,i)
    surf(squeeze(X(i,:,:)))%,squeeze(X(i,:,:)))
    axis tight
    grid off
    set(gca,'color','none')
    set(gca,'YTick',[], 'XTick', [])%, 'ZTick',[])
    set(gca,'view', [-37.5, 30])
    colormap parula
    shading flat% interp
    title(sprintf('Sample %i', i))
    xlabel('Excitation')
    ylabel('Emission')
    
    v = get(gca,'view');
    
    xh = get(gca,'XLabel'); % Handle of the x label
    set(xh, 'Units', 'Normalized')
    set(xh, 'Position',[0.57, 0.03,0],'Rotation',v(2)+5)
    
    yh = get(gca,'YLabel'); % Handle of the y label
    set(yh, 'Units', 'Normalized')
    set(yh, 'Position',[0.35,0.03, 0],'Rotation',v(1)-9)


    
end


end

function plotDataCleanDataAndDiff(X, X_clean)%, Xoverlay)
figure('units', 'normalized', 'position', [0.05,0.05,0.7,0.25*3])
N=size(X,1);
for j = 1:3
    if j == 1
        A = X;
    elseif j==2
        A = X_clean;
    elseif j==3
        A = X-X_clean;
    end

    for i = 1:N
        %subplot(2,ceil(N/2),i)
        subplot(3,N,i+N*(j-1))
        surf(squeeze(A(i,:,:)))%,squeeze(X(i,:,:)))
        axis tight
        grid off
        set(gca,'color','none')
        set(gca,'YTick',[], 'XTick', [])%, 'ZTick',[])
        set(gca,'view', [-37.5, 30])
        colormap parula
        shading flat% interp
        if j == 1
            title(sprintf('Sample %i', i))
        end
        xlabel('Excitation')
        ylabel('Emission')
        

        v = get(gca,'view');

        xh = get(gca,'XLabel'); % Handle of the x label
        set(xh, 'Units', 'Normalized')
        pos = get(xh, 'Position');
        %set(xh, 'Position',pos.*[1,1,1],'Rotation',v(2)+8)
        set(xh, 'Position',[0.57, 0.03,0],'Rotation',v(2)-5)
        yh = get(gca,'YLabel'); % Handle of the y label
        set(yh, 'Units', 'Normalized')
        set(yh, 'Position',[0.375,-0.01, 0],'Rotation',v(1))
        
        if i == 1
            zlabel('Intensity')
            zh = get(gca,'ZLabel'); % Handle of the z label
%             set(zh, 'Units', 'Normalized')
%             set(zh, 'Position',[0.35,0.03, 0])%,'Rotation',v(1)-9)
        end
        
        if j>1
            pos = get(gca, 'Position');
            if j == 2
                pos(2) = 0.49;
                set(gca,'Position',pos)%.*[1,1.2,1,1])
            elseif j==3
                pos(2) = 0.25;
                set(gca,'Position',pos)%.*[1,1.4,1,1])
            end
        end
    end
end
%savefig('./eusipco/results/aminoacid','-deps3')
print('./eusipco/results/aminoacid','-dpng')
print('./eusipco/results/aminoacid','-depsc')
end


function X_recon=getWorstSolution(X, n_repeats)
    %Worst solution
    D_est = [1, 3, 5, 1];
    perm_idx = [2,3,1];
    reverse_perm_idx = [3,1,2];

    worst_model = [];
    final_elbo = zeros(n_repeats,1)*inf;

    for rep = 1:n_repeats
        try
            [G,S,V, tau, elbo] = tt_prob_tensor(permute(X,perm_idx), [], D_est,'conv_crit',1e-6, 'verbose','no');
        catch
    %         [G,S,V, tau, elbo] = tt_prob_tensor(X, [], D_est,'conv_crit',1e-6);
        end
        final_elbo(rep) = elbo(end);

        [~, idx] = min(final_elbo);
        if idx == rep
            worst_model = G;
        end
    end
    G = worst_model;
    X_recon = permute(constructTensorTrain(G), reverse_perm_idx);
    
end

function plotDataCleanDataAndDiffBestWorst(X, X_best, X_worst)%, Xoverlay)
figure('units', 'normalized', 'position', [0.05,0.05,0.7,0.25*3])
N=size(X,1);
for j = 1:5
    if j == 1
        A = X;
    elseif j==2
        A = X_best;
    elseif j==3
        A = X-X_best;
    elseif j==4
        A = X_worst;
    elseif j==5
        A = X-X_worst;
    end

    for i = 1:N
        %subplot(2,ceil(N/2),i)
        subplot(5,N,i+N*(j-1))
        surf(squeeze(A(i,:,:)))%,squeeze(X(i,:,:)))
        axis tight
        grid off
        set(gca,'color','none')
        set(gca,'YTick',[], 'XTick', [])%, 'ZTick',[])
        set(gca,'view', [-37.5, 30])
        colormap parula
        shading flat% interp
        if j == 1
            title(sprintf('Sample %i', i))
        end
        
        if j == 5
            xlabel('Excitation')
            ylabel('Emmission')


            v = get(gca,'view');

            xh = get(gca,'XLabel'); % Handle of the x label
            set(xh, 'Units', 'Normalized')
            pos = get(xh, 'Position');
            set(xh, 'Position',[0.53, 0.01,0],'Rotation',v(2)-14)
            yh = get(gca,'YLabel'); % Handle of the y label
            set(yh, 'Units', 'Normalized')
            set(yh, 'Position',[0.425,-0.03, 0],'Rotation',v(1)+14)
        end
        
        if i == 1 
            if any(j == [1,2,4])
                zlabel('Intensity')
                zh = get(gca,'ZLabel'); % Handle of the z label
    %             set(zh, 'Units', 'Normalized')
    %             set(zh, 'Position',[0.35,0.03, 0])%,'Rotation',v(1)-9)
            else 
                zlabel('\Delta(Intensity)')
            end
        end
        
        if j>1
            no_label_offset = 0.02;
            pos = get(gca, 'Position');
            if j == 2
                pos(2) = 0.67+no_label_offset;%59;
            elseif j==3
                pos(2) = 0.53+no_label_offset*2;
            elseif j==4
                pos(2) = 0.39+no_label_offset*3;
            elseif j==5
                pos(2) = 0.25+no_label_offset*4;
            end
            set(gca,'Position',pos)
        end
    end
end
%savefig('./eusipco/results/aminoacid','-deps3')
print('./eusipco/results/aminoacid_bestworst','-dpng')
print('./eusipco/results/aminoacid_bestworst','-depsc')
end