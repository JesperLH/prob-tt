%% Visualize EUSIPCO results

%% Amino Acid Model Order
clear
load('../Data/Amino-Acid/claus.mat');
N = size(X);

max_elbo = -inf; max_elbo_ord = nan;
min_elbo = inf; min_elbo_ord = nan;
max_cp = 0; max_cp_ord = nan;
min_rmse = inf; min_rmse_ord = nan;
max_tau = -inf; max_tau_ord = nan;

all_permutes = perms(1:3); all_permutes = all_permutes(end:-1:1, :);
for p = 1:size(all_permutes,1)
    order_idx = all_permutes(p,:);
    load(sprintf('./eusipco/amino_acid_modelorder%i_%i_%i.mat',order_idx))
    
    [max_elbo, id] = max([max_elbo, nanmax(final_elbo(:))]);
    if id > 1
        max_elbo_ord = p;
    end
    [min_elbo, id] = min([min_elbo, nanmin(final_elbo(:))]);
    if id > 1
        min_elbo_ord = p;
    end
    [min_rmse, id] = min([min_rmse, nanmin(final_rmse(:))]);
    if id > 1
        min_rmse_ord = p;
    end
    
    [max_cp, id] = max([max_cp, nanmax(final_cp(:))]);
    if id > 1
        max_cp_ord = p;
    end
    
    [max_tau, id] = max([max_tau, nanmax(final_tau(:))]);
    if id > 1
        max_tau_ord = p;
    end
    
end


%%

figure('units','normalized','position',[0.1,0.0,0.8,0.9]);
fID = fopen('./eusipco/results/amino_table.tex','w');
% fprintf(fID, 'Mode Order\t & ELBO\t & D\\textsubscript{ELBO} \t & RMSE\t & D\\textsubscript{RMSE} \t & Avg. corr & $|\\cG|$\\\\\n');
fprintf(fID, 'Mode Order\t & ELBO\t & D\\textsubscript{ELBO} \t& $|\\cG_{ELBO}|$ & RMSE\t & D\\textsubscript{RMSE} \t &$|\\cG_{RMSE}|$\n')
for p = 1:size(all_permutes,1)
    order_idx = all_permutes(p,:);
    load(sprintf('./eusipco/amino_acid_modelorder%i_%i_%i.mat',order_idx))
    Nperm = N(order_idx);
    current_elbo_ord = [];
    current_elbo_corr = [];
    current_rmse_ord = [];
    current_rmse_corr = [];
    for i = 1:4
        
        if i == 1
            img = final_elbo/abs(max_elbo);
            s_title = 'ELBO';
        elseif i == 2
            img = final_rmse;
            s_title = 'RMSE';
        elseif i == 3
            img = final_cp;
            s_title = 'Avg. corr';
        elseif i == 4
            img = final_tau;
            s_title = 'Noise precision';
        end
        
        
        if i == 2
            best_img = nanmin(img,[],3);
            [x_,y_] = ind2sub(size(best_img), find(best_img == min(best_img(:))));
            current_rmse_ord = [1, D_1mode(x_), D_2mode(y_),1];
            current_rmse_corr = max(squeeze(final_cp(x_, y_,:)));
        else
            best_img = nanmax(img,[],3);
            [x_,y_] = ind2sub(size(best_img), find(best_img == max(best_img(:))));
            if i == 1
                current_elbo_ord = [1, D_1mode(x_), D_2mode(y_),1];
                current_elbo_corr = max(squeeze(final_cp(x_, y_,:)));
            end
        end
        
        if length(D_1mode) > 10
            if length(D_1mode) > 64
                skip_ = 40;
            else
                skip_ = 10;
            end
            y_tick = D_1mode(1:skip_:end);
            y_tick_pos = 1:skip_:max(D_1mode);
        else
            y_tick = D_1mode;
            y_tick_pos = 1:length(D_1mode);
        end
        
        if length(D_2mode) > 10
            if length(D_2mode) > 64
                skip_ = 20;
            else
                skip_ = 5;
            end
            x_tick = D_2mode(1:skip_:end);
            x_tick_pos = 1:skip_:max(D_2mode);
        else
            x_tick = D_2mode;
            x_tick_pos = 1:length(D_2mode);
        end
        
        subplot(6,4,i+(p-1)*4);
        imagesc(best_img); colorbar
        title(s_title);
        axis tight
        text(y_,x_,'*', 'FontSize',15)
        set(gca, 'XTick', x_tick_pos, 'YTick', y_tick_pos)
        set(gca, 'XTickLabel', x_tick, 'YTickLabel', y_tick)
        
%         subplot(2,4,i+4)
%         imagesc(nanvar(img,[],3)); colorbar
%         text(x_,y_,'*', 'FontSize',15)
%         set(gca, 'XTick', x_tick, 'YTick', y_tick)
%         set(gca, 'XTickLabel', x_tick, 'YTickLabel', y_tick)
        xlabel(sprintf('D=(1,%i,%i,1)',D_1mode(x_), D_2mode(y_)))
        if i == 1
            ylabel(sprintf('The original modes are %s',num2str(order_idx)),'Rotation',0);
        end
        
    end
    %suptitle(sprintf('The original modes are %s',num2str(order_idx)))
    %% Write to file
    %fprintf(fID, 'Mode Order\t & ELBO\t & D\textsubscript{ELBO} \t & RMSE\t & D\textsubscript{RMSE} \t & Avg. corr & $|\cG_{ELBO}|,|\cG_{RMSE}|$\n')
%     fprintf(fID, '(%i,%i,%i)\t & %6.4f\t & (%i,%i,%i,%i) \t & %6.4f\t & (%i,%i,%i,%i)\t & %6.4f,%6.4f & (%4i,%4i)\\\\\n',...
%         order_idx,...
%         max(final_elbo(:))/max_elbo,...
%         current_elbo_ord,...
%         min(final_rmse(:)),...
%         current_rmse_ord,...
%         [current_elbo_corr,current_rmse_corr],...
%         [sum(current_elbo_ord(2:3).*Nperm(1:2))+current_elbo_ord(3)*Nperm(3),...
%         sum(current_rmse_ord(2:3).*Nperm(1:2))+current_rmse_ord(3)*Nperm(3)]);
%fprintf(fID, 'Mode Order\t & ELBO\t & D\textsubscript{ELBO} \t $|\cG_{ELBO}|$ & RMSE\t & D\textsubscript{RMSE} \t $\cG_{RMSE}|$\n')
        fprintf(fID, '(%i,%i,%i)\t & %5.3f\t & (%i,%i) \t &\t%i & %5.3f\t & (%i,%i)\t & %i\\\\\n',...
        order_idx,...
        max(final_elbo(:))/max_elbo,...
        current_elbo_ord(2:3),...
        sum(current_elbo_ord(2:3).*Nperm(1:2))+current_elbo_ord(3)*Nperm(3),...
        min(final_rmse(:)),...
        current_rmse_ord(2:3),...
        sum(current_rmse_ord(2:3).*Nperm(1:2))+current_rmse_ord(3)*Nperm(3));
        
        
    
end
fclose(fID)

%% Find the worst solution.. (minimum elbo)
fprintf('\n\nThe worst solution is...\n')
min_elbo/abs(max_elbo)
load(sprintf('./eusipco/amino_acid_modelorder%i_%i_%i.mat',all_permutes(min_elbo_ord,:)))
worst_img = nanmin(final_elbo,[],3);
[x_,y_] = ind2sub(size(worst_img), find(worst_img == max(worst_img(:))));

fprintf('Mode order (%i,%i,%i)\n',all_permutes(min_elbo_ord,:))

fprintf('D=[%i,%i,%i,%i]\n',[1, D_1mode(x_), D_2mode(y_),1])