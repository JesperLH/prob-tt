%% Visualize EUSIPCO results

%% Amino Acid Model Order
clear
max_elbo = -inf; max_elbo_ord = nan;
max_cp = 0; max_cp_ord = nan;
min_rmse = inf; min_rmse_ord = nan;
max_tau = -inf; max_tau_ord = nan;

all_permutes = perms(1:3); all_permutes(end:-1:1, :);
for p = 1:size(all_permutes,1)
    order_idx = all_permutes(p,:);
    load(sprintf('./eusipco/amino_acid_modelorder%i_%i_%i.mat',order_idx))
    
    [max_elbo, id] = max([max_elbo, nanmax(final_elbo(:))]);
    if id > 1
        max_elbo_ord = p;
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
for p = 1:size(all_permutes,1)
    order_idx = all_permutes(p,:);
    load(sprintf('./eusipco/amino_acid_modelorder%i_%i_%i.mat',order_idx))
    
    for i = 1:4
        
        if i == 1
            img = final_elbo/max_elbo;
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
        else
            best_img = nanmax(img,[],3);
            [x_,y_] = ind2sub(size(best_img), find(best_img == max(best_img(:))));
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
end