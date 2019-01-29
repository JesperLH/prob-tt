%% Visualize synthetic experiment..
% known has 121 perms, unknown has 120..
clear
is_known = false;
is_rmse = true;

for is_known = [true, false]
    for is_rmse = [true, false]
        
        if is_known
            load('C:\Users\jehi\Coding\prob-tt\eusipco\synthetic_knownD.mat')
        else
            load('C:\Users\jehi\Coding\prob-tt\eusipco\synthetic_unknownD.mat')
        end
        
        
        [n_perms, n_rep, n_snr] = size(final_elbo);
        %% ELBO over perms..
        idx_sort = 1:n_perms;
        %snr_colors = distinguishable_colors(9,{'w'});
        snr_colors = ones(9,3).*(fliplr(linspace(0.1,0.7,9)))';
        figure('units','normalized','position',[0.1,0.5,0.3,0.3])
        for i_snr = 1:1:9
            
            if is_rmse
                x = final_rrmse(idx_sort,:,i_snr);
                %         x_best = min(x,[],2);
            else
                x = final_elbo(idx_sort,:,i_snr);
                %         x_best = max(x,[],2);
            end
            mu = nanmean(x,2);
            sd = sqrt(nanvar(x,[],2));
            ci = [mu-10*sd, mu+10*sd]; % Note, there is almost no variation over repeats...
            
            %plotCI(1:n_perms, x_best, ci, snr_colors(i_snr,:), snr_colors(i_snr,:))
            plotCI(1:n_perms, mu, ci, snr_colors(i_snr,:), snr_colors(i_snr,:))
        end
        % Shrink figure and add x-axis
        pos = get(gca,'position');
        set(gca,'Position',pos.*[1,2.8,1,0.7])
        xlim([0,n_perms+0.5])
        xticks([1:5:n_perms])
        str_perms = strcat('(',strrep(string(num2str(all_permuations)),'  ',','),')');
        xticklabels(str_perms(1:5:n_perms))
        xtickangle(90)
        xlabel('All possible permutations')
        
        %Setup colorbar
        cbh = colorbar; colormap(snr_colors);
        c_min = min(snr_colors(:,1));
        c_max = max(snr_colors(:,1));
        caxis([c_min,c_max])
        set(cbh, 'Ticks', [c_min, c_max], 'TickLabels',{'Low SNR', 'High SNR'})
        
        % Fix ylabel and title
        if is_rmse
            ylabel('Relative RMSE')
            save_ext = 'rmse';
        else
            ylabel('Evidence Lowerbound (ELBO)')
            save_ext = 'elbo';
        end
        
        if is_known
            title('Known model order')
            save_pre = 'known';
        else
            title('Unknown model order')
            save_pre = 'unknown';
        end
        
        print(strcat('./eusipco/results/synthtensor_',save_pre,'_',save_ext), '-dpng')
        print(strcat('./eusipco/results/synthtensor_',save_pre,'_',save_ext), '-depsc')
        
    end
end

%% What about when the model order is known ? Look at TT-compare..
% well ELBO and RMSE are already plotted, but we can look at the average
% correlation between G_est and G_true
load('C:\Users\jehi\Coding\prob-tt\eusipco\synthetic_knownD.mat')
disp('Initialized with G_true')
disp(nanmean(squeeze(final_tt_comp(1,:,:)),1)) % 
disp('Initialized with random orthognal matrices (standard approach). The mean and sd, per snr level, are:')
disp(nanmean(squeeze(final_tt_comp(2,:,:)),1))
disp(sqrt(nanvar(squeeze(final_tt_comp(2,:,:)),[],1)))

disp('Mean and variance of tt_comp for random init')
disp([nanmean(nanmean(squeeze(final_tt_comp(2,:,:)),1)), sqrt(nanvar(nanmean(squeeze(final_tt_comp(2,:,:)),1),[],2))])