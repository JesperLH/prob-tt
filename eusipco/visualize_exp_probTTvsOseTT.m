%%
load('C:\Users\jehi\Coding\prob-tt\eusipco\synthetic_tensortrain_vs_probtt.mat')

%%
line_colors = ones(length(name_methods),3).*(1-(0.3+0.7*(length(name_methods):-1:1)/length(name_methods)))';
line_colors(end-1:end,:) = [0.9,0,0; 0, 0.1, 1];
%
for i = 1:3
    
    figure('Position',[100*(i)^3,500,500,300])
    h_plot = cell(length(name_methods),1);
    hold on
    for j = 1:length(name_methods)
        if i == 1
            y = all_rmse_noise;
            title('Reconstruction error on X_{noise}')
        elseif i==2
            y = all_rmse_true;
            title('Reconstruction error compared to noiseless data')
        elseif i==3
            y = all_numel;
            title('TensorTrain Size vs Noise Level')
        end
        if j < length(name_methods)-1
            h_plot{j} = plot(y(:,j), 'Color', line_colors(j,:));
        elseif j == length(name_methods)-1
            h_plot{j} = plot(y(:,j),'--o', 'Color', line_colors(j,:));
        else
            h_plot{j} = plot(y(:,j),'--x', 'Color', line_colors(j,:));
        end
            
            
    end
    
    hold off
    if i == 3
        %pass
        ylabel('Elements in the TensorTrain')
        ylim([0,nanmax(all_numel(:))])
    else
        set(gca,'YScale','log')
        ylabel('log(RMSE)')
    end
    set(gca,'XTick',1:2:length(snr_list), 'XTickLabel',snr_list(1:2:end))
    xlabel('Signal to noise ratio (SNR) in dB')
    if any(i == [2,3])
        i_leg = [1,16,17,18];
        legend([h_plot{i_leg}], name_methods(i_leg),'Location','southwest')
    end
    
    axis tight
    if i == 1
       print('./eusipco/results/probTTvsOseTT_noise','-dpng')
       print('./eusipco/results/probTTvsOseTT_noise','-depsc') 
    elseif i == 2
       print('./eusipco/results/probTTvsOseTT_nonoise','-dpng')
       print('./eusipco/results/probTTvsOseTT_nonoise','-depsc')
    elseif i == 3
       print('./eusipco/results/probTTvsOseTT_numel','-dpng')
       print('./eusipco/results/probTTvsOseTT_numel','-depsc')
    end
end

%%
return

figure, plot(1:25, all_tt_comp(:,end-1), '-or', 1:25, all_tt_comp(:,end), '-xb')
legend(name_methods(end-1:end))

%%
figure; plot(all_numel)
set(gca,'XTick',1:2:length(snr_list), 'XTickLabel',snr_list(1:2:end))
xlabel('Signal to noise ratio (SNR) in dB')