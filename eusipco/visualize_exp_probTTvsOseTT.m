%%
load('C:\Users\jehi\Coding\prob-tt\eusipco\synthetic_tensortrain_vs_probtt.mat')
font_size=12

name_methods{end-3} = 'TTD(fixed, D=6)';
name_methods{end-2} = 'TTD(fixed, D=12)';
name_methods{end-1} = 'PTTD(D=6)';
name_methods{end} = 'PTTD(D=12)';
%%
line_colors = ones(length(name_methods),3).*(1-(0.3+0.7*(length(name_methods):-1:1)/length(name_methods)))';
line_colors(end-3:end,:) = [0.9,0,0; 0.9,0,0; 0, 0.1, 1; 0, 0.1, 1];
%
for i = 1:3
    
    figure('Position',[100*(i)^2,500,500,300])
    h_plot = cell(length(name_methods),1);
    hold on
    for j = 1:length(name_methods)
        if i == 1
            y = all_rmse_noise;
            title('Reconstruction error on X_{noise}','Fontsize',font_size+2)
        elseif i==2
            y = all_rmse_true;
            title('Reconstruction error','Fontsize',font_size+2)
        elseif i==3
            y = all_numel;
            title('TensorTrain Size vs Noise Level','Fontsize',font_size+2)
        end
        if j < length(name_methods)-3
            h_plot{j} = plot(y(:,j), 'Color', line_colors(j,:));
        else
            line = '-x';
            switch j
                case length(name_methods)-3
                    line = '-x';
                case length(name_methods)-2
                    line = '-o';
                case length(name_methods)-1
                    line = '-s';
                case length(name_methods)
                    line = '-*';
            end
            h_plot{j} = plot(y(:,j),line, 'Color', line_colors(j,:));
        end
            
            
    end
    
    hold off
    if i == 3
        %pass
        set(gca,'YScale','log')
        ylabel('Elements in G','Fontsize',font_size)
        ylim([0,nanmax(all_numel(:))*2])
    else
        set(gca,'YScale','log')
        ylabel('RMSE','Fontsize',font_size)
    end
    set(gca,'XTick',1:2:length(snr_list), 'XTickLabel',snr_list(1:2:end))
    xlabel('Signal to noise ratio (SNR) in dB','Fontsize',font_size)
    if any(i == [2,3])
        i_leg = [1,length(name_methods)-4:length(name_methods)];
        if i == 2
            legend([h_plot{i_leg}], name_methods(i_leg),'Location','southwest','Fontsize',font_size-2)
        else
            legend([h_plot{i_leg}], name_methods(i_leg),'Location','northwest','Fontsize',font_size-2)
        end
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