% Re-plot figures directly from the structure larvae saved in data.m of each experiment result folder. 
% It's available to select the certain larvae. 
% Plot each larva at a time, so the size of the graph is tunable.

%% Get cell figlocation_list
basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3#T_Bl_Sq_2to7P_15_1_3'
    };
x_cell = {
    1:21
%     [1, 6, 7, 10]
%     [25, 26, 27, 28, 29, 30, 31, 32, 33, 34]  % load the x-th set of data from the first basedir to analyze
    };
larvaNum = {[1], [7], [2], [24]};  % larva index in each expt, could be [1 2 4]



% create figlocation_list, a cell containing all folders of results
nexp = length([x_cell{:}]);  % number of experiments
figlocation_list = cell(nexp, 1);  % initialize a nexpt-by-1 empty cell
j = 0;  % the j-th experiment
for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    basedir = basedir_cell{folder_index};
    x = x_cell{folder_index};
    d = dir(fullfile(basedir, 'matfiles', '*.mat'));
    for k = 1 : length(x)  % loop for each expt in the eset
        j = j + 1;
        figlocation_list{j} = fullfile(basedir, ['results', d(x(k)).name(end-16:end-4)]);
    end
end

%% plot info of individual larva to one graph
download = true;
plot_pturn = false;
plot_speed_individual = true;
plot_turnStart = false;
plot_turnStart_turnOrder = false;
plot_turnStartTime_turnOrder = false;
plot_cdf_individual =false;
tbin = 3;  
tperiod = 15;
stepsize = 0.1;
stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start
color_pad = {'blue', 'red', 'black'};
edges = 0:tbin: tperiod;
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
%     for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
    for i = larvaNum{j}  % select certain larva to plot
        pturn = getfield(larvae, ['larva', num2str(i)], 'pturn');
        pturn_error = getfield(larvae, ['larva', num2str(i)], 'pturn_error');
        turnStart_mean = getfield(larvae, ['larva', num2str(i)], 'turnStart_mean');  % the mean turnStart time in each bin
        turnStart_std = getfield(larvae, ['larva', num2str(i)], 'turnStart_std');  % the std turnStart time in each bin
        response = [getfield(larvae, ['larva', num2str(i)], 'response', 'blue'), getfield(larvae, ['larva', num2str(i)], 'response', 'red'), getfield(larvae, ['larva', num2str(i)], 'response', 'bluered')];
        expt_info = getfield(larvae, ['larva', num2str(i)], 'expt_info');
        trackNum = getfield(larvae, ['larva', num2str(i)], 'trackNum');
        speed = getfield(larvae, ['larva', num2str(i)], 'speed');
        speed_std = getfield(larvae, ['larva', num2str(i)], 'speed_std');
        turnStart = getfield(larvae, ['larva', num2str(i)], 'turnStart');
        turnStartTime = getfield(larvae, ['larva', num2str(i)], 'turnStartTime');
        if plot_pturn
            figure;
            for k = 1 : length(fieldnames(pturn))  % loop for each stimulation
                ax(k) = subplot(length(fieldnames(pturn)), 1, k);
                pturn_temp = pturn.(stim_color{k});
                pturn_error_temp = pturn_error.(stim_color{k});
                turnStart_mean_temp = turnStart_mean.(stim_color{k});
                turnStart_std_temp = turnStart_std.(stim_color{k});
                bar(xbar, pturn_temp, 1); hold on; % the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
                errorbar(turnStart_mean_temp, pturn_temp, pturn_error_temp / 2, 'k.'); % the length of the error bar is the error.
                errorbar(turnStart_mean_temp, pturn_temp, turnStart_std_temp/2, 'horizontal', 'k.'); hold off;
                xticks(edges); ylim([0, 1]); xline(6, '--');
                pos = get(gca, 'Position');
                pos(3) = pos(4) * 1.5;  % width = height * 1.5
                set(gca, 'Position', pos)
            end
            sgtitle(response);
            if download
                savename = strcat(pwd, '\results_fig', ['\', response, '_', 'pturn_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
                savefig(gcf, savename);
            end  % end download
        end  % end plot_pturn

        if plot_speed_individual
            figure;
            x_ton = (stepsize/2):stepsize:(tperiod-stepsize/2);
            for k = 1 : length(fieldnames(speed))  % loop for each stimulation
%                 subplot(length(fieldnames(pturn)), 1, k);
                hold on;  % to put v at different stimulation together
                y_v = speed.(stim_color{k});
                std = speed_std.(stim_color{k});
                uppercurve = y_v + 0.5*std;
                lowercurve = y_v - 0.5*std;
                x_tofill = [x_ton, fliplr(x_ton)];  % the x axis of the ploygon to fill
                y_tofill = [lowercurve, fliplr(uppercurve)];
                pathObj = fill(x_tofill, y_tofill, color_pad{k}, 'FaceAlpha', 0.2, 'LineStyle', 'none'); % hold on;  % no edges for the patch
                plot(x_ton, y_v, color_pad{k}, 'LineWidth', 2); xline(6, '--'); % hold off;
                ax = gca; ax.FontSize = 16; 
                xlabel('ton (s)', 'FontSize', 18); ylabel('Speed (cm/min)', 'FontSize', 18);
                pos = get(gca, 'Position');
%                 pos(3) = pos(4) * 1.5;  % width = height * 1.5
%                 set(gca, 'Position', pos)
                ylim([0, 2.2]);
            end  % end looping stimulation
            sgtitle(response, 'FontSize', 18); hold off;
            if download
                savename = strcat(pwd, '\results_fig', ['\', response, '_', 'speed_combined_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
                savefig(gcf, savename);
                close;
            end
        end  % end plot_speed_individual

        if plot_turnStart
            figure;
            for k = 1 : length(fieldnames(turnStart))  % loop for each stimulation
                plot(turnStart.(stim_color{k}), k, [color_pad{k}, '.']); hold on;
            end
            yticks(1:k); yticklabels(stim_color); ylim([0.5, k+0.5]);
            xlim([0, tperiod]); xlabel('Turn start time in period');
            title(response); hold off;
            if download
                savename = strcat(pwd, '\results_fig', ['\', response, '_', 'turnStart_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
                savefig(gcf, savename);
            end  % end download
        end  % end plot_turnStart

        if plot_turnStart_turnOrder
            figure;
            clear nturn; nturn(1) = 1;
            for k = 1 : length(fieldnames(turnStart))  % loop for each stimulation
                turnStart_temp = turnStart.(stim_color{k});
                nturn(k+1) = length(turnStart_temp);  % nturn = [1, nturn(stim_color)]
                plot(sum(nturn(1:k)) : (sum(nturn(1:(k+1))) - 1), turnStart_temp, [color_pad{k}, 'o']); hold on;
            end
            ylim([0, tperiod]); hold off;
            xlabel('The order of turn'); ylabel('Turn start time in period (s)');
            title(response);
            if download
                savename = strcat(pwd, '\results_fig', ['\', response, '_', 'turnStart_turnOrder_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
                savefig(gcf, savename);
            end  % end download
        end  % end plot_turnStart_turnOrder

        if plot_turnStartTime_turnOrder
            figure;
            plot(turnStartTime, '.'); hold on;
            yline([600, 1200]);
            ylim([0, 1800]); hold off;
            xlabel('The order of turn'); ylabel('Turn start time (s)');
            title(response);
            if download
                savename = strcat(pwd, '\results_fig', ['\', response, '_', 'turnStartTime_turnOrder_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
                savefig(gcf, savename);
            end  % end download
        end


        if plot_cdf_individual
            figure;
            for k = 1 : length(fieldnames(turnStart))  % loop for each stimulation
                turnStart_temp = turnStart.(stim_color{k});
                turnStart_toff = ton_to_toff(turnStart_temp, 6, 9);
                h_temp = cdfplot(turnStart_toff);
                h_temp.Color = color_pad{k}; h_temp.LineWidth = 2;
                xlim([0, tperiod]); xline(9, '--');
                hold on;
            end
            hold off; grid off;
            ax = gca; ax.FontSize = 16; 
            xlabel('toff (s)'); ylabel('Pturn(t <= toff)');
            title(response);
            if download
                savename = strcat(pwd, '\results_fig', ['\', response, '_', 'cdf_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
                savefig(gcf, savename);
            end  % end download
        end  % end plot_cdf_individual
                
    end  % end loop each desired larva in the result folder
    
end  % end loop all result folders

%% plot info of a group of larvae to one graph
tperiod = 15; download = false;
color_pad = {'blue', 'red', 'black'};
stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start
cdf9 = {[], [], []};  % cdf at toff = 9 for each stim_color respectively;
for k = 1 : length(fieldnames(turnStart))  % loop for each stimulation
    figure;
    for j = 1 : length(figlocation_list)
        savename = strcat(figlocation_list{j}, '\data.mat');
        load(savename, 'larvae');
        for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
%         for i = larvaNum{j}  % select certain larva to plot
            turnStart = getfield(larvae, ['larva', num2str(i)], 'turnStart');

            turnStart_temp = turnStart.(stim_color{k});
            if ~isempty(turnStart_temp)
                turnStart_toff = ton_to_toff(turnStart_temp, 6, 9);
                [h_temp, stats_temp] = cdfplot(turnStart_toff);
                h_temp.Color = color_pad{k};
                F9_index = find(h_temp.XData < 9, 1, 'last');
                F9 = h_temp.YData(F9_index);
                cdf9{k} = [cdf9{k}, F9];
            end
            xlim([0, tperiod]);
            hold on;
        end  % end looping every larva in one experiment
    end  % end looping every experiment
    ax = gca; ax.FontSize = 20; 
    xline(9, 'w-', 'LineWidth', 3); hold off;
    xlabel('toff (s)', 'FontSize', 20); ylabel('Pturn(t <= toff)', 'FontSize', 20); 
    title(stim_color{k}, 'FontSize', 24);
    if download
        savename = strcat(basedir_cell{1}, '\results', ['\cdf_every_', stim_color{k}]);
        savefig(gcf, savename);
    end
end  % end looping every stimulationn condition


figure; hold on;
for k = 1 : length(cdf9)
    histogram(cdf9{k}, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'FaceColor', color_pad{k})
end
ax = gca; ax.FontSize = 20; 
hold off;
xlabel('CDF of pturn at toff = 9s', 'FontSize', 20); ylabel('Proportion of larvae', 'FontSize', 20);