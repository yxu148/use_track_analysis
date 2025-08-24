% Re-plot figures directly from the structure larvae saved in data.m of each experiment result folder. 
% It's available to select the certain larvae. 
% Plot each larva at a time, so the size of the graph is tunable.

%% Get cell figlocation_list
basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data_new\variability_new_try_extracted\Or42a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3_#T_Bl_Sq_2to7P_15_1_3_'
    };
x_cell = {
%     [6, 8, 6, 4]
%     [1, 6, 7, 10]
    1:22
    };
% larvaNum = {[1], [7], [2], [24]};  % larva index in each expt, could be [1 2 4]
% larvaNum = {[5], [2], [7], [4]};  % larva index in each expt, could be [1 2 4]
pause('off');


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
plot_pturn = true;
plot_speed_individual = true;
plot_turnStart = false;
plot_turnStart_turnOrder = false;
plot_turnStartTime_turnOrder = false;
plot_cdf_individual =true;
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
        speed_run = getfield(larvae, ['larva', num2str(i)], 'speed_run');
        speed_run_std = getfield(larvae, ['larva', num2str(i)], 'speed_run_std');
        turnStart = getfield(larvae, ['larva', num2str(i)], 'turnStart');
        turnStartTime = getfield(larvae, ['larva', num2str(i)], 'turnStartTime');
        nperiod = getfield(larvae, ['larva', num2str(i)], 'nperiod');
        if plot_pturn
            figure;
            for k = 1 : length(fieldnames(pturn))  % loop for each stimulation
                ax(k) = subplot(length(fieldnames(pturn)), 1, k);
                pturn_temp = pturn.(stim_color{k});
                pturn_error_temp = pturn_error.(stim_color{k});
                turnStart_mean_temp = turnStart_mean.(stim_color{k});
                turnStart_std_temp = turnStart_std.(stim_color{k});
                nperiod_temp = nperiod.(stim_color{k});
                turnStart_temp = turnStart.(stim_color{k});
                bar(xbar, pturn_temp, 1); hold on; % the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
                errorbar(turnStart_mean_temp, pturn_temp, pturn_error_temp / 2, 'k.'); % the length of the error bar is the error.
                errorbar(turnStart_mean_temp, pturn_temp, turnStart_std_temp/2, 'horizontal', 'k.'); hold off;
                xticks(edges); ylim([0, 1]); xline(6, '--');
                pos = get(gca, 'Position');
                pos(3) = pos(4) * 1.5;  % width = height * 1.5
                set(gca, 'Position', pos)
                % title([num2str(length(turnStart_temp)), ' turns in ', num2str(nperiod_temp), ' periods'])
                ax(k).FontSize = 16;
                % xlabel('ton (s)', 'FontSize', 18); ylabel('Turn probability', 'FontSize', 18);
            end
            sgtitle(response);
            if download
                savename = strcat(basedir_cell{1}, '\results_new', '\results_fig', ['\', response, '_', 'pturn_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
                savefig(gcf, savename); close;
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
                ylim([0, 2.5]);
            end  % end looping stimulation
            sgtitle(response, 'FontSize', 18); hold off;
            if download
                savename = strcat(basedir_cell{1}, '\results_new', '\results_fig', ['\', response, '_', 'speed_combined_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
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
                savename = strcat(basedir_cell{1}, '\results_new', '\results_fig', ['\', response, '_', 'turnStart_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
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
                savename = strcat(basedir_cell{1}, '\results_new', '\results_fig', ['\', response, '_', 'turnStart_turnOrder_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
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
                savename = strcat(basedir_cell{1}, '\results_new', '\results_fig', ['\', response, '_', 'turnStartTime_turnOrder_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
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
                savename = strcat(basedir_cell{1}, '\results_new', '\results_fig', ['\', response, '_', 'cdf_track_', num2str(trackNum), '_', expt_info(end-11: end)]);
                savefig(gcf, savename); close;
            end  % end download
        end  % end plot_cdf_individual
                
    end  % end loop each desired larva in the result folder
    
end  % end loop all result folders

%% re-plot figures that are done in load_multi_data.m
download = true;
plot_speed_run = true;
tbin = 3;
tperiod = 15;
stepsize = 0.1;
stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
color_pad = {'blue', 'red', 'm', 'black'};
for expt_index = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{expt_index}, '\data.mat');
    load(savename, 'larvae');

    if plot_speed_run
        nLarvae = length(fieldnames(larvae));
        cols = ceil(sqrt(nLarvae));         % Number of columns, plot subplots close to square
        rows = ceil(nLarvae / cols);        % Number of rows
        figure;
        for j = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
    %     for j = larvaNum{expt_index}  % select certain larva to plot
            subplot(rows, cols, j); hold on;
            expt_info = getfield(larvae, ['larva', num2str(j)], 'expt_info');
            trackNum = getfield(larvae, ['larva', num2str(j)], 'trackNum');
            speed_run = getfield(larvae, ['larva', num2str(j)], 'speed_run');
            speed_run_std = getfield(larvae, ['larva', num2str(j)], 'speed_run_std');
            nperiod = getfield(larvae, ['larva', num2str(j)], 'nperiod');
            response = strjoin(struct2cell(larvae.larva1.response)', '');  % str
            x_ton = (stepsize/2):stepsize:(tperiod-stepsize/2);
            for k = 1 : length(fieldnames(speed_run))  % loop for each stimulation
                y_v = speed_run.(stim_color{k});
                y_v(y_v==0) = NaN;  % to plot zeros as broken
                std = speed_run_std.(stim_color{k});
                uppercurve = y_v + 0.5*std;
                lowercurve = y_v - 0.5*std;
                x_tofill = [x_ton, fliplr(x_ton)];  % the x axis of the ploygon to fill
                y_tofill = [lowercurve, fliplr(uppercurve)];
                pathObj = fill(x_tofill, y_tofill, color_pad{k}, 'FaceAlpha', 0.2, 'LineStyle', 'none'); % hold on;  % no edges for the patch
                plot(x_ton, y_v, color_pad{k}); 
            end  % end looping stimulation
            xline(6, '--'); ylim([0, 3]); hold off;
            xlabel('ton (s)'); ylabel('Run speed (cm/min)');
            title(['(', response, ', ', num2str(trackNum), ')']);
%             pos = get(gca, 'Position');
%             pos(3) = pos(4) * 1.5;  % width = height * 1.5
%             set(gca, 'Position', pos)
        end  % end looping each larva j
        sgtitle('(response, trackNum)');
        if download
            savename = strcat(basedir_cell{1}, ['\results_', expt_info(end-11: end)], '\vRun_ton_individual');
            savefig(gcf, savename);
            close;
        end
    end  % end plot_speed_run
end  % end looping each expt

%% plot info (CDF) of a group of larvae to one graph, or other plots about CDF
tperiod = 15; 
download = true; 
pause on;
color_pad = {'b', 'r', 'k', 'b', 'r', 'k', 'k'};
stim_color = {'blue1', 'red1', 'bluered1', 'blue2', 'red2', 'bluered2', 'dark'};  % descripiton of the t_stim_start
cdf9 = {[], [], [], [], [], [], []};  % cdf at toff = 9 for each stim_color respectively;
tx = 0 : 0.1 : tperiod;  % the time to calculate the average cdf and std
cdf_tx = {[], [], [], [], [], [], []};
plot_cdf9 = false;
plot_cdf_mean = true;
plot_cdf_trying = false;
plot_pturn_cdf = true;


for k = 1 : 7  % loop for each stimulation
    figure;
    for j = 1 : length(figlocation_list)
        savename = strcat(figlocation_list{j}, '\data.mat');
        load(savename, 'larvae');
        for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
    %         for i = larvaNum{j}  % select certain larva to plot
            if getfield(larvae, ['larva', num2str(i)], 'valid')
                turnStart = getfield(larvae, ['larva', num2str(i)], 'turnStart');
    
                turnStart_temp = turnStart.(stim_color{k});
                if ~isempty(turnStart_temp)
                    turnStart_toff = ton_to_toff(turnStart_temp, 6, 9);
    
                    [F, X] = ecdf(turnStart_toff);
                    stairs(X, F, 'Color', color_pad{k});
                    X = X(2 : end); F = F(2 : end);  % delete the duplicated initial point
                    if X(end) < tperiod  % add an end one to the cdf
                        X = [X ; tperiod]; F = [F; 1];
                    end
                    if X(1) > 0  % add initial zero to the cdf
                        X = [0; X]; F = [0; F];
                    end
    
                    % find the cdf value at 9
                    F9 = interp1(X, F, 9, 'previous');
                    cdf9{k} = [cdf9{k}, F9];
    
                    % get the cdf at desired time for each larva
                    cdf_tx{k} = [cdf_tx{k}; interp1(X, F, tx, 'previous')];
                end  % end if turnStart_temp is not empty
                xlim([0, tperiod]);
                hold on;
            end  % end if valid
        end  % end looping every larva in one experiment
    end  % end looping every experiment
    ax = gca; ax.FontSize = 20; 
    xline(9, 'w-', 'LineWidth', 3); hold off;
    xlabel('toff (s)', 'FontSize', 20); ylabel('Pturn(t <= toff)', 'FontSize', 20); 
    title(stim_color{k}, 'FontSize', 24);
    pause;
    if download
        savename = strcat(basedir_cell{1}, '\results', '\results_fig', ['\cdf_every_', stim_color{k}]);
        savefig(gcf, savename);
    end 
    close;
end  % end looping every stimulationn condition


% plot the histogram of cdf at 9 seconds
if plot_cdf9
    figure; hold on;
    for k = 1 : length(cdf9)
        histogram(cdf9{k}, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'FaceColor', color_pad{k})
    end
    ax = gca; ax.FontSize = 20; 
    hold off;
    xlabel('CDF of pturn at toff = 9s', 'FontSize', 20); ylabel('Proportion of larvae', 'FontSize', 20);
    pause;
    if download
        savename = strcat(basedir_cell{1}, '\results', '\results_fig', ['\cdf9_hist_', stim_color{k}]);
        savefig(gcf, savename);
    end
    close;
end


% plot the average cdf at each stimulation
color_pad = {'b', 'r', 'k', 'b--', 'r--', 'k--', 'k:'};
if plot_cdf_mean
    figure;
    for k = 1 : length(cdf_tx)
        cdf_temp = cdf_tx{k};
        errorbar(tx, mean(cdf_temp, 1), 0.5 * std(cdf_temp, 0, 1), color_pad{k});
        hold on;
    end
    hold off; ax = gca; ax.FontSize = 20; 
    xlim([0, tperiod]); ylim([0, 1]); legend(stim_color)
    xlabel('toff (s)', 'FontSize', 20); ylabel('Average CDF', 'FontSize', 20);
    pause;
    if download
        savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\cdf_average');
        savefig(gcf, savename);
    end
    close;
end


% trying to make sense of conditional probability by trying different cdf
if plot_cdf_trying
    figure; hold on;
    plot(tx, mean(cdf_tx{1}, 1) .* mean(cdf_tx{2}, 1), 'LineWidth', 2);
    plot(tx, mean(cdf_tx{3}, 1), 'LineWidth', 2);
    plot(tx, mean(cdf_tx{1}, 1) .* mean(cdf_tx{2}, 1) ./ mean(cdf_tx{3}, 1), 'LineWidth', 2)    
    plot(tx, 0.5 * (mean(cdf_tx{1}, 1) + mean(cdf_tx{2}, 1)), 'LineWidth', 2)
    plot(tx, 1/3 * (mean(cdf_tx{1}, 1) + mean(cdf_tx{2}, 1) + mean(cdf_tx{3}, 1)), 'LineWidth', 2);
    plot(tx, mean(cdf_tx{1}, 1), 'LineWidth', 2);
    plot(tx, mean(cdf_tx{2}, 1), 'LineWidth', 2);
    xline(9, 'k--');
    hold off;
    lg = legend('P(turn|blue) * P(turn|red)', 'P(turn|both)', 'P(turn) = P(turn|blue) * P(turn|red) / P(turn|both)', ...
        '1/2 * (P(turn|blue) + P(turn|red))', '1/3 * (P(turn|blue) + P(turn|red) + P(turn|both))', 'P(turn|blue)', 'P(turn|red)');
    lg.Location = 'best';
    ax = gca; ax.FontSize = 20;
    title('Different plots about average cdf P', 'FontSize', 24)
    xlabel('toff (s)', 'FontSize', 20);
    pause;
    if download
        savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\thinking about Bayesian');
        savefig(gcf, savename);
    end
    close;
end

% trying to make sense of Bayes theorem by plotting pturn from cdf
if plot_pturn_cdf
    figure; hold on;
    for k = 1 : length(cdf_tx)
        cdf_mean_temp = mean(cdf_tx{k}, 1);
        pturn_mean.(stim_color{k}) = (cdf_mean_temp(2:end) - cdf_mean_temp(1:(end-1))) ./ diff(tx);
        plot(tx(1: (end-1)), pturn_mean.(stim_color{k}), color_pad{k});
    end
    plot(tx(1: (end-1)), pturn_mean.(stim_color{1}) .* pturn_mean.(stim_color{2}) ./ pturn_mean.(stim_color{3}));
    xline(9, 'k--');
    hold off;
    lg = legend('P(turn|blue)', 'P(turn|red)', 'P(turn|bluered)', 'P(turn) = P(turn|blue) * P(turn|red) / P(turn|bluered)');
    lg.Location = 'best';
    ax = gca; ax.FontSize = 20;
    xlabel('toff (s)', 'FontSize', 20); ylabel('Average Pturn', 'FontSize', 20);
    pause;
    if download
        savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\thinking about Bayesian from pturn');
        savefig(gcf, savename);
    end
    close;
end


%% average turn rate of multiple larvae across different experiments
tperiod = 15; 
download = true; 
pause on;
color_pad = {'b', 'r', 'm', 'k'};
stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
stepsize = 0.1; binsize = 0.5;  % seconds

figure;
for k = 1 : length(color_pad)  % loop for each stimulation
    nperiod_all = 0;
    turnStart_all = [];
    for j = 1 : length(figlocation_list)
        savename = strcat(figlocation_list{j}, '\data.mat');
        load(savename, 'larvae');
        for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
%         for i = larvaNum{j}  % select certain larva to plot
            if getfield(larvae, ['larva', num2str(i)], 'valid')
                turnStart = getfield(larvae, ['larva', num2str(i)], 'turnStart');
                turnStartTime = getfield(larvae, ['larva', num2str(i)], 'turnStartTime');
                nperiod = getfield(larvae, ['larva', num2str(i)], 'nperiod');
                turnStart_all = [turnStart_all, turnStart.(stim_color{k})];
                nperiod_all = nperiod_all + nperiod.(stim_color{k});
            end  % end if valid larva
        end  % end looping larvae
    end  % end looping experiments
    
    turnStart_all_toff = ton_to_toff(turnStart_all, 6, 9);

    [turnrate_all, std] = rate_from_time(turnStart_all_toff, tperiod, stepsize, binsize);
    turnrate_all = turnrate_all./ double(nperiod_all) * 60;
    std = std ./ double(nperiod_all) * 60;
    time_timestep = [0 : fix(tperiod/stepsize)] * stepsize;

    ax(k) = subplot(length(color_pad),1,k);
    plot(time_timestep, turnrate_all); xline(9, 'k--');
    xlabel('ton (s)'); ylabel('Reorientation Rate (per min)'); 
    title([num2str(length(turnStart_all)), ' turns, in ', num2str(nperiod_all), ' periods, ', num2str(k), '-th intensity of stimulation']);

    hold on;  % to add standard deviation
    uppercurve = turnrate_all + 0.5*std;
    lowercurve = turnrate_all - 0.5*std;
    x_tofill = [time_timestep, fliplr(time_timestep)];  % the x axis of the ploygon to fill
    y_tofill = [lowercurve, fliplr(uppercurve)];
    pathObj = fill(x_tofill, y_tofill, color_pad{k}, 'FaceAlpha', 0.2, 'LineStyle', 'none');  % no edges for the patch
    hold off;

end  % end looping stimluation conditions
sgtitle(['Step size = ', num2str(stepsize), ', bin size = ', num2str(binsize)]); linkaxes(ax, 'y');  % align subplots with y axis
pause;
if download
    savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\turnrate');
    savefig(gcf, savename);
end
close;

%% average turn rate of multiple larvae across different experiments, plot together
tperiod = 15; 
download = true; 
pause on;
color_pad = {'b', 'r', 'k', 'k--'};
stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
handles = gobjects(1, length(color_pad)); % Preallocate array for handles
stepsize = 0.1; binsize = 0.5;  % seconds
figure; hold on;
turns_count = 0;
for k = 1 : length(color_pad)  % loop for each stimulation
    nperiod_all = 0;
    turnStart_all = [];
    for j = 1 : length(figlocation_list)
        savename = strcat(figlocation_list{j}, '\data.mat');
        load(savename, 'larvae');  % use larvae/tracks to plot
        for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
%         for i = larvaNum{j}  % select certain larva to plot
            if getfield(larvae, ['larva', num2str(i)], 'valid')
                turnStart = getfield(larvae, ['larva', num2str(i)], 'turnStart');
                turnStartTime = getfield(larvae, ['larva', num2str(i)], 'turnStartTime');
                nperiod = getfield(larvae, ['larva', num2str(i)], 'nperiod');
                turnStart_all = [turnStart_all, turnStart.(stim_color{k})];
                nperiod_all = nperiod_all + nperiod.(stim_color{k});
            end  % end if valid larva
        end  % end looping larvae
    end  % end looping experiments

    turns_count = turns_count + length(turnStart_all);  % length(turnStart_all) is the turns_count for this stimulation

    turnStart_all_toff = ton_to_toff(turnStart_all, 6, 9);

    [turnrate_all, std] = rate_from_time(turnStart_all_toff, tperiod, stepsize, binsize);
    turnrate_all = turnrate_all./ double(nperiod_all) * 60;
    std = std ./ double(nperiod_all) * 60;
    time_timestep = [0 : fix(tperiod/stepsize)] * stepsize;

    handles(k) = plot(time_timestep, turnrate_all, color_pad{k}); xline(9, 'k--');

    uppercurve = turnrate_all + 0.5*std;
    lowercurve = turnrate_all - 0.5*std;
    x_tofill = [time_timestep, fliplr(time_timestep)];  % the x axis of the ploygon to fill
    y_tofill = [lowercurve, fliplr(uppercurve)];
    pathObj = fill(x_tofill, y_tofill, color_pad{k}, 'FaceAlpha', 0.2, 'LineStyle', 'none');  % no edges for the patch

end  % end looping stimluation conditions
hold off;
legend(handles, stim_color)
xlabel('toff (s)'); ylabel('Reorientation Rate (per min)'); 
disp([num2str(turns_count), ' turns']); 
pause;
if download
    mkdir(fullfile(basedir_cell{1}, '\results', 'results_fig'));
    savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\turnrate_combined_valid_larvae');
    savefig(gcf, savename);
end
close;
exportgraphics(gcf, [savename, '.pdf']);
print(gcf, savename, '-dpdf');  % saves as 'myplot.pdf'
%% average turn rate of all tracks across different experiments verses texpt
texpt = 2400;  % second, experiment time
download = true; 
pause on;
stepsize = 0.1; binsize = 0.5;  % seconds

nperiod_all = 0;
turnStartTime_all = [];
ntracks_frame_all = 0;
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'tracks', 'GQled1Val', 'GQled2Val', 'ntracks_frame');
    for i = 1 : length(fieldnames(tracks))  % loop every track in the experiment
%         for i = larvaNum{j}  % select certain larva to plot
        turnStartTime = getfield(tracks, ['track', num2str(i)], 'turnStartTime');
        turnStartTime_all = [turnStartTime_all, turnStartTime];
    end  % end looping larvae
    ntracks_frame_downsampling = ntracks_frame(round(linspace(1, length(ntracks_frame), texpt/stepsize + 1)));
    ntracks_frame_all = ntracks_frame_all + ntracks_frame_downsampling;
end  % end looping experiments

[turnrate_all, turnrate_error] = rate_from_time(turnStartTime_all, texpt, stepsize, binsize);
turnrate = turnrate_all./  ntracks_frame_all * 60;
turnrate_error = turnrate_error ./ ntracks_frame_all * 60;  % error scales linearly with the value, so also need to divide
time_timestep = [0 : fix(texpt/stepsize)] * stepsize;
uppercurve = turnrate + 0.5*turnrate_error;
lowercurve = turnrate - 0.5*turnrate_error;
x_tofill = [time_timestep, fliplr(time_timestep)];  % the x axis of the ploygon to fill
y_tofill = [lowercurve, fliplr(uppercurve)];
y_tofill(isnan(y_tofill)) = 0;  % replace the NA to 0
figure; hold on;
yyaxis left
fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none');
plot(time_timestep, turnrate, '-');
xlabel('Reorientation Start time (s)'); ylabel('Reorientation Rate (per min)'); 
title([num2str(length(turnStartTime_all)), ' turns, ', num2str(max(ntracks_frame_all)), ' larvae']);
yyaxis right
plot(GQled1Val.xData, GQled1Val.yData, 'r-')
plot(GQled2Val.xData, GQled2Val.yData, 'b-')
ylabel('ledVal (PWM)'); 
hold off;
pause;
if download
    savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\turnrate_texpt_all_tracks');
    savefig(gcf, savename);
end
close;

%% create a new response with new criteria

save_data = true;
stim_color = {'blue1', 'blue2', 'blue3'};  % descripiton of the t_stim_start
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
%     for i = larvaNum{j}  % select certain larva to plot
        larva_index = ['larva', num2str(i)];
        pturn = getfield(larvae, larva_index, 'pturn');
        for k = 1 : length(fieldnames(pturn))  % loop for each stimulation
            if pturn.(stim_color{k})(1) - mean(pturn.(stim_color{k})(3:end)) > 0.2  % if the first bin of pturn is much larger than low, call it response
                larvae.(larva_index).response.(stim_color{k}) = '1';
            elseif mean(pturn.(stim_color{k})(1:2)) - mean(pturn.(stim_color{k})(3:end)) > 0.2
                larvae.(larva_index).response.(stim_color{k}) = '1';
            else
                larvae.(larva_index).response.(stim_color{k}) = '0';
            end  % end defining response2
        end  % end looping for every stimulation
    end  % end looping for every larva

    if save_data
        savename = strcat(figlocation_list{j}, '\data.mat');
        if isfile(savename)
            save(savename, 'larvae', '-append')
        else
            save(savename, 'larvae')
        end
    end  % end saving data

end  % end looping for every experiments 


%% create a label 'valid' for each larva based on if the speed is ever bigger than 0.8 cm/min

save_data = true;
stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
        larva_index = ['larva', num2str(i)];
        speed = getfield(larvae, larva_index, 'speed');
        %  if the maximum speed under any stimulation condition is less
        %  than 0.8 cm/min, label the larva invalid.
        if (max(speed.(stim_color{1})) < 0.8) || (max(speed.(stim_color{2})) < 0.8) || (max(speed.(stim_color{3})) < 0.8)
            larvae.(larva_index).valid = 0;
        else
            larvae.(larva_index).valid = 1;
        end  % end if criteria
    end  % end looping all larva in the larvae

    if save_data
        if isfile(savename)
            save(savename, 'larvae', '-append')
        else
            save(savename, 'larvae')
        end
    end  % end saving data
end  % end looping each experiment


%% distribution of max speed under different stimulation

stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
max_speed = [];
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
        larva_index = ['larva', num2str(i)];
        speed = getfield(larvae, larva_index, 'speed');
        max_speed = [max_speed, max(speed.(stim_color{1})), max(speed.(stim_color{2})), max(speed.(stim_color{3}))];
    end
end

% fit with raw data
pd = fitdist(init_speed{i}', 'Normal');  mu = pd.mu; sigma = pd.sigma;% normal distribution
h = histogram(data, 0:0.2:4); hold on;
binWidth = h.BinWidth; nSamples = length(data);
x = linspace(min(data), max(data), 100);
y = nSamples * binWidth * normpdf(x, mu, sigma);
plot(x, y)
legend('Data Histogram', sprintf('Gaussian Fit (\\mu = %.2f, \\sigma = %.2f)', mu, sigma), 'Location', 'eastoutside')


figure;
% fit raw data with normal distribution
pd = fitdist(max_speed', 'Normal');  mu = pd.mu; sigma = pd.sigma;
h = histogram(max_speed); hold on;
binWidth = h.BinWidth; nSamples = length(max_speed);
x = linspace(min(max_speed), max(max_speed), 100);
y = nSamples * binWidth * normpdf(x, mu, sigma);
plot(x, y); hold off;
legend('Data Histogram', sprintf('Gaussian Fit \n (\\mu = %.2f, \\sigma = %.2f)', mu, sigma), 'FontSize',9)
xlabel('Maximum speed under one stimulation (cm/min)'); ylabel('Count');
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\hist_max_speed');
savefig(gcf, savename);

%% Mean run speed of each early start track verses the duration of the track
latest_start = 120;  % seconds, select the tracks that start earlier than latest_start time. 80% length
frame_rate = 20;  % number of frames per second
v_run_mean = [];  % cm/min
duration = [];  % min
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'tracks');
    for i = 1 : length(fieldnames(tracks))  % loop every track in the experiment
        track_index = ['track', num2str(i)];
        if getfield(tracks, track_index, 'startFrame') <= latest_start * frame_rate
            v_run_mean_temp = getfield(tracks, track_index, 'speed_run_mean');
            endtime_temp = getfield(tracks, track_index, 'endFrame') / frame_rate /60;  % min
            v_run_mean = [v_run_mean, v_run_mean_temp];
            duration = [duration, endtime_temp];  % end time of early start track is called duration here
        end  % end if early track
    end  % end looping tracks
end  % end looping expts

figure; plot(duration, v_run_mean, 'o');
xlabel('Duration of larvae (min)'); ylabel('Mean run speed (cm/min)');
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\speed_mean_run-duration_early_start_tracks');
savefig(gcf, savename);
%% create a label 'init_speed' for each larva,
% based on the mean speed within 1 second before stimulation on
% with unit of cm/min, and 1 float number for each stimulation
save_data = true;
stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
stepsize = 0.1;  % in second, when generate the speed vectors
dt = 1;  % in second, get mean speed of 1 second before stimulation on
init_speed = [];
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
        larva_index = ['larva', num2str(i)];
        speed = getfield(larvae, larva_index, 'speed');
        for s = 1 : length(stim_color)  % -1 if exclude the last stim_color, 'dark'
            speed_temp = getfield(speed, stim_color{s});
            larvae.(larva_index).init_speed.(stim_color{s}) = mean(speed_temp(end-dt/stepsize + 1 : end));
        end  % end looping stim_color
    end  % end looping larvae

    if save_data
        if isfile(savename)
            save(savename, 'larvae', '-append')
        else
            save(savename, 'larvae')
        end
    end  % end saving data
end

%% distribution of initial speed (within 1 s before stimulation) for each stimulation

stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
init_speed = [];  % num_larvae-by-num_stim
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
        v_init = getfield(larvae, ['larva', num2str(i)], 'init_speed');
        init_speed = [init_speed; getfield(v_init, stim_color{1}), getfield(v_init, stim_color{2}), getfield(v_init, stim_color{3}), getfield(v_init, stim_color{4})];
    end
end

figure; 
for i = 1: size(init_speed, 2)
    ax(i) = subplot(size(init_speed, 2), 1, i); histogram(init_speed(:, i)); title(stim_color{i}); xlim([0, 4]);
end
sgtitle('Mean speed within 1 second before stimulation on of all larvae (cm/min)');
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\hist_1smean_speed_stim');
savefig(gcf, savename);

%% distribution of initial speed (within 1 s before stimulation) for each pulse for each stimulation

stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
title_name = {'Visual stimulation', 'Or42a stimulation', 'Visual and Or42a stimulation', 'No stimulation'};
init_speed = cell(1, length(stim_color));  % num_larvae-by-num_stim
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
        v_init_pulse = getfield(larvae, ['larva', num2str(i)], 'init_speed_each');
        for s = 1 : length(stim_color)
            init_speed{s} = [init_speed{s}, getfield(v_init_pulse, stim_color{s})];
        end  % end looping different stimulation
    end  % end looping different larva
end  % end looping different expt

figure; 
for i = 1: size(init_speed, 2)
    data = init_speed{i};
    ax(i) = subplot(size(init_speed, 2), 1, i); 
    % fit with raw data
    pd = fitdist(init_speed{i}', 'Normal');  mu = pd.mu; sigma = pd.sigma;% normal distribution
    h = histogram(data, 0:0.2:4); hold on;
    binWidth = h.BinWidth; nSamples = length(data);
    x = linspace(min(data), max(data), 100);
    y = nSamples * binWidth * normpdf(x, mu, sigma);
    plot(x, y)
    legend('Data Histogram', sprintf('Gaussian Fit (\\mu = %.2f, \\sigma = %.2f)', mu, sigma), 'Location', 'eastoutside')
    title(title_name{i}); xlim([0, 4]);
end
sgtitle({'Mean speed within 1 second before stimulation on', 'of each pulse of each larvae (cm/min)'});  % multiple lines
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\hist_1smean_speed_stim_pulse_GaussianFit');
savefig(gcf, savename);

%% edit the larvae.larva1.longturn_time
long_turn = 6;  % second, turn time longer than long_turn is considered as long turn, default is tperiod
each_stim = true;  % if calculate longturn_time for each stimulation period
inner_plot = true;  % if add detailed histogram as inner plot
stim_color = {'blue', 'red', 'bluered'};  % description of the t_stim_start
t_stim_start = [0, 600, 1200];  % start time (s) of each intensity of stimulation
t_stim_end = [600, 1200, 1800];
longturn_time_larvae_all = [];  % all expts
longturn_time_stim_larvae_all = [];
for k = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{k}, '\data.mat');
    load(savename, 'larvae');
    longturn_time_larvae = zeros(1, length(fieldnames(larvae)));
    longturn_time_stim_larvae = zeros(length(t_stim_end), length(fieldnames(larvae)));
    for j = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
        larva_index = ['larva', num2str(j)];
        turnStartTime = larvae.(larva_index).turnStartTime;
        turnEndTime = larvae.(larva_index).turnEndTime;
        if each_stim
            for i = 1 : length(t_stim_end)  % ith intensity of stimulation---------
                overlap = max(0, min(t_stim_end(i), turnEndTime) - max(t_stim_start(i), turnStartTime));
                longturn_time_stim_larvae(i, j) = sum(overlap(overlap > long_turn));
            end  % end looping i-th stimulation
        else
            duration_turn = turnEndTime - turnStartTime;  % duration of all turns
            longturn_time = sum(duration_turn(duration_turn > long_turn));  % float, second
            longturn_time_larvae(j) = longturn_time;
        end  % end if each_stim
    end  % end looping larva j
    longturn_time_larvae_all = [longturn_time_larvae_all, longturn_time_larvae];
    longturn_time_stim_larvae_all = [longturn_time_stim_larvae_all, longturn_time_stim_larvae];
end  % end looping experiment k

figure;
for i = 1 : size(longturn_time_stim_larvae_all, 1)
    ax(i) = subplot(size(longturn_time_stim_larvae_all, 1), 1, i);
    histogram(longturn_time_stim_larvae_all(i, :)/60, 0:1:10);
    xlabel('Long turn time (min)'); ylabel('Number of larvae'); title(stim_color{i})
    if inner_plot
        mainPos = get(ax(i), 'Position');
        insetPos = [mainPos(1)+(1-0.6)*mainPos(3), mainPos(2)+(1-0.6)*mainPos(4), 0.4*mainPos(3), 0.4*mainPos(4)];  % [left bottom width height]
        insetAxes = axes('Position', insetPos);
        temp = longturn_time_stim_larvae_all(i, :);
        histogram(temp(temp<=60)); box on;
        xlabel('Second'); % title('Detailed histogram of the first 1 min');
    end  % end plotting inner_plot
end
% sgtitle([num2str(size(longturn_time_stim_larvae_all, 2)), ' larvae, long turn is longer than ', num2str(long_turn), ' s']);
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\hist_long_turn_6s_stim_inset');
savefig(gcf, savename);


%% histogram of nperiod of all long tracks (larvae), % of long turn time.

stim_color = {'blue', 'red', 'bluered', 'dark'};  % description of the t_stim_start
time_larvae = [];
longTurnPercentLarvae = [];
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop every larva in the experiment
        larva_index = ['larva', num2str(i)];

        if larvae.(larva_index).valid
            nperiod = getfield(larvae, larva_index, 'nperiod');
            time = 0;  % time of larva under stim_color
            for s = 1 : length(stim_color)
                time = time + nperiod.(stim_color{s}) * 15 / 60;  % min
            end
            time_larvae = [time_larvae, time];
            longturn_time = getfield(larvae, larva_index, 'longturn_time');
            time = double(getfield(larvae, larva_index, 'endFrame') - getfield(larvae, larva_index, 'startFrame')) / 20;  % seconds
            longTurnPercentLarvae = [longTurnPercentLarvae, longturn_time / time * 100];
        end
    end
end

figure; histogram(time_larvae, 26:2:40);  %, 'Normalization', 'probability');
xlabel('Length of long tracks (min)'); ylabel('Count');
disp(['Total ', num2str(length(time_larvae)), ' valid larvae']);
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\hist_long_track_length_full');
savefig(gcf, savename);

% hist_longTurnTime_long_track without seperating different stimulation
figure; histogram(longTurnPercentLarvae, 'Normalization', 'probability');
xlabel('Percentage of long turn time of each larvae (%)'); ylabel('Probability');
disp(['Total ', num2str(length(longTurnPercentLarvae)), ' valid larvae']);
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\hist_longTurnTime_long_track');
savefig(gcf, savename);

%% histogram of start time of all tracks
startFrame_all = [];
endFrame_all = [];
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'tracks');
    for i = 1 : length(fieldnames(tracks))  % loop every track in the experiment
        track_index = ['track', num2str(i)];
        startFrame = getfield(tracks, track_index, 'startFrame');
        startFrame_all = [startFrame_all, startFrame];
        endFrame = getfield(tracks, track_index, 'endFrame');
        endFrame_all = [endFrame_all, endFrame];
    end
end
figure; histogram(startFrame_all/1200, 30);  % nbins
xlabel('Start time of all tracks (min)'); ylabel('Count');
disp([num2str(length(startFrame_all)), ' tracks in total'])
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\hist_track_startFrame');
savefig(gcf, savename);

figure; histogram(endFrame_all/1200, 30);  % nbins
xlabel('End time of all tracks (min)'); ylabel('Count');
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\hist_track_endFrame');
savefig(gcf, savename);

figure; h = histogram(endFrame_all(startFrame_all < 2400)/1200, 0:2:40);  % nbins, 'Normalization', 'probability'
xlabel('Duration of larval tracks (min)'); ylabel('Count of larvae'); % pbaspect([3, 2, 1]);
disp([num2str(sum(h.Values)), ' early start tracks in total'])
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\hist_start_early_track_endFrame');
savefig(gcf, savename); 
% raw data is at h.Values
%% save a copy
savename = strcat(basedir_cell{1}, '\results');
copyfile('plot_from_data.m', savename)  % copy the file to the location