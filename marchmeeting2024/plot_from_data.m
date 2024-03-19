basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3#T_Bl_Sq_2to7P_15_1_3'
    };
x_cell = {
    [1, 6, 7, 10]
%     [25, 26, 27, 28, 29, 30, 31, 32, 33, 34]  % load the x-th set of data from the first basedir to analyze
    };
larvaNum = {[1], [7], [2], [24]};  % larva index in each expt, could be [1 2 4]
plot_pturn = false;
plot_speed_individual = true;


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


tbin = 3;  
tperiod = 15;
stepsize = 0.1;
stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start
edges = 0:tbin: tperiod;
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
%     for i = 1 : length(fieldnames(larvae))  % loop for each larva in larvae
    for i = larvaNum{j}  % which larva to plot
        pturn = getfield(larvae, ['larva', num2str(i)], 'pturn');
        pturn_error = getfield(larvae, ['larva', num2str(i)], 'pturn_error');
        turnStart_mean = getfield(larvae, ['larva', num2str(i)], 'turnStart_mean');
        turnStart_std = getfield(larvae, ['larva', num2str(i)], 'turnStart_std');
        response = [getfield(larvae, ['larva', num2str(i)], 'response', 'blue'), getfield(larvae, ['larva', num2str(i)], 'response', 'red'), getfield(larvae, ['larva', num2str(i)], 'response', 'bluered')];
        expt_info = getfield(larvae, ['larva', num2str(i)], 'expt_info');
        trackNum = getfield(larvae, ['larva', num2str(i)], 'trackNum');
        speed = getfield(larvae, ['larva', num2str(i)], 'speed');
        speed_std = getfield(larvae, ['larva', num2str(i)], 'speed_std');
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
            savename = strcat(response, '_', 'pturn_track_', num2str(trackNum), '_', expt_info(end-11: end));
            savefig(gcf, savename);
        end

        if plot_speed_individual
            figure;
            x_ton = (stepsize/2):stepsize:(tperiod-stepsize/2);
            for k = 1 : length(fieldnames(speed))  % loop for each stimulation
                ax(k) = subplot(length(fieldnames(pturn)), 1, k);
                y_v = speed.(stim_color{k});
                std = speed_std.(stim_color{k});
                uppercurve = y_v + 0.5*std;
                lowercurve = y_v - 0.5*std;
                x_tofill = [x_ton, fliplr(x_ton)];  % the x axis of the ploygon to fill
                y_tofill = [lowercurve, fliplr(uppercurve)];
                pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none'); hold on;  % no edges for the patch
                plot(x_ton, y_v, 'Color', 0.2*[1 1 1]); xline(6, '--'); hold off;
                xlabel('ton (s)'); ylabel('Speed (cm/min)');
                pos = get(gca, 'Position');
                pos(3) = pos(4) * 1.5;  % width = height * 1.5
                set(gca, 'Position', pos)
                ylim([0, 2.2]);
            end  % end looping stimulation
            sgtitle(response);
            savename = strcat(response, '_', 'speed_track_', num2str(trackNum), '_', expt_info(end-11: end));
            savefig(gcf, savename);
        end  % end plotting speed_individual
    end  % end plotting all larvae in the expt
end  % end plotting all desired expts
