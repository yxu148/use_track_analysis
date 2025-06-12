% Plot figures directly from the data.mat of each experiment result folder. 

%% Get cell figlocation_list
basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Tim@Chrimson(3)\T_Re_Sq_44to1098P_10_every_44P#C_Bl_3P'
    };
x_cell = {
%     [6, 8, 6, 4]
%     [1, 6, 7, 10]
    1:4
    };
% larvaNum = {[1], [7], [2], [24]};  % larva index in each expt, could be [1 2 4]
% larvaNum = {[5], [2], [7], [4]};  % larva index in each expt, could be [1 2 4]
% pause('on');


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

%% plot
download = true;  % false, true
texpt = 250;  % second, experiment time 
frame_rate = 20;  % number of frames per second
plot_turnrate_texpt = true;

if plot_turnrate_texpt
    stepsize = 0.5; binsize = 1;  % seconds
    % initialize some values
    turnStartTime_comb=[];
    ntracks_frame_comb = 0;
    for j = 1 : length(figlocation_list)
        savename = strcat(figlocation_list{j}, '\data.mat');
        load(savename, 'turnStartTime_total', 'ntracks_frame', 'GQled1Val');
        turnStartTime_comb = [turnStartTime_comb, turnStartTime_total];
        ntracks_frame_downsampling = ntracks_frame(round(linspace(1, length(ntracks_frame), texpt/stepsize + 1)));
        ntracks_frame_comb = ntracks_frame_comb + ntracks_frame_downsampling;
    end
    [turnrate, turnrate_error] = rate_from_time(turnStartTime_comb, texpt, stepsize, binsize);
    turnrate = turnrate ./ ntracks_frame_comb * 60;
    time_timestep = (0 : fix(texpt/stepsize)) * stepsize;
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
    title([num2str(length(turnStartTime_comb)), ' turns, ', num2str(max(ntracks_frame_comb)), ' larvae']);
    yyaxis right
    plot(GQled1Val.xData, GQled1Val.yData, 'r-')
    plot(GQled2Val.xData, GQled2Val.yData, 'b-')
    ylabel('ledVal (PWM)'); 
    hold off;
    
    pause;
    if download
        mkdir(fullfile(basedir, 'results'));
        savename = strcat(basedir, '\results', '\rate_turn_texpt');
        savefig(gcf, savename); 
    end
    close;
end  % plot_turnrate_texpt

%% save a copy
savename = strcat(basedir_cell{1}, '\results');
copyfile('plot_from_data_a27h.m', savename)  % copy the file to the location