basedir_cell = {
   'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Tim@Chrimson(3)\T_Re_Sq_44to1098P_10_every_44P#C_Bl_3P'
    };
x_cell = {
    1:4
    %[25, 26, 27, 28, 29, 30, 31, 32, 33, 34]  % load the x-th set of data from the first basedir to analyze.
    };
pause('on');  % 'on'---ask the user to press any key to save the figure, and continue; 'off'--directly save without asking

%% load data, multiple eset, multiple expt, into esets
for folder_index = 1 : length(x_cell)  % loop for each basedir folder

    basedir = basedir_cell{folder_index};
    x = x_cell{folder_index};
    d = dir(fullfile(basedir, 'matfiles', '*.mat'));
    eset_name = ['eset', num2str(folder_index)];  % dynamic structure name
    % load .mat files into esets
    esets.(eset_name) = ExperimentSet.fromMatFiles(fullfile(basedir, 'matfiles', {d(x).name}));  % create a structure esets to hold eset
    d_tracks = dir(fullfile(basedir, 'matfiles', '*tracks'));
    for k = 1 : length(x)  % loop for each expt in the eset
        disp('Loading tracks ...');
        tracks_mat = dir(fullfile(basedir, 'matfiles', d_tracks(x(k)).name, '*.mat'));
        clear s;
        for i = 1 : length(tracks_mat)
            load(fullfile(tracks_mat(i).folder, tracks_mat(i).name));
            s(i) = track;  % a structure to hold track from each .mat
        end
        esets.(eset_name).expt(k).track = s;  % save the track to x(k) or 1
    end
    
    disp('Sorting tracks into run, reorientation, and headswing...');
    esets.(eset_name).executeTrackFunction('segmentTrack');
    
end


%% Draw basic figures, like led1Val, track_length, turn_start_time, num_moving_maggots
tperiod = 10;  % depend on the name of .mat file loaded, '_18_', or from the plot led2Val-Time-----------------
frame_rate = 20;  % number of frames per second
download = true;  % if true, plot and save figures including led1Val, Track length, Num of maggots; otherwise not plotting, but all valuables are prepared
save_data = true;

for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    basedir = basedir_cell{folder_index};
    eset_name = ['eset', num2str(folder_index)];
    x = x_cell{folder_index};
    d = dir(fullfile(basedir, 'matfiles', '*.mat'));
    
    for k = 1 : length(x)  % loop for each expt in the eset

        % eset.expt.globalQuantity.fieldname to check what field is, the field doesn't depend on which track
        field_name = {esets.(eset_name).expt(k).globalQuantity.fieldname};
        GQled1Val = esets.(eset_name).expt(k).globalQuantity(strcmp({esets.(eset_name).expt(k).globalQuantity.fieldname}, {'led1Val'}));
        GQled2Val = esets.(eset_name).expt(k).globalQuantity(strcmp({esets.(eset_name).expt(k).globalQuantity.fieldname}, {'led2Val'}));
        
        % create a new globalquantity which is square wave across the whole experiment time
        esets.(eset_name).expt(k).addTonToff('led1Val', 'square', 'period', tperiod);  % create time on/off field based a global quantity fieldname 'led1Val'
        GQton = esets.(eset_name).expt(k).globalQuantity(strcmp({esets.(eset_name).expt(k).globalQuantity.fieldname}, {'led1Val_ton'}));
        GQtoff = esets.(eset_name).expt(k).globalQuantity(strcmp({esets.(eset_name).expt(k).globalQuantity.fieldname}, {'led1Val_toff'}));

        % correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
        mkdir(fullfile(basedir, ['results', d(x(k)).name(end-16:end-4)]));
        if download
            % correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
            figure;
            plot(GQled2Val.xData/60, GQled1Val.yData, 'r'); hold on;
            plot(GQled2Val.xData/60, GQled2Val.yData, 'b');
            plot(GQled2Val.xData/60, GQtoff.yData, 'k');
            xlabel('Time (min)'); legend('led1Val', 'led2Val', 'toff'); hold off;
            pause;
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\led1Val');
            savefig(gcf,savename);
            close;
        end

        
        % plot the number of recognized maggots verses frame. Maggots in collision,
        % out of ROI, or discarded by many different tests in processBIN won't be recognized.
        ntracks_frame = zeros(1, length(esets.(eset_name).expt(k).elapsedTime));
        start_frame_tracks = [esets.(eset_name).expt(k).track.startFrame] + 1;  % min is 1
        end_frame_tracks = [esets.(eset_name).expt(k).track.endFrame] + 1;
        for j = 1:length(esets.(eset_name).expt(k).track)
            ntracks_frame(start_frame_tracks(j) : end) = ntracks_frame(start_frame_tracks(j) : end) + 1;
            ntracks_frame(end_frame_tracks(j) : end) = ntracks_frame(end_frame_tracks(j) : end) - 1;
        end
        if download
            figure; plot(ntracks_frame);
            xlabel('Frame number'); ylabel('Number of recognized maggots'); pause;
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\num_maggots_recognized');
            savefig(gcf, savename); close;
        end
    
        % if save the useful basic variables, like ntracks_frame to data
        if save_data
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
            if isfile(savename)
                save(savename, 'ntracks_frame', 'GQled1Val', '-append')
            else
                save(savename, 'ntracks_frame', 'GQled1Val')
            end
        end

    end  % end looping different experiment expt(k)
    
end  % end looping different folder


%% Plot behavior info for each experiment
download = true;  % if download true, save figures; if false, don't save.
texpt = 250;  % second, experiment time 
plot_mean_speed_all_t = true;  % mean speed of all tracks verses experiment time
plot_mean_speed_all_t_dots = false;  % get mean speed and std in a wide bin size.
plot_turnrate_texpt = true;
% plot across different expts
plot_mean_speed_all_expt_t_dots = true;  % get mean speed and std in a wide bin size.
plot_mean_speed_all_expt_t = true;
save_data = true;
for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    
    basedir = basedir_cell{folder_index};
    eset_name = ['eset', num2str(folder_index)];
    x = x_cell{folder_index};
    
    for k = 1 : length(x)  % loop for each expt in the eset
        
        if plot_mean_speed_all_t
            stepsize = 0.1; % second
            v_all = esets.(eset_name).expt(k).gatherField('speed') * 60;  % cm/min, [v11, v12, ..., v21, v22, ..., ], with vij the speed of track i at frame j
            t_all = esets.(eset_name).expt(k).gatherField('eti');  % interpolated time (s) for each frame of each track, 1-by-(sum of all tracks' frame) 
            [x_ton,y_v, ~, stdev] = meanyvsx(t_all, v_all, 0:stepsize:texpt);
            uppercurve = y_v + 0.5*stdev;
            lowercurve = y_v - 0.5*stdev;
            x_tofill = [x_ton, fliplr(x_ton)];  % the x axis of the ploygon to fill
            y_tofill = [lowercurve, fliplr(uppercurve)];
            figure; hold on;
            yyaxis left
            pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none');
            plot(x_ton, y_v, '-')
            xlabel('Time (s)'); ylabel('Speed (cm/min)');
            yyaxis right
            plot(GQled1Val.xData, GQled1Val.yData, 'r')
            ylabel('led1Val (PWM)'); 
            hold off;
            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\speed_mean_all_t');
                savefig(gcf, savename);
            end
            close;
        end


        if plot_mean_speed_all_t_dots
            stepsize = 1; % second
            v_all = esets.(eset_name).expt(k).gatherField('speed') * 60;  % cm/min, [v11, v12, ..., v21, v22, ..., ], with vij the speed of track i at frame j
            t_all = esets.(eset_name).expt(k).gatherField('eti');  % interpolated time (s) for each frame of each track, 1-by-(sum of all tracks' frame) 
            [x_ton,y_v, ~, stdev] = meanyvsx(t_all, v_all, 0:stepsize:texpt);
            figure; hold on;
            yyaxis left
%             plot(x_ton, y_v, 'ro');
            errorbar(x_ton, y_v, stdev/2, 'ok--');
            xlabel('Time (s)'); ylabel('Speed (cm/min)'); 
            yyaxis right
            plot(GQled1Val.xData, GQled1Val.yData, 'r')
            ylabel('led1Val (PWM)'); 
            hold off;
            pause; 
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\speed_mean_all_t_sparse');
                savefig(gcf, savename);
            end
            close;
        end


        if plot_turnrate_texpt
            % turn rate of all tracks at a certain period of time
            t = esets.(eset_name).expt(k).track;
            turnStartTime_total=[];
            stepsize = 0.5; binsize = 1;  % seconds
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
            load(savename, 'ntracks_frame');
            for j = 1: length(t)
                turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');
                turnStartTime_total = [turnStartTime_total turnStartTime];
            end
            [turnrate, turnrate_error] = rate_from_time(turnStartTime_total, texpt, stepsize, binsize);
            ntracks_frame_downsampling = ntracks_frame(round(linspace(1, length(ntracks_frame), length(turnrate))));
            turnrate = turnrate ./ ntracks_frame_downsampling * 60;
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
            title([num2str(length(turnStartTime_total)), ' turns, ', num2str(max(ntracks_frame)), ' larvae']);
            yyaxis right
            plot(GQled1Val.xData, GQled1Val.yData, 'r')
            ylabel('led1Val (PWM)');
            hold off;

            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\rate_turn_texpt');
                savefig(gcf, savename); 
            end
            close;
        end % end plot_pturn
        
        % if save the useful variables, like ntracks_frame, GQled1Val, turnStartTime_total to data
        if save_data
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
            if isfile(savename)
                save(savename, 'turnStartTime_total', '-append')
            else
                save(savename, 'turnStartTime_total')
            end
        end

        disp([d(x(k)).name(end-15:end-4), ' is done']);

    end  % end expt loop
    
    % plot with the combined data of all expts
    if plot_mean_speed_all_expt_t_dots
        stepsize = 1; % second
        v_all = esets.(eset_name).gatherField('speed') * 60;  % cm/min, [v11, v12, ..., v21, v22, ..., ], with vij the speed of track i at frame j
        t_all = esets.(eset_name).gatherField('eti');  % interpolated time (s) for each frame of each track, 1-by-(sum of all tracks' frame) 
        [x_ton,y_v, ~, stdev] = meanyvsx(t_all, v_all, 0:stepsize:texpt);
        figure; hold on;
        yyaxis left
%             plot(x_ton, y_v, 'ro');
        errorbar(x_ton, y_v, stdev/2, 'ok--');
        xlabel('Time (s)'); ylabel('Speed (cm/min)'); 
        yyaxis right
        plot(GQled1Val.xData, GQled1Val.yData, 'r')
        ylabel('led1Val (PWM)'); 
        hold off;
        pause; 
        if download
            mkdir(fullfile(basedir, 'results'));
            savename = strcat(basedir,'\results', '\speed_mean_t_allExpt_sparse');
            savefig(gcf, savename);
        end
        close;
    end  % end plot_mean_speed_all_expt_t_dots
    
    if plot_mean_speed_all_expt_t
        stepsize = 0.1; % second
        v_all = esets.(eset_name).gatherField('speed') * 60;  % cm/min, [v11, v12, ..., v21, v22, ..., ], with vij the speed of track i at frame j
        t_all = esets.(eset_name).gatherField('eti');  % interpolated time (s) for each frame of each track, 1-by-(sum of all tracks' frame) 
        [x_ton,y_v, ~, stdev] = meanyvsx(t_all, v_all, 0:stepsize:texpt);
        uppercurve = y_v + 0.5*stdev;
        lowercurve = y_v - 0.5*stdev;
        x_tofill = [x_ton, fliplr(x_ton)];  % the x axis of the ploygon to fill
        y_tofill = [lowercurve, fliplr(uppercurve)];
        figure; hold on;
        yyaxis left
        pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none');
        plot(x_ton, y_v, '-')
        xlabel('Time (s)'); ylabel('Speed (cm/min)'); 
        yyaxis right
        plot(GQled1Val.xData, GQled1Val.yData, 'r')
        ylabel('led1Val (PWM)'); 
        hold off;
        pause; 
        if download
            mkdir(fullfile(basedir, 'results'));
            savename = strcat(basedir,'\results', '\speed_mean_t_allExpt');
            savefig(gcf, savename);
        end
        close;
    end  % end plot_mean_speed_all_expt_t_dots

end  % end base folder loop


%% copy this file to all results folders that used it
for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    basedir = basedir_cell{folder_index};
    x = x_cell{folder_index};
    d = dir(fullfile(basedir, 'matfiles', '*.mat'));
    for k = 1 : length(x)  % loop for each expt in the eset
        savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)]);
        copyfile('load_multi_data_a27h.m', savename);
    end
end