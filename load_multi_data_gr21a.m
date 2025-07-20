basedir_cell = {
   'G:\AS-Filer\PHY\mmihovil\Shared\data_variability\variability_extracted\T_Re_Sq_219to436P_15_2_3#T_Bl_Sq_2to7P_15_1_3'
    };
x_cell = {
    1:21
    %[25, 26, 27, 28, 29, 30, 31, 32, 33, 34]  % load the x-th set of data from the first basedir to analyze.
    };
pause('off');  % 'on'---ask the user to press any key to save the figure, and continue; 'off'--directly save without asking

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
tperiod = 15;  % in seconds, period of stimulation time
latest_start = 120;  % seconds, select the tracks that start earlier than latest_start time. 80% length
earliest_end = 1680;  % seconds, select the tracks that end later than earliest_end time.
t_stim_start = [0, 600, 1200];  % start time (s) of each intensity of stimulation
t_stim_end = [600, 1200, 1800];
stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start
frame_rate = 20;  % number of frames per second
download = true;  % if true, plot and save figures including led1Val, Track length, Num of maggots; otherwise not plotting, but all valuables are prepared
create_larvae =true;  % create/edit structure larvae
save_data = true;  % save larvae to data.mat
% add led12Val, and select long tracks
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
        xdata = GQled2Val.xData;  % xdata of all Global Quantity ledVals should be the same
        ydata = GQled1Val.yData + GQled2Val.yData;  % LED intensity in PWM
        index = find(GQled1Val.yData==0, 1, 'last');  % find the index of the last zero element of led1Val, this index may be different from the initial setting because of hardware noise
        ydata(1:index) = ydata(1:index) + 60;  % this is to make the square wave fluctruate around a center quantity, use index to make the combined square wave cleaner

        % add ton toff
        esets.(eset_name).expt(k).addGlobalQuantity('eti', 'led12Val', xdata, ydata);  % create a man-made global field to add ton/toff
        esets.(eset_name).expt(k).addTonToff('led12Val', 'square', 'period', tperiod);  % create time on/off field based a global quantity fieldname 'led2Val'
        GQton = esets.(eset_name).expt(k).globalQuantity(strcmp({esets.(eset_name).expt(k).globalQuantity.fieldname}, {'led12Val_ton'}));
        GQtoff = esets.(eset_name).expt(k).globalQuantity(strcmp({esets.(eset_name).expt(k).globalQuantity.fieldname}, {'led12Val_toff'}));
%         figure; plot(xdata, ydata)

        mkdir(fullfile(basedir, ['results', d(x(k)).name(end-16:end-4)]));
        if download
            % correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
            figure;
            plot(GQled2Val.xData/60, GQled1Val.yData, 'r'); hold on;
            plot(GQled2Val.xData/60, GQled2Val.yData, 'b');
            plot(GQled2Val.xData/60, GQtoff.yData, 'k');
            plot(GQled2Val.xData/60, GQtoff.yData+60, 'k'); hold off;
            xlabel('Time (min)'); legend('led1Val', 'led2Val', 'toff');
            pause;
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\led12Val_toff');
            savefig(gcf,savename);
            close;
        end
        
        if download
            figure;
            histogram(round(esets.(eset_name).expt(k).elapsedTime([esets.(eset_name).expt(k).track.npts])/tperiod), 0:20:(t_stim_end(end)/tperiod+1));  % round 60.001 to 60
            xlabel('Number of Periods'); ylabel('Number of Tracks'); title(['Histogram of The Length of All ', num2str(length(esets.(eset_name).expt(k).track)), ' Tracks']);
            pause;
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\track_length');
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
        savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
        if isfile(savename)
            save(savename, 'ntracks_frame', 'GQled1Val', 'GQled2Val', '-append')
        else
            save(savename, 'ntracks_frame', 'GQled1Val', 'GQled2Val')
        end
        if download
            figure; plot(ntracks_frame);
            xlabel('Frame number'); ylabel('Number of recognized maggots'); pause;
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\num_maggots_recognized');
            savefig(gcf, savename); close;
        end


        % plot start time and end time of all tracks
        if download
            figure;
            for j = 1: length(esets.(eset_name).expt(k).track)
                    plot(esets.(eset_name).expt(k).track(j).startFrame + 1, j, 'bo'); hold on;
                    plot(esets.(eset_name).expt(k).track(j).endFrame, j, 'rx'); hold on;
            end
            xlabel('Frame number'); ylabel('Index of tracks'); hold off; pause;
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\track_start_end');
            savefig(gcf, savename); close;
        end

        % create structure larvae, add basic info
        if create_larvae
            t = esets.(eset_name).expt(k).track;
            t = t(([t.startFrame] <= latest_start * frame_rate) & ([t.endFrame] >= earliest_end * frame_rate));
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
            clear larvae  % if firstly run, to clear the larvae of last ex pt
            larvae = struct;  % initialize an empty structure in case no long track in an expt.
            load(savename, 'larvae');  % if edit larvae
            for j = 1 : length(t)
                larva_index = ['larva', num2str(j)];
                % create a new structure larvae, and add basic info
                larvae.(larva_index).trackNum = t(j).trackNum;  % integer
                larvae.(larva_index).expt_info = d(x(k)).name(1:end-4);  % string
                larvae.(larva_index).startFrame = t(j).startFrame;
                larvae.(larva_index).endFrame = t(j).endFrame;
                % add longturn_time (s), unlabeled_time (s), startLabelTime (s)
                t_running = t(j).getSubFieldDQ('run', 'eti') ;  % second, time when larva is running
                t_turning = t(j).getSubFieldDQ('reorientation', 'eti') ;  % second, time when larva is turning
                % figure; plot(t_running, ones(1, length(t_running)), '.k', t_turning, 2 *ones(1, length(t_turning)), '.r' ); ylim([-10, 10])
                turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % time (s) not in period
                turnEndTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'end');  % time (s) not in period
                duration_turn = turnEndTime - turnStartTime;  % duration of all turns
                larvae.(larva_index).turnStartTime = turnStartTime;
                larvae.(larva_index).turnEndTime = turnEndTime;
                larvae.(larva_index).longturn_time = sum(duration_turn(duration_turn > 1*tperiod));  % second, double, total pause (labeled as turn but too long,  longer than a period) time of a track
                larvae.(larva_index).unlabeled_time = (t(j).npts - length(t_running) - length(t_turning)) / frame_rate;  % second, double, total time not labeled as either run or turn, invalid, mostly in the beginning
                larvae.(larva_index).startLabelTime = min(t_running(1), t_turning(1));  % replace this with track start time to get nperiod in the future
                larvae.(larva_index).bodyLength = median(t(j).dq.spineLength);  % cm, float, median of spine length is body length
            end  % end looping long tracks
        end  % end create_larvae

        % if save the structure larvae to data
        if save_data
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
            if isfile(savename)
                save(savename, 'larvae', '-append')
            else
                save(savename, 'larvae')
            end
        end

        % display the track length info
        t = esets.(eset_name).expt(k).track;
        t = t(([t.startFrame] <= latest_start * frame_rate) & ([t.endFrame] >= earliest_end * frame_rate));
        disp([num2str(length(t)), ' long tracks out of ', num2str(length(esets.(eset_name).expt(k).track)), ' tracks, from ', num2str(max(ntracks_frame)), ' moving maggots']);

    end  % end looping different experiment expt(k)
    
end  % end looping different folder


%% Plot behavior info for each experiment
download = true;  % if download true, save figures; if false, don't save.
save_data = true;  % save larvae into data.m
plot_pturn = true;
tbin = 3;  edges = 0:tbin: tperiod;
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
plot_turnrate = true;
plot_turnrate_individual = true;
plot_speed_individual = true;
plot_speed_run_individual = true;
plot_angular_speed_individual = true;
speed_each_larva_each_stim = true;
track_pturn =  true;  % save pturn of all tracks into structure tracks
track_info = true;  % save the temporal and spatial info about start/end of all tracks
stepsize = 0.1; binsize = 0.5;  % seconds
texpt = t_stim_end(end);  % used to plot v-texpt
for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    
    basedir = basedir_cell{folder_index};
    eset_name = ['eset', num2str(folder_index)];
    x = x_cell{folder_index};
    
    for k = setdiff(1 : length(x), 4)  % loop for each expt in the eset if 4 is []
        
        if plot_pturn
            % save long tracks to t
            t = esets.(eset_name).expt(k).track;
            t = t(([t.startFrame] <= latest_start * frame_rate) & ([t.endFrame] >= earliest_end * frame_rate));
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
            clear larvae;
            if isfile(savename)
                load(savename, 'larvae');
            else
                clear larvae;
            end
            turnStart_mean = zeros(1, length(xbar));  % the average turnTime that falls to one bin
            turnStart_std = zeros(1, length(xbar));  % standard error = standard deviation / sqrt(N)
            for i = 1 : length(t_stim_start)  % ith intensity of stimulation--------------------
                for j = 1 : length(t)
                    % 'led12Val_ton' fails some time if stimulation isn't strict square wave
                    turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_ton', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % turn start time in period, ton means period starts with light on
                    turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % time (s) not in period
                    turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i))); %only keep the reorientation whose start time falls into the i-th intensity of stimulation
                    % number of period of stimulation that falls into certain intensity
                    t_start = max([t_stim_start(i), esets.(eset_name).expt(k).elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
                    t_end = min([t_stim_end(i), esets.(eset_name).expt(k).elapsedTime(t(j).endFrame-2)]);  % temporal - 2
                    nperiod = (t_end - t_start) / tperiod;
                    index_subplot = j + (i-1)*length(t);
                    ax(index_subplot) = subplot(length(t_stim_start), length(t), index_subplot);
                    [N, ~] = histcounts(turnStart, edges);  % assume larvae can only turn one time within tbin
                    pturn = N/nperiod;
                    pturn_error = sqrt(N)/nperiod;
                    for r = 1 : (length(edges) - 1)
                        turnStart_bin = turnStart((turnStart > edges(r)) & (turnStart < edges(r+1)));
                        turnStart_mean(r) = mean(turnStart_bin);  % the mean turn start time in period within each bin
                        turnStart_std(r) = std(turnStart_bin, 1);
                    end
                    bar(xbar, pturn, 1);  hold on;% the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
                    errorbar(turnStart_mean, pturn, pturn_error / 2, 'k.'); % the length of the error bar is sqrt(N)/nperiod
                    errorbar(turnStart_mean, pturn, turnStart_std/2, 'horizontal', 'k.'); hold off;
                    xline(6, 'k--');
                    xticks(edges);    ylim([0, 1]);  % comment ylim first, change to ymax at the second run ----------------

                    larva_index = ['larva', num2str(j)];
                    % create a new structure larvae
                    larvae.(larva_index).nperiod.(stim_color{i}) = nperiod;
                    larvae.(larva_index).pturn.(stim_color{i}) = N/nperiod;
                    larvae.(larva_index).pturn_error.(stim_color{i}) = sqrt(N)/nperiod;
                    larvae.(larva_index).turnStart.(stim_color{i}) = turnStart;
                    larvae.(larva_index).turnStartTime = turnStartTime;
                    larvae.(larva_index).turnStart_mean.(stim_color{i}) = turnStart_mean;
                    larvae.(larva_index).turnStart_std.(stim_color{i}) = turnStart_std;
                    if N(1)/nperiod - mean(N(3:end)/nperiod) >= 0.2  % if the first bin of pturn is much larger than the low intensity bins, call it response
                        larvae.(larva_index).response.(stim_color{i}) = '1';
                    elseif mean(N(1:2)/nperiod) - mean(N(3:end)/nperiod) >= 0.2  % if the first 2 bins of pturn are much larger than the rest
                        larvae.(larva_index).response.(stim_color{i}) = '1';
                    else
                        larvae.(larva_index).response.(stim_color{i}) = '0';
                    end

                    title([stim_color{i}, ' (', num2str(t(j).trackNum), ', ', num2str(round(nperiod,1)), ', ', num2str(sum(N)), ', ', larvae.(larva_index).response.(stim_color{i}), ')']);

                end
            end
            xlabel(ax(1), 'Time in period (s)'); ylabel(ax(1), 'Probability of Starting to Turn'); 
            sgtitle('(Track number, Number of Stimulation, Number of Turn, Response)');
            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\p_turn_ton');
                savefig(gcf, savename); 
            end
            close;
        end % end plot_pturn
        
        
        if plot_turnrate
            % turn rate of all tracks at a certain period of time
            t = esets.(eset_name).expt(k).track;
            figure;
            % only keep the reorientation whose start time falls into the i-th intensity of stimulation
            for i = 1: length(t_stim_start)
                nperiod = 0;
                turnStart_total=[];
                for j = 1: length(t)

                    turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');
                    turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_toff', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');
                    turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i)));
                    turnStart_total = [turnStart_total turnStart];
                    if (t(j).startFrame < t_stim_end(i)*frame_rate) && (t(j).endFrame > t_stim_start(i)*frame_rate)
                        nperiod = nperiod + (min(t(j).endFrame, t_stim_end(i)*frame_rate) - max(t(j).startFrame, t_stim_start(i)*frame_rate)) / frame_rate / tperiod;
                    end
                end
                turnrate = rate_from_time(turnStart_total, tperiod, stepsize, binsize) ./ double(nperiod) * 60;
                time_timestep = (0 : fix(tperiod/stepsize)) * stepsize;
                ax(i) = subplot(length(t_stim_start),1,i);
                plot(time_timestep, turnrate);
                xline(9, 'k--');  % stimulus change time at -------------
                xlabel('Reorientation Start Time in Period (s)'); ylabel('Reorientation Rate (per min)'); 
                title([stim_color{i}, ', ', num2str(length(turnStart_total)), ' turns, in ', num2str(nperiod), ' periods']);
            end
            sgtitle(['Step size = ', num2str(stepsize), ', bin size = ', num2str(binsize)]);
            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\rate_turn_toff');
                savefig(gcf, savename); 
            end
            close;
        end % end plot_pturn

        if plot_turnrate_individual
            % Variability, each track in each column, each stimulation period in each row
            % save info into a existing structure larvae
            t = esets.(eset_name).expt(k).track;
            t = t(([t.startFrame] <= latest_start * frame_rate) & ([t.endFrame] >= earliest_end * frame_rate));
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
%             if isfile(savename)  % uncomment when only run this block
%                 load(savename, 'larvae');
%             end
            figure;
            for i = 1 : length(t_stim_start)  % ith intensity of stimulation--------------------
                for j = 1 : length(t)  % j-th track
                    turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_ton', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % turn start time in period, ton means period starts with light on
                    turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % time (s) not in period
                    turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i))); %only keep the reorientation whose start time falls into the i-th intensity of stimulation
                    %number of period of stimulation that falls into certain intensity
                    t_start = max([t_stim_start(i), esets.(eset_name).expt(k).elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th intensity of stimulation 
                    t_end = min([t_stim_end(i), esets.(eset_name).expt(k).elapsedTime(t(j).endFrame-2)]);  % temporal - 2
                    nperiod = (t_end - t_start) / tperiod;
            %         index_subplot = j*length(t_stim_start) - (length(t_stim_start)-i);  % to plot by column
            %         ax(index_subplot) = subplot(length(t), length(t_stim_start), index_subplot);
                    index_subplot = j + (i-1)*length(t);  % to put each track in each column, and each stimulation period in each row
                    ax(index_subplot) = subplot(length(t_stim_start), length(t), index_subplot);
                    turnrate = rate_from_time(turnStart, tperiod, stepsize, binsize) ./ double(nperiod) * 60;
                    time_timestep = (0 : fix(tperiod/stepsize)) * stepsize;
                    plot(time_timestep, turnrate); xline(6, 'k--');
                    larva_index = ['larva', num2str(j)];
                    title([stim_color{i}, ' (', num2str(t(j).trackNum), ', ', num2str(round(nperiod,1)), ', ', num2str(length(turnStart)), ', ', larvae.(larva_index).response.(stim_color{i}),  ')']);
                    
                    if save_data
                        larvae.(larva_index).turnrate.(stim_color{i}) = turnrate;
                    end
                end  % end looping each track j in the expt
            end  % end looping each stimulation condition i
            linkaxes(ax, 'y');  % align subplots with y axis
            xlabel(ax(1), 'ton (s)'); ylabel(ax(1), 'Turn rate (min^{-1})'); sgtitle(['Stepsize = ', num2str(stepsize), ', binsize = ', num2str(binsize), ' (track number, number of period, number of turn, response)']);
            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\rate_turn_ton_individual');
                savefig(gcf, savename); 
            end
            close;
        end  % end plot turn_rate_Individual

        
        if plot_speed_individual
            % Variability, each track in each column, each stimulation period in each row
            % save speed into a existing structure larvae
            t = esets.(eset_name).expt(k).track;
            t = t(([t.startFrame] <= latest_start * frame_rate) & ([t.endFrame] >= earliest_end * frame_rate));
%             savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
%             if isfile(savename)  % uncomment when only run this block
%                 load(savename, 'larvae');
%             end
            figure;  % to get one average speed for stepsize seconds
            % only keep the speed whose start time falls into the i-th intensity of stimulation
            for i = 1: length(t_stim_start)
                for j = 1: length(t)
                    larva_index = ['larva', num2str(j)];
                    turnStart = larvae.(larva_index).turnStart.(stim_color{i});
                    bodyL = larvae.(larva_index).bodyLength;  % cm

                    v_frame = t(j).dq.speed * 60;  % cm/min
                    ton_frame = t(j).dq.led12Val_ton;  % interpolated time (s) for each frame of track j, 1-by-(number of the track's frame) 
                    t_frame = t(j).dq.eti;
                    % select the start time in period of turns that happen between t_stim_start(i) and t_stim_start(i)
                    v_frame = v_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
                    ton_frame = ton_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
                    %number of period of stimulation that falls into certain intensity
                    t_start = max([t_stim_start(i), esets.(eset_name).expt(k).elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
                    t_end = min([t_stim_end(i), esets.(eset_name).expt(k).elapsedTime(t(j).endFrame-2)]);  % temporal - 2
                    nperiod = (t_end - t_start) / tperiod;
                    index_subplot = j + (i-1)*length(t);
                    ax(index_subplot) = subplot(length(t_stim_start), length(t), index_subplot);
                    [x_ton,y_v, ~, stdev] = meanyvsx(ton_frame, v_frame, 0:stepsize:tperiod);
                    uppercurve = y_v + 0.5*stdev;
                    lowercurve = y_v - 0.5*stdev;
                    x_tofill = [x_ton, fliplr(x_ton)];  % the x axis of the polygon to fill
                    y_tofill = [lowercurve, fliplr(uppercurve)];
                    pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none'); hold on;  % no edges for the patch. use y_tofill/bodyL if...
                    plot(x_ton, y_v, 'Color', 0.2*[1 1 1]); xline(6, '--'); hold off;  % use y_v /bodyL if use bodyL as unit instead of cm
                    xlabel('ton (s)'); ylabel('Speed (cm/min)');

                    if save_data  % add new contents to the structure larvae
                        larvae.(larva_index).speed.(stim_color{i}) = y_v;
                        larvae.(larva_index).speed_std.(stim_color{i}) = stdev;
                    end
            
                    title([stim_color{i}, ' (', num2str(t(j).trackNum), ', ', num2str(round(nperiod,1)), ', ', num2str(length(turnStart)), ', ', larvae.(larva_index).response.(stim_color{i}) ')']);
                end
            end
            sgtitle(['Step size = ', num2str(stepsize), ' s, (Track number, Number of Stimulation, Number of Turn, Response)']);  linkaxes(ax, 'y');  % align subplots with y axis
            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\v_ton_stim_individual');
                savefig(gcf, savename); 
            end
            close;
        end  % end plot_speed_individual


        if plot_speed_run_individual
            % Variability, each track in each column, each stimulation period in each row
            % save speed of each long track at a certain period of time of stimulation into a existing structure larvae
            t = esets.(eset_name).expt(k).track;
            t = t(([t.startFrame] <= latest_start * frame_rate) & ([t.endFrame] >= earliest_end * frame_rate));
%             savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
%             if isfile(savename)  % uncomment when only run this block
%                 load(savename, 'larvae');
%             end
            figure;
            % only keep the speed whose start time falls into the i-th intensity of stimulation
            for i = 1: length(t_stim_start)
                for j = 1: length(t)
                    v_frame = t(j).getSubFieldDQ('run', 'speed') * 60 ;  % cm/min, speed of every frame when larva is running
                    ton_frame = t(j).getSubFieldDQ('run', 'led12Val_ton');  % second
                    t_frame = t(j).getSubFieldDQ('run', 'eti');
                    
                    % select the start time in period of turns that happen between t_stim_start(i) and t_stim_start(i)
                    if i < length(t_stim_start)
                        v_frame = v_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
                        ton_frame = ton_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
                    else
                        v_frame = v_frame((t_stim_start(i) <= t_frame) & (t_frame <= t_stim_end(i)));
                        ton_frame = ton_frame((t_stim_start(i) <= t_frame) & (t_frame <= t_stim_end(i)));
                    end
            
                    % number of period of stimulation that fa  to certain intensity
                    if (t(j).startFrame < t_stim_end(i)*frame_rate) && (t(j).endFrame > t_stim_start(i)*frame_rate)
                        t_start = max([t_stim_start(i), esets.(eset_name).expt(k).elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
                        t_end = min([t_stim_end(i), esets.(eset_name).expt(k).elapsedTime(t(j).endFrame-2)]);  % temporal - 2
                        nperiod = (t_end - t_start) / tperiod;
                    else 
                        nperiod = 0;
                    end
                    index_subplot = j + (i-1)*length(t);
                    ax(index_subplot) = subplot(length(t_stim_start), length(t), index_subplot);
                    [x_ton,y_v, ~, stdev] = meanyvsx(ton_frame, v_frame, 0:stepsize:tperiod);  % std
                    uppercurve = y_v + 0.5*stdev;
                    lowercurve = y_v - 0.5*stdev;
                    x_tofill = [x_ton, fliplr(x_ton)];  % the x axis of the ploygon to fill
                    y_tofill = [lowercurve, fliplr(uppercurve)];
                    pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none'); hold on;  % no edges for the patch
                    plot(x_ton, y_v, 'Color', 0.2*[1 1 1]); xline(6, '--'); hold off;
                    xlabel('ton (s)'); ylabel('Run speed (cm/min)');
            
                    % add new contents to the structure larvae, or read info from
                    % larvae
                    larva_index = ['larva', num2str(j)];
                    turnStart = larvae.(larva_index).turnStart.(stim_color{i});
                    if save_data  % add new contents to the structure larvae
                        larvae.(larva_index).speed_run.(stim_color{i}) = y_v;
                        larvae.(larva_index).speed_run_std.(stim_color{i}) = stdev;
                    end
            
                    title([stim_color{i}, ' (', num2str(t(j).trackNum), ', ', num2str(round(nperiod,1)), ', ', num2str(length(turnStart)), ', ', larvae.(larva_index).response.(stim_color{i}) ')']);
            
                end
            end
            sgtitle(['Step size = ', num2str(stepsize), ' s, (track index, Nperiod, Nturn, response)']);  linkaxes(ax, 'y');  % align subplots with y axis
            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\vRun_ton_stim_individual');
                savefig(gcf, savename); 
            end
            close;
        end  % end plotting speed_run_individual

        if plot_angular_speed_individual
            % Variability, each track in each column, each stimulation period in each row
            % save angular speed into a existing structure larvae
            t = esets.(eset_name).expt(k).track;
            t = t(([t.startFrame] <= latest_start * frame_rate) & ([t.endFrame] >= earliest_end * frame_rate));
%             savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
%             if isfile(savename)  % uncomment when only run this block
%                 load(savename, 'larvae');
%             end
            figure;  % to get one average speed for stepsize seconds
            % only keep the speed whose start time falls into the i-th intensity of stimulation
            for i = 1: length(t_stim_start)
                for j = 1: length(t)
                    larva_index = ['larva', num2str(j)];
                    turnStart = larvae.(larva_index).turnStart.(stim_color{i});
                    bodyL = larvae.(larva_index).bodyLength;  % cm

                    w_frame = t(j).dq.deltatheta;  % omega, angular speed rad/s, angle [-pi, pi]
%                     w_frame = abs(w_frame);  % absolute value, [0, pi]
                    ton_frame = t(j).dq.led12Val_ton;  % interpolated time (s) for each frame of track j, 1-by-(number of the track's frame) 
                    t_frame = t(j).dq.eti;
                    % select the start time in period of turns that happen between t_stim_start(i) and t_stim_start(i)
                    w_frame = w_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
                    ton_frame = ton_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
                    %number of period of stimulation that falls into certain intensity
                    t_start = max([t_stim_start(i), esets.(eset_name).expt(k).elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
                    t_end = min([t_stim_end(i), esets.(eset_name).expt(k).elapsedTime(t(j).endFrame-2)]);  % temporal - 2
                    nperiod = (t_end - t_start) / tperiod;
                    index_subplot = j + (i-1)*length(t);
                    ax(index_subplot) = subplot(length(t_stim_start), length(t), index_subplot);
                    [x_ton,y_w, ~, stdev] = meanyvsx(ton_frame, w_frame, 0:stepsize:tperiod);
                    uppercurve = y_w + 0.5*stdev;
                    lowercurve = y_w - 0.5*stdev;
                    x_tofill = [x_ton, fliplr(x_ton)];  % the x axis of the polygon to fill
                    y_tofill = [lowercurve, fliplr(uppercurve)];
                    pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none'); hold on;  % no edges for the patch. use y_tofill/bodyL if...
                    plot(x_ton, y_w, 'Color', 0.2*[1 1 1]); xline(6, '--'); hold off;  % use y_v /bodyL if use bodyL as unit instead of cm
                    xlabel('ton (s)'); ylabel('Angular Speed (rad/s)');

                    if save_data  % add new contents to the structure larvae
                        larvae.(larva_index).angular_speed.(stim_color{i}) = y_w;
                        larvae.(larva_index).angular_speed_std.(stim_color{i}) = stdev;
                    end
            
                    title([stim_color{i}, ' (', num2str(t(j).trackNum), ', ', num2str(round(nperiod,1)), ', ', num2str(length(turnStart)), ', ', larvae.(larva_index).response.(stim_color{i}) ')']);
                end
            end
            sgtitle(['Step size = ', num2str(stepsize), ' s, (Track number, Number of Stimulation, Number of Turn, Response)']);  linkaxes(ax, 'y');  % align subplots with y axis
            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\w_ton_stim_individual');
                savefig(gcf, savename); 
            end
            close;
        end  % end plot_angular_speed_individual


        if speed_each_larva_each_stim
%             savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');  % uncomment when only run this block
%             if isfile(savename)
%                 load(savename, 'larvae');
%             end
            t = esets.(eset_name).expt(k).track;
            t = t(([t.startFrame] <= latest_start * frame_rate) & ([t.endFrame] >= earliest_end * frame_rate));
            for i = 1: length(t_stim_start)  % i-th stimulation
                for j = 1: length(t)  % track j
                    larva_index = ['larva', num2str(j)];
                    v_frame = t(j).dq.speed * 60;  % cm/min
                    ton_frame = t(j).dq.led12Val_ton;  % interpolated time (s) for each frame of track j, 1-by-(number of the track's frame) 
                    t_frame = t(j).dq.eti;
                    % select the v_frame and ton_frame that fall between t_stim_start(i) and t_stim_end(i)
                    v_frame = v_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
                    ton_frame = ton_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
                    interest = (14 < ton_frame) & (ton_frame < 15);  % logic array, pick all within 1s before stimulation on
                    % the length of the array is dependent on the number of periods that the track has spanned on
                    larvae.(larva_index).init_speed_each.(stim_color{i}) = getClusterAvg(v_frame, interest);
                end  % end looping different tracks t
            end  % end looping different stimulation i
        end  % end save speed_each_larva_each_stim to larvae


        % create a structure tracks that have all tracks' info saved,
        % including short tracks
        if track_pturn
            t = esets.(eset_name).expt(k).track;
%             savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
            clear tracks
%             if isfile(savename)
%                 load(savename, 'tracks');
%             end
            turnStart_mean = zeros(1, length(xbar));  % the average turnTime that falls to one bin
            turnStart_std = zeros(1, length(xbar));  % standard error = standard deviation / sqrt(N)
            for i = 1 : length(t_stim_start)  % ith intensity of stimulation--------------------
                for j = 1 : length(t)
                    % 'led12Val_ton' fails some time if stimulation isn't strict square wave
                    turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_ton', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % turn start time in period, ton means period starts with light on
                    turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % time (s) not in period
                    turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i))); %only keep the reorientation whose start time falls into the i-th intensity of stimulation
                    % number of period of stimulation that falls into certain intensity
                    t_start = max([t_stim_start(i), esets.(eset_name).expt(k).elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
                    t_end = min([t_stim_end(i), esets.(eset_name).expt(k).elapsedTime(t(j).endFrame-2)]);  % temporal - 2
                    nperiod = (t_end - t_start) / tperiod;
                    [N, ~] = histcounts(turnStart, edges);  % assume larvae can only turn one time within tbin
                    pturn = N/nperiod;
                    pturn_error = sqrt(N)/nperiod;
                    for r = 1 : (length(edges) - 1)
                        turnStart_bin = turnStart((turnStart > edges(r)) & (turnStart < edges(r+1)));
                        turnStart_mean(r) = mean(turnStart_bin);  % the mean turn start time in period within each bin
                        turnStart_std(r) = std(turnStart_bin, 1);
                    end

                    track_index = ['track', num2str(esets.(eset_name).expt(k).track(j).trackNum)];
                    % create a new structure tracks
                    tracks.(track_index).startFrame = esets.(eset_name).expt(k).track(j).startFrame;
                    tracks.(track_index).endFrame = esets.(eset_name).expt(k).track(j).endFrame;
                    tracks.(track_index).expt_info = d(x(k)).name(1:end-4);  % string
                    tracks.(track_index).nperiod.(stim_color{i}) = nperiod;
                    tracks.(track_index).pturn.(stim_color{i}) = N/nperiod;
                    tracks.(track_index).pturn_error.(stim_color{i}) = sqrt(N)/nperiod;
                    tracks.(track_index).turnStart.(stim_color{i}) = turnStart;
                    tracks.(track_index).turnStartTime = turnStartTime;
                    tracks.(track_index).turnStart_mean.(stim_color{i}) = turnStart_mean;
                    tracks.(track_index).turnStart_std.(stim_color{i}) = turnStart_std;
                    tracks.(track_index).speed_run_mean = mean(t(j).getSubFieldDQ('run', 'speed', 'position', 'mean') * 60);  % cm/min, float

                    if N(1)/nperiod - mean(N(3:end)/nperiod) >= 0.2  % if the first bin of pturn is much larger than the low intensity bins, call it response
                        tracks.(track_index).response.(stim_color{i}) = '1';
                    elseif mean(N(1:2)/nperiod) - mean(N(3:end)/nperiod) >= 0.2  % if the first 2 bins of pturn are much larger than the rest
                        tracks.(track_index).response.(stim_color{i}) = '1';
                    else
                        tracks.(track_index).response.(stim_color{i}) = '0';
                    end
                end  % end looping tracks
            end  % end looping different stimulation condition
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
            if isfile(savename)
                save(savename, 'tracks', '-append')
            else
                save(savename, 'tracks')
            end
        end % end track_info

        if track_info
            % Track info [track_ID, startFrame, startX(pixel), endY(pixel), endFrame,
            % endX(pixel), endY(pixel)]
            trackInfo = [];
            t = esets.(eset_name).expt(k).track;
            for j = 1: length(t)
                startFrame = t(j).startFrame;
                startLoc = transpose(camPtsFromRealPts(esets.(eset_name).expt(k).camcalinfo, t(j).pt(1).loc));  % [x, y]
                endFrame = t(j).endFrame;
                endLoc = transpose(camPtsFromRealPts(esets.(eset_name).expt(k).camcalinfo, t(j).pt(end).loc));  % [x, y]
                trackInfo = [trackInfo; j, startFrame, startLoc, endFrame, endLoc];
            end
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], ['\trackInfo', d(x(k)).name(end-16:end-4), '.csv']);
            writematrix( trackInfo, savename)
        end
 
        % if save the structure larvae to data
        if save_data
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\data.mat');
            if isfile(savename)
                save(savename, 'larvae', '-append')
            else
                save(savename, 'larvae')
            end
        end


            


        disp([d(x(k)).name(end-15:end-4), ' is done']);
    end  % end expt loop
    
    % analyze across all experiments together

end  % end base folder loop




%% copy this file to all results folders that used it
for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    basedir = basedir_cell{folder_index};
    x = x_cell{folder_index};
    d = dir(fullfile(basedir, 'matfiles', '*.mat'));
    for k = 1 : length(x)  % loop for each expt in the eset
        savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)]);
        copyfile('load_multi_data.m', savename);
    end
end