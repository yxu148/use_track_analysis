basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_1_3#T_Bl_Sq_2,5to6P_20_2_3'
    };
x_cell = {
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    %[25, 26, 27, 28, 29, 30, 31, 32, 33, 34]  % load the x-th set of data from the first basedir to analyze
    };

% load data, multiple eset, multiple expt, into esets
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


tperiod = 20;
nperiods = 69;  % select tracks that have [nperiods, Nperiods] length
Nperiods = 91;
download = true;  % if true, plot and save figures including led1Val, Track length, Num of maggots; otherwise not plotting, but all valuables are prepared
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
        esets.(eset_name).expt(k).addGlobalQuantity('eti', 'led12Val', xdata, ydata);  % create a man-made global field to add ton/toff
        esets.(eset_name).expt(k).addTonToff('led12Val', 'square');  % create time on/off field based a global quantity fieldname 'led2Val'
        GQton = esets.(eset_name).expt(k).globalQuantity(strcmp({esets.(eset_name).expt(k).globalQuantity.fieldname}, {'led12Val_ton'}));
        GQtoff = esets.(eset_name).expt(k).globalQuantity(strcmp({esets.(eset_name).expt(k).globalQuantity.fieldname}, {'led12Val_toff'}));
        
        % correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
        mkdir(fullfile(basedir, ['results', d(x(k)).name(end-16:end-4)]));
        if download
            % correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
            figure;
            plot(GQled2Val.xData/60, GQled1Val.yData, 'r'); hold on;
            plot(GQled2Val.xData/60, GQled2Val.yData, 'b'); hold on;
            plot(GQled2Val.xData/60, GQtoff.yData, 'k'); hold on;
            plot(GQled2Val.xData/60, GQtoff.yData+60, 'k'); hold off;
            xlabel('Time (min)'); legend('led1Val','led2Val', 'toff');
            pause;
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\led12Val_toff');
            savefig(gcf,savename);
            close;
        end
        
        Ntracks = size(esets.(eset_name).expt(k).track);  % 1-by-Number_of_Tracks ( number of maggots)
        if download
            figure;
            histogram(round(esets.(eset_name).expt(k).elapsedTime([esets.(eset_name).expt(k).track.npts])/tperiod), [0:10:90]);  % round 60.001 to 60
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
        if download
            figure; plot(ntracks_frame);
            xlabel('Frame number'); ylabel('Number of recognized maggots'); pause;
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\num_maggots_recognized');
            savefig(gcf, savename); close;
        end
        

    end
    
end


download = true;  % always plot, if download true, save; if false, don't save.
t_stim_start = [0, 600, 948];  % start time (s) of each intensity of stimulation
t_stim_end = [600, 948, 1548];
frame_rate = 20;  % number of frames per second
plot_pturn = true;  
tbin = 3;  edges = [0:tbin: tperiod];
%edges = [0,4,7,10,13,16,20];
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
plot_turnrate = true;
stepsize = 0.1; binsize = 0.5;
for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    
    basedir = basedir_cell{folder_index};
    eset_name = ['eset', num2str(folder_index)];
    x = x_cell{folder_index};
    
    for k = 1 : length(x)  % loop for each expt in the eset
        
        if plot_pturn
            % save long tracks to t
            t = esets.(eset_name).expt(k).track;
            minNpoints = nperiods * tperiod / (esets.(eset_name).expt(k).elapsedTime(end)/length(esets.(eset_name).expt(k).elapsedTime));
            maxNpoints = Nperiods * tperiod / (esets.(eset_name).expt(k).elapsedTime(end)/length(esets.(eset_name).expt(k).elapsedTime));
            disp(['Long tracks / tracks: ', num2str(nnz((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints))), ' / ', num2str(length(t))]);  % nnz (number of nonzero elements)
            t = t((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints));  % select tracks longer than requirement
            ymax = 0;
            for i = 1 : length(t_stim_start)  % ith intensity of stimulation--------------------
                for j = 1 : length(t)
                    % 'led2Val_ton' fails some time if stimulation isn't strict square wave
                    turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_ton', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % turn start time in period, ton means period starts with light on
                    turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % time (s) not in period
                    turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i))); %only keep the reorientation whose start time falls into the i-th intensity of stimulation
                    %number of period of stimulation that falls into certain intensity
                    t_start = max([t_stim_start(i), esets.(eset_name).expt(k).elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
                    t_end = min([t_stim_end(i), esets.(eset_name).expt(k).elapsedTime(t(j).endFrame-2)]);  % temporal - 2
                    nperiod = (t_end - t_start) / tperiod;
                    index_subplot = j + (i-1)*length(t);
                    ax(index_subplot) = subplot(length(t_stim_start), length(t), index_subplot);
                    [N, e] = histcounts(turnStart, edges);  % make sure larvae can only turn one time within tbin
                    bar(xbar, N/nperiod, 1);  % the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
                    if max(N/nperiod) > ymax  % check to make this work!!!!!!!!!!!!
                        ymax = max(N/nperiod);
                    end
                    xticks(edges);    ylim([0, 1]);  % comment ylim first, change to ymax at the second run ----------------
                    title(['Track ', num2str(t(j).trackNum), ' (', num2str(nperiod), ', ', num2str(sum(N)), ')']);
                end
            end
            xlabel(ax(1), 'Time in period (s)'); ylabel(ax(1), 'Probability of Starting to Turn'); 
            sgtitle('(Number of Stimulation, Number of Turn)');
            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\p_turn_no_pause_all');
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
            %         turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'position', 'start');  % time (s) not in period
            %         turnStartTime = turnStartTime((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i)));
            %         turnStart = mod(turnStartTime, tperiod);  % in period
            %         turnStart =  t(j).getSubFieldDQ('reorientation', 'led2Val_toff', 'position', 'start');  % turn start time in period, ton means period starts with light on
                    turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');
                    turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_toff', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');
                    turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i)));
                    turnStart_total = [turnStart_total turnStart];
                    if (t(j).startFrame < t_stim_end(i)*frame_rate) && (t(j).endFrame > t_stim_start(i)*frame_rate)
                        nperiod = nperiod + (min(t(j).endFrame, t_stim_end(i)*frame_rate) - max(t(j).startFrame, t_stim_start(i)*frame_rate)) / frame_rate / tperiod;
                    end
                end
                turnrate = rate_from_time(turnStart_total, tperiod, stepsize, binsize) ./ double(nperiod) * 60;
                time_timestep = [0 : fix(tperiod/stepsize)] * stepsize;
                ax(i) = subplot(length(t_stim_start),1,i);
                plot(time_timestep, turnrate);
                xlabel('Reorientation Start Time in Period (s)'); ylabel('Reorientation Rate (per min)'); 
                title([num2str(length(turnStart_total)), ' turns, in ', num2str(nperiod), ' periods, ', num2str(i), '-th intensity of stimulation']);
            end
            sgtitle(['Step size = ', num2str(stepsize), ', bin size = ', num2str(binsize)]);
            pause;
            if download
                savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\rate_turn_period_no_pause');
                savefig(gcf, savename); 
            end
            close;
        end % end plot_pturn
        
        
    end  % end expt loop
end  % end base folder loop

        
