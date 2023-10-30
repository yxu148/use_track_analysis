basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\stitch_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3',
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3',
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\stitch_extracted\Gr21a@Chrimson(3)_s\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3'
    };
x_cell = {
    [1, 2, 3, 4];  % load the x-th set of data from the first basedir to analyze
    [20, 21];
    [12]
    };

% load data, multiple eset, multiple expt, into esets
for folder_index = 1 : length(x_cell)  % loop for each basedir folder

    basedir = basedir_cell{folder_index};
    x = x_cell{folder_index};
    d = dir(fullfile(basedir, 'matfiles', '*.mat'));
    disp('Loading data...');
    eset_name = ['eset', num2str(folder_index)];  % dynamic structure name
    esets.(eset_name) = ExperimentSet.fromMatFiles(fullfile(basedir, 'matfiles', {d(x(k)).name}));  % create a structure esets to hold eset

    % load .mat files containing track information into eset
    disp('Loading tracks ...');
    d_tracks = dir(fullfile(basedir, 'matfiles', '*tracks'));
    for k = 1 : length(x)  % loop for each expt in the eset
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
download = true;  % true if plot figures
% add led12Val, and select long tracks
for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    basedir = basedir_cell{folder_index};
    eset_name = ['eset', num2str(folder_index)];
    x = x_cell{folder_index};
    d = dir(fullfile(basedir, 'matfiles', '*.mat'));
    
    for k = 1 : length(x)  % loop for each expt in the eset
        
        led1Val = esets.(eset_name).expt(k).gatherField('led1Val');
        led2Val = esets.(eset_name).expt(k).gatherField('led2Val');
        inds_led1Val = find(strcmpi('led1Val', {esets.(eset_name).expt(k).globalQuantity.fieldname}));
        inds_led2Val = find(strcmpi('led2Val', {esets.(eset_name).expt(k).globalQuantity.fieldname}));
        xdata = esets.(eset_name).expt(k).globalQuantity(inds_led1Val).xData; 
        ydata = esets.(eset_name).expt(k).globalQuantity(inds_led1Val).yData + esets.(eset_name).expt(k).globalQuantity(inds_led2Val).yData;
        ydata(1:12e3) = ydata(1:12e3) + 100;  % this is to make the square wave fluctruate around a center quantity
        esets.(eset_name).expt(k).addGlobalQuantity('eti', 'led12Val', xdata, ydata)
        led12Val = esets.(eset_name).expt(k).gatherField('led12Val');
        esets.(eset_name).expt(k).addTonToff('led12Val', 'square');
        
        ton = esets.(eset_name).expt(k).gatherField('led12Val_ton');  % get all values of 'led2Val_on' for all track in ExperimentSet eset
        % return a k-N array of values, where N is total number of points, k is the dimension of the values of fieldname 'led2Val_off'
        toff= esets.(eset_name).expt(k).gatherField('led12Val_toff'); 
        % correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
        if download
            figure; t_end = 36e3-10;
            plot(esets.(eset_name).expt(k).elapsedTime(1 : t_end)/60, led1Val(1:t_end), 'r'); hold on;
            plot(esets.(eset_name).expt(k).elapsedTime(1 : t_end)/60, led2Val(1:t_end), 'b'); hold on;
            plot(esets.(eset_name).expt(k).elapsedTime(1:t_end)/60, toff(1:t_end), 'k'); hold on;
            plot(esets.(eset_name).expt(k).elapsedTime(1:t_end)/60, toff(1:t_end)+100, 'k'); hold off;
            xlabel('Time (min)'); legend('led1Val','led2Val', 'toff', 'toff');
            pause;
            savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\led12Val_toff');
            savefig(gcf,savename);
            close;
        else
            continue;
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
        else
            continue;
        end

        
        t = esets.(eset_name).expt(k).track;
        minNpoints = nperiods * tperiod / (esets.(eset_name).expt(k).elapsedTime(end)/length(esets.(eset_name).expt(k).elapsedTime));
        maxNpoints = Nperiods * tperiod / (esets.(eset_name).expt(k).elapsedTime(end)/length(esets.(eset_name).expt(k).elapsedTime));
        disp(['Long tracks / tracks: ', num2str(nnz((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints))), ' / ', num2str(length(t))]);  % nnz (number of nonzero elements)
        esets.(eset_name).expt(k).track = t((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints));  % select tracks longer than requirement

    end
    
end


t_stim_start = [0, 600, 1200];  % start time (s) of each intensity of stimulation
t_stim_end = [600, 1200, 1800];
frame_rate = 20;  % number of frames per second
tbin = 3;  edges = [0:tbin: tperiod];
%edges = [0,4,7,10,13,16,20];
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    
    basedir = basedir_cell{folder_index};
    eset_name = ['eset', num2str(folder_index)];
    x = x_cell{folder_index};
    
    for k = 1 : length(x)  % loop for each expt in the eset
        t = esets.(eset_name).expt(k).track;  % save each specific track in t
        for i = 1 : length(t_stim_start)  % ith intensity of stimulation--------------------
            for j = 1 : length(t)
                % 'led2Val_ton' fails some time if stimulation isn't strict square wave
                turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_ton', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % turn start time in period, ton means period starts with light on
                turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % time (s) not in period
                turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i))); %only keep the reorientation whose start time falls into the i-th intensity of stimulation
                %number of period of stimulation that falls into certain intensity
                t_start = max([t_stim_start(i), eset.expt.elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
                t_end = min([t_stim_end(i), eset.expt.elapsedTime(t(j).endFrame-2)]);  % temporal - 2
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
        savename = strcat(basedir,['\results', d(x(k)).name(end-16:end-4)], '\p_turn_no_pause_all');
        savefig(gcf, savename); 
        close;
    end
end

        
