basedir = 'G:\AS-Filer\PHY\mmihovil\Shared\Isabel\chemotaxis_extracted\wt@na\5dpa_P_food';  % ---------------------------
d = dir(fullfile(basedir, 'matfiles', '*.mat'));
% reload experiment from mat files, called experiment set (eset), belong to @ExperimentSet object
disp('Loading data...');
x = [1];  % load the x-th set of data to analyze, x is a list [1], or [1, 2, 5], or delete (x) below for all--------------------
eset = ExperimentSet.fromMatFiles(fullfile(basedir, 'matfiles', {d(x).name}));  % d(x) or d
pause('off');  % 'on'---ask the user to press any key to save the figure, and continue; 'off'--directly save without asking


% load .mat files containing track information into eset
disp('Loading tracks ...');
d_tracks = dir(fullfile(basedir, 'matfiles', '*tracks'));
for k = 1 : length(x)
    tracks_mat = dir(fullfile(basedir, 'matfiles', d_tracks(x(k)).name, '*.mat'));
    clear s;
    for i = 1 : length(tracks_mat)
        load(fullfile(tracks_mat(i).folder, tracks_mat(i).name));
        s(i) = track;  % a structure to hold track from each .mat
    end
    eset.expt(1).track = s;  % save the track to x(k) or 1
end
 
disp('Sorting tracks into run, reorientation, and headswing...');
eset.executeTrackFunction('segmentTrack')
v_mean = 60 * mean([eset.expt.track.getSubFieldDQ('run', 'speed', 'position', 'mean')]);
disp(['Mean run speed of this group is ', num2str(v_mean), ' cm/min'])
mkdir(fullfile(basedir, ['results', d(x).name(end-16:end-4)]));  % auto name after date of expt, e.g. results_202309141250

% eset.expt.globalQuantity.fieldname  % to check what field is

%  led1Val = eset.gatherField('led1Val');  % Red Light intensity of LED at each frame
% plot(eset.expt.elapsedTime(1:15e3), led1Val(1:15e3)); xlabel('Time (s)'); ylabel('led1Val');
% savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\led1Val_time');
% savefig(gcf,savename);

led2Val = eset.gatherField('led2Val');  % Blue Light intensity (PWM) of LED at each frame
% eset.expt goes to @Experiment object; waveform: square.

% % create a new globalquantity which is square wave across the whole experiment time
% inds_led1Val = find(strcmpi('led1Val', {eset.expt.globalQuantity.fieldname}));
% inds_led2Val = find(strcmpi('led2Val', {eset.expt.globalQuantity.fieldname}));
% xdata = eset.expt(1).globalQuantity(inds_led1Val).xData; 
% ydata = eset.expt(1).globalQuantity(inds_led1Val).yData + eset.expt(1).globalQuantity(inds_led2Val).yData;
% ydata(12e3:24e3) = ydata(12e3:24e3) + 100;  % this is to make the square wave fluctruate around a center quantity
% eset.expt(1).addGlobalQuantity('eti', 'led12Val', xdata, ydata)
% led12Val = eset.gatherField('led12Val');

% eset.expt(1).addTonToff('led1Val', 'square');  % create time on/off field based a global quantity fieldname 'led2Val'
% ton = eset.gatherField('led1Val_ton');  % get all values of 'led2Val_on' for all track in ExperimentSet eset
% % return a k-N array of values, where N is total number of points, k is the dimension of the values of fieldname 'led2Val_off'
% toff= eset.gatherField('led1Val_toff'); 
% correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
% figure; t_end = 24e3-10;
% plot(eset.expt(1).elapsedTime(1 : t_end)/60, led1Val(1:t_end), 'r'); hold on;
% plot(eset.expt(1).elapsedTime(1 : t_end)/60, led2Val(1:t_end), 'b'); hold on;
% plot(eset.expt(1).elapsedTime(1:t_end)/60, toff(1:t_end), 'k'); hold on;
% plot(eset.expt(1).elapsedTime(1:t_end)/60, toff(1:t_end)+100, 'k'); hold off;
% xlabel('Time (min)'); legend('led1Val','led2Val', 'toff', 'toff');
% pause;
% savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\led12Val_toff');
% savefig(gcf,savename);
% close;
% plot(eset.expt.elapsedTime(1:5e3),toff(1:5000), 'b'); xlabel('Time (s)'); hold off;  

texpt = 1200; % time of experiment in seconds
tperiod = 20;  % depend on the name of .mat file loaded, '_18_', or from the plot led2Val-Time
disp(['Time for one frame is ', num2str(eset.expt.elapsedTime(2)), ' s']);
Ntracks = size(eset.expt(1).track);  % 1-by-Number_of_Tracks ( number of maggots)
figure;
histogram(eset.expt.elapsedTime([eset.expt.track.npts]), [0:texpt/10:texpt]);  % round 60.001 to 60
xlabel('Temporal length of the track (second)'); ylabel('Number of Tracks'); title(['Histogram of The Temporal Length of All ', num2str(length(eset.expt(1).track)), ' Tracks']);
pause;
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\track_length');
savefig(gcf,savename);
% close;


% nperiods = 69;  % select tracks that have [nperiods, Nperiods] length
% Nperiods = 91;  %expt time is 20 min, i.e. 20s periods at most for a 60 cycles, use 61 to include 60.001
% disp(['After filtering out tracks within [', num2str(nperiods), ', ', num2str(Nperiods), '] periods']);
% t = eset.expt.track;
% minNpoints = nperiods * tperiod / (eset.expt.elapsedTime(end)/length(eset.expt.elapsedTime));
% maxNpoints = Nperiods * tperiod / (eset.expt.elapsedTime(end)/length(eset.expt.elapsedTime));
% disp(['There are still ', num2str(nnz((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints))), ' tracks left']);  % nnz (number of nonzero elements)
% t = t((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints));  % select tracks longer than requirement
% disp('The filtered tracks are stored in t');


% plot start time and end time of all tracks
figure;
for j = 1: length(eset.expt.track)
        plot(eset.expt.track(j).startFrame + 1, j, 'bo'); hold on;
        plot(eset.expt.track(j).endFrame, j, 'rx'); hold on;
end
xlabel('Frame number'); ylabel('Index of tracks'); hold off;
pause;
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\track_start_end');
savefig(gcf, savename); 


% play video
% create the video outside of playMovie. Recover the playMovie.m.
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\video2');
videoObject = VideoWriter(savename);
open(videoObject);
figure; eset.expt.track(8).playMovie('frameRate', 100, 'startTime', 100, 'stopTime', 700, 'vidObj', videoObject)
close(videoObject);


% Plot one or multiple tracks to determine which to stitch
width = 2048;  % in pixel, in x-axis, geometry of Region of Interest (ROI) read from Image Recorder front panel
height = 2048;  % in pixel, in y-axis
len_pixel = realUnitsPerPixel(eset.expt.camcalinfo);  % how many cm per pixel
color_pad = ['r', 'g', 'b', 'k', 'c', 'm', 'y'];
for j = [1, 8, 9, 14, 15]
    track_path = [1, 8, 9, 14, 15];  % index of track to plot, could be [1], or [1, 3, 8], [j] -----------------------------
    figure; 
    eset.expt.track.plotPath('sloc', 'color', 0.8*[1, 1, 1], 'LineWidth', 2); hold on;  % the larger the whiter
    for i = 1 : length(track_path)  %  The index of track to plotPath, should be shorter than color_pad
        eset.expt.track(track_path(i)).plotPath('sloc', color_pad(i), 'LineWidth', 2);   hold on;  % 
        xy_s = eset.expt.track(track_path(i)).getDerivedQuantity('sloc');
        plot(xy_s(1, 1), xy_s(2, 1), append('o', color_pad(i)));  % o marks start
        plot(xy_s(1, end), xy_s(2, end), append('x', color_pad(i)));  % x marks end
    end
    rectangle('Position', [0, 0, width * len_pixel, height * len_pixel]); 
    axis equal;  % use the same length for data unit
    title(['Path of tracks ', num2str(track_path), ', with color ', color_pad(1: length(track_path)) ]); xlabel('x (cm)'); ylabel('y (cm)'); hold off;
    savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], ['\path_of_track', num2str(track_path(1))]);
    savefig(gcf, savename); 
    % saveas(gcf, savename, 'png');
end

track_path = [1, 8, 9, 14, 15];  % index of track to plot, could be [1], or [1, 3, 8]
for j = track_path  % the one track to plot
    t_single = eset.expt.track(j).dq.eti;
    v_single = eset.expt.track(j).dq.speed * 60;
    figure;
    plot(t_single, v_single);
    xlabel('Time (s)'); ylabel('Speed of Single Track (cm/min)');
    title(['Track ', num2str(j)]);
    savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], ['\vsingle_t_track', num2str(j)]);
    savefig(gcf,savename);
end


