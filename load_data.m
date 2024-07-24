basedir = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3#T_Bl_Sq_2to7P_15_1_3';  % ---------------------------
d = dir(fullfile(basedir, 'matfiles', '*.mat'));
% reload experiment from mat files, called experiment set (eset), belong to @ExperimentSet object
disp('Loading data...');
x = [19];  % load the x-th set of data to analyze, x is a list [1], or [1, 2, 5], or delete (x) below for all--------------------
eset = ExperimentSet.fromMatFiles(fullfile(basedir, 'matfiles', {d(x).name}));  % d(x) or d


pause('on');  % 'on'---ask the user to press any key to save the figure, and continue; 'off'--directly save without asking
download = true;  % if true, plot and save basic figures including led1Val, Track length, Num of maggots; otherwise not plotting, but all valuables are prepared
tperiod = 15;  % depend on the name of .mat file loaded, '_18_', or from the plot led2Val-Time-----------------
% paramerters when generating BIN files of stimulation --------------------
t_stim_start = [0, 600, 1200];  % start time (s) of each intensity of stimulation
t_stim_end = [600, 1200, 1800];
stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start
frame_rate = 20;  % number of frames per second
% To save longer tracks in t
nperiods = 79;  % select tracks that have [nperiods, Nperiods] length, minimum 23 min
Nperiods = 121;  %expt time is 20 min, i.e. 20s periods at most for a 60 cycles, use 61 to include 60.001



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

% eset.expt.globalQuantity.fieldname to check what field is, the field doesn't depend on which track
field_name = {eset.expt.globalQuantity.fieldname};
GQled1Val = eset.expt.globalQuantity(strcmp({eset.expt.globalQuantity.fieldname}, {'led1Val'}));
GQled2Val = eset.expt.globalQuantity(strcmp({eset.expt.globalQuantity.fieldname}, {'led2Val'}));

% create a new globalquantity which is square wave across the whole experiment time
xdata = GQled2Val.xData;  % xdata of all Global Quantity ledVals should be the same
ydata = GQled1Val.yData + GQled2Val.yData;  % LED intensity in PWM
index = find(GQled1Val.yData==0, 1, 'last');  % find the index of the last zero element of led1Val, this index may be different from the initial setting because of hardware noise
ydata(1:index) = ydata(1:index) + 60;  % this is to make the square wave fluctruate around a center quantity, use index to make the combined square wave cleaner
eset.expt(1).addGlobalQuantity('eti', 'led12Val', xdata, ydata);  % create a man-made global field to add ton/toff
eset.expt(1).addTonToff('led12Val', 'square');  % create time on/off field based a global quantity fieldname 'led2Val'
GQton = eset.expt.globalQuantity(strcmp({eset.expt.globalQuantity.fieldname}, {'led12Val_ton'}));
GQtoff = eset.expt.globalQuantity(strcmp({eset.expt.globalQuantity.fieldname}, {'led12Val_toff'}));

% figure; plot(xdata, ydata);
% figure; plot(GQled2Val.xData, GQled2Val.yData);
% hold on;
% plot(GQtoff.xData, GQtoff.yData);
% hold off;
%% plot basic graphs
if download
    % correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
    figure;
%     plot(GQled2Val.xData/60, GQled1Val.yData, 'r'); hold on;
    plot(GQled2Val.xData/60, GQled2Val.yData, 'b'); hold on;
    plot(GQled2Val.xData/60, GQtoff.yData, 'k'); hold off;
%     plot(GQled2Val.xData/60, GQtoff.yData + 60, 'k'); hold off;
    xlabel('Time (min)'); legend('led2Val', 'toff');
    pause;
    savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\led12Val_toff');
    savefig(gcf,savename);
    close;
end


disp(['Time for one frame is ', num2str(eset.expt.elapsedTime(2)), ' s']);
Ntracks = size(eset.expt(1).track);  % 1-by-Number_of_Tracks ( number of maggots)
if download
    figure;
    histogram(round(eset.expt.elapsedTime([eset.expt.track.npts])/tperiod), 0:10:(Nperiods+1));  % round 60.001 to 60
    xlabel('Number of Periods'); ylabel('Number of Tracks'); title(['Histogram of The Length of All ', num2str(length(eset.expt(1).track)), ' Tracks']);
    pause;
    savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\track_length');
    savefig(gcf,savename);
    close;
end


% plot start time and end time of all tracks
if download
    figure;
    for j = 1: length(eset.expt.track)
            plot(eset.expt.track(j).startFrame + 1, j, 'bo'); hold on;
            plot(eset.expt.track(j).endFrame, j, 'rx'); hold on;
    end
    xlabel('Frame number'); ylabel('Index of tracks'); hold off; pause;
    savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\track_start_end');
    savefig(gcf, savename); close;
end


% plot the number of recognized maggots verses frame. Maggots in collision,
% out of ROI, or discarded by many different tests in processBIN won't be recognized.
ntracks_frame = zeros(1, length(eset.expt.elapsedTime));
start_frame_tracks = [eset.expt.track.startFrame] + 1;  % min is 1
end_frame_tracks = [eset.expt.track.endFrame] + 1;
for j = 1:length(eset.expt.track)
    ntracks_frame(start_frame_tracks(j) : end) = ntracks_frame(start_frame_tracks(j) : end) + 1;
    ntracks_frame(end_frame_tracks(j) : end) = ntracks_frame(end_frame_tracks(j) : end) - 1;
end
if download
    figure; plot(ntracks_frame);
    xlabel('Frame number'); ylabel('Number of recognized maggots'); pause;
    savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\num_maggots_recognized');
    savefig(gcf, savename); close;
end


t = eset.expt.track;
minNpoints = nperiods * tperiod / (eset.expt.elapsedTime(end)/length(eset.expt.elapsedTime));
t = t([t.npts] >= minNpoints);  % select tracks longer than requirement
disp([num2str(length(t)), ' long tracks out of ', num2str(length(eset.expt(1).track)), ' tracks, from ', num2str(max(ntracks_frame)), ' moving maggots']);


%% copy this file to the results folder
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)]);  % folder location to save file
copyfile('load_data.m', savename)  % copy the file to the location