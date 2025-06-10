basedir = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\demo\test_extracted\CS@NA\T_Bl_Sq_0to32uWcm2_20';  % ---------------------------
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

% to plot the light intensity
led2Val = eset.gatherField('led2Val');  % Blue Light intensity (PWM) of LED at each frame
eset.expt(1).addTonToff('led2Val', 'square');  % create time on/off field based a global quantity fieldname 'led2Val'
ton = eset.gatherField('led2Val_ton');  % get all values of 'led2Val_on' for all track in ExperimentSet eset
% return a k-N array of values, where N is total number of points, k is the dimension of the values of fieldname 'led2Val_off'
toff= eset.gatherField('led2Val_toff'); 
figure; t_end = 24e3-10;
plot(eset.expt(1).elapsedTime(1 : t_end)/60, led2Val(1:t_end), 'b'); hold on;
plot(eset.expt(1).elapsedTime(1:t_end)/60, toff(1:t_end), 'k'); 
xlabel('Time (min)'); legend('led2Val', 'toff'); hold off;
pause;
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\led2Val_toff');
savefig(gcf,savename);
close;

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
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\video1');
videoObject = VideoWriter(savename);
open(videoObject);
figure; eset.expt.track(1).playMovie('frameRate', 100, 'startTime', 500, 'stopTime', 510, 'vidObj', videoObject)
close(videoObject);


track_path = [1, 8];  % index of track to plot, could be [1], or [1, 3, 8]
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


figure;
tbin = 3;  edges = [0:tbin: tperiod];
%edges = [0,4,7,10,13,16,20];
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
for j = 1 :  10
    turnStart =  t(j).getSubFieldDQ('reorientation', 'led2Val_ton', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % turn start time in period, toff means period starts with light off
    %number of cycle for all points in the track, start with light on
    t_start = eset.expt.elapsedTime(t(j).startFrame + 1);  % time (s) of start for track j under i-th  intensity of stimulation 
    t_end = eset.expt.elapsedTime(t(j).endFrame);
    nperiod = (t_end - t_start) / tperiod;
    ax(j) = subplot(2,5,j);
    [N, e] = histcounts(turnStart, edges);  % make sure larvae can only turn one time within tbin
    bar(xbar, N/nperiod, 1);  % the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
    xticks(edges);  ylim([0, 1]);  % comment ylim first, change to ymax at the second run --------------
    title(['Track ', num2str(t(j).trackNum), ' (', num2str(nperiod), ', ', num2str(sum(N)), ')']);
end
xlabel(ax(1), 'Time in period (s)'); ylabel(ax(1), 'Probability of Start Turning'); 
sgtitle('(Number of Stimulation, Number of Turn)');
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\p_turn_ton');
savefig(gcf, savename); 



