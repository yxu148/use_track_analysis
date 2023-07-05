basedir = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Or42a@NA\T_Bl_Sq_0to100uWcm2_20';  % ---------------------------
d = dir(fullfile(basedir, 'matfiles', '*.mat'));
% reload experiment from mat files, called experiment set (eset), belong to @ExperimentSet object
disp('Loading data...');
eset = ExperimentSet.fromMatFiles(fullfile(basedir, 'matfiles', {d(1).name})); % load the first .mat if there are more ------------------


% load .mat files containing track information into eset
disp('Loading tracks ...');
tracks_mat = dir(fullfile(basedir, 'matfiles', 'Or42a@NA_202306231456 - tracks', '*.mat'));  %------------------------------
num_of_tracks = size(tracks_mat);
for i = 1 : num_of_tracks(1)
    load(fullfile(tracks_mat(i).folder, tracks_mat(i).name));
    s(i) = track;  % a structure to hold track from each .mat
end
eset.expt(1).track = s; 

 
disp('Sorting tracks into run, reorientation, and headswing...');
eset.executeTrackFunction('segmentTrack')  % input segmentation algorithm 'segmentTrack' ???????????
mkdir(fullfile(basedir, 'results1'));

% eset.expt.globalQuantity.fieldname  % to check what field is

%  led1Val = eset.gatherField('led1Val');  % Red Light intensity of LED at each frame -------------
% plot(eset.expt.elapsedTime(1:15e3), led1Val(1:15e3)); xlabel('Time (s)'); ylabel('led1Val');
% savename = strcat(basedir,'\results1', '\led1Val_time');
% savefig(gcf,savename);

led2Val = eset.gatherField('led2Val');  % Blue Light intensity (PWM) of LED at each frame  -----------------------------------------
% eset.expt goes to @Experiment object; waveform: square.

% % create a new globalquantity which is square wave across the whole experiment time
% inds_led1Val = find(strcmpi('led1Val', {eset.expt.globalQuantity.fieldname}));
% inds_led2Val = find(strcmpi('led2Val', {eset.expt.globalQuantity.fieldname}));
% xdata = eset.expt(1).globalQuantity(inds_led1Val).xData(1: end-7200); 
% ydata = eset.expt(1).globalQuantity(inds_led1Val).yData(7201: end) + eset.expt(1).globalQuantity(inds_led2Val).yData(1:end-7200);
% eset.expt(1).addGlobalQuantity('eti', 'led12Val', xdata, ydata)
% led12Val = eset.gatherField('led12Val');

eset.expt(1).addTonToff('led2Val', 'square');  % create time on/off field based a global quantity fieldname 'led2Val'----------------------
ton = eset.gatherField('led2Val_ton');  % get all values of 'led2Val_on' for all track in ExperimentSet eset----------------------
% return a k-N array of values, where N is total number of points, k is the dimension of the values of fieldname 'led2Val_off'
toff= eset.gatherField('led2Val_toff');  %-------------------------------------------
% correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
figure; t_end = 24e3-6;
% plot(eset.expt(1).elapsedTime(1 : t_end)/60, led1Val(1:t_end), 'r'); hold on;
plot(eset.expt(1).elapsedTime(1 : t_end)/60, led2Val(1:t_end), 'b'); hold on;
plot(eset.expt(1).elapsedTime(1:t_end)/60, toff(1:t_end), 'k'); hold off;
xlabel('Time (min)'); legend('led2Val', 'toff');
savename = strcat(basedir,'\results1', '\led2Val_toff');
savefig(gcf,savename);
% plot(eset.expt.elapsedTime(1:5e3),toff(1:5000), 'b'); xlabel('Time (s)'); hold off;  


tperiod = 20;  % depend on the name of .mat file loaded, '_18_', or from the plot led2Val-Time -----------------------------------------
disp(['Time for one frame is ', num2str(eset.expt.elapsedTime(2)), ' s']);
Ntracks = size(eset.expt(1).track);  % 1-by-Number_of_Tracks ( number of maggots)
figure;
histogram(round(eset.expt.elapsedTime([eset.expt.track.npts])/tperiod), [0:10:60]);  % round 60.001 to 60 ----------------------
xlabel('Number of Periods'); ylabel('Number of Tracks'); title(['Histogram of The Length of All ', num2str(length(eset.expt(1).track)), ' Tracks']);
savename = strcat(basedir,'\results1', '\track_length');
savefig(gcf,savename);


nperiods = 54;  % select tracks that have [nperiods, Nperiods] length
Nperiods = 91;  %expt time is 20 min, i.e. 20s periods at most for a 60 cycles, use 61 to include 60.001------------------------------------
disp(['After filtering out tracks within [', num2str(nperiods), ', ', num2str(Nperiods), '] periods']);
t = eset.expt.track;
minNpoints = nperiods * tperiod / (eset.expt.elapsedTime(end)/length(eset.expt.elapsedTime));
maxNpoints = Nperiods * tperiod / (eset.expt.elapsedTime(end)/length(eset.expt.elapsedTime));
disp(['There are still ', num2str(nnz((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints))), ' tracks left']);  % nnz (number of nonzero elements)
t = t((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints));  % select tracks longer than requirement
disp('The filtered tracks are stored in t');

