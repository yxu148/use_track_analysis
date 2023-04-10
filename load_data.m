basedir = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\test_extracted_small\NA@NA\T_Bl_Sq_0to10_0'  % ---------------------------
d = dir(fullfile(basedir, 'matfiles', '*.mat'));
% reload experiment from mat files, called experiment set (eset), belong to @ExperimentSet object
disp('Loading data...');
eset = ExperimentSet.fromMatFiles(fullfile(basedir, 'matfiles', {d(1).name})); % load the first .mat if there are more


% load .mat files containing track information into eset
disp('Loading tracks ...');
tracks_mat = dir(fullfile(basedir, 'matfiles', 'NA@NA_202302241751 - tracks', '*.mat'));  %------------------------------
num_of_tracks = size(tracks_mat);
for i = 1 : num_of_tracks(1)
    load(fullfile(tracks_mat(i).folder, tracks_mat(i).name));
    s(i) = track;  % a structure to hold track from each .mat
end
eset.expt.track = s; 

 
disp('Sorting tracks into run, reorientation, and headswing...');
eset.executeTrackFunction('segmentTrack')  % input segmentation algorithm 'segmentTrack' ???????????
mkdir(fullfile(basedir, 'results'));

% eset.expt.globalQuantity.fieldname  % to check what field is

%  led1Val = eset.gatherField('led1Val');  % Red Light intensity of LED at each frame -------------
% plot(eset.expt.elapsedTime(1:15e3), led1Val(1:15e3)); xlabel('Time (s)'); ylabel('led1Val');
% savename = strcat(basedir,'\results', '\led1Val_time');
% savefig(gcf,savename);

led2Val = eset.gatherField('led2Val');  % Blue Light intensity of LED at each frame  -----------------------------------------
plot(eset.expt.elapsedTime(1:5e3), led2Val(1:5e3)); xlabel('Time (s)'); ylabel('led2Val');
savename = strcat(basedir,'\results', '\led2Val');
savefig(gcf,savename);
% eset.expt goes to @Experiment object; waveform: square.
eset.expt.addTonToff('led2Val', 'square');  % create time on/off field based a global quantity fieldname 'led2Val'----------------------
ton = eset.gatherField('led2Val_ton');  % get all values of 'led2Val_on' for all track in ExperimentSet eset----------------------
% return a k-N array of values, where N is total number of points, k is the dimension of the values of fieldname 'led2Val_off'
toff= eset.gatherField('led2Val_toff');  %-------------------------------------------
% correspond frame number to time, seems like the stimulation happens at 10 s of the 20 s period
% plot(eset.expt.elapsedTime(1:5e3),ton(1:5000), 'r'); hold on; 
% plot(eset.expt.elapsedTime(1:5e3),toff(1:5000), 'b'); xlabel('Time (s)'); hold off;  


tperiod = 6;  % depend on the name of .mat file loaded, '_18_', or from the plot led2Val-Time -----------------------------------------
disp(['Time for one frame is ', num2str(eset.expt.elapsedTime(2)), ' s']);
Ntracks = size(eset.expt(1).track);  % 1-by-Number_of_Tracks ( number of maggots)
figure;
histogram(round(eset.expt.elapsedTime([eset.expt.track.npts])/tperiod), [0:10:200]);  % round 60.001 to 60
xlabel('Number of Periods'); ylabel('Number of Tracks'); title(['Histogram of The Length of All ', num2str(length(eset.expt(1).track)), ' Tracks']);
savename = strcat(basedir,'\results', '\track_length');
savefig(gcf,savename);


nperiods = 100;  % select tracks that have [nperiods, Nperiods] length
Nperiods = 201;  %expt time is 20 min, i.e. 20s periods at most for a 60 cycles, use 61 to include 60.001------------------------------------
disp(['After filtering out tracks within [', num2str(nperiods), ', ', num2str(Nperiods), '] periods']);
t = eset.expt.track;
minNpoints = nperiods * tperiod / (eset.expt.elapsedTime(end)/length(eset.expt.elapsedTime));
maxNpoints = Nperiods * tperiod / (eset.expt.elapsedTime(end)/length(eset.expt.elapsedTime));
disp(['There are still ', num2str(nnz((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints))), ' tracks left']);  % nnz (number of nonzero elements)
t = t((maxNpoints >= [t.npts]) & ([t.npts] >= minNpoints));  % select tracks longer than requirement
disp('The filtered tracks are stored in t');

