% Plot all tracks moving together with the raw video (.avi),
% adapted from Isabel

basedir = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data_new\variability_new_try_extracted\Or42a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3_#T_Bl_Sq_2to7P_15_1_3_';  % ---------------------------
d = dir(fullfile(basedir, 'matfiles', '*.mat'));
% reload experiment from mat files, called experiment set (eset), belong to @ExperimentSet object
disp('Loading data...');
x = [1];  % load the x-th set of data to analyze, x is a list [1], or [1, 2, 5], or delete (x) below for all--------------------
eset = ExperimentSet.fromMatFiles(fullfile(basedir, 'matfiles', {d(x).name}));  % d(x) or d
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
eset.executeTrackFunction('segmentTrack')
mkdir(fullfile(basedir, ['results', d(x).name(end-16:end-4)])); 
results_dir = strcat(basedir,['\results', d(x).name(end-16:end-4)]);


savename = fullfile(results_dir, 'collision3_tracks_trace_label.avi');
videoObject = VideoWriter(savename);
videoObject.FrameRate = 20;  % adjust as desired
open(videoObject);


%% 2) Figure / Axes Setup
fig = figure('Units','pixels','Position',[0,0,1980,1980],'Color','w');
ax  = axes('Parent',fig,'Color','k','XColor','k','YColor','k');
hold(ax,'on');
axis equal; 
set(ax,'YDir','normal');


% Suppose we want to plot a 25×25 cm region
xlim_plot = [6, 9];  % 'auto'; [0, 28];  [12, 15]; % in cm
ylim_plot = [2, 5];  % 'auto'; % [0, 28]; [10, 13]; 

xlim(ax,xlim_plot);
ylim(ax,ylim_plot);

xlabel(ax,'X (cm)','Color','k');
ylabel(ax,'Y (cm)','Color','k');


%% 3) Frame Range & Skip
frameStart  = 25500;       % e.g. 1
frameEnd    = 27000;    % e.g. 9000, 48000
skipFrames  =10;      % e.g. only show every 20th frame, 1 for not skipping
max_history = 1200;  % frames, maximum passing 200 frames of track will show

%% 4) Data / Track Info
eset.expt.openDataFile;
fid = eset.expt.fid;


all_tracks = eset.expt.track;
num_tracks = length(all_tracks); % length(all_tracks); 1;  


% Precompute each track's sloc + startFrame
sloc_cells   = cell(1,num_tracks);
start_frames = zeros(1,num_tracks);
for tIdx = 1:num_tracks
    sloc_cells{tIdx}  = all_tracks(tIdx).getDerivedQuantity('sloc');
    start_frames(tIdx) = all_tracks(tIdx).startFrame;
end


%% 5) Use distinguishable_colors for unique track colors
% Make sure you have 'distinguishable_colors.m' on your path.
% Optionally specify a background color or multiple backgrounds, e.g.:
% colorMap = distinguishable_colors(num_tracks, {'w','k'}); 
% for ensuring colors differ from white or black
colorMap = distinguishable_colors(num_tracks,[0 0 0.5; 0 0 0]); %so none are too close to black. can change to {'w','k'} to include white also
for i = 1:num_tracks
    hsvVal = rgb2hsv(colorMap(i,:));
    if hsvVal(2) < 0.1
        hsvVal(2) = 0.1;  % e.g. min saturation
    end
    if hsvVal(3) < 0.7
        hsvVal(3) = 0.7;  % e.g. min brightness
    end
    colorMap(i,:) = hsv2rgb(hsvVal);
end
%idx = randperm(num_tracks);
%colorMap = colorMap(idx, :);
%% 7) Animation Loop (draw images first, then path & labels)
for frameIdx = frameStart : skipFrames : frameEnd
    cla(ax); 
    hold(ax,'on');

    % -------- PHASE 1: DRAW ALL WORM IMAGES ----------
    for tIdx = 1:num_tracks
        track = all_tracks(tIdx);
        xy_s  = sloc_cells{tIdx};
        npts  = track.npts;  % size(xy_s,2); since plot images need track.pt which isn't smoothed
        sf    = start_frames(tIdx);


        if frameIdx < sf
            continue;  % track not started
        end
        localFrame = frameIdx - sf + 1;
        drawIdx = min(localFrame, npts);


        % Draw worm image for this track/frame
        track.pt(drawIdx).drawTrackImage(eset.expt.camcalinfo, 'fid', fid, ...
            'Axes', ax, 'pretty', true, ...
            'drawSpine', false, 'drawHeadArrow', false, ...
            'contourColor','y','drawContour',true,'contourWidth',1);
    end  % loop for each frame and draw images


    % -------- PHASE 2: DRAW ALL TRACK PATHS & LABELS ----------
    for tIdx = 1:num_tracks
        track = all_tracks(tIdx);
        xy_s  = sloc_cells{tIdx};
        npts  = size(xy_s,2);
        sf    = start_frames(tIdx);
        cTrack = colorMap(tIdx,:);  % color from distinguishable_colors


        if frameIdx < sf
            continue;  
        end
        localFrame = frameIdx - sf + 1;
        drawIdx = min(localFrame, npts);
        track_start_draw = max(1, drawIdx-max_history);


        % (B1) Plot the path in cTrack
        plot(ax, xy_s(1,track_start_draw:drawIdx), xy_s(2,track_start_draw:drawIdx), '-', ...
            'LineWidth',1.5,'Color',cTrack);


        % (B2) Yellow dot at track's start (once)
        if npts >= 1
            plot(ax, xy_s(1,1), xy_s(2,1),'yo','MarkerSize',4,'Color',cTrack);
        end


        % (B3) Mark the wormos current position with a white dot
        plot(ax, xy_s(1,drawIdx), xy_s(2,drawIdx), 'w.',...
            'MarkerSize',2);


        % (B4) Label with track #, in track color
        text(ax, xy_s(1,drawIdx), xy_s(2,drawIdx),...
            sprintf('%d', tIdx),...
            'Color', cTrack,'FontSize',18,'FontWeight','bold',...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
    end  % loop for each track and draw 


    % Keep axis at 25×25 cm 
    xlim(ax, xlim_plot);
    ylim(ax, ylim_plot);
    title(ax,sprintf('Frame %d/%d (skip=%d)', frameIdx, frameEnd, skipFrames),'Color','k');


    % Capture video frame
    frame = getframe(fig);
    writeVideo(videoObject, frame.cdata);


    pause(0.001);
end  % loop for each frame and draw tracks & labels


%% 8) Cleanup
close(videoObject);
close(fig);


disp('Multi-track animation with distinguishable_colors, star at start, path above images done!');