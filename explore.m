% fieldname
led2Val = eset.gatherField('led2Val');  % Blue Light intensity (PWM) of LED at each frame
led2ValDeriv = eset.gatherField('led2ValDeriv');  % Derivative of Blue Light intensity (PWM) of LED at each frame
led2ValDiff = eset.gatherField('led2ValDiff');  % Differiente of Blue Light intensity (PWM) of LED at each frame
led2Val_high = eset.gatherField('led2Val_high');  % 1 for high, 0 for low of Blue Light intensity (PWM) of LED at each frame
led2Val_ton = eset.gatherField('led2Val_ton');  % 1 for high, 0 for low of Blue Light intensity (PWM) of LED at each frame

% To read BIN file:
fileID = fopen('G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_1_3#T_Bl_Sq_2,5to6P_20_2_3\202310251024\Gr21a@Chrimson(3)_T_Re_Sq_318to532P_20_1_3#T_Bl_Sq_2,5to6P_20_2_3_202310251024 led1 values.bin');
A = fread(fileID);
fclose(fileID);
fileID2 = fopen('G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_1_3#T_Bl_Sq_2,5to6P_20_2_3\202310251024\Gr21a@Chrimson(3)_T_Re_Sq_318to532P_20_1_3#T_Bl_Sq_2,5to6P_20_2_3_202310251024 led2 values.bin');
B = fread(fileID2);
fclose(fileID2);
figure; plot(A, 'r'); hold on; plot(B, 'b'); hold off;
C = A + B;
C(12001:24e3) = C(12001:24e3) + 100;
figure; plot(C, 'k');


% coordinate property of maggots
eset.expt.track(11).validDQName;  % will return all valid derived quantity names
% the derived quantities are saved in eset.expt.track(11).dq after calling getDerivedQuantity for the first time
% the derived quantities have the full length, which means some points are made up if some information is lost.
xy_s = eset.expt.track(11).getDerivedQuantity('sloc');  % smoothed location, [x1, x2, ... ; y1, y2, ...] in cm, use a lowpass Gaussian filter to smooth
xy_i = eset.expt.track(11).getDerivedQuantity('iloc');  % interpolated location, [x1, x2, ...; y1, y2, ...] in cm, use linear interp function to insert missing points.
vh = eset.expt(1).track(1).getDerivedQuantity('vhead');  % velocity of the head, vector
vt = eset.expt(1).track(1).getDerivedQuantity('vtail');  % velocity of the tail
tx = eset.expt(1).track(1).getDerivedQuantity('eti');  % time
ih = eset.expt(1).track(1).getDerivedQuantity('ihead');  % x and y pos of the head, in cm
it = eset.expt(1).track(1).getDerivedQuantity('itail');  % x and y pos of the tail
iiscollede = eset.expt.track(2).getDerivedQuantity('iiscollision');  % interpolated iscollision
% tderiv = 0.9;  % defines derivative time
% vh = deriv(ih, tderiv/eset.expt(1).dr.interpTime); 
% vt = deriv(it, tderiv/eset.expt(1).dr.interpTime);  % new calculations of the head and tail velocities
tmdir = eset.expt(1).track(1).getDerivedQuantity('itmdir');  % tail mid direction, vector
mhdir = eset.expt(1).track(1).getDerivedQuantity('imhdir');  % head mid direction
% plot (tx, dot(vh,tmdir)); xlabel('Time'); ylabel('Velocity of Tail');
[xc,lags] = xcov(dot(vh,mhdir), dot(vt,tmdir), ceil(4/eset.expt(1).dr.interpTime),'coeff');  % cross-covariance and lags ???????????????
plot (lags*eset.expt(1).dr.interpTime, xc); xlabel('Lags'); ylabel('Cross-covariance of Velocities of Head and Tail'); title('Track 1');
savename = strcat(basedir,'\results', '\cov_vh_vt');
savefig(gcf,savename);  % get current figure
% mean speed of each run
figure; plot(60 * eset.expt.track(1).getSubFieldDQ('run', 'speed', 'position', 'mean'));
xlabel('The number-th of Run'); ylabel('Mean speed (cm/min)');


%%%%%%%%%%%%%%%           turn probability per larva            %%%%%%%%%%%%%%%%%%
figure;
reorientation_rate = eset.makeReorientationHistogram('led2Val_toff', 0:0.1:20);  % per minute ?????????????????????????
plot(reorientation_rate); xlabel('led2Val\_toff'); ylabel('Reorientation Rate (per min)');
savename = strcat(basedir,'\results', '\reorientation_rate_power');
savefig(gcf,savename);


% ?????????????????????????????
for j = 1:7  % 4 tracks in total
    tbefore(j) = [sum(t(j).run.getDerivedQuantity('led2Val_ton') > 17)]*t(j).dr.interpTime;
    tafter(j) = [sum(t(j).run.getDerivedQuantity('led2Val_ton') < 3 & t(j).run.getDerivedQuantity('led2Val_ton') > 0)]*t(j).dr.interpTime;
    nturnafter(j) = [sum(t(j).reorientation.getDerivedQuantity('led2Val_ton','position', 'start') < 3 & t(j).reorientation.getDerivedQuantity('led2Val_ton','position', 'start') > 0)];
    nturnbefore(j) = [sum(t(j).reorientation.getDerivedQuantity('led2Val_ton','position', 'start') > 17)];
end

rbefore = nturnbefore./tbefore;  % rate of turn
rafter = nturnafter./tafter;
xx = [zeros(size(rbefore)); ones(size(rafter))];
plot (xx, [rbefore;rafter], 'bo-'); xlim([-3 3]); ylabel('Rate of Turn'); xlabel('Before 0 and after 1');\


% plot video, is running the function in @MaggotTrack
% time in second, frameRate default to be 20 Hz,
figure; eset.expt.track(1).playMovie('frameRate', 50, 'startTime', 0, 'stopTime', 10)


% Check collision time for each track, this is collision without breaking
% track (seperate automatically by software)
for j = 1:length(eset.expt.track)  % track j
    frame_number_collision = find(eset.expt.track(j).iscollision);
    if frame_number_collision
        time_collision = [eset.expt.track(j).pt(frame_number_collision).et];  % collision is a process
        % treat collisions seperated shorter than 0.5 s as the same collision
        if diff(time_collision) < 0.5  % true when all elements in diff(time_collision) are less than 0.5
            time_collision_start = time_collision(1); % select the beginning time of collision
        else  % treat collisions seperated longer than 0.5 s as different collisions
            time_collision_start = time_collision([false, diff(time_collision) >= 0.5]);
        end
        disp(['Track ', num2str(j), ' collides at ', num2str(time_collision_start), ' s']);
    else
        disp(['Track ', num2str(j), ' does not collide']);
    end
end


% plot start time and end time of each track, to figure out how to stitch
% them up. t = eset.expt.track
figure;
for j = 1: length(t)
        plot(eset.expt(1).elapsedTime(t(j).startFrame + 1), j, 'bo'); hold on;  % frame number starts with 0
        plot(eset.expt(1).elapsedTime(t(j).endFrame), j, 'rx'); hold on;
end
xlabel('Time (s)'); ylabel('Index of tracks'); hold off;
savename = strcat(basedir,'\results2', '\track_time');
savefig(gcf, savename); 

% plot the number of recognized maggots verses frame. Maggots in collision,
% out of ROI, or discarded by many different tests in processBIN won't be recognized.
ntracks_frame = zeros(1, length(eset.expt.elapsedTime));
start_frame_tracks = [eset.expt.track.startFrame] + 1;  % min is 1
end_frame_tracks = [eset.expt.track.endFrame] + 1;
for j = 1:length(eset.expt.track)
    ntracks_frame(start_frame_tracks(j) : end) = ntracks_frame(start_frame_tracks(j) : end) + 1;
    ntracks_frame(end_frame_tracks(j) : end) = ntracks_frame(end_frame_tracks(j) : end) - 1;
end
figure; plot(ntracks_frame);
xlabel('Frame number'); ylabel('Number of recognized maggots');
savename = strcat(basedir,'\results1', '\num_maggots_recognized');
savefig(gcf, savename); 


% size of Region of Interest (ROI) read from Image Recorder front panel
width = 2048;  % in pixel, in x-axis
height = 2048;  % in pixel, in y-axis
% scale info
len_pixel = realUnitsPerPixel(eset.expt.camcalinfo);  % how many cm per pixel
% correspond track to larva in video
eset.expt.track(j).pt(1).loc  % location (cm) of the j-th track at the first frame, 2*1 coloumn vector
% change the location in cm to location in pixels, so you can find it from MMF opened in ImageJ
camPtsFromRealPts(eset.expt.camcalinfo, eset.expt.track(j).pt(1).loc)  
figure;
eset.expt.track.plotPath(); hold on;
rectangle('Position', [0, 0, width * len_pixel, height * len_pixel]); hold off;
xlabel('x (cm)'); ylabel('y (cm)');


% geometry of Region of Interest (ROI) read from Image Recorder front panel
width = 1920;  % in pixel, in x-axis
height = 1920;  % in pixel, in y-axis
len_pixel = realUnitsPerPixel(eset.expt.camcalinfo);  % how many cm per pixel
color_pad = ['r', 'g', 'b', 'y', 'k', 'c', 'm'];
track_path = [3, 17, 18];  % index of track to plot, could be [1], or [1, 3, 8]-----------------------------
figure;
for i = 1 : length(track_path)  %  The index of track to plotPath, should be shorter than color_pad
    eset.expt.track(track_path(i)).plotPath('sloc', color_pad(i));   hold on;  % 
    plot(eset.expt.track(track_path(i)).pt(1).loc(1), eset.expt.track(track_path(i)).pt(1).loc(2), append('o', color_pad(i)));  % o marks start
    plot(eset.expt.track(track_path(i)).pt(end).loc(1), eset.expt.track(track_path(i)).pt(end).loc(2), append('x', color_pad(i)));  % x marks end
end
xlim([0, width * len_pixel]); ylim([0, height * len_pixel]); hold off;  % limit graph size by ROI
title('Path'); xlabel('x (cm)'); ylabel('y (cm)');
savename = strcat(basedir,'\results10', '\stitch_paths');
savefig(gcf, savename); 


figure;
eset.expt.track(2).plotPath('sloc', 'highlightinds', 'iscollision')





