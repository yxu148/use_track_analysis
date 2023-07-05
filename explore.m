% fieldname
led2Val = eset.gatherField('led2Val');  % Blue Light intensity (PWM) of LED at each frame
led2ValDeriv = eset.gatherField('led2ValDeriv');  % Derivative of Blue Light intensity (PWM) of LED at each frame
led2ValDiff = eset.gatherField('led2ValDiff');  % Differiente of Blue Light intensity (PWM) of LED at each frame
led2Val_high = eset.gatherField('led2Val_high');  % 1 for high, 0 for low of Blue Light intensity (PWM) of LED at each frame
led2Val_ton = eset.gatherField('led2Val_ton');  % 1 for high, 0 for low of Blue Light intensity (PWM) of LED at each frame

% To read BIN file:
fileID = fopen('CS@NA_T_Bl_Sq_8to15PWM_20_202306021600 led2 values.bin');
A = fread(fileID);
fclose(fileID);


% scale info
len_pixel = realUnitsPerPixel(eset.expt.camcalinfo);  % how many cm per pixel

% correspond track to larva in video
eset.expt.track(1).pt(1).loc  % location (cm) of the first track at the first frame, 2*1 coloumn vector
camPtsFromRealPts(eset.expt.camcalinfo, eset.expt.track(1).pt(1).loc)  % change the location in cm to location in pixels, so you can find it from MMF opened in ImageJ
figure; eset.expt.track(1).plotPath  % plot the path of singe path
figure; eset.expt.track.plotPath  % plot the path of all paths

% coordinate property of maggots
vh = eset.expt(1).track(1).getDerivedQuantity('vhead');  % velocity of the head, vector
vt = eset.expt(1).track(1).getDerivedQuantity('vtail');  % velocity of the tail
tx = eset.expt(1).track(1).getDerivedQuantity('eti');  % time
ih = eset.expt(1).track(1).getDerivedQuantity('ihead');  % x and y pos of the head, in cm
it = eset.expt(1).track(1).getDerivedQuantity('itail');  % x and y pos of the tail
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
figure; eset.expt.track(2).playMovie('frameRate', 50, 'startTime', 1300, 'stopTime', 1700)


% Check collision situation
j = 6;  % track j
frame_number_collision = find(eset.expt.track(j).iscollision);
time_collision = [eset.expt.track(j).pt(frame_number_collision).et];


% plot start time and end time of each track, to figure out how to stitch
% them up. t = eset.expt.track
figure;
for j = 1: length(t)
        plot(eset.expt(1).elapsedTime(t(j).startFrame + 1), j, 'bo'); hold on;
        plot(eset.expt(1).elapsedTime(t(j).endFrame), j, 'rx'); hold on;
end
xlabel('Time (s)'); ylabel('Index of tracks'); hold off;
savename = strcat(basedir,'\results1', '\track_time');
savefig(gcf, savename); 





