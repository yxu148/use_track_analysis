ntracks = 11;  % number of tracks to plot ---------------------------
n = 4;  % --------------------------
m = fix(ntracks/n) + 1;  % m-by-n subplot ---------------
tbin = 1;  % bin size of [0, tperiod], for histogram
stepsize = 0.2; binsize = 2;  % for rate

% paramerters when generating BIN files of stimulation --------------------
t_stim_start = [0, 600, 1200];  % start time (s) of each intensity of stimulation
t_stim_end = [600, 1200, 1800];
frame_rate = 20;  % number of frames per second

% head swing rate in period
figure;
for j = 1 : ntracks
    % use led1Val if red lights are used, use off if start with off in one period--------------------
    nperiod = max(t(j).getDerivedQuantity('led2Val_cyclenum_on')) - min(t(j).getDerivedQuantity('led2Val_cyclenum_on')) + 1;
    turnStart =  t(j).getSubFieldDQ('headSwing', 'led2Val_toff', 'position', 'start');
    HSrate = rate_from_time(turnStart, tperiod, stepsize, binsize) / nperiod * 60;  % per min if t is in second
    time_timestep = [0 : fix(tperiod/stepsize)] * stepsize;
    ax(j) = subplot(m,n,j);
    plot(time_timestep, HSrate); ax(j).XAxis.TickValues = [0, tperiod/2, tperiod]; xline(tperiod/2, '--'); 
    ylim([0, 25]);
    title(['Track ', num2str(t(j).trackNum)]);
end
xlabel(ax(1),'Time in Period (s)'); ylabel(ax(1), 'Head Swing Rate (per min)');
sgtitle(['Step size is ', num2str(stepsize), ', bin size is ', num2str(binsize)]);
savename = strcat(basedir,'\results', '\rate_headSwing');
savefig(gcf, savename); 


% histogram of start head swing time in period in each track
figure;
for j=1:ntracks
    %the track divided by the period)
    % eset.expt.track(j) calls for the j-th track (maggot)
    HSstarttime =  t(j).getSubFieldDQ('headSwing', 'eti', 'position', 'start');  % get the beginning of reorientation
    HSstart = toff(round(HSstarttime/eset.expt.elapsedTime(2)));
    
    ax(j) = subplot(m,n,j);
    Turns = histogram(HSstart, 0:tbin:tperiod); title(['Track ', num2str(t(j).trackNum)]);  % change the size of bin according to period
    %ylim([0, 12]);  % unitify the y axis range to compare easily
end
xlabel(ax(1),'Head Swing Start Time in Period (s)'); ylabel(ax(1), 'Count');  % label the first plot
savename = strcat(basedir,'\results', '\num_headSwing');
savefig(gcf, savename);  % save figNumTurn as num_turn.fig in folder \results
% saveas(gcf, savename, 'png');  % run after insert texts


% histogram of start time in period of rejected head swing
figure;
for j = 1 : ntracks
    % 
    HSstarttime =  t(j).getSubFieldDQ('headSwing', 'eti', 'position', 'start');
    HSrejectedStartTime = HSstarttime([t(j).headSwing.accepted] == 0);  % select the start time of rejected swing
    HSrejStart = toff(round(HSrejectedStartTime/eset.expt.elapsedTime(2)));  % don't use mod, because tperiod isn't exact
    ax(j) = subplot(m,n,j);
    histogram(HSrejStart, 0:tbin:tperiod); title(['Track ', num2str(t(j).trackNum)]);
end
xlabel(ax(1),'Rejected Head Swing Start Time in Period (s)'); ylabel(ax(1), 'Count');
savename = strcat(basedir,'\results', '\num_headSwingRej');
savefig(gcf, savename); 


% reorientation rate in period
figure;
for j = 1 : ntracks
    % return the numbe-th of ton
    nperiod = max(t(j).getDerivedQuantity('led2Val_cyclenum_on')) - min(t(j).getDerivedQuantity('led2Val_cyclenum_on'));
    turnStart =  t(j).getSubFieldDQ('reorientation', 'led2Val_ton', 'position', 'start');  % time in period
    % turnStart = turnStart([t(j).reorientation.startInd] < t_work*frame_rate);  % only keey the reorientation whose start frame number is less than some value
    turnrate = rate_from_time(turnStart, tperiod, stepsize, binsize) / nperiod * 60;  % per min if t is in second
    Nturn = nnz(turnStart>=0 & turnStart<=tperiod) ;
    time_timestep = [0 : fix(tperiod/stepsize)] * stepsize;
    ax(j) = subplot(m,n,j);
    plot(time_timestep, turnrate); ax(j).XAxis.TickValues = [0, tperiod/4, tperiod/2, 3*tperiod/4, tperiod]; xline(tperiod/4, '--');
    % ylim([0, 20]);
    title(['Track ', num2str(t(j).trackNum), ' (', num2str(nperiod), ', ', num2str(Nturn), ')']);
end
xlabel(ax(1),'Time in Period (s)'); ylabel(ax(1), 'Reorientation Rate (per min)');
sgtitle(['Step size is ', num2str(stepsize), ', bin size is ', num2str(binsize), ', (Number of Stimulation, Number of Turn)']);
savename = strcat(basedir,'\results', '\rate_turn_ton');
savefig(gcf, savename); 
% saveas(gcf, savename, 'png');


% Probability of start to turn in one period for each track
figure;
tbin = 3;  edges = [0:tbin: tperiod];
%edges = [0,4,7,10,13,16,20];
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
ymax=0;
for j = 1 : ntracks
    turnStart =  t(j).getSubFieldDQ('reorientation', 'led2Val_ton', 'position', 'start');  % turn start time in period, toff means period starts with light off
    %number of cycle for all points in the track, start with light on
    nperiod = max(t(j).getDerivedQuantity('led2Val_cyclenum_on')) - min(t(j).getDerivedQuantity('led2Val_cyclenum_on'));
    ax(j) = subplot(m,n,j);
    [N, e] = histcounts(turnStart, edges);  % make sure larvae can only turn one time within tbin
    bar(xbar, N/nperiod, 1);  % the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
    if N/nperiod > ymax
        ymax = N/nperiod;
    end
    xticks(edges);  %ylim([0, 0.4]);  % comment ylim first, change to ymax at the second run --------------
    title(['Track ', num2str(t(j).trackNum), ' (', num2str(nperiod), ', ', num2str(sum(N)), ')']);
end
xlabel(ax(1), 'Time in period (s)'); ylabel(ax(1), 'Probability of Start Turning'); 
sgtitle('(Number of Stimulation, Number of Turn)');
savename = strcat(basedir,'\results', '\p_turn_ton');
savefig(gcf, savename); 


% Probability of start to turn in one period for each track
% Variability
figure;
tbin = 3;  edges = [0:tbin: tperiod];
%edges = [0,4,7,10,13,16,20];
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
ymax=0;
i = 3;  % ith intensity of stimulation
for j = 1 : ntracks
    % 'led2Val_ton' fails some time if stimulation isn't strict square wave
    turnStart =  t(j).getSubFieldDQ('reorientation', 'led2Val_ton', 'position', 'start');  % turn start time in period, ton means period starts with light on
%     turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'position', 'start');  % time (s) not in period
% 	turnStart = mod(turnStartTime, tperiod);  % in period
    % only keep the reorientation whose start time falls into the i-th
    % intensity of stimulation
    turnStart = turnStart((t_stim_start(i) * frame_rate <= [t(j).reorientation.startInd]) & ([t(j).reorientation.startInd] < t_stim_end(i) * frame_rate));
    %number of period of stimulation that falls into certain intensity
    t_start = max([t_stim_start(i), eset.expt.elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
    t_end = min([t_stim_end(i), eset.expt.elapsedTime(t(j).endFrame)]);
    nperiod = (t_end - t_start) / tperiod;
    ax(j) = subplot(m,n,j);
    [N, e] = histcounts(turnStart, edges);  % make sure larvae can only turn one time within tbin
    bar(xbar, N/nperiod, 1);  % the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
    if N/nperiod > ymax
        ymax = N/nperiod;
    end
    xticks(edges);    ylim([0, 0.5]);  % comment ylim first, change to ymax at the second run ----------------
    title(['Track ', num2str(t(j).trackNum), ' (', num2str(nperiod), ', ', num2str(sum(N)), ')']);
end
xlabel(ax(1), 'Time in period (s)'); ylabel(ax(1), 'Probability of Starting to Turn'); 
sgtitle([num2str(i), '-th intensity of stimulation (Number of Stimulation, Number of Turn)']);
savename = strcat(basedir,'\results', '\p_turn_', num2str(i));
savefig(gcf, savename); 


% histogram of start turn time in period in each track
figure;
for j=1:ntracks
    %the track divided by the period)
    % eset.expt.track(j) calls for the j-th track (maggot)
    turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'position', 'start');  % get the beginning of reorientation not in period
    turnStart = toff(round(turnStartTime/eset.expt.elapsedTime(2)));  %the problem is the slope of time-frame number isn't eset.expt.elapsedtime(2)
    ax(j) = subplot(m,n,j);
    histogram(turnStart, 0:tbin:tperiod); title(['Track ', num2str(t(j).trackNum)]);  % change the size of bin according to period
    %ylim([0, 12]);  % unitify the y axis range to compare easily
end
xlabel(ax(1),'Reorientation Start Time in Period (s)'); ylabel(ax(1), 'Count');  % label the first plot
savename = strcat(basedir,'\results', '\num_turn');
savefig(gcf, savename);  % save figNumTurn as num_turn.fig in folder \results
% saveas(gcf, savename, 'png');  % run after insert texts

% Histogram of probability difference before and after stimulation
tbin = 3;  edges = [0:tbin: tperiod];
% edges = [0,3,6,10,14,17,20];
%edges = [0,4,8,12,16,20];
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
increase_p = NaN(1, ntracks);
abs_p = NaN(1, ntracks);
for j = 1 : ntracks
    turnStart =  t(j).getSubFieldDQ('reorientation', 'led2Val_ton', 'position', 'start');  % turn start time in period, toff means period starts with light off
    nperiod =max(t(j).getDerivedQuantity('led2Val_cyclenum_on')) - min(t(j).getDerivedQuantity('led2Val_cyclenum_on'));  %number of cycle for all points in the track, start with light on
    [N, e] = histcounts(turnStart, edges);  % make sure larvae can only turn one time within tbin
    increase_p(j) = (N(1) - (N(5) + N(6))/2) / nperiod;  % increasement of turn probability is the first bin after stimulation minus the average of [12s, 18s].
    abs_p(j) = N(1) / nperiod;
end
savename = strcat(basedir,'\results', '\data');  % save in file data.mat
save(savename, 'increase_p');  % if data.mat not exist, create one, otherwise overwrite all.
save(savename, 'abs_p', '-append');  % data.mat has to exist, append new variable abs_p, or overwrite abs_p.
figure;
histogram(increase_p, [-0.5: 0.1: 1]); xlabel('Probability Increasement After Stimulation'); ylabel('Number of Tracks');
title([num2str(ntracks), ' Tracks Longer Than ', num2str(nperiods),  ' Periods with Bin Edges ', num2str(edges)]);
savename = strcat(basedir,'\results', '\diff_p_turn_1');  % no . in name of file
savefig(gcf, savename);
figure;
histogram(abs_p, [0: 0.1: 1]); xlabel('Probability of Turn within the first 3 s of Stimulation'); ylabel('Number of Tracks');
title([num2str(ntracks), ' Tracks Longer Than ', num2str(nperiods),  ' Periods']);
savename = strcat(basedir,'\results', '\abs_p_turn_1');
savefig(gcf, savename); 


% not make much sense, because track selection is random.
% it will make sense if cutting tracks into begining at 0 s, and end at 60
% s in period
turnStart =  eset.expt.track.getSubFieldDQ('reorientation', 'led2Val_ton', 'position', 'start');
nperiod = sum(eset.expt.elapsedTime(eset.gatherField('npts'))) / tperiod;  % total number of period for all tracks, assume all period is uniform
figure;
histogram(turnStart, 0:1:tperiod);
xlabel('Reorientation Start Time in Period (s)'); ylabel('Count'); title('All tracks');
savename = strcat(basedir,'\results', '\num_turn_all');
savefig(gcf, savename); 


% turn rate of all tracks at a certain period of time
figure;
stepsize = 0.1; binsize = 0.5;
% only keep the reorientation whose start time falls into the i-th intensity of stimulation
for i = 1: length(t_stim_start)
    nperiod = 0;
    turnStart_total=[];
    for j = 1: ntracks
%         turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'position', 'start');  % time (s) not in period
%         turnStart = mod(turnStartTime, tperiod);  % in period
        turnStart =  t(j).getSubFieldDQ('reorientation', 'led2Val_ton', 'position', 'start');  % turn start time in period, ton means period starts with light on
        turnStart = turnStart((t_stim_start(i) * frame_rate <= [t(j).reorientation.startInd]) & ([t(j).reorientation.startInd] < t_stim_end(i) * frame_rate));
        turnStart_total = [turnStart_total turnStart];
        if (t(j).startFrame < t_stim_end(i)*frame_rate) && (t(j).endFrame > t_stim_start(i)*frame_rate)
            nperiod = nperiod + (min(t(j).endFrame, t_stim_end(i)*frame_rate) - max(t(j).startFrame, t_stim_start(i)*frame_rate)) / frame_rate / tperiod;
        end
    end
    display(nperiod)
    turnrate = rate_from_time(turnStart_total, tperiod, stepsize, binsize) ./ double(nperiod) * 60;
    time_timestep = [0 : fix(tperiod/stepsize)] * stepsize;
    ax(i) = subplot(length(t_stim_start),1,i);
    plot(time_timestep, turnrate);
    xlabel('Reorientation Start Time in Period (s)'); ylabel('Reorientation Rate (per min)'); 
    title([num2str(i), '-th intensity of stimulation']);
end
sgtitle(['Step size = ', num2str(stepsize), ', bin size = ', num2str(binsize)]);
savename = strcat(basedir,'\results', '\rate_turn_period');
savefig(gcf, savename); 

% Changing of average turn rate
t_period = 300;  % calculate one rate for every 5 min (300 s)
t_rate_start = [0:t_period:1500];  % in second
t_rate_end = [t_period:t_period:1800];  % calculate one rate every 5 min
turnrate_array = zeros(1, length(t_rate_start));  % array to hold average rate
for i = 1: length(t_rate_start)  % loop through different time intervals
    nperiod = 0;  % initialize number of 
    turnStart_total=[];
    for j = 1: ntracks
        turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'position', 'start');  % time (s) not in period
        turnStart = mod(turnStartTime, t_period);
        turnStart = turnStart((t_rate_start(i) * frame_rate <= [t(j).reorientation.startInd]) & ([t(j).reorientation.startInd] < t_rate_end(i) * frame_rate));
        turnStart_total = [turnStart_total turnStart];
        if (t(j).startFrame < t_rate_end(i)*frame_rate) && (t(j).endFrame > t_rate_start(i)*frame_rate)
            nperiod = nperiod + (min(t(j).endFrame, t_rate_end(i)*frame_rate) - max(t(j).startFrame, t_rate_start(i)*frame_rate)) / frame_rate / t_period;
        end
    end
    display(nperiod)
    turnrate = rate_from_time(turnStart_total, t_period, stepsize, binsize) ./ double(nperiod) * 60;
    turnrate_array(i) = mean(turnrate);
end
figure
plot((t_rate_start +  t_rate_end)/2/60, turnrate_array, 'o-');
xlabel('Time (min)'); ylabel('Reorientation Rate (per min)'); 
title(['Step size = ', num2str(stepsize), ', bin size = ', num2str(binsize)]);
savename = strcat(basedir,'\results', '\rate_turn_changing');
savefig(gcf, savename); 


% turn rate of all tracks
figure;
stepsize = 0.1; binsize = 0.8;
% turnStartTime =  eset.expt.track.getSubFieldDQ('reorientation', 'eti', 'position', 'start');  % time (s) not in period
% turnStart = mod(turnStartTime, tperiod);  % in period
turnStart =  t.getSubFieldDQ('reorientation', 'led2Val_ton', 'position', 'start');
nperiod = sum(eset.expt.elapsedTime(eset.gatherField('npts'))) / tperiod;  % total number of period for all tracks, assume all period is uniform
turnrate = rate_from_time(turnStart, tperiod, stepsize, binsize) / double(nperiod) * 60;
time_timestep = [0 : fix(tperiod/stepsize)] * stepsize;
plot(time_timestep, turnrate);
xlabel('Reorientation Start Time in Period (s)'); ylabel('Reorientation Rate (per min)'); 
title(['All Tracks, ', 'Step size = ', num2str(stepsize), ', bin size = ', num2str(binsize)]);
savename = strcat(basedir,'\results', '\rate_turn_all');
savefig(gcf, savename); 

% to save some variables into a file, so that data of multiple files can be
% plotted together when load the data
turnStart_4 = turnStart;  % change the saving name _1 for the next one
nperiod_4 = nperiod;
savename = strcat(basedir,'\results', '\data_turn_rate_all_tracks.mat');
save(savename, 'turnStart_3', 'nperiod_3')  % run this at the first time to create new saving file
save(savename, 'turnStart_4', 'nperiod_4', '-append')  % run this next time to add new data to the same file
load(savename)
turnStart = [turnStart_1 turnStart_2 turnStart_3 turnStart_4];
nperiod = nperiod_1 + nperiod_2 + nperiod_3 + nperiod_4;

isrun = eset.gatherField('isrun');  % for all tracks in the expt, 1 for run, 0 for not
[tx,fracinrun] = meanyvsx(ton, isrun, 0:0.5:tperiod);  % toff transfer frame index to time between 0 and 60s
figure;
plot (tx, fracinrun, 'r'); % xline(1, '--'); 
xlabel('Time (s)');
ylabel('Fraction in run for all tracks');
savename = strcat(basedir,'\results', '\frac_run');
savefig(gcf,savename);




% raw data of first 5 tracks
for j = 1:5
    turnStartTime = t(j).getSubFieldDQ('reorientation', 'eti', 'position', 'start');  % does it need to calibrate by toff?????????????
    indstart = t(j).startFrame + 1;  % avoid 0
    indend = t(j).endFrame + 1;
    subplot(5, 1, j);
    plot(eset.expt.elapsedTime(indstart:indend), led2Val(indstart:indend)); hold on;
    plot(turnStartTime,12, 'rx'); hold off;  %-----------------------------------------------
    xlabel('Time (s)'); ylim([5, 18]); title(['Track ', num2str(t(j).trackNum)]);  %--------------------------------
end
sgtitle('Raw Data of Reorientation Start Time and Stimulation');
savename = strcat(basedir,'\results', '\raw_turn');
savefig(gcf, savename); 
% saveas(gcf, savename, 'png');