ntracks =52;  % number of tracks to plot ---------------------------
n = 4;  % --------------------------row
m = fix(ntracks/n);  % m-by-n subplot ---------------
tbin = 1;  % bin size of [0, tperiod], for histogram
stepsize = 0.2; binsize = 2;  % for rate

%%  Head Swing
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
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\rate_headSwing');
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
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\num_headSwing');
savefig(gcf, savename);  % save figNumTurn as num_turn.fig in folder \results_202308121650
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
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\num_headSwingRej');
savefig(gcf, savename); 


%% Probability of Rurn

% Probability of start to turn in one period for each track
figure;
tbin = 3;  edges = [0:tbin: tperiod];
%edges = [0,4,7,10,13,16,20];
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
ymax=0;
for j = 1 : ntracks
    turnStart =  t(j).getSubFieldDQ('reorientation', 'led2Val_ton', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % turn start time in period, toff means period starts with light off
    %number of cycle for all points in the track, start with light on
    nperiod = max(t(j).getDerivedQuantity('led2Val_cyclenum_on')) - min(t(j).getDerivedQuantity('led2Val_cyclenum_on'));
    ax(j) = subplot(m,n,j);
    [N, e] = histcounts(turnStart, edges);  % make sure larvae can only turn one time within tbin
    bar(xbar, N/nperiod, 1);  % the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
    if N/nperiod > ymax
        ymax = N/nperiod;
    end
    xticks(edges);  ylim([0, 0.4]);  % comment ylim first, change to ymax at the second run --------------
    title(['Track ', num2str(t(j).trackNum), ' (', num2str(nperiod), ', ', num2str(sum(N)), ')']);
end
xlabel(ax(1), 'Time in period (s)'); ylabel(ax(1), 'Probability of Start Turning'); 
sgtitle('(Number of Stimulation, Number of Turn)');
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\p_turn_no_pause_ton');
savefig(gcf, savename); 


% Probability of start to turn in one period for each track
% Variability
figure;
tbin = 3;  edges = [0:tbin: tperiod];
%edges = [0,4,7,10,13,16,20];
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
ymax=0;
for j = 1 : length(t)
    % 'led2Val_ton' fails some time if stimulation isn't strict square wave
    turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_ton', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % turn start time in period, ton means period starts with light on
    turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % time (s) not in period
    turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i))); %only keep the reorientation whose start time falls into the i-th intensity of stimulation
% 	turnStart = mod(turnStartTime, tperiod);  % in period
    %number of period of stimulation that falls into certain intensity
    t_start = max([t_stim_start(i), eset.expt.elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
    t_end = min([t_stim_end(i), eset.expt.elapsedTime(t(j).endFrame)]);
    nperiod = (t_end - t_start) / tperiod;
    ax(j) = subplot(m, n, j);
    [N, e] = histcounts(turnStart, edges);  % make sure larvae can only turn one time within tbin
    bar(xbar, N/nperiod, 1);  % the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
    if N/nperiod > ymax
        ymax = N/nperiod;
    end
    xticks(edges);    %ylim([0, 0.85]);  % comment ylim first, change to ymax at the second run ----------------
    title(['Track ', num2str(t(j).trackNum), ' (', num2str(nperiod), ', ', num2str(sum(N)), ')']);
end
xlabel(ax(1), 'Time in period (s)'); ylabel(ax(1), 'Probability of Starting to Turn'); 
sgtitle([num2str(i), '-th intensity of stimulation (Number of Stimulation, Number of Turn)']);
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\p_turn_ton_no_pause_', num2str(i));
savefig(gcf, savename); 


% Probability of start to turn in one period for each track
% Variability, each track in each column, each period in each row
figure;
tbin = 3;  edges = 0:tbin: tperiod;
%edges = [0,4,7,10,13,16,20];
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
turnStart_mean = zeros(1, length(xbar));  % the average turnTime that falls to one bin
turnStart_str = zeros(1, length(xbar));  % standard error = standard deviation / sqrt(N)
load(strcat(basedir,['\results', d(x).name(end-16:end-4)], '\data.mat'), 'larvae');
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
%         index_subplot = j*length(t_stim_start) - (length(t_stim_start)-i);  % to plot by column
%         ax(index_subplot) = subplot(length(t), length(t_stim_start), index_subplot);
        index_subplot = j + (i-1)*length(t);
        ax(index_subplot) = subplot(length(t_stim_start), length(t), index_subplot);
        
%         set(ax(index_subplot), 'Units', 'inches');
%         pos = get(ax(index_subplot), 'Position');
%         pos(3:4) =  [2, 1];  % width and height
%         set(ax(index_subplot), 'Position', pos);

        [N, ~] = histcounts(turnStart, edges);  % make sure larvae can only turn one time within tbin
        for r = 1 : (length(edges) - 1)
            turnStart_bin = turnStart((turnStart > edges(r)) & (turnStart < edges(r+1)));
            turnStart_mean(r) = mean(turnStart_bin);
%             turnStart_str(r) = std(turnStart_bin, 1) / sqrt(N(r));  % standard error
            turnStart_str(r) = std(turnStart_bin, 1);
        end
        bar(xbar, N/nperiod, 1); hold on; % the value at [10, 13], describe the  possibility of turning within 2 seconds after stimulation
        errorbar(turnStart_mean, N/nperiod, sqrt(N)/nperiod / 2, 'b.'); % the length of the error bar is the error.
        errorbar(turnStart_mean, N/nperiod, turnStart_str/2, 'horizontal', 'k.'); hold off;
        xticks(edges); ylim([0, 1]); xline(6, '--');
        
        
%         % save info into a new structure larvae
        larva_index = ['larva', num2str(j)];
%         larvae.(larva_index).trackNum = t(j).trackNum;  % integer
%         larvae.(larva_index).expt_info = d(x).name(1:end-4);  % string
%         larvae.(larva_index).nperiod.(stim_color{i}) = nperiod;
%         larvae.(larva_index).pturn.(stim_color{i}) = N/nperiod;
%         larvae.(larva_index).c.(stim_color{i}) = sqrt(N)/nperiod;
%         larvae.(larva_index).turnStart.(stim_color{i}) = turnStart;
%         if N(1)/nperiod - mean(N(2:end)/nperiod) > 0.2  % if the first bin of pturn is much larger than the rest, call it response
%             larvae.(larva_index).response.(stim_color{i}) = '1';
%         elseif mean(N(1:2)/nperiod) - mean(N(3:end)/nperiod) > 0.2  % if the first 2 bins of pturn are much larger than the rest
%             larvae.(larva_index).response.(stim_color{i}) = '1';
%         else
%             larvae.(larva_index).response.(stim_color{i}) = '0';
%         end

        title(['(', num2str(t(j).trackNum), ', ', num2str(nperiod), ', ', num2str(sum(N)), ', ', larvae.(larva_index).response.(stim_color{i}), ')']);
    end
end
xlabel(ax(1), 'Time in period (s)'); ylabel(ax(1), 'Probability of Starting to Turn'); 
sgtitle('(Track number, Number of Stimulation, Number of Turn, Response)');
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\p_turn_no_pause_all_errorbar');
savefig(gcf, savename); 

% Histogram of probability difference before and after stimulation
tbin = 3;  edges = 0:tbin: tperiod;
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
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\data');  % save in file data.mat
save(savename, 'increase_p');  % if data.mat not exist, create one, otherwise overwrite all.
save(savename, 'abs_p', '-append');  % data.mat has to exist, append new variable abs_p, or overwrite abs_p.
figure;
histogram(increase_p, [-0.5: 0.1: 1]); xlabel('Probability Increasement After Stimulation'); ylabel('Number of Tracks');
title([num2str(ntracks), ' Tracks Longer Than ', num2str(nperiods),  ' Periods with Bin Edges ', num2str(edges)]);
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\diff_p_turn_1');  % no . in name of file
savefig(gcf, savename);
figure;
histogram(abs_p, [0: 0.1: 1]); xlabel('Probability of Turn within the first 3 s of Stimulation'); ylabel('Number of Tracks');
title([num2str(ntracks), ' Tracks Longer Than ', num2str(nperiods),  ' Periods']);
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\abs_p_turn_1');
savefig(gcf, savename); 

%% Turn Time
% histogram of start turn time in period in each track
figure;
for j=1:ntracks
    %the track divided by the period)
    % eset.expt.track(j) calls for the j-th track (maggot)
    turnStart =  t(j).getSubFieldDQ('reorientation', 'led2Val_toff', 'position', 'start');  % turn start time in period, ton means period starts with light off
%     turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'position', 'start');  % get the beginning of reorientation not in period
%     turnStart = toff(round(turnStartTime/eset.expt.elapsedTime(2)));  %the problem is the slope of time-frame number isn't eset.expt.elapsedtime(2)
    ax(j) = subplot(m,n,j);
    histogram(turnStart, 0:tbin:tperiod); title(['Track ', num2str(t(j).trackNum)]);  % change the size of bin according to period
    %ylim([0, 12]);  % unitify the y axis range to compare easily
end
xlabel(ax(1),'Reorientation Start Time in Period (s)'); ylabel(ax(1), 'Count');  % label the first plot
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\num_turn');
savefig(gcf, savename);  % save figNumTurn as num_turn.fig in folder \results_202308121650
% saveas(gcf, savename, 'png');  % run after insert texts


% not make much sense, because track selection is random.
% it will make sense if cutting tracks into begining at 0 s, and end at 60
% s in period
turnStart =  eset.expt.track.getSubFieldDQ('reorientation', 'led2Val_ton', 'position', 'start');
nperiod = sum(eset.expt.elapsedTime(eset.gatherField('npts'))) / tperiod;  % total number of period for all tracks, assume all period is uniform
figure;
histogram(turnStart, 0:1:tperiod);
xlabel('Reorientation Start Time in Period (s)'); ylabel('Count'); title('All tracks');
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\num_turn_all');
savefig(gcf, savename); 


%% Turn Rate

% turn rate in period for each track
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
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\rate_turn_ton');
savefig(gcf, savename); 
% saveas(gcf, savename, 'png');


% turn rate of all tracks at a certain period of time
figure; t = eset.expt.track;
stepsize = 0.1; binsize = 0.5;
% only keep the reorientation whose start time falls into the i-th intensity of stimulation
for i = 1: length(t_stim_start)
    nperiod = 0;  % initialize the number of periods for each stimulation condition
    turnStart_total=[];  % initialize the turn start time of all tracks for each stimulation condition
    for j = 1: length(t)
        % turn start time (second) of track j, excluding pause
        turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');
        % turn start time in period (second) of track j, excluding pause, start from intensity low. 
        turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_toff', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');
        % select the start time in period of turns that happen between t_stim_start(i) and t_stim_start(i)
        turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i)));
        turnStart_total = [turnStart_total turnStart];  % to contain more turnStart after each loop of j
        if (t(j).startFrame < t_stim_end(i)*frame_rate) && (t(j).endFrame > t_stim_start(i)*frame_rate)
            % add to nperiod after each loop of j
            nperiod = nperiod + (min(t(j).endFrame, t_stim_end(i)*frame_rate) - max(t(j).startFrame, t_stim_start(i)*frame_rate)) / frame_rate / tperiod;
        end
    end
    turnrate = rate_from_time(turnStart_total, tperiod, stepsize, binsize) ./ double(nperiod) * 60;
    time_timestep = [0 : fix(tperiod/stepsize)] * stepsize;
    ax(i) = subplot(length(t_stim_start),1,i);
    plot(time_timestep, turnrate); xline(9, 'k--');
    xlabel('Reorientation Start Time in Period (s)'); ylabel('Reorientation Rate (per min)'); 
    title([num2str(length(turnStart_total)), ' turns, in ', num2str(nperiod), ' periods, ', num2str(i), '-th intensity of stimulation']);
end
sgtitle(['Step size = ', num2str(stepsize), ', bin size = ', num2str(binsize)]);
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\rate_turn_period_no_pause');
savefig(gcf, savename); 



% Changing of average turn rate (including pause)
stepsize = 0.1; binsize = 0.5; frame_rate = 20;
t_period = 300;  % calculate one rate for every 5 min (300 s)
t_rate_start = [0:t_period:3300];  % in second
t_rate_end = [t_period:t_period:3600];  % calculate one rate every 5 min
turnrate_array = zeros(1, length(t_rate_start));  % array to hold average rate
nperiod_array = zeros(1, length(t_rate_start));  % arrray to hold the nperiod for each period
for i = 1: length(t_rate_start)  % loop through different time intervals
    nperiod = 0;  % initialize number of period
    turnStart_total=[];
    for j = 1: ntracks
        turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'position', 'start');  % time (s) not in period
        turnStartTime = turnStartTime((t_rate_start(i) <= turnStartTime) & (turnStartTime <= t_rate_end(i)));
        turnStart = mod(turnStartTime, tperiod);
        % turnStart = turnStart((t_rate_start(i) * frame_rate <= [t(j).reorientation.startInd]) & ([t(j).reorientation.startInd] < t_rate_end(i) * frame_rate));
        turnStart_total = [turnStart_total turnStart];
        if (t(j).startFrame < t_rate_end(i)*frame_rate) && (t(j).endFrame > t_rate_start(i)*frame_rate)
            nperiod = nperiod + (min(t(j).endFrame, t_rate_end(i)*frame_rate) - max(t(j).startFrame, t_rate_start(i)*frame_rate)) / frame_rate / tperiod;
        end
    end
    turnrate = rate_from_time(turnStart_total, tperiod, stepsize, binsize) ./ double(nperiod) * 60;
    turnrate_array(i) = mean(turnrate);
    nperiod_array(i) = nperiod;
end
figure
plot((t_rate_start +  t_rate_end)/2/60, turnrate_array, 'o-');
text((t_rate_start +  t_rate_end)/2/60 + 1, turnrate_array, string(nperiod_array));
xlabel('Time (min)'); ylabel('Reorientation Rate (per min)'); 
title(['Step size = ', num2str(stepsize), ', bin size = ', num2str(binsize)]);
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\rate_turn_changing_1h_includ_pause');
savefig(gcf, savename); 


% turn rate of all tracks
figure; eset.makeReorientationHistogram('led12Val_ton', 0:stepsize:tperiod);

% turn rate of all tracks
figure;  t = eset.expt.track;
stepsize = 0.1; binsize = 0.5;
turnStart = [];
for i = 1: length(t)  % when include indsExpression in function getSubFieldDQ, cannot get quantities for all tracks
    turnStart_single =  t(i).getSubFieldDQ('reorientation', 'led1Val_toff', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');
    turnStart = [turnStart turnStart_single];
end
nperiod = sum(eset.expt.elapsedTime(eset.gatherField('npts'))) / tperiod;  % total number of period for all tracks, assume all period is uniform
turnrate = rate_from_time(turnStart, tperiod, stepsize, binsize) / double(nperiod) * 60;
time_timestep = [0 : fix(tperiod/stepsize)] * stepsize;
plot(time_timestep, turnrate);
xlabel('Reorientation Start Time in Period (s)'); ylabel('Reorientation Rate (per min)'); 
title([num2str(fix(nperiod)), ' nperiods, ',num2str(length(turnStart)),  ' turns, Step size = ', num2str(stepsize), ', bin size = ', num2str(binsize)]);
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\rate_turn_all_toff_no_pause');
savefig(gcf, savename); 


% Turn rate in period for each long track
% Variability, each track in each column, each stimulation period in each row
% save info into a existing structure larvae
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\data.mat');
load(savename, 'larvae');
figure;
stepsize = 0.1; binsize = 0.5;
for i = 1 : length(t_stim_start)  % ith intensity of stimulation--------------------
    for j = 1 : length(t)  % j-th track
        turnStart =  t(j).getSubFieldDQ('reorientation', 'led12Val_ton', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % turn start time in period, ton means period starts with light on
        turnStartTime =  t(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % time (s) not in period
        turnStart = turnStart((t_stim_start(i) <= turnStartTime) & (turnStartTime < t_stim_end(i))); %only keep the reorientation whose start time falls into the i-th intensity of stimulation
        %number of period of stimulation that falls into certain intensity
        t_start = max([t_stim_start(i), eset.expt.elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th intensity of stimulation 
        t_end = min([t_stim_end(i), eset.expt.elapsedTime(t(j).endFrame-2)]);  % temporal - 2
        nperiod = (t_end - t_start) / tperiod;
%         index_subplot = j*length(t_stim_start) - (length(t_stim_start)-i);  % to plot by column
%         ax(index_subplot) = subplot(length(t), length(t_stim_start), index_subplot);
        index_subplot = j + (i-1)*length(t);  % to put each track in each column, and each stimulation period in each row
        ax(index_subplot) = subplot(length(t_stim_start), length(t), index_subplot);
        turnrate = rate_from_time(turnStart, tperiod, stepsize, binsize) ./ double(nperiod) * 60;
        time_timestep = (0 : fix(tperiod/stepsize)) * stepsize;
        plot(time_timestep, turnrate); xline(6, 'k--');
        title(['(', num2str(t(j).trackNum), ', ', num2str(nperiod), ', ', num2str(length(turnStart)), ')']);
        
        % add new contents to the structure larvae
        larva_index = ['larva', num2str(j)];
        larvae.(larva_index).turnrate.(stim_color{i}) = turnrate;
    end
end
linkaxes(ax, 'y');  % align subplots with y axis
xlabel(ax(1), 'ton (s)'); ylabel(ax(1), 'Turn rate (min^{-1})'); sgtitle(['Stepsize = ', num2str(stepsize), ', binsize = ', num2str(binsize), ' (track number, number of period, number of turn)']);
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\rate_turn_ton_individual');
savefig(gcf, savename); 


%% Speed
% speed of run track part verses the index of run
figure; v_mean_run = t(1).getSubFieldDQ('run', 'speed', 'position', 'mean') * 60 ; plot(v_mean_run);


% speed of single larvae verses experimental time
j = 1;  % track number
v_frame = eset.expt.track(j).dq.speed * 60;  % cm/min
t_frame = eset.expt.track(j).dq.eti;  % interpolated time (s) for each frame of track j, 1-by-(number of the track's frame) 
turnStartTime =  eset.expt.track(j).getSubFieldDQ('reorientation', 'eti', 'indsExpression', '[track.reorientation.numHS] >= 1', 'position', 'start');  % time (s) not in period
figure;
plot(t_frame, v_frame); hold on;
plot(turnStartTime, 0.3, 'ok'); hold off; legend('Speed', 'Turn start')
xlabel('Time (s)'); ylabel('Speed (cm/min)'); title(['Track ', num2str(j)])






% speed of single larvae verses time in period
j = 2;  % track number
v_frame = eset.expt.track(j).dq.speed * 60;  % cm/min
toff_frame = eset.expt.track(j).dq.led12Val_toff;  % interpolated time (s) for each frame of track j, 1-by-(number of the track's frame) 
stepsize = 0.1;  % to get one average speed for stepsize seconds
[x_toff,y_v, stderror] = meanyvsx(toff_frame, v_frame, 0:stepsize:tperiod);
uppercurve = y_v + 0.5*stderror;
lowercurve = y_v - 0.5*stderror;
x_tofill = [x_toff, fliplr(x_toff)];  % the x axis of the ploygon to fill
y_tofill = [lowercurve, fliplr(uppercurve)];
figure;
pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none'); hold on;  % no edges for the patch
plot(x_toff, y_v, 'Color', 0.2*[1 1 1]); 
xline(10, '--'); hold off;
xlabel('Time toff in period (s)'); ylabel(['Speed of track ', num2str(j), ' (cm/min)']);

% average speed of all tracks verses time in period at certain period of time/stimulation condition
v_all = eset.gatherField('speed') * 60;  % cm/min
t_all = eset.gatherField('eti');  % interpolated time (s) for each frame of each track, 1-by-(sum of all tracks' frame) 
toff_all = eset.gatherField('led12Val_toff');  % interpolated time (s) for each frame of track j, 1-by-(number of the track's frame) 


% Changing of mean speed of each run for one track
j=3;  % track index
figure;
plot(60 * eset.expt.track(j).getSubFieldDQ('run', 'speed', 'position', 'mean'));
xlabel('The number-th Run'); ylabel('Mean speed (cm/min)'); title(['Track ', num2str(j)]);


% speed of all tracks at a certain period of time of stimulation
figure; t = eset.expt.track;
% only keep the speed whose start time falls into the i-th intensity of stimulation
for i = 1: length(t_stim_start)
    nperiod = 0;  % initialize the number of periods for each stimulation condition
    toff_frame_total= [];  % initialize the turn start time of all tracks for each stimulation condition
    v_frame_total = [];
    for j = 1: length(t)
        v_frame = eset.expt.track(j).dq.speed * 60;  % cm/min
        toff_frame = eset.expt.track(j).dq.led12Val_toff;  % interpolated time (s) for each frame of track j, 1-by-(number of the track's frame) 
        t_frame = eset.expt.track(j).dq.eti;
        % select the start time in period of turns that happen between t_stim_start(i) and t_stim_start(i)
        v_frame = v_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
        toff_frame = toff_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
        toff_frame_total = [toff_frame_total, toff_frame];
        v_frame_total = [v_frame_total, v_frame];
        if (t(j).startFrame < t_stim_end(i)*frame_rate) && (t(j).endFrame > t_stim_start(i)*frame_rate)
            % add to nperiod after each loop of j
            nperiod = nperiod + (min(t(j).endFrame, t_stim_end(i)*frame_rate) - max(t(j).startFrame, t_stim_start(i)*frame_rate)) / frame_rate / tperiod;
        end
    end
    stepsize = 0.1;  % to get one average speed for stepsize seconds
    [x_toff,y_v, stderror] = meanyvsx(toff_frame_total, v_frame_total, 0:stepsize:tperiod);
    uppercurve = y_v + 0.5*stderror;
    lowercurve = y_v - 0.5*stderror;
    x_tofill = [x_toff, fliplr(x_toff)];  % the x axis of the ploygon to fill
    y_tofill = [lowercurve, fliplr(uppercurve)];
    ax(i) = subplot(length(t_stim_start),1,i);
    pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none'); hold on;  % no edges for the patch
    plot(x_toff, y_v, 'Color', 0.2*[1 1 1]); xline(9, '--'); hold off;
    xlabel('Time toff in period (s)'); ylabel('Average speed of all tracks (cm/min)');
    title([num2str(nperiod), ' periods, ', num2str(i), '-th intensity of stimulation']);
end
sgtitle(['Step size = ', num2str(stepsize), ' s']); linkaxes(ax, 'y');  % align subplots with y axis
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\v_tperiod_stim_group');
savefig(gcf, savename); 


% speed of each long track at a certain period of time of stimulation
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\data.mat');
if isfile(savename)
    load(savename, 'larvae');
else
    clear larvae;
end
figure; stepsize = 0.1;  % to get one average speed for stepsize seconds
% only keep the speed whose start time falls into the i-th intensity of stimulation
for i = 1: length(t_stim_start)
    for j = 1: length(t)
        v_frame = t(j).dq.speed * 60;  % cm/min
        ton_frame = t(j).dq.led12Val_ton;  % interpolated time (s) for each frame of track j, 1-by-(number of the track's frame) 
        t_frame = t(j).dq.eti;
        % select the start time in period of turns that happen between t_stim_start(i) and t_stim_start(i)
        v_frame = v_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
        ton_frame = ton_frame((t_stim_start(i) <= t_frame) & (t_frame < t_stim_end(i)));
        %number of period of stimulation that falls into certain intensity
        t_start = max([t_stim_start(i), eset.expt.elapsedTime(t(j).startFrame + 1)]);  % time (s) of start for track j under i-th  intensity of stimulation 
        t_end = min([t_stim_end(i), eset.expt.elapsedTime(t(j).endFrame-2)]);  % temporal - 2
        nperiod = (t_end - t_start) / tperiod;
        index_subplot = j + (i-1)*length(t);
        ax(index_subplot) = subplot(length(t_stim_start), length(t), index_subplot);
        [x_ton,y_v, ~, stderror] = meanyvsx(ton_frame, v_frame, 0:stepsize:tperiod);  % std
        uppercurve = y_v + 0.5*stderror;
        lowercurve = y_v - 0.5*stderror;
        x_tofill = [x_ton, fliplr(x_ton)];  % the x axis of the ploygon to fill
        y_tofill = [lowercurve, fliplr(uppercurve)];
        pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none'); hold on;  % no edges for the patch
        plot(x_ton, y_v, 'Color', 0.2*[1 1 1]); xline(6, '--'); hold off;
        xlabel('ton (s)'); ylabel('Speed (cm/min)');

        % add new contents to the structure larvae
        larva_index = ['larva', num2str(j)];
        larvae.(larva_index).speed.(stim_color{i}) = y_v;

        title(['(', num2str(t(j).trackNum), ', ', num2str(nperiod), ', ', larvae.(larva_index).response.(stim_color{i}) ')']);

    end
end
sgtitle(['Step size = ', num2str(stepsize), ' s, (track index, nperiod, response)']);  linkaxes(ax, 'y');  % align subplots with y axis
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\v_ton_stim_individual');
savefig(gcf, savename); 


%% Body angle
% Body angle of single larvae verses time in period
j = 1;  % track number
theta_frame = eset.expt.track(j).dq.sbodytheta * 180 / pi;  % degree 0-360
toff_frame = eset.expt.track(j).dq.led12Val_toff;  % interpolated time (s) for each frame of track j, 1-by-(number of the track's frame) 
stepsize = 0.1;  % to get one average speed for stepsize seconds
[x_toff,y_theta, stderror] = meanyvsx(toff_frame, theta_frame, 0:stepsize:tperiod);
uppercurve = y_theta + 0.5*stderror;
lowercurve = y_theta - 0.5*stderror;
x_tofill = [x_toff, fliplr(x_toff)];  % the x axis of the ploygon to fill
y_tofill = [lowercurve, fliplr(uppercurve)];
figure;
pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none'); hold on;  % no edges for the patch
plot(x_toff, y_theta, 'Color', 0.2*[1 1 1]);
xline(9, '--'); hold off;
xlabel('toff (s)'); ylabel(['Body angle of track ', num2str(j), ' (^o)']);

%% Other functions
% to save some variables into a file, so that data of multiple files can be
% plotted together when load the data
% turnStart_total.(['turnStart', num2str(d(x).name(end-15:end-4))]) = turnStart;
% nperiod_total.(['nperiod', num2str(d(x).name(end-15:end-4))]) = nperiod;
% turnrate_array_1 =turnrate_array;  % change the saving name _1 for the next one
% nperiod_array_1 = nperiod_array;
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\data.mat');

save(savename, 'turnStart', 'nperiod')
save(savename, 'turnStart', 'nperiod', '-append')
if isfile(savename)
    save(savename, 'larvae', '-append')
else
    save(savename, 'larvae')
end

load(savename)
turnrate_array = [turnrate_array_1 turnrate_array_2];
nperiod_array = [nperiod_array_1 nperiod_array_2];


isrun = eset.gatherField('isrun');  % for all tracks in the expt, 1 for run, 0 for not
[tx,fracinrun] = meanyvsx(toff, isrun, 0:0.5:tperiod);  % toff transfer frame index to time between 0 and 60s
figure;
plot (tx, fracinrun, 'r');  xline(10, '--'); 
xlabel('Time (s)');
ylabel('Fraction in run for all tracks');
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\frac_run_toff');
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
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\raw_turn');
savefig(gcf, savename); 
% saveas(gcf, savename, 'png');





% Plot one or multiple tracks to determine which to stitch
width = 2048;  % in pixel, in x-axis, geometry of Region of Interest (ROI) read from Image Recorder front panel
height = 2048;  % in pixel, in y-axis
len_pixel = realUnitsPerPixel(eset.expt.camcalinfo);  % how many cm per pixel
color_pad = ['r', 'g', 'b', 'k', 'c', 'm', 'y'];
track_path = stitched_track_index(7 : end);  % index of track to plot, could be [1], or [1, 3, 8]-----------------------------
figure; 
eset.expt.track.plotPath('sloc', 'color', [0.8, 0.8, 0.8]); hold on;  % the larger the whiter
for i = 1 : length(track_path)  %  The index of track to plotPath, should be shorter than color_pad
    eset.expt.track(track_path(i)).plotPath('sloc', color_pad(i));   hold on;  % 
    xy_s = eset.expt.track(track_path(i)).getDerivedQuantity('sloc');
    plot(xy_s(1, 1), xy_s(2, 1), append('o', color_pad(i)));  % o marks start
    plot(xy_s(1, end), xy_s(2, end), append('x', color_pad(i)));  % x marks end
    collisions_index = find(collision_events(:, 1) == track_path(i));  % find the row index of collisions in one track from collision events
    plot(collision_events(collisions_index, 5), collision_events(collisions_index, 6), append('>', color_pad(i)));
end
rectangle('Position', [0, 0, width * len_pixel, height * len_pixel]); 
axis equal;  % use the same length for data unit
title(['Path of tracks ', num2str(track_path), ', with color ', color_pad(1: length(track_path)) ]); xlabel('x (cm)'); ylabel('y (cm)'); hold off;
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\stitch_paths_2');
savefig(gcf, savename); 



% Find the collision location and time, and save collision_events [track_1, x1(pixel), y1(pixel), t1(frame), x1(cm), y1(cm); track_1, x2(pixel), y2(pixel), t1(frame), x2(cm), y2(cm); ...]
stitched_track_index = find([eset.expt.track.nt] ~= 1);  % track index in which collision happened, whose nt(number of tracks) ~= 1.
collision_events = [];  % [track_ID, x(pixel), y(pixel), t(frame), x(cm), y(cm)] in units of pixels and frames
for j = 1 : length(stitched_track_index)  % loop all stitched tracks
    frames_collision = find(eset.expt.track(stitched_track_index(j)).iscollision);  % frames that collision is in progress
    frame_collision_start = frames_collision([true, diff(frames_collision) ~= 1]);  % [collision 1 start frame, collision 2 start frame, ...]
    frame_collision_end = frames_collision([diff(frames_collision) ~= 1, true]);  % this collision end frame will be off a lot
    loc_collision_start_cm = [eset.expt.track(stitched_track_index(j)).pt(frame_collision_start).loc];  % [collision 1 start x, collision 2 start x, ... ; collision 1 start y, collision 2 start y, ...] in cm
    loc_collision_start_pixel = camPtsFromRealPts(eset.expt.camcalinfo, loc_collision_start_cm);  % [x1, x2, ...; y1, y2, ...] in pixel
    % collision event(s) for this single track, [track_1, x1(pixel), y1(pixel), t1(frame), x1(cm), y1(cm); track_1, x2(pixel), y2(pixel), t1(frame), x2(cm), y2(cm); ...]
    collision_event = [repmat(stitched_track_index(j), length(frame_collision_start), 1), transpose(loc_collision_start_pixel), transpose(frame_collision_start), transpose(loc_collision_start_cm)];
    collision_events = [collision_events; collision_event];  % combine collision events of different tracks together
end
% collision_events = fix(collision_events);  % change from double to integer
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\data.mat');
save(savename, 'collision_events')  % run this at the first time to create new saving file
save(savename, 'collision_events', '-append')  % run this at the first time to create new saving file
% save(savename, 'turnrate_array', 'nperiod_array', '-append')  % run this next time to add new data to the same file

