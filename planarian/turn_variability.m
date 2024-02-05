% Fraction in run for all tracks verses time
isrun = eset.gatherField('isrun');  % for all tracks in the expt, 1 for run, 0 for not, 1-by-(sum of all tracks' frame)
t_all = eset.gatherField('eti');  % interpolated time (s) for each frame of each track, 1-by-(sum of all tracks' frame) 
stepsize = 0.5;
[tx,fracinrun] = meanyvsx(t_all, isrun, 0: stepsize :1200);  % toff transfer frame index to time between 0 and 60s
figure;
plot (tx, fracinrun, 'r'); ylim([0, 1]);
xlabel('Time (s)'); ylabel('Fraction in run for all tracks');
title(['Stepsize ', num2str(stepsize), ' s']);
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\frac_run_toff');
savefig(gcf,savename);


% Mean speed of all planarians (cm/min) verses time
v_all = eset.gatherField('speed') * 60;  % cm/min
t_all = eset.gatherField('eti');  % interpolated time (s) for each frame of each track, 1-by-(sum of all tracks' frame) 
stepsize = 0.5;  % to get one average speed for stepsize seconds
[tx,vmean, stderror] = meanyvsx(t_all, v_all, 0:stepsize:1200);
uppercurve = vmean + 0.5*stderror;
lowercurve = vmean - 0.5*stderror;
x_tofill = [tx, fliplr(tx)];  % the x axis of the ploygon to fill
y_tofill = [lowercurve, fliplr(uppercurve)];
figure;
pathObj = fill(x_tofill, y_tofill, 0.8*[1 1 1], 'LineStyle', 'none'); hold on;  % no edges for the patch
plot (tx, vmean, 'Color', 0.2*[1 1 1]); 
ylim([0, 8]);
xlabel('Time (s)'); ylabel('Mean speed of all planarians (cm/min)'); 
title(['Stepsize ', num2str(stepsize), ' s']);
hold off;
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], '\vmean_t');
savefig(gcf,savename);

j = 1;  % the one track to plot
t_single = eset.expt.track(j).dq.eti;
v_single = eset.expt.track(j).dq.speed * 60;
figure;
plot(t_single, v_single);
xlabel('Time (s)'); ylabel('Speed of Single Track (cm/min)');
title(['Track ', num2str(j)]);
savename = strcat(basedir,['\results', d(x).name(end-16:end-4)], ['\vsingle_t_track', num2str(j)]);
savefig(gcf,savename);


