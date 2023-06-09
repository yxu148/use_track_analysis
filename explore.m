
vh = eset.expt(1).track(1).getDerivedQuantity('vhead');  % velocity of the head, vector
vt = eset.expt(1).track(1).getDerivedQuantity('vtail');  % velocity of the tail
tx = eset.expt(1).track(1).getDerivedQuantity('eti');  % time
ih = eset.expt(1).track(1).getDerivedQuantity('ihead');  % x and y pos of the head
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

reorientation_rate = eset.makeReorientationHistogram('led2Val_toff', 0:1:20);  % per minute ?????????????????????????
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
plot (xx, [rbefore;rafter], 'bo-'); xlim([-3 3]); ylabel('Rate of Turn'); xlabel('Before 0 and after 1');


