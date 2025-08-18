

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handedness Analysis, version 1.0 (11/01/2016)
%
% Combines JustLoadAndSegment.m with Misc2B.m, with some changes
%
% Idea is to save a table of information about RUNS, and
% also the trajectories, so that subsequent analysis in Igor
% can stitch tracks and determine handedness statistics
%

% PART 0: USER INPUT
fileName = input('\nEnter file name (without `.` or `\\`): ', 's');

folder=0;
while (folder ~= 1) && (folder ~= 2)
    fprintf('\nWhere do you want to save the results? \n')
    disp(['1 : current folder : ' pwd])
    disp('2 : choose a folder')
    folder=input('Type 1 or 2 : ');
end 
if folder==2
    folder_name = uigetdir;
    saveFileName = [folder_name '\' fileName];
elseif folder ==0 
    saveFileName = '';
end

saveFileName3 = [saveFileName '_RUNS.txt'];
saveFileName4 = [saveFileName '_TRAJ'];
saveFileName5 = [saveFileName '_HS.txt'];




% PART I: LOAD THE FILES AND CLEAN UP TRACKS
% (we might not do much track cleaning for this, compared to the normal
% spatial analysis, as we might want to keep short tracks for stitching
% with longer tracks)

% (A) LOAD FILE(S)
if (~exist('eset','var'))
    eset=ExperimentSet.fromFiles();
end

% (B) CLEAN THE TRACKS
%
existsAndDefault('cleanEset', 'true');
if (cleanEset)
    ecl = ESetCleaner;
    ecl.minHTValid = 0.65;
   % ecl.minDist = 50;
   
    ecl.minSpeed = 0.65;
    ecl.minPts = 1000;
    ecl.clean(eset);
    
    cleanEset = false;
end

% (C) FIX HEAD-TAIL ORIENTATION 
existsAndDefault('fixht','true');
if (fixht)
    eset.executeTrackFunction('fixHTOrientation');
    fixht = false;
end

% 
% (D) SET SEGMENTATION SPEED
existsAndDefault('autosetspeeds', true);
if (autosetspeeds)
    eset.executeTrackFunction('setSegmentSpeeds');
    autosetspeeds = false;
end

% (E) SEGMENT THE TRACKS
existsAndDefault('segment', true);
if (segment)
    eset.executeTrackFunction('segmentTrack');
    segment = false; 
end

%%%%%%%%%%%%%%%%%%%
%% PART II: SAVE RUN INFORMATION

% (A) choose set number and name of file to save
setNumber = 1;
%fileName = 'TEST.txt';

% (B) gather list of runs and their important indexes
RUNall = eset.gatherField('run');
RUNallTracks = [RUNall.track];
RUNtrackstart = [RUNallTracks.startFrame];
RUNendind = [RUNall.endInd];
RUNinds = [RUNtrackstart] + int32([RUNendind]);
RUNvalid = [RUNendind]>0;
%
RUN = RUNall(RUNvalid);
RUNtracks = [RUN.track];
RUNexpts = [RUNtracks.expt];
RUNstartind = [RUN.startInd];
RUNendind = RUNendind(RUNvalid);

% (C) Basic information about each run
RUNtheta = [RUN.meanTheta];     % mean direction of travel during run
RUNtheta0 = [RUN.startTheta];   % initial direction of travel in the run
RUNlength = [RUN.pathLength];   % path length in pixels
RUNtime = [RUN.runTime];        % i.e. duration of the run in seconds
%
RUNtrackNum = [RUNtracks.trackNum];
for i=1:length(RUN)
   RUNexptNum(i)=find(eset.expt==RUNexpts(i));
   RUNtime0(i) = RUNtracks(i).dq.eti(RUNstartind(i));
   RUNsetNum(i)=setNumber;
end
% RUNsetNum = zeros(1,length(RUN));
%[RUNexptNum] = find(eset.expt==[RUNexpts]);
%RUNtime0 = [RUNtracks.dq.eti([RUNstartind])];


% (D) Make empty arrays to hold more information
RUNxpos = zeros(1,length(RUN));
RUNreoYN = boolean(zeros(1,length(RUN)));
RUNreoHS = zeros(1,length(RUN));
RUNreotheta1 = zeros(1,length(RUN));
RUNreotheta2 = zeros(1,length(RUN));
RUNreoHS1 = zeros(1,length(RUN));
%
RUNx0 = zeros(1,length(RUN));
RUNy0 = zeros(1,length(RUN));
RUNx1 = zeros(1,length(RUN));
RUNy1 = zeros(1,length(RUN));

% (E) loop through all the runs and get what we need
for i=1:length(RUN)
    
    % start position of each run
    runStartInd = RUN(i).startInd;
    %loc = RUN(i).track.pt(runStartInd).loc;
    loc = RUN(i).track.dq.sloc(:,runStartInd);
    RUNx0(i) = loc(1);
    RUNy0(i) = loc(2);
    
    % find the end position of each run
    runEndInd = RUN(i).endInd;
   % if(runEndInd>RUN(i).track.npts)
    %   runEndInd = RUN(i).track.npts; 
   % end
    %loc = RUN(i).track.pt(runEndInd).loc;
    loc = RUN(i).track.dq.sloc(:,runEndInd);
    RUNxpos(i) = loc(1);
    %
    RUNx1(i) = loc(1);
    RUNy1(i) = loc(2);
    
    % find info. about the turn at the end of the run, if there is one
    if(~isempty(RUN(i).nextReorientation))
       %if(~isempty(RUN(i).nextReorientation.headSwing))
           %
           RUNreoYN(i)=true;
           %
           RUNreoHS(i)=RUN(i).nextReorientation.numHS;
           % RUNreoHS(i)=0;
           %
           %RUNreoHS1(i) = RUN(i).nextReorientation.headSwing(1).maxTheta;
           RUNreoHS1(i) = 0;
           %
           RUNreotheta1(i) = RUN(i).nextReorientation.prevDir;
           %
           RUNreotheta2(i) = RUN(i).nextReorientation.nextDir;        
           %
       %end
    end
     
end

% (F) Save the results in a text file
fileID = fopen(saveFileName3,'w');
fprintf(fileID,'%14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\r\n','set','expt','track','time0','reoYN','runQ','runL','runT','runX','reo#HS','reoQ1','reoQ2','reoHS1','runQ0','runX0','runY0','runX1','runY1');
A = [RUNsetNum;RUNexptNum;RUNtrackNum;RUNtime0;RUNreoYN;RUNtheta;RUNlength;RUNtime;RUNxpos;RUNreoHS;RUNreotheta1;RUNreotheta2;RUNreoHS1;RUNtheta0;RUNx0;RUNy0;RUNx1;RUNy1];
fprintf(fileID,'%14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\r\n',A);
fclose(fileID);

%%%%%%%%%%%%%%%%%%%
% PART III: SAVE HEAD SWEEP INFORMATION

% (A) choose set number 
setNumber = 1;

% (B) gather list of head sweeps and their important indexes
HSall = eset.gatherField('headSwing');
HSvalid = [HSall.valid];
%HS = HSall(HSallValid);
HS = [HSall];

HStracks = [HS.track];
HSexpts = [HStracks.expt];
HSprevRun = [HS.prevRun];
%HSreos = [HSprevRun.nextReorientation];
%
HSstartInd = [HS.startInd];
HSendInd = [HS.endInd];
HSmaxInd = [HS.maxInd];

% (C) Basic information about each head swing
HStrackNum = [HStracks.trackNum];
HSaccepted = [HS.accepted];
HSdir = [HS.sign];
HSnum = [HS.num];
HStheta1 = [HS.prevDir];
HStheta2 = [HS.nextDir];
HSsize = [HS.maxTheta];

% (D) Make empty arrays to hold more information
HSsetNum = zeros(1,length(HS));
HSexptNum = zeros(1,length(HS));
HSprevRunNum = zeros(1,length(HS));
HSt0 = zeros(1,length(HS));
HSt1 = zeros(1,length(HS));
HStmax = zeros(1,length(HS));
HSx = zeros(1,length(HS));
HSy = zeros(1,length(HS));

% (E) Loop through the head swings to get what we need
for i=1:length(HS)

    HSsetNum(i) = setNumber;
    HSexptNum(i) = find(eset.expt==HSexpts(i));
    HSprevRunNum(i) = find(HStracks(i).run==HSprevRun(i));

    HSt0(i) = HStracks(i).dq.eti(HSstartInd(i));
    HSt1(i) = HStracks(i).dq.eti(HSendInd(i));
    HStmax(i) = HStracks(i).dq.eti(HSmaxInd(i));

    HSx(i) = HStracks(i).dq.sloc(1,HSstartInd(i));
    HSy(i) = HStracks(i).dq.sloc(2,HSstartInd(i));
    
end

% (F) Save the head swing information in a text file
fileID = fopen(saveFileName5,'w');
fprintf(fileID,'%14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\r\n','set','expt','track','run','valid','num','accept','dir','size','theta1','theta2','t0','t1','tmax','x0','y0');
B = [HSsetNum;HSexptNum;HStrackNum;HSprevRunNum;HSvalid;HSnum;HSaccepted;HSdir;HSsize;HStheta1;HStheta2;HSt0;HSt1;HStmax;HSx;HSy];
fprintf(fileID,'%14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\t %14f\r\n',B);
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART IV: Save Track Trajectories
% 
% For now just save into a single Excel file, but could
% later incorporate code from TrackSaving2.m if we wanted
% to save each track as a text file
%
% Also, changing some of this so that we save as pixels and not mm


fprintf('\n\n - DATA SAVING...\n')
saveFileName = input('\nEnter trajectory file name (without `.` or `\\`) : ', 's');

folder=0;
while (folder ~= 1) && (folder ~= 2)
    fprintf('\nWhere do you want to save the results? \n')
    disp(['1 : current folder : ' pwd])
    disp('2 : choose a folder')
    folder=input('Type 1 or 2 : ');
end
 
if folder==2
    folder_name = uigetdir;
    saveFileName = [folder_name '\' saveFileName];
end

saveFileName4 = [saveFileName '_TRACKS'];

%SEXCEL FILE:
warning('off', 'MATLAB:xlswrite:AddSheet');

%TRACKS:
k = 1;
for i=1:length(eset.expt)
    for j=1:length(eset.expt(i).track)
        % Labels
        sheetName = ['Track' num2str(k) '(e' num2str(i) 't' num2str(j) ')'];
        labels = {'time','x','y','run?'};
        % Times:
        times = eset.expt(i).track(j).dq.eti;
        numPoints = length(times);
        % Positions:
        pos = eset.expt(i).track(j).getDerivedQuantity('sloc');
        xpos = pos(1,:);
       % xpos = xpos*lengthPerPixel;
        ypos = pos(2,:);
      %  ypos = ypos*lengthPerPixel;
        % Run Y/N:
        runYN = eset.expt(i).track(j).isrun;
        % Body Angle:
        %bodyangle = eset.expt(i).track(j).getDerivedQuantity('sbodytheta');
        % Combined Matrix:
        combined = [times;xpos;ypos;runYN];
        combined = transpose(combined);
        combinedrange = ['A2:D' num2str(numPoints+1)];
        % Write to Excel:
        xlswrite(saveFileName4,labels,sheetName);
        xlswrite(saveFileName4,combined,sheetName,combinedrange);
        %
        k=k+1;
    end
end
clear times pos xpos ypos runYN combined combinedrange;

fprintf('\ndone\n');



