% @Yiming

warning backtrace off
startdir1 = 'G:\AS-Filer\PHY\mmihovil\Shared\Isabel\temporal'; % Start folder, contains data from experiments, and *check*.png ---------------
% startdir2 = ...

outputdir1 = 'G:\AS-Filer\PHY\mmihovil\Shared\Isabel\temporal_extracted'; % Extracted folder, contains data after extraction, must exist ---------------


% start a log file to record
logname = "processBIN" + datestr(now, 'yyyy-mm-dd')+ ".txt";  % doesn't exist before diary()
logpath = 'G:\AS-Filer\PHY\mmihovil\Shared\Isabel\temporal_process_logs';  % this folder must exist ---------------
pathname = fullfile(logpath, logname);
diary(pathname);
disp(datestr(now));
disp('Processing bin to mat');
disp('Comments: Change the threshold to 20 and area to 10 in .bxx file')  % Comments here -----------------
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool();
    closepool = true;
else
    closepool = false;
end
verbose = true;

%{
%processLoops(startdir1,outputdir);  
%processLoops(startdir2,outputdir);


for i=1:length(startdir)
    disp(['Processing data from ' char(startdir(i))]);
    processLoops(char(startdir(i)),outputdir);    
end
%}

extractinbackground = true;

% if (exist(fullfile (outputdir, 'processingProblems.mat'), 'file'))
%     load(fullfile (outputdir, 'processingProblems.mat', 'processingProblems'));
% else
%     processingProblems = {};
% end

startdirs = {startdir1};  % include startdir2 if there is-----------------
outputdirs = {outputdir1};
for i = 1:length(startdirs)
    startdir = startdirs{i};
    outputdir = outputdirs{i};
    clear dstdirs relpath
    existsAndDefault('redobtd', false);

    if (~exist(startdir, 'file'))
        disp (['Cannot find data directory: ' startdir]);
        return;
    end
    
    cc = CalibrationCopier();
    cc.startdir = fullfile(startdir, 'calibrations');
    cc = cc.updateAll();

    basedirs = findBaseDirectories(startdir);

    for j = 1:length(basedirs)
        bd = java.io.File(basedirs{j});
        sd = java.io.File(startdir);
        relative = sd.toURI().relativize(bd.toURI()).getPath();
        relpath{j} = fullfile(char(relative),''); %#ok<SAGROW> %this is required to convert '/' to '\' on windows
        dstdirs{j} = fullfile(outputdir, char(relative)); %#ok<SAGROW>
    end

   % processingOptions = {'trimrectpixels', [10 10 2582 1934]}; 
   % CHANGED 2/7/2015 by MHG to not trim tracks
%     processingOptions = {'trimrect', [], 'trimrectpixels',[],'buffer', [],'ccInSupDataDir', true}; %no trimming!
% Default 'frameDiff' is 7, default 'maxDist' is 0.1 (cm)---------
    processingOptions = {'trimrect', [], 'trimrectpixels',[],'buffer', [],'ccInSupDataDir', true, 'frameDiff', 1000, 'maxDist', 0.3};
    for j = 1:length(dstdirs)
%         if (any(strcmpi(dstdirs{j}, processingProblems)))
%             warning ('ap2:pp', [dstdirs{j} ' is marked as a problem -- skipping']);
%         end
         try
            if (verbose)
                disp (['processing ' dstdirs{j}]); %#ok<*UNRCH>
            end
            %now looks for camera calibration in supplemental data
            %directory; should be copied over by process_mmfs_gershow, but
            %we give it a second chance here
            [np, nc, nf] = processToMatfiles_gershow(basedirs{j}, dstdirs{j},'processingOptions', processingOptions, 'redobtd', redobtd, 'verbose', verbose, 'CalCopier', cc);
            if (np <= 0 && nf <= 0)
                 if (verbose)
                    disp ('nothing to do');
                 end
                continue;
            end
             if (verbose)
                disp (dstdirs{j});
                disp ([num2str(nc) ' bin files already processed']);
                disp ([num2str(np) ' bin files newly processed']);
             end
            if (nf > 0)
                disp([dstdirs{j} num2str(nf) ' bin files failed to process']);
                %disp ([num2str(nf) ' bin files failed to process']); 
            end
         catch me
             disp([dstdirs{j} ' failed to process.']);
             disp (me.getReport());
%              processingProblems = [processingProblems dstdirs{j}]; %#ok<AGROW>
             
         end
    end
    
end
%save(fullfile (outputdir, 'processingProblems.mat', 'processingProblems'));
if (closepool)
    delete(poolobj);
end
diary off;
disp ('finished autoprocess2');