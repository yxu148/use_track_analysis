% @ Yiming
% Run it at the first time
% setupDirectories('Ruben');

startdir1 = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability'; % Start folder, contains data from experiments ---------------
% startdir2 = ...

outputdir1 = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_extracted_small'; % Extracted folder, contains data after extraction, must exist ---------------

extraparams = {{'processing_params', processing_params_paulsmall},{},{}};  % pass the paulsmall extraction setting to the first start folder

% start a log file to record
logname = "processMMF" + datestr(now, 'yyyy-mm-dd')+ ".txt";  % doesn't exist before diary()
logpath = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_process_logs';  % this folder must exist ---------------
pathname = fullfile(logpath, logname);
diary(pathname);
disp(datestr(now));
disp('Processing mmf to bin');

extractinbackground = true;
verbose = false;

startdirs = {startdir1};  % include startdir2 if there is ---------------
outputdirs = {outputdir1}; % ---------------

for i = 1:length(startdirs)
    startdir = startdirs{i};
    outputdir = outputdirs{i};
    % checkboard should be a PNG file with 'check' in the file name
%     checkerLocation = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\behavior box\calibrations';
%     %Find the most recent checkerboard
%     calibs = dir(fullfile(startdir,'calibrations'));
%     calibs = calibs([calibs.isdir] & ~cellfun(@(s) s(1) == '.', {calibs.name}));
%     datenums = cell2mat({calibs.datenum});
%     [~,ind] = max(datenums);
%     checkerLocation = fullfile(startdir,'calibrations',calibs(ind).name);

   

    if (~exist(startdir, 'file'))
        disp (['Cannot find data directory: ' startdir]);
        return;
    end

    % copy correct camera calibration files to start folder------------
    % checkboard should be a PNG file with 'check' in the file name, could
    % some subfolders
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

    for j = 1:length(basedirs)
        try
            if (verbose)
                disp (basedirs{j});
            end
            process_mmfs_gershow(basedirs{j}, dstdirs{j},cc, 'extractinbackground', extractinbackground,  'verbose', verbose, extraparams{i}{:});
        catch me
            disp (me.getReport());
        end
    end
    
end

diary off;  % disable recording log, used before opening the log

disp ('finished with auto process 1');