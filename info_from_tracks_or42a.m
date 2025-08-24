% information from pturn for short tracks

%% create pturn matrix from the structure larvae
basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data_new\variability_new_try_extracted\Or42a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3_#T_Bl_Sq_2to7P_15_1_3_'
    };
x_cell = {
    1:22
%     [25, 26, 27, 28, 29, 30, 31, 32, 33, 34]  % load the x-th set of data from the first basedir to analyze
    };

% create figlocation_list, a cell containing all folders of results
nexp = length([x_cell{:}]);  % number of experiments
figlocation_list = cell(nexp, 1);  % initialize a nexpt-by-1 empty cell
j = 0;  % the j-th experiment
for folder_index = 1 : length(x_cell)  % loop for each basedir folder
    basedir = basedir_cell{folder_index};
    x = x_cell{folder_index};
    d = dir(fullfile(basedir, 'matfiles', '*.mat'));
    mkdir(fullfile(basedir, 'results'));  % create a overall results folder to contain general data
    for k = 1 : length(x)  % loop for each expt in the eset
        j = j + 1;
        figlocation_list{j} = fullfile(basedir, ['results', d(x(k)).name(end-16:end-4)]);
    end
end


frame_rate = 20;  % number of frames per second
latest_start = 120;  % seconds, select the tracks that start earlier than latest_start time. 80% length
earliest_end = 1680;  % seconds, select the tracks that end later than earliest_end time.
stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start

pturn_all = [];
endFrame_all = [];  % to order the tracks
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'tracks');
    for i = 1 : length(fieldnames(tracks))  % loop for each larva in larvae
        % select the tracks that start early and end early.
        if (getfield(tracks, ['track', num2str(i)], 'startFrame') <= latest_start * frame_rate) && (getfield(tracks, ['track', num2str(i)], 'endFrame') < earliest_end * frame_rate)
            pturn_temp = [];  % for each track
            for s = 1 : length(stim_color)  % put pturn together
                pturn_temp = [pturn_temp, getfield(tracks, ['track', num2str(i)], 'pturn', stim_color{s})];
            end
            pturn_all = [pturn_all; pturn_temp];
            endFrame_all = [endFrame_all; getfield(tracks, ['track', num2str(i)], 'endFrame')];
        end
    end  % end looping all tracks
end  % end looping all expts

% sort the tracks based on decending endFrame
[endFrame_all_des, sortIdx] = sort(endFrame_all, 'descend');
pturn_all = pturn_all(sortIdx, :);
writematrix(pturn_all, strcat(basedir_cell{1}, '\results', '\pturn_short.xlsx'));


figure; imagesc(pturn_all, [0, 1]); 
ylabel('Index of larva'); set(gca,'XTick',[])
cbar = colorbar; cbar.Label.String = 'Turn possibility with the 3-s bin';
cbar.Ticks = 0:0.2:1;  % tick every 0.2
savename = strcat(basedir, '\results', '\results_fig', '\pturn_short');
savefig(gcf, savename);


%% save the current file to \results in the basedir
savename = strcat(basedir_cell{1}, '\results');
copyfile('info_from_tracks.m', savename)  % copy the file to the location