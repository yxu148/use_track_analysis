% information from speed

%% create speed matrix from the structure larvae
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

% create speed_all in order of categories and save in data.mat
% speed_all = [];  % nlarva-by-15 array, 3s bins
stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start
ntype = 2 ^ (length(stim_color));  % int, -1 if exclude 'dark'
speed_type = cell(1, ntype);  % 1-by-ntype cell, 1-indexed, each cell is one type, thing in cell is an array
speed_run_type = cell(1, ntype);
response_type = cell(1, ntype);
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop for each larva in larvae
        if getfield(larvae, ['larva', num2str(i)], 'valid')
%         if 1
            speed = getfield(larvae, ['larva', num2str(i)], 'speed');
            speed_run = getfield(larvae, ['larva', num2str(i)], 'speed_run');
            % to load speed in order of response
            larva_index = ['larva', num2str(i)];
            response = '';  % e.g. '001'
            for s = 1 : length(stim_color)  % -1 if exclude the last stim_color, 'dark'-----------------
                response = [response, getfield(larvae, ['larva', num2str(i)], 'response', stim_color{s})];
            end
            typeth = bin2dec(response) + 1;  % positive int
            speed_temp = [];
            speed_run_temp = [];
            for s = 1 : length(stim_color)  % put speed together
                speed_temp = [speed_temp, getfield(speed, stim_color{s})];
                speed_run_temp = [speed_run_temp, getfield(speed_run, stim_color{s})];
            end
            speed_type{typeth} = [speed_type{typeth};  speed_temp];
            speed_run_type{typeth} = [speed_run_type{typeth};  speed_run_temp];
            response_type{typeth} = [response_type{typeth}; response];

        end  % end if valid larvae
    end  % end looping each larva in one experiment
end  % end loop each experiment

% put array from different cell together to one array
[speed_all, speed_run_all, response_all] = deal([]);
for s = 1 : length(speed_type)
    speed_all = [speed_all; speed_type{s}];
    speed_run_all = [speed_run_all; speed_run_type{s}];
    response_all = [response_all; response_type{s}];
end

savename = strcat(basedir_cell{1}, '\results', '\data.mat');  % there should be only 1 basedir in basedir_cell
if isfile(savename)
    save(savename, 'speed_all', 'speed_type', 'speed_run_all', 'speed_run_type', '-append');
else
    save(savename, 'speed_all', 'speed_type', 'speed_run_all', 'speed_run_type')
end
% writematrix(speed_run_all, strcat(basedir_cell{1}, '\results_new', '\speed_run.xlsx'));
% writematrix(speed_all, 'speed.csv');

%% Plot speed from data.mat ===================================================
basedir = basedir_cell{1};
savename = strcat(basedir, '\results', '\data.mat');  % there should be only 1 basedir in basedir_cell
load(savename, 'speed_all', 'speed_run_all', 'speed_type');

stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start
ntype = 2 ^ (length(stim_color));  % int, -1 if exclude 'dark'
name_type = dec2bin(0:ntype-1);  % 1-indexed
name_type_cell = {};  % 1-by-numofBehaviorType cell
for temp = 1 : length(name_type)
    name_type_cell(temp) = {name_type(temp, :)};
end
nlarva_type = [];  % number of larvae belonging to every type in order of dec2bin(0:63)
for s = 1 : length(speed_type)
    nlarva_type = [nlarva_type, size(speed_type{s}, 1)];
end
% find index to seperate diff types, only work if nlarva_type is longer than 1
boundary_type = zeros(1, length(nlarva_type) - 1) - 1;  % initialize with -1
boundary_type(1) = nlarva_type(1);  % 1-indexed
boundary_temp = nlarva_type(1);
for temp = 2 : length(nlarva_type)
    boundary_temp = boundary_temp + nlarva_type(temp);
    boundary_type(temp) = boundary_temp;
end
% replace the name of types with less than 5 larvae to '_'
name_type_cell(nlarva_type<5) = {'\_'};

% plot the matrix speed_all
figure;
imagesc(speed_all, [0, 4]); 
xline([1, 150, 300], 'w', stim_color); yline(boundary_type + 0.5, 'w', name_type_cell);
ylabel('Index of larva'); set(gca,'XTick',[])
cbar = colorbar; cbar.Label.String = 'Speed (cm/min)'; %cbar.Limits(2) = 1;
mkdir(fullfile(basedir, '\results', 'results_fig'));
savename = strcat(basedir, '\results', '\results_fig', '\speed_all');
savefig(gcf, savename);

% plot the matrix speed_run_all
figure;
imagesc(speed_run_all, [0, 4]); 
xline([1, 150, 300], 'w', stim_color); yline(boundary_type + 0.5, 'w', name_type_cell);
ylabel('Index of larva'); set(gca,'XTick',[])
cbar = colorbar; cbar.Label.String = 'Run Speed (cm/min)';
mkdir(fullfile(basedir, '\results', 'results_fig'));
savename = strcat(basedir, '\results', '\results_fig', '\speed_run_all');
savefig(gcf, savename);

%% save the current file to \results in the basedir
savename = strcat(basedir_cell{1}, '\results');
copyfile('info_from_speed.m', savename)  % copy the file to the location