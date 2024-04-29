basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3#T_Bl_Sq_2to7P_15_1_3'
    };
x_cell = {
    1:10
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
    for k = 1 : length(x)  % loop for each expt in the eset
        j = j + 1;
        figlocation_list{j} = fullfile(basedir, ['results', d(x(k)).name(end-16:end-4)]);
    end
end


% create nlarva-by-15 matrix to contain all pturn from figures
pturn_hist_all = [];
for j = 1 : length(figlocation_list)

    % Extract the data of the figure saved in figlocation
    fig = openfig(fullfile(figlocation_list{j}, 'p_turn_no_pause_all.fig'));  % class of Figure
    nLongTracks = (length(fig.Children) - 1) / 3;
    pturn_hist = zeros(nLongTracks, 15);  % number of bin * number of hist for each maggots (5 * 3 = 15)----------
    for i = 2 : length(fig.Children)  % fig.Children is the class of Axes
        title(fig.Children(i), ['Children ', num2str(i)]);
        larva_index = mod(i, nLongTracks) + 1;  % row index
        if (i - 1) / nLongTracks > 2
            stim = 1;  % the first kind of stimulation
        elseif (i - 1) / nLongTracks <= 1
            stim = 3; % the third kind of stimulation
        else
            stim = 2;
        end
        if length(fig.Children(i).Children) == 1  % without xline
            pturn_hist(larva_index, (1 + (stim - 1) * 5) : (stim * 5)) = fig.Children(i).Children.YData;  % 5 values for 1 stimulation
        else  % with xline, the bar object is the second Child
            pturn_hist(larva_index, (1 + (stim - 1) * 5) : (stim * 5)) = fig.Children(i).Children(2).YData;
        end
    end
%     save(fullfile(figlocation_list{j}, 'data'), 'pturn_hist')
    
    pturn_hist_all = [pturn_hist_all; pturn_hist];
    
end

close all;
save('data.mat', 'pturn_hist_all');
