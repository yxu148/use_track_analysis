figlocation = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3\results_202308251718';


% Extract the data of the figure saved in figlocation, 
% This won't change the figure, so open the original plot
fig = openfig(fullfile(figlocation, 'p_turn_no_pause_all.fig'));  % class of Figure
nLongTracks = (length(fig.Children) - 1) / 3;
pturn_hist = zeros(nLongTracks, 18);  % number of bin * number of hist for each maggots (6 * 3 = 18)
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
    pturn_hist(larva_index, (1 + (stim - 1) * 6) : (stim * 6)) = fig.Children(i).Children.YData;
end
save(fullfile(figlocation, 'data'), 'pturn_hist')

% Edit the copy of the figure saved in figlocation
copyfile(fullfile(figlocation, 'p_turn_no_pause_all.fig'), fullfile(figlocation, 'p_turn_no_pause_all_cshl2023.fig'));
fig = openfig(fullfile(figlocation, 'p_turn_no_pause_all_cshl2023.fig'));  % class of Figure
set(fig, 'Units', 'inches');
for i = 2 : length(fig.Children)  % the first children is the big title, others are the Axes of subplots
    fig.Children(i).YLim = [0, 0.8];
    fig.Children(i).XLim = [-3, 21];
    fig.Children(i).Position = [fig.Children(i).Position(1:2) 0.08 0.08];  % left, bottom, width, height
    fig.Children(i).XLabel.String = '';
    fig.Children(i).YLabel.String = '';
    title(fig.Children(i), '');
    xline(fig.Children(i), 10, 'k--');
end
fig.Children(1).String = '';  % delete the big title
savefig(fig, fullfile(figlocation, 'p_turn_no_pause_all_cshl2023.fig'));


