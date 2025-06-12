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


% save multiple figs into PDF
basedir = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Or42a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3_#T_Bl_Sq_2to7P_15_1_3_';
d = dir(fullfile(basedir, 'results_2025*'));
filename = fullfile(basedir, 'results', 'vRun_combined.pdf');
for k = 1 : length(d)  % k-th experiment
    figlocation = fullfile(d(k).folder, d(k).name, 'vRun_ton_individual.fig');
    fig = openfig(figlocation);
    exportgraphics(fig, filename, 'Append', true);
    close(fig);
end


% Edit the copy of the figure saved in figlocation
figlocation = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\A27h@Chrimson(3)\T_Bl_Sq_5to125P_20_every_5P\results';
copyfile(fullfile(figlocation, 'speed_mean_t_allExpt.fig'), fullfile(figlocation, 'speed_mean_t_allExpt_rescale.fig'));  % source, destination
fig = openfig(fullfile(figlocation, 'speed_mean_t_allExpt_rescale.fig'));  % class of Figure
ax = gca;  % current axes (with two y‐axes internally)
hLeftLines = findobj(ax, 'Type','line', 'YAxis','left');  % find lines attached to the LEFT y‐axis
a = 0.6;  b =  2.6;   % source limits
c =  0;  d =  2.5;   % desired limits
% remap just those curves
for k = 1:numel(hLeftLines)
    y = hLeftLines(k).YData;
    y2 = (y - a) / (b - a) * (d - c) + c;
    hLeftLines(k).YData = y2;
end
% now fix the left axis display only:
yyaxis left
ylim([c d]);
savefig(fig, fullfile(figlocation, 'speed_mean_t_allExpt_rescale.fig'));