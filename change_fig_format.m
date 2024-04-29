% Location of the figures need to be changed
results_dir = 'D:\Apps-SU\Matlab-Track-Analysis-SkanataLab\user specific\Yiming\marchmeeting2024';
d = dir(fullfile(results_dir, '*hist_pturn_stim*.fig'));

mkdir(fullfile(results_dir, 'results_epsc'));

figs = fullfile(results_dir, {d.name});
for i = 1 : length(figs)
    [~, fig_name] = fileparts(figs{i});
    fig = openfig(figs{i});
    fig.Renderer='painter';  % use this so every line/dot is editable
    savename = fullfile(results_dir, 'results_epsc', fig_name);
    saveas(fig, savename, 'epsc');  % EPS level 3 color, .eps
    close;
end