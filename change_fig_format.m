% Location of the figures need to be changed
results_dir = 'C:\Users\yxu148\OneDrive - Syracuse University\General\Yiming\plots for Mirna';
d = dir(fullfile(results_dir, '*.fig'));

mkdir(fullfile(results_dir, 'epsc'));

figs = fullfile(results_dir, {d.name});
for i = 1 : length(figs)
    [~, fig_name] = fileparts(figs{i});
    fig = openfig(figs{i});
    fig.Renderer='painter';  % use this so every line/dot is editable
    savename = fullfile(results_dir, 'epsc', fig_name);
    saveas(fig, savename, 'epsc');  % EPS level 3 color, .eps
    close;
end