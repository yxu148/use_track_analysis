% Location of the figures need to be changed
results_dir = 'G:\AS-Filer\PHY\mmihovil\Shared\Isabel\chemotaxis_extracted\unc@na\C_Re_0_7dpa_food\results_202401181654';
d = dir(fullfile(results_dir, '*.fig'));

mkdir(fullfile(results_dir, '..', 'results_epsc'));

figs = fullfile(results_dir, {d.name});
for i = 1 : length(figs)
    [~, fig_name] = fileparts(figs{i});
    fig = openfig(figs{i});
    fig.Renderer='painter';  % use this so every line/dot is editable
    savename = fullfile(results_dir, '..', 'results_epsc', fig_name);
    saveas(fig, savename, 'epsc');  % EPS level 3 color, .eps
    close;
end