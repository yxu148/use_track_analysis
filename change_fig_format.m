% Location of the figures need to be changed
results_dir = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3#T_Bl_Sq_2to7P_15_1_3\results\results_fig';
d = dir(fullfile(results_dir, '*hist_pturn_denoise_average*.fig'));

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