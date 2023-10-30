figlocation_list = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\stitch_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3\results_202307271555_stitched',
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\stitch_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3\results_202307281451',
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\stitch_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3\results_202308021332',
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3\results_202308161356',
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3\results_202308241620',
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3\results_202308251718',
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\stitch_extracted\Gr21a@Chrimson(3)_s\T_Re_Sq_318to532P_20_2_3#T_Bl_Sq_2,5to6P_20_1_3\results_202309141250'
    };
    
pturn_hist_all = [];

for j = 1 : length(figlocation_list)
    
    % Extract the data of the figure saved in figlocation
    fig = openfig(fullfile(figlocation_list{j}, 'p_turn_no_pause_all.fig'));  % class of Figure
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
    save(fullfile(figlocation_list{j}, 'data'), 'pturn_hist')
    
    pturn_hist_all = [pturn_hist_all; pturn_hist];
    
end

close all;
save('data.mat', 'pturn_hist_all');

% ===================================================================
load('data.mat', 'pturn_hist_all');
denoise = true;

if denoise  % the p of first bin minuses the average of rest bins
    Pblue = pturn_hist_all(:,1) - mean(pturn_hist_all(:,2:6), 2);
    Pred = pturn_hist_all(:,7) - mean(pturn_hist_all(:,8:12), 2);
    Pbluered = pturn_hist_all(:,13) - mean(pturn_hist_all(:,14:18), 2);
else
    Pblue = pturn_hist_all(:,1);
    Pred = pturn_hist_all(:, 7);
    Pbluered = pturn_hist_all(:, 13);
end


figure;  % p( 0<t<3s, blue+red) as a function of p( 0<t<3s, blue) and p( 0<t<3s, red)
x = Pblue;  % p( 0<t<3s, blue)
y = Pred;
sz = 65;  % size of dots in pixels
% c = 1 - pturn_hist_all(:, 13) * [1 1 1];  % the larger p, the smaller c, the darker
c = Pbluered;
scatter(x, y, sz, c, 'filled')
axis equal;
xlabel('Pturn during the first 3 s of blue'); ylabel('Pturn during the first 3 s of red');
cbar = colorbar;
colormap parula;
cbar.Label.String = 'Pturn during the first 3 s of both blue and red';
if denoise
    xlim([-0.2, 0.6]); ylim([-0.2, 0.6]); savefig(gcf, 'pbluered-pb_pr_denoise')
else
    savefig(gcf, 'pbluered-pb_pr')
end


figure;  % p( 0<t<3s, blue+red) as a function of p( 0<t<3s, blue) + p( 0<t<3s, red)
x = Pblue + Pred;
y = Pbluered;
sz = 65;
scatter(x, y, sz, 'filled', 'DisplayName', 'Maggots')
lfit = lsline;  % least-squares line
lfit.DisplayName = ['y = ', num2str(diff(lfit.YData) / diff(lfit.XData)), 'x + ', num2str(lfit.YData(1))];
hold on; plot([0, 1], [0, 1], 'DisplayName', 'y = x'); hold off;
axis equal; xlim([0,max(x) + 0.2]); ylim([0, max(y) + 0.2]);
xlabel('Pturn during the first 3 s of blue plus that of red'); ylabel('Pturn during the first 3 s of both blue and red');
legend();
if denoise
    savefig(gcf, 'pbluered--pb+pr_denoise')
else
    savefig(gcf, 'pbluered--pb+pr')
end


figure;  % Raw data of Pturn at the first 3 s under different color LEDs
nmaggots = size(pturn_hist_all, 1);
x = transpose(1 : nmaggots);
y = [Pblue, Pred, Pbluered];  % blue, red, both
c = [0 0 1; 1 0 0; 0 0 0];  % color in RGB vector
mkr = ['_'; '|'; '+'];
for i = 1 : 3  % have to scatter multiple times if using different signs for different sets of data
    scatter(x, y(:, i), [], c(i, :), mkr(i), 'LineWidth', 2);
    hold on;
end
xline(x, 'LineWidth', 0.1);
legend('Blue', 'Red', 'Blue and Red');
xlabel('Indexes of Maggots'); ylabel('Pturn during the first 3s of specific color of LED');
hold off;
if denoise
    savefig(gcf, 'ps--eachmaggot_denoise')
else
    savefig(gcf, 'ps--eachmaggot')
end


figure;  % 3-d plot
x = Pblue;
y = Pred;
z = Pbluered;
plot3(x, y, z, 'o');
xlabel('Pturn_{blue}'); ylabel('Pturn_{red}'); zlabel('Pturn_{blue+red}');
axis equal; grid on;


figure;  % 3-d plane fit
x = Pblue;
y = Pred;
z = Pbluered;
[sf, goodness_of_fit, fitting_algorithm] = fit([x, y], z, 'poly11');  % plane fit
plot(sf, [x, y], z); axis equal;
xlabel('Pturn_{blue}'); ylabel('Pturn_{red}'); zlabel('Pturn_{blue+red}');
title(['Pturn_{blue + red} = ', num2str(sf.p00), ' + ', num2str(sf.p10), ' * Pturn_{blue} ', ' + ', num2str(sf.p01), ' * Pturn_{red} ']);
if denoise    
    savefig(gcf, '3d_plane_fit_denoise')
else
    savefig(gcf, '3d_plane_fit')
end


% histogram of pturn for different light with Gaussian fit
% determine the edges to plot
edge_min = min([Pblue; Pred; Pbluered]);
if fix(10 * edge_min) - 1 == 0
    edge_min = 0;  % 0.15 or 0.1 give 0
else
    edge_min = (fix(10 * edge_min) - 1) / 10;  % -0.354 or -0.3 give -0.4
end
edge_max = max([Pblue; Pred; Pbluered]);
if fix(10 * edge_max) + 1 == 0
    edge_max = 0;
else
    edge_max = (fix(10 * edge_max) + 1) / 10;
end
edges = edge_min : 0.1 : edge_max;
% get the counts
nblue = histcounts(Pblue, edges, 'Normalization', 'probability');
nred = histcounts(Pred, edges, 'Normalization', 'probability');
nbluered = histcounts(Pbluered, edges, 'Normalization', 'probability');
centers = edges(1 : end-1) + diff(edges) / 2; 
fblue = fit(centers.', nblue.', 'gauss1');  % .' is transpose, fit model for blue
fred = fit(centers.', nred.', 'gauss1');
fbluered = fit(centers.', nbluered.', 'gauss1');
% to plot with Gaussian fit
figure;
lblue = plot(fblue, centers, nblue, 'o');  % return the line objects, raw data and fitted line
lblue(1).Color = [0, 0, 1]; lblue(1).MarkerFaceColor = [0, 0, 1];  lblue(1).DisplayName = 'Blue'; % properties of dot data
lblue(2).Color = [0, 0, 1]; lblue(2).DisplayName = ['\mu = ', num2str(fblue.b1), ', \sigma = ', num2str(fblue.c1/sqrt(2))];  % properties of fit line
hold on;
lred = plot(fred, centers, nred, '*');  % return the line object
lred(1).Color = [1, 0, 0];  lred(1).DisplayName = 'Red'; % properties of dot data
lred(2).Color = [1, 0, 0]; lred(2).DisplayName = ['\mu = ', num2str(fred.b1), ', \sigma = ', num2str(fred.c1/sqrt(2))];  % properties of fit line
lbluered = plot(fbluered, centers, nbluered, 'o');  % return the line object
lbluered(1).Color = [0, 0, 0]; lbluered(1).DisplayName = 'Blue and red'; % properties of dot data
lbluered(2).Color = [0, 0, 0]; lbluered(2).DisplayName = ['\mu = ', num2str(fbluered.b1), ', \sigma = ', num2str(fbluered.c1/sqrt(2))];  % properties of fit line
hold off;
title('Fit with Gaussian distribution'); xlabel('Turn probability'); ylabel('Proportion of maggots');
if denoise
    savefig(gcf, 'hist_all_fit_denoise');
else
    savefig(gcf, 'hist_all_fit');
end

% to plot without fit
figure;
plot(centers, nblue, '-ob', centers, nred, '-or', centers, nbluered, '-ok');
xlabel('Turn probability'); ylabel('Proportion of maggots');
legend('Blue', 'Red', 'Blue and Red');
if denoise
    savefig(gcf, 'hist_all_denoise');
else
    savefig(gcf, 'hist_all');
end


figure;  % Superlinear if p(blue+red) > pblue + pred, vice verse sublinear
x = Pblue;  % p( 0<t<3s, blue)
y = Pred;
sz = 65;  % size of dots in pixels
% c = 1 - pturn_hist_all(:, 13) * [1 1 1];  % the larger p, the smaller c, the darker
c = (x + y) >= Pbluered;
scatter(x, y, sz, c, 'filled')
axis equal;
xlabel('Pturn during the first 3 s of blue'); ylabel('Pturn during the first 3 s of red');
cbar = colorbar;
colormap parula;
cbar.Label.String = 'If pblue + pred >= p(blue+red)';
cbar.Ticks = [0, 1];
cbar.TickLabels = {'Superlinear'; 'Sublinear'};
hold on; plot([0, 1], [1, 0]); legend('', 'x + y = 1'); hold off;
if denoise
    xlim([-0.2, 0.6]); ylim([-0.2, 0.6]); savefig(gcf, 'super_sub-lineaer_denoise');
else
    savefig(gcf, 'super_sub-lineaer');
end


figure;  % Linearity-pblue_pred
x = Pblue;  % p( 0<t<3s, blue)
y = Pred;
sz = 65;  % size of dots in pixels
% c = 1 - pturn_hist_all(:, 13) * [1 1 1];  % the larger p, the smaller c, the darker
c = x + y - Pbluered;
scatter(x, y, sz, c, 'filled')
axis equal;
xlabel('Pturn during the first 3 s of blue'); ylabel('Pturn during the first 3 s of red');
cbar = colorbar;
colormap parula;
cbar.Label.String = 'pblue + pred - p(blue+red)';
cbar.Limits = [-0.6, 0.6];
hold on; plot([0, 1], [1, 0]); legend('', 'x + y = 1'); hold off;
if denoise
    xlim([-0.2, 0.6]); ylim([-0.2, 0.6]); savefig(gcf, 'linearity_denoise');
else
    savefig(gcf, 'linearity');
end
