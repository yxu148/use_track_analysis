%% create pturn matrix from the structure larvae
basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\variability_new_extracted\Gr21a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3#T_Bl_Sq_2to7P_15_1_3'
    };
x_cell = {
    1:21
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

% create pturn_all in order of categories and save in data.mat
% pturn_all = [];  % nlarva-by-15 array, 3s bins
pturn_000 = []; pturn_001 = []; pturn_010 = []; pturn_011 = [];
pturn_100 = []; pturn_101 = []; pturn_110 = []; pturn_111 = [];
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop for each larva in larvae
        pturn = getfield(larvae, ['larva', num2str(i)], 'pturn');
        % pturn_all = [pturn_all; pturn.blue, pturn.red pturn.bluered];
        % to load pturn in order of response
        response = [getfield(larvae, ['larva', num2str(i)], 'response', 'blue'), getfield(larvae, ['larva', num2str(i)], 'response', 'red'), getfield(larvae, ['larva', num2str(i)], 'response', 'bluered')];
        switch response
            case '000'
                pturn_000 = [pturn_000; pturn.blue, pturn.red pturn.bluered];
            case '001'
                pturn_001 = [pturn_001; pturn.blue, pturn.red pturn.bluered];
            case '010'
                pturn_010 = [pturn_010; pturn.blue, pturn.red pturn.bluered];
            case '011'
                pturn_011 = [pturn_011; pturn.blue, pturn.red pturn.bluered];
            case '100'
                pturn_100 = [pturn_100; pturn.blue, pturn.red pturn.bluered];
            case '101'
                pturn_101 = [pturn_101; pturn.blue, pturn.red pturn.bluered];
            case '110'
                pturn_110 = [pturn_110; pturn.blue, pturn.red pturn.bluered];
            case '111'
                pturn_111 = [pturn_111; pturn.blue, pturn.red pturn.bluered];
        end  % end switch
    end  % end looping each larva in one experiment
end  % end loop each experiment
pturn_all = [pturn_000; pturn_001; pturn_010; pturn_011; pturn_100; pturn_101; pturn_110; pturn_111];

savename = strcat(pwd, '\data.mat');
if isfile(savename)
    save(savename, 'pturn_all', 'pturn_000', 'pturn_001', 'pturn_010', 'pturn_011', 'pturn_100', 'pturn_101', 'pturn_110', 'pturn_111', '-append');
else
    save(savename, 'pturn_all', 'pturn_000', 'pturn_001', 'pturn_010', 'pturn_011', 'pturn_100', 'pturn_101', 'pturn_110', 'pturn_111')
end
% writematrix(pturn_all, strcat(pwd, '\pturn.xlsx'));


% % remove certain fields from a structure
% for j = 1 : length(figlocation_list)  % loop for each expt
%     savename = strcat(figlocation_list{j}, '\data.mat');
%     load(savename, 'larvae');
%     rm_list = [];  % the larva index to be removed
%     for i = 1 : length(fieldnames(larvae))  % loop for each larva in larvae
%         if length(fieldnames(larvae.(['larva', num2str(i)]))) == 1  % select the larva index who contain only 1 field
%             rm_list = [rm_list, i];
%         end
%     end
%     for k = 1 : length(rm_list)
%         larvae = rmfield(larvae, ['larva', num2str(rm_list(k))]);  % remove the certain larva indexes
%     end
%     save(savename, 'larvae', '-append')
%     disp(['Larva ', num2str(rm_list), ' is removed from expt ', num2str(j)]);
% end

%% Plot pturn from data.mat ===================================================
load('data.mat');



denoise =false;
% Initialize some values for plotting below
pturn = pturn_111;  % what pturn to use, could be pturn_111 or pturn_all, etc.
prefix = '111_';  % use '001_' or 'all_' etc
if denoise  % the p of first bin minuses the average of rest bins
    Pbaseline = (pturn(:, 5) + pturn(:, 10) + pturn(:, 15)) / 3;
    Pblue = pturn(:,1) - Pbaseline;
    Pred = pturn(:,6) - Pbaseline;
    Pbluered = pturn(:,11) - Pbaseline;
else
    Pblue = pturn(:,1);
    Pred = pturn(:, 6);
    Pbluered = pturn(:, 11);
end


% plot the matrix pturn_all
imagesc(pturn_all, [0, 1]); 
xlabel('3-second bins'); ylabel('Index of larva');
% number of larva for each type, in order of 000, 001, 010, 011, 100, 101, 110, 111
nlarva_type = [size(pturn_000, 1), size(pturn_001, 1), size(pturn_010, 1), size(pturn_011, 1),size(pturn_100, 1), size(pturn_101, 1), size(pturn_110, 1), size(pturn_111, 1)];
yline([nlarva_type(1), sum(nlarva_type(1:2)), sum(nlarva_type(1:3)), sum(nlarva_type(1:4)), sum(nlarva_type(1:5)), sum(nlarva_type(1:6)), ...
    sum(nlarva_type(1:7))] + 0.5, 'w-');  % draw white line to seperate different response type
title(['1-', num2str(nlarva_type(1)), ' larva-000, ', ...
    num2str(nlarva_type(1) + 1), '-', num2str(sum(nlarva_type(1:2))), ' larva-001, ', ...,
    num2str(sum(nlarva_type(1:2)) + 1), '-', num2str(sum(nlarva_type(1:3))), ' larva-010, ', ...,
    num2str(sum(nlarva_type(1:3)) + 1), '-', num2str(sum(nlarva_type(1:4))), ' larva-011, ', ...,
    num2str(sum(nlarva_type(1:4)) + 1), '-', num2str(sum(nlarva_type(1:5))), ' larva-100, ', ...,
    num2str(sum(nlarva_type(1:5)) + 1), '-', num2str(sum(nlarva_type(1:6))), ' larva-101, ', ...,
    num2str(sum(nlarva_type(1:6)) + 1), '-', num2str(sum(nlarva_type(1:7))), ' larva-110, ', ...,
    num2str(sum(nlarva_type(1:7)) + 1), '-', num2str(sum(nlarva_type(1:8))), ' larva-111, '
    ])
cbar = colorbar; cbar.Label.String = 'Turn possibility with the 3-s bin';
savename = strcat(pwd, '\results_fig', '\pturn_all');
savefig(gcf, savename);



figure;  % p( 0<t<3s, blue+red) as a function of p( 0<t<3s, blue) and p( 0<t<3s, red)
x = Pblue;  % p( 0<t<3s, blue)
y = Pred;
sz = 40;  % size of dots in pixels
% c = 1 - pturn_hist_all(:, 13) * [1 1 1];  % the larger p, the smaller c, the darker
c = Pbluered;
scatter(x, y, sz, c, 'filled')
axis equal;
xlabel('Pturn during the first 3 s of blue'); ylabel('Pturn during the first 3 s of red');
cbar = colorbar;
colormap parula;
cbar.Label.String = 'Pturn during the first 3 s of both blue and red';
if denoise
    xlim([-0.2, 0.6]); ylim([-0.2, 0.6]); savefig(gcf, [prefix, 'pbluered-pb_pr_denoise'])
else
    xlim([0, 0.85]); ylim([0, 0.85]); savefig(gcf, [prefix, 'pbluered-pb_pr'])
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
    savefig(gcf, [prefix, 'pbluered--pb+pr_denoise'])
else
    savefig(gcf, [prefix, 'pbluered--pb+pr'])
end


figure;  % 3-d plane fit
x = Pblue;
y = Pred;
z = Pbluered;
[sf, goodness_of_fit, fitting_algorithm] = fit([x, y], z, 'poly11');  % plane fit
plot(sf, [x, y], z); axis equal;
xlabel('Pturn_{blue}'); ylabel('Pturn_{red}'); zlabel('Pturn_{blue+red}');
title(['Pturn_{blue + red} = ', num2str(sf.p00), ' + ', num2str(sf.p10), ' * Pturn_{blue} ', ' + ', num2str(sf.p01), ' * Pturn_{red} ']);
if denoise    
    savefig(gcf, [prefix, '3d_plane_fit_denoise'])
else
    savefig(gcf, [prefix, '3d_plane_fit'])
end


download = true;
for denoise = [true, false]
    % Initialize some values for plotting below
    pturns = {pturn_all, pturn_001, pturn_011, pturn_101, pturn_111};
    prefixes = {'all_', '001_', '011_', '101_', '111_'};
    for i = 1 : length(pturns)
        pturn = pturns{i};  % what pturn to use, could be pturn_111 or pturn_all, etc.
        prefix = prefixes{i};  % use '001_' or 'all_' etc
        if denoise  % the p of first bin minuses the average of rest bins
            Pbaseline = (pturn(:, 5) + pturn(:, 10) + pturn(:, 15)) / 3;
            Pblue = pturn(:,1) - Pbaseline;
            Pred = pturn(:,6) - Pbaseline;
            Pbluered = pturn(:,11) - Pbaseline;
        else
            Pblue = pturn(:,1);
            Pred = pturn(:, 6);
            Pbluered = pturn(:, 11);
        end
        
        % histogram of pturn for different light with different fit
        % % use the same size of bins for each histogram
        fit_func = 'Gaussian';  % 'Gaussian', or 'Poisson'
        edge_min = min([Pblue; Pred; Pbluered]);
        edge_min = floor(edge_min * 10) / 10;  % find it's left boundary with 1 decimal
        edge_max = max([Pblue; Pred; Pbluered]);
        edge_max = ceil(edge_max * 10) / 10;  % find it's right boundary with 1 decimal
        edges = edge_min : 0.05 : edge_max;
        % get the counts
        nblue = histcounts(Pblue, edges, 'Normalization', 'probability');
        nred = histcounts(Pred, edges, 'Normalization', 'probability');
        nbluered = histcounts(Pbluered, edges, 'Normalization', 'probability');
        centers = edges(1 : end-1) + diff(edges) / 2; 
        centers_blue = centers; centers_red = centers; centers_bluered = centers;
        % % use the same number of bins (nbins) for each histogram but different size of bin, not good
        % nbins = 6;
        % [nblue, edges_blue] = histcounts(Pblue, nbins, 'Normalization', 'probability');
        % [nred, edges_red] = histcounts(Pred, nbins, 'Normalization', 'probability');
        % [nbluered, edges_bluered] = histcounts(Pbluered, nbins, 'Normalization', 'probability');
        % centers_blue = edges_blue(1 : end-1) + diff(edges_blue) / 2; 
        % centers_red = edges_red(1 : end-1) + diff(edges_red) / 2; 
        % centers_bluered = edges_bluered(1 : end-1) + diff(edges_bluered) / 2; 
        switch fit_func 
            case 'Gaussian'
                fblue = fit(centers_blue(nblue~=0).', nblue(nblue~=0).', 'gauss1');  % .' is transpose, fit model for blue
                fred = fit(centers_red(nred~=0).', nred(nred~=0).', 'gauss1');
                fbluered = fit(centers_bluered(nbluered~=0).', nbluered(nbluered~=0).', 'gauss1');
                % to plot with Gaussian fit
                figure; hold on;
                lblue = plot(fblue, centers_blue(nblue~=0), nblue(nblue~=0), 'o');  % return the line objects, raw data and fitted line
                bblue = bar(centers_blue(nblue~=0), nblue(nblue~=0), 1, 'FaceColor', [0, 0, 1]);
                lblue(1).Color = [0, 0, 1]; lblue(1).MarkerFaceColor = [0, 0, 1]; lblue(1).DisplayName = 'Blue'; % properties of dot data
                lblue(2).Color = [0, 0, 1]; lblue(2).DisplayName = ['\mu = ', num2str(fblue.b1), ', \sigma = ', num2str(fblue.c1/sqrt(2))];  % properties of fit line
                bblue.FaceAlpha = 0.5; bblue.DisplayName = 'Blue'; bblue.BarWidth = 1;
        
                lred = plot(fred, centers_red(nred~=0), nred(nred~=0), 'o');  % return the line object
                bred = bar(centers_red(nred~=0), nred(nred~=0), 1, 'FaceColor', [1, 0, 0]);
                lred(1).Color = [1, 0, 0];  lred(1).MarkerFaceColor = [1, 0, 0]; lred(1).DisplayName = 'Red'; % properties of dot data
                lred(2).Color = [1, 0, 0]; lred(2).DisplayName = ['\mu = ', num2str(fred.b1), ', \sigma = ', num2str(fred.c1/sqrt(2))];  % properties of fit line
                bred.FaceAlpha = 0.5; bred.DisplayName = 'Red'; bred.BarWidth = 0.8;
        
                lbluered = plot(fbluered, centers_bluered(nbluered~=0), nbluered(nbluered~=0), 'o');  % return the line object
                bbluered = bar(centers_bluered(nbluered~=0), nbluered(nbluered~=0), 1, 'FaceColor', [0, 0, 0]);
                lbluered(1).Color = [0, 0, 0]; lbluered(1).MarkerFaceColor = [0, 0, 0]; lbluered(1).DisplayName = 'Blue and red'; % properties of dot data
                lbluered(2).Color = [0, 0, 0]; lbluered(2).DisplayName = ['\mu = ', num2str(fbluered.b1), ', \sigma = ', num2str(fbluered.c1/sqrt(2))];  % properties of fit line
                bbluered.FaceAlpha = 0.5;  bbluered.DisplayName = 'Blue and red'; bbluered.BarWidth = 0.6;
                hold off; xticks(edges); 
                legend('Location', 'eastoutside');
                xlabel('Pturn during the first 3 s of stimulation high'); ylabel('Proportion of maggots');
            case 'Poisson'
                edges_continuous = 0 : 40;
                centers_continuous = edges_continuous(1 : end-1) + diff(edges_continuous) / 2; 
                lambda_blue = poissfit(40*centers_blue(nblue~=0).', nblue(nblue~=0).');
                lambda_red = poissfit(40*centers_red(nred~=0).', nred(nred~=0).');
                lambda_bluered = poissfit(40*centers_bluered(nbluered~=0).', nbluered(nbluered~=0).');
                figure;
                dblue = plot(40*centers_blue(nblue~=0), nblue(nblue~=0), 'o'); hold on;
                dblue.Color = [0, 0, 1]; dblue.MarkerFaceColor = [0, 0, 1]; dblue.DisplayName = 'Blue'; % properties of dot data
                lblue = plot(edges_continuous, poisspdf(edges_continuous, lambda_blue), '-');
                lblue.Color = [0, 0, 1]; lblue.DisplayName = ['\lambda_{blue} = ', num2str(lambda_blue)];
                dred = plot(40*centers_red(nred~=0), nred(nred~=0), '*');
                dred.Color = [1, 0, 0]; dred.DisplayName = 'Red'; % properties of dot data
                lred = plot(edges_continuous, poisspdf(edges_continuous, lambda_red), '-');
                lred.Color = [1, 0, 0]; lred.DisplayName = ['\lambda_{red} = ', num2str(lambda_red)];
                dbluered = plot(40*centers_bluered(nbluered~=0), nbluered(nbluered~=0), 'o');
                dbluered.Color = [0, 0, 0]; dbluered.DisplayName = 'Blue and red'; % properties of dot data
                lbluered = plot(edges_continuous, poisspdf(edges_continuous, lambda_bluered), '-');
                lbluered.Color = [0, 0, 0]; lbluered.DisplayName = ['\lambda_{bluered} = ', num2str(lambda_bluered)];
                legend('Location', 'eastoutside'); hold off
                xlabel('Average Nturns for a larva during the first 3 s of stimulation high'); ylabel('Proportion of maggots');
        end
        title([prefix(1:end-1), ' ', fit_func]);
        if denoise
            title([prefix(1:end-1), ' ', fit_func, ' Denoise']);
            if download
                savename = strcat(pwd, '\results_fig', '\', [prefix, 'hist_pturn_stim_bin0p05_denoise_', fit_func]);
                savefig(gcf, savename);
            end
        else
            if download
                savename = strcat(pwd, '\results_fig', '\', [prefix, 'hist_pturn_stim_bin0p05_', fit_func]);
                savefig(gcf, savename);
            end
        end

    end  % end loop for each stimulation type
end  % end loop for denoise


% histogram of pturn for different light with polynomial fit
% determine the edges to plot
edge_min = min([Pblue; Pred; Pbluered]);
edge_min = floor(edge_min * 10) / 10;  % find it's left boundary with 1 decimal
edge_max = max([Pblue; Pred; Pbluered]);
edge_max = ceil(edge_max * 10) / 10;  % find it's right boundary with 1 decimal
edges = edge_min : 0.1 : edge_max;
% get the counts
nblue = histcounts(Pblue, edges, 'Normalization', 'probability');
nred = histcounts(Pred, edges, 'Normalization', 'probability');
nbluered = histcounts(Pbluered, edges, 'Normalization', 'probability');
centers = edges(1 : end-1) + diff(edges) / 2; 
% only fit with non-zero points
fblue = fit(centers.', nblue.', 'poly3');  % .' is transpose, fit model for blue
fred = fit(centers.', nred.', 'poly5');
fbluered = fit(centers.', nbluered.', 'poly7');
% to plot with fit
figure;
lblue = plot(fblue, centers, nblue, 'o');  % return the line objects, raw data and fitted line
lblue(1).Color = [0, 0, 1]; lblue(1).MarkerFaceColor = [0, 0, 1];  lblue(1).DisplayName = 'Blue'; % properties of dot data
lblue(2).Color = [0, 0, 1]; % lblue(2).DisplayName = ['\mu = ', num2str(fblue.b1), ', \sigma = ', num2str(fblue.c1/sqrt(2))];  % properties of fit line
hold on;
lred = plot(fred, centers, nred, '*');  % return the line object
lred(1).Color = [1, 0, 0];  lred(1).DisplayName = 'Red'; % properties of dot data
lred(2).Color = [1, 0, 0]; % lred(2).DisplayName = ['\mu = ', num2str(fred.b1), ', \sigma = ', num2str(fred.c1/sqrt(2))];  % properties of fit line
lbluered = plot(fbluered, centers, nbluered, 'o');  % return the line object
lbluered(1).Color = [0, 0, 0]; lbluered(1).DisplayName = 'Blue and red'; % properties of dot data
lbluered(2).Color = [0, 0, 0]; % lbluered(2).DisplayName = ['\mu = ', num2str(fbluered.b1), ', \sigma = ', num2str(fbluered.c1/sqrt(2))];  % properties of fit line
hold off;
title([prefix, 'Fit with polynoimal']); xlabel('Pturn during the first 3 s of stimulation high'); ylabel('Proportion of maggots');
if denoise
    savefig(gcf, [prefix, 'hist_pturn_stim_gaussian_denoise']);
else
    savefig(gcf, [prefix, 'hist_pturn_stim_gaussian_fit']);
end


figure;  % Superlinear if p(blue+red) > pblue + pred, vice verse sublinear
x = Pblue;  % p( 0<t<3s, blue)
y = Pred;
sz = 30;  % size of dots in pixels
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
    xlim([-0.2, 0.6]); ylim([-0.2, 0.6]); savefig(gcf, [prefix, 'super_sub-lineaer_denoise']);
else
    xlim([0, 1]); ylim([0, 1]); savefig(gcf, [prefix, 'super_sub-lineaer']);
end


figure;  % Linearity-pblue_pred
x = Pblue;  % p( 0<t<3s, blue)
y = Pred;
sz = 40;  % size of dots in pixels
% c = 1 - pturn_hist_all(:, 13) * [1 1 1];  % the larger p, the smaller c, the darker
c = x + y - Pbluered;
scatter(x, y, sz, c, 'filled')
axis equal;
xlabel('Pturn during the first 3 s of blue'); ylabel('Pturn during the first 3 s of red');
cbar = colorbar;
colormap turbo;
cbar.Label.String = 'pblue + pred - p(blue+red)';
cbar.Limits = [-0.6, 0.6];
hold on; plot([0, 1], [1, 0]); legend('', 'x + y = 1'); hold off;
if denoise
    xlim([-0.2, 0.6]); ylim([-0.2, 0.6]); savefig(gcf, [prefix, 'linearity_denoise']);
else
    xlim([0, 1]); ylim([0, 1]); savefig(gcf, [prefix, 'linearity']);
end



% Kullback–Leibler divergence
kld_blue = zeros(size(pturn_all, 1), 1);
kld_red = zeros(size(pturn_all, 1), 1);
kld_blue_terms = zeros(size(pturn_all, 1), 5);
kld_red_terms = zeros(size(pturn_all, 1), 5);
for j = 1 : size(pturn_all, 1)  % the j-th row (larva) of the pturn
    P_blue = pturn_all(j, 1:5);  % + 0.025 to eliminate 0 probability
    P_red = pturn_all(j, 6:10);  % 1/40, add 1 turn to every bin
    Q = pturn_all(j, 11:15);
    [kld_blue(j), kld_blue_terms(j, :)] = KL_divergence(Q, P_blue);
    [kld_red(j), kld_red_terms(j, :)] = KL_divergence(Q, P_red);
end
% number of larva for each type, in order of 000, 001, 010, 011, 100, 101, 110, 111
nlarva_type = [size(pturn_000, 1), size(pturn_001, 1), size(pturn_010, 1), size(pturn_011, 1),size(pturn_100, 1), size(pturn_101, 1), size(pturn_110, 1), size(pturn_111, 1)];
boundary_type = 0.5 + [nlarva_type(1), sum(nlarva_type(1:2)), sum(nlarva_type(1:3)), sum(nlarva_type(1:4)), sum(nlarva_type(1:5)), sum(nlarva_type(1:6)), sum(nlarva_type(1:7))];
response_type = {'000', '001', '010', '011', '100', '101', '110', '111'};

figure;  % all 3s bins
ax(1) = subplot(2, 1, 1);
p1 = plot(kld_blue, 'b-');
xline(boundary_type, 'k--');
text([0.5 boundary_type] + nlarva_type/2, 0.1 + max(kld_blue)*ones(1, length(response_type)), response_type);
xlabel('Larva Index'); ylabel('D_{KL}(Pbluered||Pblue)');
ax(2) = subplot(2, 1, 2);
p2 = plot(kld_red, 'r-'); 
xlabel('Larva Index'); ylabel('D_{KL}(Pbluered||Pred)')
xline(boundary_type, 'k--');  % draw black dash line to seperate different response type
text([0.5 boundary_type] + nlarva_type/2, 0.1 + max(kld_red)*ones(1, length(response_type)), response_type);
% legend([p1, p2], {'D_{KL}(Pblue||Pbluered)', 'D_{KL}(Pred||Pbluered)'});
ylim([ax(1), ax(2)], [-1, 4]);
sgtitle('Raw pturn');
savename = strcat(pwd, '\results_fig', '\Dkl_larvaIndex_reverse');
savefig(gcf, savename);

figure;  % a 3s bin
ax(1) = subplot(2, 3, [1, 2]);
p1 = plot(kld_blue_terms(:, 1), 'b-');
xline([nlarva_type(1), sum(nlarva_type(1:2)), sum(nlarva_type(1:3)), sum(nlarva_type(1:4)), sum(nlarva_type(1:5)), sum(nlarva_type(1:6)), ...
    sum(nlarva_type(1:7))] + 0.5, 'k--');
xlabel('Larva Index in order of 000, 001, 010, 011, 100, 101, 110, 111'); ylabel('D_{KL}(Pblue||Pbluered)');
ax(2) = subplot(2, 3, 3);
histogram(kld_blue_terms(:, 1), 'Orientation', 'horizontal'); xlabel('Count of larvae');
ax(3) = subplot(2, 3, [4, 5]);
p2 = plot(kld_red_terms(:, 1), 'r-'); 
xlabel('Larva Index in order of 000, 001, 010, 011, 100, 101, 110, 111'); ylabel('D_{KL}(Pred||Pbluered)')
xline([nlarva_type(1), sum(nlarva_type(1:2)), sum(nlarva_type(1:3)), sum(nlarva_type(1:4)), sum(nlarva_type(1:5)), sum(nlarva_type(1:6)), ...
    sum(nlarva_type(1:7))] + 0.5, 'k--');  % draw black dash line to seperate different response type
ax(4) = subplot(2, 3, 6);
histogram(kld_red_terms(:, 1), 'Orientation', 'horizontal'); xlabel('Count of larvae');
% legend([p1, p2], {'D_{KL}(Pblue||Pbluered) of 2nd bin', 'D_{KL}(Pred||Pbluered) of 2nd bin'});
ylim([ax(1), ax(3), ax(2), ax(4)], [-1, 4]);
sgtitle('Add 1 turn to each bin');
savename = strcat(pwd, '\results_fig', '\Dkl_larvaIndex_add1turn_1stbin_reverse');
savefig(gcf, savename);



% KL divergence test, plot pturn with DKL title
index_larva = [16];
edges = 0:3:15;
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
for j = 1:length(index_larva)  % plot the j-th larva
    P_blue = pturn_all(index_larva(j), 1:5);  % 1-by-5 array
    P_red = pturn_all(index_larva(j), 6:10);
    Q = pturn_all(index_larva(j), 11:15);
    [kld_blue, kld_blue_terms] = KL_divergence(Q, P_blue);
    [kld_red, kld_red_terms] = KL_divergence(Q, P_red);
    subplot(3, length(index_larva), j)
    bar(xbar, P_blue);
    title([num2str(kld_blue,'%.2f '), ' = sum(', num2str(kld_blue_terms,'%.2f '), ')'])
    xticks(edges); ylim([0, 1]); xline(6, '--');
    pos = get(gca, 'Position');
    pos(3) = pos(4) * 1.5;  % width = height * 1.5
    set(gca, 'Position', pos)
    subplot(3, length(index_larva), j+length(index_larva))
    bar(xbar, P_red); ylim([0, 1]);
    title([num2str(kld_red,'%.2f '), ' = sum(', num2str(kld_red_terms,'%.2f '), ')'])
    xticks(edges); ylim([0, 1]); xline(6, '--');
    pos = get(gca, 'Position');
    pos(3) = pos(4) * 1.5;  % width = height * 1.5
    set(gca, 'Position', pos)
    subplot(3, length(index_larva), j+length(index_larva)*2)
    bar(xbar, Q);
    title(['larva ', num2str(index_larva(j))])
    xticks(edges); ylim([0, 1]); xline(6, '--');
    pos = get(gca, 'Position');
    pos(3) = pos(4) * 1.5;  % width = height * 1.5
    set(gca, 'Position', pos)
end


% to find the larvae with DKL blue/red very close but belong to different
% response type. By looking at DKLblue-DKLred plot, I have found there are 4 different larvae fall into 
% DKLblue (0.59, 0.64), and DKLred(0.51, 0.58)
index_larva = find((kld_blue > 0.59) & (kld_blue < 0.64) & (kld_red > 0.51) & (kld_red < 0.58));
% pturn test, by plotting certain pturn column with DKL info
index = find((Pblue >0.2) & (Pblue < 0.3));  % pick the index of larva in pturn (not pturn_all) whose pblue(1) < 0.2
index_larva = index(1:min([4, length(index)]));  % only select the first 4 of them
edges = 0:3:15;
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
for j = 1:length(index_larva)  % plot the j-th larva
    P_blue = pturn(index_larva(j), 1:5) + 0.025;  % 1-by-5 array
    P_red = pturn(index_larva(j), 6:10) + 0.025;
    Q = pturn(index_larva(j), 11:15) + 0.025;
    [kld_blue, kld_blue_terms] = KL_divergence(Q, P_blue);
    [kld_red, kld_red_terms] = KL_divergence(Q, P_red);
    subplot(3, length(index_larva), j)
    bar(xbar, P_blue); ylim([0, 1]);
    title([num2str(kld_blue,'%.2f '), ' = sum(', num2str(kld_blue_terms,'%.2f '), ')'])
    subplot(3, length(index_larva), j+length(index_larva))
    bar(xbar, P_red); ylim([0, 1]);
    title([num2str(kld_red,'%.2f '), ' = sum(', num2str(kld_red_terms,'%.2f '), ')'])
    subplot(3, length(index_larva), j+length(index_larva)*2)
    bar(xbar, Q); ylim([0, 1]);
    index_in_pturn_all = find(ismember(pturn_all, pturn(index_larva(j), :), 'rows' ));
    title(['larva ', num2str(index_in_pturn_all)])
end


% DKL_blue vs DKL red
% number of larva for each type, in order of 000, 001, 010, 011, 100, 101, 110, 111
nlarva_type = [size(pturn_000, 1), size(pturn_001, 1), size(pturn_010, 1), size(pturn_011, 1),size(pturn_100, 1), size(pturn_101, 1), size(pturn_110, 1), size(pturn_111, 1)];
boundary_type = 0.5 + [nlarva_type(1), sum(nlarva_type(1:2)), sum(nlarva_type(1:3)), sum(nlarva_type(1:4)), sum(nlarva_type(1:5)), sum(nlarva_type(1:6)), sum(nlarva_type(1:7))];
response_type = {'000', '001', '010', '011', '100', '101', '110', '111'};
response_marker = {'.', 'o', '.', '|', '.', '_', '.', '+'};
response_color = {'c', 'g', 'c', 'r', 'c', 'b', 'c', 'k'};
pturns = {pturn_000, pturn_001, pturn_010, pturn_011, pturn_100, pturn_101, pturn_110, pturn_111};
figure; hold on;
for m = 1 : length(pturns)  % m types of larva
    pturn = pturns{m};
    kld_blue = zeros(size(pturn, 1), 1);
    kld_red = zeros(size(pturn, 1), 1);
    kld_blue_terms = zeros(size(pturn, 1), 5);
    kld_red_terms = zeros(size(pturn, 1), 5);
    for j = 1 : size(pturn, 1)  % the j-th row (larva) of the pturn
        P_blue = pturn(j, 1:5) + 0.025;  % + 0.025 to eliminate 0 probability
        P_red = pturn(j, 6:10) + 0.025;  % 1/40, add 1 turn to every bin
        Q = pturn(j, 11:15) + 0.025;
        [kld_blue(j), kld_blue_terms(j, :)] = KL_divergence(Q, P_blue);
        [kld_red(j), kld_red_terms(j, :)] = KL_divergence(Q, P_red);
    end
    ax(m) = plot(kld_blue, kld_red, [response_marker{m}, response_color{m}]);
end  % end looping different types of larvae
hold off;
xlabel('DKL blue'); ylabel('DKL red');
title('pturn = pturn + 0.025, D_{KL}(single||both color)');
legend(ax([2, 4, 6, 8, 1]), {response_type{[2, 4, 6, 8]}, 'Others'}, 'Location', 'best');
savename = strcat(pwd, '\results_fig', '\DKLblue-DKLred_add1turn');
savefig(gcf, savename);


% DKL_blue histogram of each response type
% number of larva for each type, in order of 000, 001, 010, 011, 100, 101, 110, 111
nlarva_type = [size(pturn_000, 1), size(pturn_001, 1), size(pturn_010, 1), size(pturn_011, 1),size(pturn_100, 1), size(pturn_101, 1), size(pturn_110, 1), size(pturn_111, 1)];
boundary_type = 0.5 + [nlarva_type(1), sum(nlarva_type(1:2)), sum(nlarva_type(1:3)), sum(nlarva_type(1:4)), sum(nlarva_type(1:5)), sum(nlarva_type(1:6)), sum(nlarva_type(1:7))];
response_type = {'000', '001', '010', '011', '100', '101', '110', '111'};
response_type_major = {'001', '011', '101', '111'};
response_marker = {'.', 'o', '.', '|', '.', '_', '.', '+'};
response_color = {'c', 'g', 'c', 'r', 'c', 'b', 'c', 'k'};
pturns = {pturn_000, pturn_001, pturn_010, pturn_011, pturn_100, pturn_101, pturn_110, pturn_111};
pturns_major = {pturn_001, pturn_011, pturn_101, pturn_111};
response_color_major = {'g', 'r', 'b', 'k'};
figure; hold on;
for m = 1 : length(pturns_major)  % m types of larva
    pturn = pturns_major{m};
    kld_blue = zeros(size(pturn, 1), 1);
    kld_red = zeros(size(pturn, 1), 1);
    kld_blue_terms = zeros(size(pturn, 1), 5);
    kld_red_terms = zeros(size(pturn, 1), 5);
    for j = 1 : size(pturn, 1)  % the j-th row (larva) of the pturn
        P_blue = pturn(j, 1:5) + 0.025;  % + 0.025 to eliminate 0 probability
        P_red = pturn(j, 6:10) + 0.025;  % 1/40, add 1 turn to every bin
        Q = pturn(j, 11:15) + 0.025;
        [kld_blue(j), kld_blue_terms(j, :)] = KL_divergence(P_blue, Q);
        [kld_red(j), kld_red_terms(j, :)] = KL_divergence(P_red, Q);
    end
    ax(m) = histogram(kld_blue, 'FaceColor', response_color_major{m}, 'FaceAlpha', 0.2);
end  % end looping different types of larvae
hold off;
xlabel('DKL blue'); ylabel('Count');
title('pturn = pturn + 0.025, D_{KL}(single||both color)');
legend(ax([1, 2, 3, 4]), response_type_major{[1, 2, 3, 4]}, 'Location', 'best');
savename = strcat(pwd, '\results_fig', '\DKLblue_hist_add1turn');
savefig(gcf, savename);
