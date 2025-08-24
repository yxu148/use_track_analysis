% information from pturn, including pturn, hist(pturn), DKL

%% create pturn matrix from the structure larvae
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

% create pturn_all in order of categories and save in data.mat
% pturn_all = [];  % nlarva-by-15 array, 3s bins
stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
ntype = 2 ^ (length(stim_color) - 1);  % int, -1 in the power if exclude dark in the end ------------------
pturn_type = cell(1, ntype);  % 1-by-ntype cell, 1-indexed, each element is an array
v_init_type = cell(1, ntype);
trackNum_type = cell(1, ntype);  % the info of each type is saved in one cell of the *_type cell
expt_type = cell(1, ntype);
larva_type = cell(1, ntype);  % the larva index in expt
response_type = cell(1, ntype);
for j = 1 : length(figlocation_list)
    savename = strcat(figlocation_list{j}, '\data.mat');
    load(savename, 'larvae');
    for i = 1 : length(fieldnames(larvae))  % loop for each larva in larvae
        if getfield(larvae, ['larva', num2str(i)], 'valid')
%         if 1
            pturn = getfield(larvae, ['larva', num2str(i)], 'pturn');
            v_init = getfield(larvae, ['larva', num2str(i)], 'init_speed');
            % pturn_all = [pturn_all; pturn.blue, pturn.red pturn.bluered];
            % to load pturn in order of response
            larva_index = ['larva', num2str(i)];
            response = '';  % e.g. '001001', initialize for every larva
            for s = 1 : length(stim_color) - 1  % -1 if exclude the last stim_color, 'dark' -------------
                response = [response, getfield(larvae, ['larva', num2str(i)], 'response', stim_color{s})];
            end
            typeth = bin2dec(response) + 1;  % positive int
            pturn_temp = [];
            v_init_temp = [];
            for s = 1 : length(stim_color)  % put pturn together, including 'dark' if there is
                pturn_temp = [pturn_temp, getfield(pturn, stim_color{s})];
            end
            for s = 1 : length(stim_color)  % -1 to exclude the last stim_color, 'dark'
                v_init_temp = [v_init_temp, getfield(v_init, stim_color{s})];
            end
            pturn_type{typeth} = [pturn_type{typeth};  pturn_temp];
            v_init_type{typeth} = [v_init_type{typeth};  v_init_temp];
            trackNum_type{typeth} = [trackNum_type{typeth}; getfield(larvae, ['larva', num2str(i)], 'trackNum')];
            expt_type{typeth} = [expt_type{typeth}; getfield(larvae, ['larva', num2str(i)], 'expt_info')];
            larva_type{typeth} = [larva_type{typeth}; i];
            response_type{typeth} = [response_type{typeth}; response];

        end  % end if valid larvae
    end  % end looping each larva in one experiment
end  % end loop each experiment
% pturn_all = [pturn_000; pturn_001; pturn_010; pturn_011; pturn_100; pturn_101; pturn_110; pturn_111];
% trackNum_all = [trackNum_000; trackNum_001; trackNum_010; trackNum_011; trackNum_100; trackNum_101; trackNum_110; trackNum_111];
% expt_all = [expt_000; expt_001; expt_010; expt_011; expt_100; expt_101; expt_110; expt_111];
% larva_all = [larva_000; larva_001; larva_010; larva_011; larva_100; larva_101; larva_110; larva_111];
% 

% put array from different cell together to one array
[pturn_all, v_init_all, trackNum_all, expt_all, larva_all, response_all] = deal([]);  % all variables have the same defination
for s = 1 : length(pturn_type)
    pturn_all = [pturn_all; pturn_type{s}];
    v_init_all = [v_init_all; v_init_type{s}];
    trackNum_all = [trackNum_all; trackNum_type{s}];
    expt_all = [expt_all; expt_type{s}];
    larva_all = [larva_all; larva_type{s}];
    response_all = [response_all; response_type{s}];
end

larva2track = table([1 : length(pturn_all)].', larva_all, trackNum_all, response_all, expt_all);

savename = strcat(basedir_cell{1}, '\results', '\data.mat');  % there should be only 1 basedir in basedir_cell
if isfile(savename)
    save(savename, 'pturn_all', 'pturn_type', '-append');
    save(savename, 'larva2track', '-append')
    save(savename, 'v_init_all', '-append')
else
    save(savename, 'pturn_all', 'pturn_type')
    save(savename, 'larva2track')
end
writematrix(pturn_all, strcat(basedir_cell{1}, '\results', '\pturn+dark.xlsx'));
writetable(larva2track, strcat(basedir_cell{1}, '\results', '\larva2track.xlsx'))
% writematrix(pturn_all, 'pturn.csv');

% proportion-lize pturn_all, so that each 5-size pturn is a probability distribution (sum is 1)
for row = 1 : size(pturn_all, 1)
    for s = 1 : length(stim_color)
        pturn_all(row, s*5-4:s*5) = proplize(pturn_all(row, s*5-4:s*5));
    end
end
pturn_all_prop = pturn_all;
save(savename, 'pturn_all_prop', '-append');
writematrix(pturn_all_prop, strcat(basedir_cell{1}, '\results', '\pturn_prop+dark.xlsx'));
clear;

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

% % remove all data.mat
% for j = 1 : length(figlocation_list)
%     savename = strcat(figlocation_list{j}, '\data.mat');
%     delete(savename);
% end

%% Plot pturn from data.mat ===================================================
basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data_new\variability_new_try_extracted\Or42a@Chrimson(3)\T_Re_Sq_219to436P_15_2_3_#T_Bl_Sq_2to7P_15_1_3_'
    };
basedir = basedir_cell{1};
savename = strcat(basedir, '\results', '\data.mat');  % there should be only 1 basedir in basedir_cell
load(savename);

% proportion-lize pturn_all, so that each 5-size pturn is a probability distribution (sum is 1)
stim_color = {'blue', 'red', 'bluered', 'dark'};  % descripiton of the t_stim_start
for row = 1 : size(pturn_all, 1)
    for s = 1 : length(stim_color)
        pturn_all(row, s*5-4:s*5) = proplize(pturn_all(row, s*5-4:s*5));
    end
end
pturn_all_prop = pturn_all;

% order types by whose numbers of larvae in descending order
name_tyep = dec2bin(0:7);  % 1-indexed ---------
nlarva_type = [];  % number of larvae belonging to every type in order of dec2bin(0:63)
for s = 1 : length(pturn_type)
    nlarva_type = [nlarva_type, size(pturn_type{s}, 1)];
end
[nlarva_type_des, sortIdx] = sort(nlarva_type, 'descend');
table_response_nlarva_des = table(name_tyep(sortIdx, :), nlarva_type_des.');
writetable(table_response_nlarva_des, strcat(basedir_cell{1}, '\results', '\response_nlarva_des.xlsx'))

denoise = false;
% Initialize some values for plotting below
pturn = pturn_all;  % what pturn to use, could be pturn_111 or pturn_all, etc.
prefix = 'all_';  % use '001_' or 'all_' etc
if denoise  % the p of first bin minuses the average of rest bins
%     Pbaseline = (pturn(:, 5) + pturn(:, 10) + pturn(:, 15)) / 3;
    Pbaseline = mean(pturn(:, end-4:end), 2);
    Pblue = pturn(:,1) - Pbaseline;
    Pred = pturn(:,6) - Pbaseline;
    Pbluered = pturn(:,11) - Pbaseline;
else
    Pblue = pturn(:,1);
    Pred = pturn(:, 6);
    Pbluered = pturn(:, 11);
end

% order the pturn matrix by one column
column_index = 1;
pturn_all_ordered = order_matrix_rows_from_column(pturn_all, column_index);
imagesc(pturn_all_ordered, [0, 1]); 
% order the pturn matrix by one column by type
pturn_all_ordered_bytype = [order_matrix_rows_from_column(pturn_000, column_index);
    order_matrix_rows_from_column(pturn_001, column_index);
    order_matrix_rows_from_column(pturn_010, column_index);
    order_matrix_rows_from_column(pturn_011, column_index);
    order_matrix_rows_from_column(pturn_100, column_index);
    order_matrix_rows_from_column(pturn_101, column_index);
    order_matrix_rows_from_column(pturn_110, column_index);
    order_matrix_rows_from_column(pturn_111, column_index)];
imagesc(pturn_all_ordered_bytype, [0, 1]); 


% automatical naming
stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start
ntype = 2 ^ (length(stim_color));  % int, -1 if exclude 'dark'
name_type = dec2bin(0:ntype-1);  % 1-indexed
name_type_cell = {};  % 1-by-numofBehaviorType cell
for temp = 1 : length(name_type)
    name_type_cell(temp) = {name_type(temp, :)};
end
nlarva_type = [];  % number of larvae belonging to every type in order of dec2bin(0:63)
for s = 1 : length(pturn_type)
    nlarva_type = [nlarva_type, size(pturn_type{s}, 1)];
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

% plot the matrix pturn_all
stim_color = {'blue', 'red', 'bluered'};  %, 'dark'
figure;
imagesc(pturn_all(:, 1:15), [0, 1]); 
xline([0.5, 5.5, 10.5], 'w', stim_color); yline(boundary_type + 0.5, 'w', name_type_cell);  %---- add 15.5 if include dark, name_type_cell, stim_color
ylabel('Index of larva'); set(gca,'XTick',[])
% cmap = makeColormap([0 0 0.5; 0 0 1; 0 0.5 0; 0 1 0; 0.5 0.5 0; 1 1 0], [4; 0; 2; 0; 4]);
% cmap = makeColormap([0 0 204; 153 204 255; 0 153 76; 0 204 102; 255 128 0; 255 255 204]/255, [4; 0; 2; 0; 4]);
colormap(parula);  % colormap(parula(4)); get the downsampled version of parula colormap that has 4 colors. Empty for gradient.
cbar = colorbar; cbar.Label.String = 'Turn possibility within 3-s bins';
cbar.Ticks = 0:0.2:1;  % tick every 0.25
mkdir(fullfile(basedir, '\results', 'results_fig'));
savename = strcat(basedir, '\results', '\results_fig', '\pturn_all');  % '\pturn_all_discrete+dark'
savefig(gcf, savename);

% plot the matrix v_init_all
figure; imagesc(v_init_all); colormap(); cbar = colorbar; cbar.Label.String = 'Mean speed of 1s before stimulation on (cm/min)';
set(gca,'XTick', 1:4); set(gca,'XTickLabel', stim_color); pbaspect([1, 4, 1]);  % ratio of axis x, y, z is [1, 4, 1]
yline(boundary_type + 0.5, 'w', name_type_cell);
savename = strcat(basedir, '\results', '\results_fig', '\init_speed_all');
savefig(gcf, savename);

% colormap of pturn of both-blue, both-red
figure;
notes = {'(Pbluered - Pblue) / Pbluered', '(Pbluered - Pred) / Pbluered'};  % to label what is it
pturn_change = [(Pbluered - Pblue) ./ Pbluered, (Pbluered - Pred) ./ Pbluered];
zmin = min(pturn_change, [], 'all');  % min, max value of the matrix
clim = [zmin, 1];  % color limit
im = imagesc([colum1, colum2], clim);
xline([0.5, 1.5], 'w', notes); yline(boundary_type + 0.5, 'w', name_type_cell);
ylabel('Index of larva'); set(gca,'XTick',[])
grey = [0.5, 0.5, 0.5];  % grey in RGB
colormap([repmat(grey, 4, 1); parula(4)]);  % keep the 2 colormaps the same size
caxis([-1, 1]);  % assign the color axis to [-1, 1], so 0 is in the middle, divide grey and parula
cbar = colorbar; cbar.Label.String = 'Turn possibility relative change';  % to show colorbar
cbar.Limits = [zmin, 1];  % only show the colorbar from min to 1
cbar.Ticks = [0:0.25:1];
set(gcf, 'Position', gcf().Position .* [1, 1, 0.5, 1])  % shink the width by half
savename = strcat(basedir, '\results', '\results_fig', '\pturn_all_both-each_relalative_discrete');
savefig(gcf, savename);

% hist of pturn of both - pturn of each
figure;
histogram(Pbluered-Pblue)
xlabel('Pbluered - Pblue'); ylabel('Larvae Count')
savename = strcat(basedir, '\results', '\results_fig', '\hist_pboth-pblue');
savefig(gcf, savename);


% colormap of pturn of (both-blue)/both, (both-red)/both
figure;
notes = {'(Pbluered - Pblue) / Pbluered', '(Pbluered - Pred) / Pbluered'};  % to label what is it
imagesc([(Pbluered - Pblue) ./ Pbluered, (Pbluered - Pred) ./ Pbluered]);
xline([0.5, 1.5], 'w', notes); yline(boundary_type + 0.5, 'w', name_type_cell);
ylabel('Index of larva'); set(gca,'XTick',[])
cbar = colorbar; cbar.Label.String = 'Turn possibility relative change';
set(gcf, 'Position', gcf().Position .* [1, 1, 0.5, 1])  % shink the width by half
savename = strcat(basedir, '\results', '\results_fig', '\pturn_all_both-each_relalative');
savefig(gcf, savename);

% hist of pturn in dark and fit
pdark = pturn_all(:, 16:20);
pdark = mean(pdark, 2);  % average each row, nlarvae-by-1
pd = fitdist(pdark, 'Normal');  mu = pd.mu; sigma = pd.sigma;% normal distribution
figure;
h = histogram(pdark, 7); hold on;
binWidth = h.BinWidth; nSamples = length(pdark);
x = linspace(min(pdark), max(pdark), 100);
y = nSamples * binWidth * normpdf(x, mu, sigma);
plot(x, y, 'LineWidth', 2); hold off;
legend('Data Histogram', sprintf('Gaussian Fit (\\mu = %.2f, \\sigma = %.2f)', mu, sigma)); %, 'Location', 'eastoutside')
xlabel('Average pturn in dark for each larva'); ylabel('Count');
savename = strcat(basedir, '\results', '\results_fig', '\hist_pdark_fit');
savefig(gcf, savename);

% study the background activity, which one to choose, p(end) or p_dark
Pturn_dark = mean(pturn_all(:, end-4:end), 2);  % mean of pturn in dark (the last 5 values here), nlarva-by-1
Pbaseline = (pturn_all(:, 5) + pturn_all(:, 10) + pturn_all(:, 15)) / 3;
Pbaseline = [Pbaseline, pturn_all(:, 5)];
Pbaseline = [Pbaseline, pturn_all(:, 10)];
Pbaseline = [Pbaseline, pturn_all(:, 15)];
Pbaseline = [Pbaseline, (pturn_all(:, 2) + pturn_all(:, 7) + pturn_all(:, 12)) / 3;];
figure; hold on;
for f = 1:size(Pbaseline, 2)
    subplot(2, 3, f); hold on;
    plot(Pturn_dark, Pbaseline(:, f), '.'); axis equal; 
    plot([0, 1], [0, 1], 'DisplayName', 'y = x');
    xlim([0, max([Pturn_dark, Pbaseline], [], 'all')]); ylim([0, max([Pturn_dark, Pbaseline], [], 'all')]);
    xlabel('Pturn in dark'); ylabel(num2str(f)); hold off;
end
figure;
histogram(Pturn_dark);
xlabel('Avereage Pturn in dark for each larva'); ylabel('Count')
Pturn_dark_bins = pturn_all(:, end-4:end);
Pturn_dark_bins = Pturn_dark_bins(:);
histogram(Pturn_dark_bins, 10);
xlabel('Pturn in dark for each larva of each 3-s bin'); ylabel('Count')


% study the threshold of deciding if the larva turns
criteria = 'pturn1-pturn345';
switch criteria
    case 'pturn1-pturn_rest'
        criteria_blue = pturn_all(:, 1) - mean(pturn_all(:, 2:5), 2);
        criteriat_red = pturn_all(:, 6) - mean(pturn_all(:, 7:10), 2);
        criteria_bluered = pturn_all(:, 11) - mean(pturn_all(:, 12:15), 2);
    case 'pturn12-pturn_rest'
        criteria_blue = mean(pturn_all(:, 1:2), 2) - mean(pturn_all(:, 3:5), 2);
        criteria_red = mean(pturn_all(:, 6:7), 2) - mean(pturn_all(:, 8:10), 2);
        criteria_bluered = mean(pturn_all(:, 11:12), 2) - mean(pturn_all(:, 13:15), 2);
    case 'pturn1-pturn345'
        criteria_blue = pturn_all(:, 1) - mean(pturn_all(:, 3:5), 2);
        criteria_red = pturn_all(:, 6) - mean(pturn_all(:, 8:10), 2);
        criteria_bluered = pturn_all(:, 11) - mean(pturn_all(:, 13:15), 2);
end
figure; hold on;
histogram(criteria_blue, 'BinWidth', 0.05, 'FaceColor', 'b', 'FaceAlpha', 0.2);
histogram(criteria_red, 'BinWidth', 0.05, 'FaceColor', 'r', 'FaceAlpha', 0.2);
histogram(criteria_bluered, 'BinWidth', 0.05, 'FaceColor', 'k', 'FaceAlpha', 0.2);
ax = gca; ax.FontSize = 20;
hold off;
switch criteria
    case 'pturn1-pturn_rest'
        xlabel('pturn(0<t<3s) - p(turn(3<t<15s)', 'FontSize', 20); ylabel('Count', 'FontSize', 20);
        savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\turn_criteria_1');
        savefig(gcf, savename);
    case 'pturn12-pturn_rest'
        xlabel('pturn(0<t<6s) - p(turn(6<t<15s)', 'FontSize', 20); ylabel('Count', 'FontSize', 20);
        savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\turn_criteria_2');
        savefig(gcf, savename);
    case 'pturn1-pturn345'
        xlabel('pturn(0<t<3s) - p(turn(6<t<15s)', 'FontSize', 20); ylabel('Count', 'FontSize', 20);
        savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\turn_criteria_3');
        savefig(gcf, savename);
end



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
    savename = strcat(basedir, '\results', '\results_fig', ['\', prefix, 'pbluered-pb_pr_denoise']);
    xlim([-0.1, 0.8]); ylim([-0.1, 0.8]); title('Denoise'); savefig(gcf, savename)
else
    savename = strcat(basedir, '\results', '\results_fig', ['\', prefix, 'pbluered-pb_pr']);
    xlim([0, 0.8]); ylim([0, 0.8]); savefig(gcf, savename)
end


figure;  % p( 0<t<3s, blue+red) as a function of p( 0<t<3s, blue) + p( 0<t<3s, red)
x = Pblue + Pred;
y = Pbluered;
sz = 65;
scatter(x, y, sz, 'filled', 'DisplayName', 'Maggots')
lfit = lsline;  % least-squares line
lfit.DisplayName = ['y = ', num2str(diff(lfit.YData) / diff(lfit.XData)), 'x + ', num2str(lfit.YData(1))];
hold on; plot([0, 1], [0, 1], 'DisplayName', 'y = x'); hold off;
axis equal; xlim([0,max(max(x), 1) + 0.2]); ylim([0, 1]);
xlabel('Pturn during the first 3 s of blue plus that of red'); ylabel('Pturn during the first 3 s of both blue and red');
legend();
if denoise
    title('Denoise');
    savename = strcat(basedir, '\results', '\results_fig', ['\', prefix, 'pbluered-pb+pr_denoise']);
    savefig(gcf, savename)
else
    savename = strcat(basedir, '\results', '\results_fig', ['\', prefix, 'pbluered-pb+pr']);
    savefig(gcf, savename)
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
    savename = strcat(basedir, '\results', '\results_fig', ['\', prefix, '3d_plane_fit_denoise']);
    savefig(gcf, savename)
else
    savename = strcat(basedir, '\results', '\results_fig', ['\', prefix, '3d_plane_fit']);
    savefig(gcf, savename)
end


% Bayes theorem check p(turn|blue+red)(1st 3s) = p(turn|blue)(1st 3s) *
% p(turn|red)(1st 3s) / p(turn), where p(turn) = sum(p(turn|i-th 3s), i)
% with i = 1:15
pturn = mean(pturn_all, 2);
plot(pturn_all(:, 11), pturn_all(:, 1) .* pturn_all(:, 6), '.');
xlabel('p(turn|blue+red)[0, 3s]'); ylabel('Predicted by Bayes Theorem');

%%  histogram of some pturn
% mean(pturn_all(:, 3:5), 'all')
% std(pturn_all(:, 3:5), 1, "all")
download = true;
check_off = false;  % check the pturn when light is off, multiple columns.
positive_pturn = false;  % if force negative pturn to be zero
for denoise = [false]  % [true, false] for both plots
    % Initialize some values for plotting below
    pturns = {pturn_all, pturn_type{4}, pturn_type{8}};
    prefixes = {'all_', '011_', '111_'};
%     pturns = {[pturn_001; pturn_011], [pturn_001; pturn_101], [pturn_101; pturn_111], [pturn_011; pturn_111]};
%     prefixes = {'0X1_', 'X01_', '1X1_', 'X11_'};  % X = 0 or 1
    figure;
    for i = 1 : length(pturns)
        pturn = pturns{i};  % what pturn to use, could be pturn_111 or pturn_all, etc.
        prefix = prefixes{i};  % use '001_' or 'all_' etc
        if denoise  % the p of first bin minuses the average of rest bins
            Pbaseline = mean(pturn(:, end-4:end), 2);
            Pblue = pturn(:,1) - Pbaseline;
            Pred = pturn(:,6) - Pbaseline;
            Pbluered = pturn(:,11) - Pbaseline;
        else
            Pblue = pturn(:,1);
            Pred = pturn(:, 6);
            Pbluered = pturn(:, 11);
            if check_off
                Pblue = pturn(:,3:5);
                Pred = pturn(:, 8:10);
                Pbluered = pturn(:, 13:15);  % get the pturn matrix when off, then transfer to array
                Pblue = Pblue(:);
                Pred = Pred(:);
                Pbluered = Pbluered(:);
            end
        end  % end if denoise

        if positive_pturn
            Pblue(Pblue<0) = 0;
            Pred(Pred<0) = 0;
            Pbluered(Pbluered<0) = 0;
        end
        
        % histogram of pturn for different light with different fit
        % % use the same size of bins for each histogram
        fit_func = 'Gaussian';  % 'Gaussian', or 'Poisson'
        edge_min = min([Pblue; Pred; Pbluered]);
        edge_min = floor(edge_min * 10) / 10;  % find it's left boundary with 1 decimal
        edge_max = max([Pblue; Pred; Pbluered]);
        edge_max = ceil(edge_max * 10) / 10;  % find it's right boundary with 1 decimal
        edges = edge_min : 0.1 : edge_max;  % use 0.05 before
        % get the counts
        nblue = histcounts(Pblue, edges);  % 'Normalization', 'probability'
        nred = histcounts(Pred, edges);
        nbluered = histcounts(Pbluered, edges);
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
                subplot(length(pturns), 1, i); hold on;
                lblue = plot(fblue, centers_blue(nblue~=0), nblue(nblue~=0), 'o');  % return the line objects, raw data and fitted line
                bblue = bar(centers_blue(nblue~=0), nblue(nblue~=0), 1, 'FaceColor', [0, 0, 1]);
                lblue(1).Color = [0, 0, 1]; lblue(1).MarkerFaceColor = [0, 0, 1]; lblue(1).DisplayName = ['Blue, mean ', num2str(mean(Pblue)), ', std ', num2str(std(Pblue))]; % properties of dot data
                lblue(2).Color = [0, 0, 1]; lblue(2).DisplayName = ['\mu = ', num2str(fblue.b1), ', \sigma = ', num2str(fblue.c1/sqrt(2))];  % properties of fit line
                bblue.FaceAlpha = 0.5; bblue.DisplayName = 'Blue'; bblue.BarWidth = 1;
        
                lred = plot(fred, centers_red(nred~=0), nred(nred~=0), 'o');  % return the line object
                bred = bar(centers_red(nred~=0), nred(nred~=0), 1, 'FaceColor', [1, 0, 0]);
                lred(1).Color = [1, 0, 0];  lred(1).MarkerFaceColor = [1, 0, 0]; lred(1).DisplayName = ['Red, mean ', num2str(mean(Pred)), ', std ', num2str(std(Pred))]; % properties of dot data
                lred(2).Color = [1, 0, 0]; lred(2).DisplayName = ['\mu = ', num2str(fred.b1), ', \sigma = ', num2str(fred.c1/sqrt(2))];  % properties of fit line
                bred.FaceAlpha = 0.5; bred.DisplayName = 'Red'; bred.BarWidth = 0.8;
        
                lbluered = plot(fbluered, centers_bluered(nbluered~=0), nbluered(nbluered~=0), 'o');  % return the line object
                bbluered = bar(centers_bluered(nbluered~=0), nbluered(nbluered~=0), 1, 'FaceColor', [0, 0, 0]);
                lbluered(1).Color = [0, 0, 0]; lbluered(1).MarkerFaceColor = [0, 0, 0]; lbluered(1).DisplayName = ['Blue and red, mean ', num2str(mean(Pbluered)), ', std ', num2str(std(Pbluered))]; % properties of dot data
                lbluered(2).Color = [0, 0, 0]; lbluered(2).DisplayName = ['\mu = ', num2str(fbluered.b1), ', \sigma = ', num2str(fbluered.c1/sqrt(2))];  % properties of fit line
                bbluered.FaceAlpha = 0.5;  bbluered.DisplayName = 'Blue and red'; bbluered.BarWidth = 0.6;
                hold off; xticks(edges); xlim([0, 1]);
                legend('Location', 'eastoutside');
                xlabel('Pturn during the first 3 s of stimulation high'); ylabel('Number of maggots');
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
        end  % end swith fit function
        title([prefix(1:end-1), ' ', fit_func]);

    end  % end loop for each stimulation type
    if denoise
        title([prefix(1:end-1), ' ', fit_func, ' Denoise']);
        if download
            pause;
            savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', ['hist_pturn_denoise_', fit_func, '_count']);
            savefig(gcf, savename); close;
        end
    else
        if download
            pause;
            savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', ['hist_pturn_', fit_func, '_count']);
            savefig(gcf, savename); close
        end
    end  % end if denoise

end  % end loop for denoise

%% histogram of pturn for different light with polynomial fit
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
%% Try different scatters of pturn

figure;  % Superlinear if p(blue+red) > pblue + pred, vice verse sublinear
x = Pblue;  % p( 0<t<3s, blue)
y = Pred;
sz = 30;  % size of dots in pixels
% c = 1 - pturn_hist_all(:, 13) * [1 1 1];  % the larger p, the smaller c, the darker
c = (x + y) >= Pbluered;
scatter(x, y, sz, c, 'filled'); hold on; 
axis equal;
xlabel('Pturn during the first 3 s of blue'); ylabel('Pturn during the first 3 s of red');
cbar = colorbar;
colormap parula;
cbar.Label.String = 'If pblue + pred >= p(blue+red)';
cbar.Ticks = [0, 1];
cbar.TickLabels = {'Superlinear'; 'Sublinear'};
plot([0, 1], [1, 0]); legend('', 'x + y = 1'); hold off;
if denoise
    xlim([-0.1, 0.8]); ylim([-0.1, 0.8]); title('Denoise');
    savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', 'super_sub-lineaer_denoise');
    savefig(gcf, savename);
else
    xlim([0, 0.8]); ylim([0, 0.8]); 
    savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', 'super_sub-lineaer');
    savefig(gcf, savename);
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
    xlim([-0.1, 0.8]); ylim([-0.1, 0.8]); title('Denoise');
    savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', [prefix, 'linearity_denoise']);
    savefig(gcf, savename);
else
    xlim([0, 0.8]); ylim([0, 0.8]); 
    savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', 'linearity');
    savefig(gcf, savename);
end


% percentage change pred-pblue-percentage change
% x = (pturn_all(:, 11) - pturn_all(:, 1)) ./ pturn_all(:, 11);
% y = (pturn_all(:, 11) - pturn_all(:, 6)) ./ pturn_all(:, 11);
x = pturn_all(:, 1); y = pturn_all(:, 6); 
sz = 40;  % size of dots in pixels
c = (pturn_all(:, 11) - pturn_all(:, 1)) ./ pturn_all(:, 1);
scatter(x, y, sz, c, 'filled')
cbar = colorbar;
colormap parula;
lims = clim;
clim([lims(1), 5]);  % to change the lim of colorbar, and the rearrange the color
% xlabel('(Pboth - Pblue)/Pboth'); ylabel('(Pboth - Pred)/Pboth'); cbar.Label.String = 'Pboth';
xlabel('Pblue'); ylabel('Pred'); cbar.Label.String = '(Pboth - Pblue)/Pblue';
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', ['all_', 'pred-pblue-pbluered_relative_increase_blue']);
savefig(gcf, savename);


% percentage change percentage_change-pblue
% x = (pturn_all(:, 11) - pturn_all(:, 1)) ./ pturn_all(:, 11); 
x = pturn_all(:, 1); 
y = pturn_all(:, 11);
% y = pturn_all(:, 11) - pturn_all(:, 6);
% y = (pturn_all(:, 11) - pturn_all(:, 6)) ./ pturn_all(:, 11);
plot(x, y, 'bo'); % hold on;
% plot([-1, 1], [-1, 1], 'k--'); hold off;
axis equal;
xlabel('Pblue'); ylabel('Pboth');
xlim([0, 1]); ylim([0, 1]);
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', ['all_', 'pboth-pblue']);
savefig(gcf, savename);


% pboth-pred vs pboth-pblue, different types with different shapes
response_marker = {'.', '.', '.', 'o', '.', '.', '.', '+'}; % {'.', 'o', '.', '|', '.', '_', '.', '+'};  % adjust based on needs------
response_color = {'c', 'c', 'c', 'r', 'c', 'c', 'c', 'k'};
figure; hold on;
for m = 1 : length(pturn_type)  % m types of larva
    pturn = pturn_type{m};
    if isempty(pturn)
        x = pturn_type{1}(:, 11) - pturn_type{1}(:, 1);  % if empty, assign it as '000'. Require '000' isn't empty
        y = pturn_type{1}(:, 11) - pturn_type{1}(:, 6);
    else
        x = pturn(:, 11) - pturn(:, 1);  % pboth-pblue
        y = pturn(:, 11) - pturn(:, 6);
    end
    ax(m) = plot(x, y, [response_marker{m}, response_color{m}]);
end 
hold off;
legend(ax([4, 8, 1]), {name_type_cell{[4, 8]}, 'Others'}, 'Location', 'best');  
%legend(ax([2, 4, 6, 8, 1]), {name_type_cell{[2, 4, 6, 8]}, 'Others'}, 'Location', 'best');  %----------
xlabel('Pboth-Pblue'); ylabel('Pboth-Pred');
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', 'pboth-predVSpboth-pblue');
savefig(gcf, savename);
%% Kullback–Leibler divergence
pturn_all = pturn_all_prop;
kld_blue = zeros(size(pturn_all, 1), 1);
kld_red = zeros(size(pturn_all, 1), 1);
kld_bluered = zeros(size(pturn_all, 1), 1);
kld_blue_terms = zeros(size(pturn_all, 1), 5);
kld_red_terms = zeros(size(pturn_all, 1), 5);
kld_bluered_terms = zeros(size(pturn_all, 1), 5);
for j = 1 : size(pturn_all, 1)  % the j-th row (larva) of the pturn
    P_blue = pturn_all(j, 1:5);  % + 0.025 to eliminate 0 probability
    P_red = pturn_all(j, 6:10);  % 1/40, add 1 turn to every bin
    P_bluered = pturn_all(j, 11:15);
    Q = pturn_all(j, 11:15);  % Q = [0.2, 0.2, 0.2, 0.2, 0.2];  %Q = pturn_all(j, 11:15);
    [kld_blue(j), kld_blue_terms(j, :)] = KL_divergence(Q, P_blue);
    [kld_red(j), kld_red_terms(j, :)] = KL_divergence(Q, P_red);
    [kld_bluered(j), kld_bluered_terms(j, :)] = KL_divergence(Q, P_bluered);
end
% % number of larva for each type, in order of 000, 001, 010, 011, 100, 101, 110, 111, done by automatical naming before
% nlarva_type = [size(pturn_000, 1), size(pturn_001, 1), size(pturn_010, 1), size(pturn_011, 1),size(pturn_100, 1), size(pturn_101, 1), size(pturn_110, 1), size(pturn_111, 1)];
% boundary_type = 0.5 + [nlarva_type(1), sum(nlarva_type(1:2)), sum(nlarva_type(1:3)), sum(nlarva_type(1:4)), sum(nlarva_type(1:5)), sum(nlarva_type(1:6)), sum(nlarva_type(1:7))];
% name_type_cell = {'000', '001', '010', '011', '100', '101', '110', '111'};

% DKL verses larva index
figure;  % all 3s bins
ax(1) = subplot(2, 1, 1);
p1 = plot(kld_blue, 'b-');
xline(boundary_type, 'k--');
text(boundary_type - 0.5*nlarva_type, 0.1 + max(kld_blue)*ones(1, length(name_type_cell)), name_type_cell);
xlabel('Larva Index'); ylabel('D_{KL}(Pbluered||Pblue)');
ax(2) = subplot(2, 1, 2);
p2 = plot(kld_red, 'r-'); 
xlabel('Larva Index'); ylabel('D_{KL}(Pbluered||Pred)')
xline(boundary_type, 'k--');  % draw black dash line to seperate different response type
% text(boundary_type - 0.5*nlarva_type, 0.3 + max(kld_red)*ones(1, length(name_type_cell)), name_type_cell);
% legend([p1, p2], {'D_{KL}(Pblue||Pbluered)', 'D_{KL}(Pred||Pbluered)'});
ylim([ax(1), ax(2)], [-0.1, 0.8]);
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\Dkl_both_single_larvaIndex_prop');
savefig(gcf, savename);

% DKL verses larva index with histogram of DKL on the side
figure;  % a 3s bin
ax(1) = subplot(2, 3, [1, 2]);
p1 = plot(kld_blue, 'b-');
xline(boundary_type, 'k--');
xlabel('Larva Index'); ylabel('D_{KL}(Pblue||Pdark)');
ax(2) = subplot(2, 3, 3);
histogram(kld_blue, 'Orientation', 'horizontal'); xlabel('Count of larvae');
ax(3) = subplot(2, 3, [4, 5]);
p2 = plot(kld_red, 'r-'); 
xlabel('Larva Index'); ylabel('D_{KL}(Pred||Pdark)')
xline(boundary_type, 'k--');  % draw black dash line to seperate different response type
ax(4) = subplot(2, 3, 6);
histogram(kld_red, 'Orientation', 'horizontal'); xlabel('Count of larvae');
% legend([p1, p2], {'D_{KL}(Pblue||Pbluered) of 2nd bin', 'D_{KL}(Pred||Pbluered) of 2nd bin'});
ylim([ax(1), ax(3), ax(2), ax(4)], [0, 1.2]);
% sgtitle('Add 1 turn to each bin');
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\Dkl_prop_single_dark_larvaIndex');
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
response_marker = {'.', 'o', '.', '|', '.', '_', '.', '+'};
response_color = {'c', 'g', 'c', 'r', 'c', 'b', 'c', 'k'};
figure; hold on;
for m = 1 : length(pturn_type)  % m types of larva
    pturn = pturn_type{m};
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
legend(ax([2, 4, 6, 8, 1]), {name_type_cell{[2, 4, 6, 8]}, 'Others'}, 'Location', 'best');
savename = strcat(pwd, '\results_fig', '\DKLblue-DKLred_add1turn');
savefig(gcf, savename);


% DKL_blue histogram of each response type
% number of larva for each type, in order of 000, 001, 010, 011, 100, 101, 110, 111
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


%% Pearson correlation coefficient
[R, P, RL, RU] = corrcoef(transpose(pturn_all(:, 6:10)));
figure; imagesc(R, [-1, 1]);
cbar = colorbar; cbar.Label.String = 'Pearson correlation coefficient';
axis equal;
axis xy;
xlabel('Larva index'); ylabel('Larva index');
% number of larva for each type, in order of 000, 001, 010, 011, 100, 101, 110, 111
yline(boundary_type, 'r-', 'LineWidth', 2);  % draw white line to seperate different response type
xline(boundary_type, 'r-', 'LineWidth', 2);  % draw white line to seperate different response type
text(boundary_type - nlarva_type/2, 4 + size(pturn_all, 1)*ones(1, length(boundary_type)), name_type_cell, ...
    'FontSize', 12, 'Color', 'red');
text(2 + size(pturn_all, 1)*ones(1, length(boundary_type)), boundary_type - nlarva_type/2, name_type_cell, ...
    'FontSize', 12, 'Color', 'red');
% title('PCC of pturn in bluered');
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\pcc_pturn_red');
savefig(gcf, savename);


%% plot the pturn from larva index 
index_larva = 95:104;
edges = 0:3:15;
xbar = edges(1: numel(edges)-1) + diff(edges)/2;
for j = 1:length(index_larva)  % plot the j-th larva
    P_blue = pturn_all(index_larva(j), 1:5);
    P_red = pturn_all(index_larva(j), 6:10);
    P_bluered = pturn_all(index_larva(j), 11:15);
    subplot(3, length(index_larva), j)
    bar(xbar, P_blue); ylim([0, 1]);
    subplot(3, length(index_larva), j+length(index_larva))
    bar(xbar, P_red); ylim([0, 1]);
    subplot(3, length(index_larva), j+length(index_larva)*2)
    bar(xbar, P_bluered); ylim([0, 1]);
    title(['larva ', num2str(index_larva(j))]);
end
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\pturn_111');
savefig(gcf, savename);

%% save the current file to \results in the basedir
savename = strcat(basedir_cell{1}, '\results');
copyfile('info_from_pturn.m', savename)  % copy the file to the location
