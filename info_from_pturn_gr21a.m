% information from pturn, including pturn, hist(pturn), DKL

%% create pturn matrix from the structure larvae
basedir_cell = {
    'G:\AS-Filer\PHY\mmihovil\Shared\data_variability\variability_extracted\T_Re_Sq_219to436P_15_2_3#T_Bl_Sq_2to7P_15_1_3'
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
    mkdir(fullfile(basedir, 'results'));  % create a overall results folder to contain general data
    for k = 1 : length(x)  % loop for each expt in the eset
        j = j + 1;
        figlocation_list{j} = fullfile(basedir, ['results', d(x(k)).name(end-16:end-4)]);
    end
end

% create pturn_all in order of categories and save in data.mat
% pturn_all = [];  % nlarva-by-15 array, 3s bins
stim_color = {'blue', 'red', 'bluered'};  % descripiton of the t_stim_start
ntype = 2 ^ (length(stim_color));  % int, -1 in the power if exclude dark in the end ------------------
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

%% Plot pturn from data.mat ===================================================
basedir = basedir_cell{1};
savename = strcat(basedir, '\results', '\data.mat');  % there should be only 1 basedir in basedir_cell
load(savename);

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
colormap(parula);  % colormap(parula(4)); get the downsampled version of parula colormap that has 4 colors. Empty for gradient.
cbar = colorbar; cbar.Label.String = 'Turn possibility within 3-s bins';
cbar.Ticks = 0:0.2:1;  % tick every 0.25
mkdir(fullfile(basedir, '\results', 'results_fig'));
savename = strcat(basedir, '\results', '\results_fig', '\pturn_all');  % '\pturn_all_discrete+dark'
savefig(gcf, savename);

%%  histogram of some pturn
download = true;
check_off = false;  % check the pturn when light is off, multiple columns.
positive_pturn = false;  % if force negative pturn to be zero
for denoise = [false]  % [true, false] for both plots
    % Initialize some values for plotting below
    pturns = {pturn_all, pturn_type{2}, pturn_type{4}, pturn_type{6}, pturn_type{8}};
    prefixes = {'all_', '001_', '011_', '101_', '111_'};
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


%% Try different scatters of pturn

% pboth-pred vs pboth-pblue, different types with different shapes
response_marker = {'.', 'o', '.', '|', '.', '_', '.', '+'} ; % {'.', '.', '.', 'o', '.', '.', '.', '+'};  % adjust based on needs------
response_color = {'c', 'c', 'c', 'r', 'c', 'b', 'c', 'k'};
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
% legend(ax([4, 8, 1]), {name_type_cell{[4, 8]}, 'Others'}, 'Location', 'best');  
legend(ax([2, 4, 6, 8, 1]), {name_type_cell{[2, 4, 6, 8]}, 'Others'}, 'Location', 'best');  %----------
xlabel('Pboth-Pblue'); ylabel('Pboth-Pred');
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\', 'pboth-predVSpboth-pblue');
savefig(gcf, savename);
%% Kullback–Leibler divergence
kld_blue = zeros(size(pturn_all, 1), 1);
kld_red = zeros(size(pturn_all, 1), 1);
kld_blue_terms = zeros(size(pturn_all, 1), 5);
kld_red_terms = zeros(size(pturn_all, 1), 5);
for j = 1 : size(pturn_all, 1)  % the j-th row (larva) of the pturn
    P_blue = pturn_all(j, 1:5);
    P_red = pturn_all(j, 6:10);
    Q = pturn_all(j, 11:15);
    [kld_blue(j), kld_blue_terms(j, :)] = KL_divergence(P_blue, Q);
    [kld_red(j), kld_red_terms(j, :)] = KL_divergence(P_red, Q);
end

% DKL verses larva index
figure;  % all 3s bins
ax(1) = subplot(2, 1, 1);
p1 = plot(kld_blue, 'b-');
xline(boundary_type, 'k--');
text(boundary_type - 0.5*nlarva_type, 0.1 + max(kld_blue)*ones(1, length(name_type_cell)), name_type_cell);
xlabel('Larva Index'); ylabel('D_{KL}(Pblue||Pbluered)');
ax(2) = subplot(2, 1, 2);
p2 = plot(kld_red, 'r-'); 
xlabel('Larva Index'); ylabel('D_{KL}(Pred||Pbluered)')
xline(boundary_type, 'k--');  % draw black dash line to seperate different response type
text(boundary_type - 0.5*nlarva_type, 0.1 + max(kld_red)*ones(1, length(name_type_cell)), name_type_cell);
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\Dkl_both_single_larvaIndex');
savefig(gcf, savename);


%% Pearson correlation coefficient
[R, P, RL, RU] = corrcoef(transpose(pturn_all(:, 1:5)));
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
savename = strcat(basedir_cell{1}, '\results', '\results_fig', '\pcc_pturn_blue');
savefig(gcf, savename);


%% save the current file to \results in the basedir
savename = strcat(basedir_cell{1}, '\results');
copyfile('info_from_pturn_gr21a.m', savename)  % copy the file to the location
