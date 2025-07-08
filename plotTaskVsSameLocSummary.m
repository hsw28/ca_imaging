function plotTaskVsSameLocSummary_updated
% Summary plot using getTaskAndSameLocSpikes to extract spike counts
% from task-defined and non-task-defined hulls for An-2, An-1, An



ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};
win = [0 2];
colors = lines(numel(ratNames));

nRats = numel(ratNames);
rate_task_taskhull = nan(nRats,1);
rate_nontask_taskhull = nan(nRats,1);
rate_task_poshull = nan(nRats,1);
rate_nontask_poshull = nan(nRats,1);
sem_task_taskhull = nan(nRats,1);
sem_nontask_taskhull = nan(nRats,1);
sem_task_poshull = nan(nRats,1);
sem_nontask_poshull = nan(nRats,1);
pvals_taskhull = nan(nRats,1);
pvals_poshull = nan(nRats,1);

figure('Color','w','Position',[200 300 1100 450]);

for r = 1:nRats
  ratNames{r}
    rat = evalin('base', ratNames{r});
    dateList = autoDateList(rat);
    idx = find(strcmp(dateList, rat.An));
    theseDays = dateList(idx-2:idx);

    rt1_all = []; rt2_all = []; rp1_all = []; rp2_all = [];

    for d = 1:3
        dateStr = theseDays{d};
        spikeMat = rat.Ca_peaks.(['CA_peaks_' dateStr]);
        ts = rat.pos.(['pos_' dateStr])(:,1);
        nNeurons = size(spikeMat,1);

        for ni = 1:nNeurons
            spikes = spikeMat(ni,:);
            spikes = spikes(~isnan(spikes));
            if numel(spikes) < 3
                continue
            end

            [rt1, rt2, rp1, rp2] = getTaskAndSameLocSpikes(rat, ni, dateStr, win, false);
            if any(isnan([rt1, rt2, rp1, rp2]))
                continue
            end
            rt1_all(end+1) = rt1;
            rt2_all(end+1) = rt2;
            rp1_all(end+1) = rp1;
            rp2_all(end+1) = rp2;
        end
    end

    % Store means and SEMs
    rate_task_taskhull(r) = mean(rt1_all);
    rate_nontask_taskhull(r) = mean(rt2_all);
    sem_task_taskhull(r) = std(rt1_all)/sqrt(numel(rt1_all));
    sem_nontask_taskhull(r) = std(rt2_all)/sqrt(numel(rt2_all));

    rate_task_poshull(r) = mean(rp1_all);
    rate_nontask_poshull(r) = mean(rp2_all);
    sem_task_poshull(r) = std(rp1_all)/sqrt(numel(rp1_all));
    sem_nontask_poshull(r) = std(rp2_all)/sqrt(numel(rp2_all));

    % Paired t-tests
    if numel(rt1_all) >= 3 && numel(rt2_all) >= 3
      fprintf('task hull')
        [~, p, ~, stats] = ttest(rt1_all, rt2_all)
        pvals_taskhull(r) = p;
          [~, pvals_taskhull2, ~, stats] = ttest(rt1_all, rt2_all)
    end
    if numel(rp1_all) >= 3 && numel(rp2_all) >= 3
            fprintf('pos hull')
        [~, p, ~, stats] = ttest(rp1_all, rp2_all)
        pvals_poshull(r) = p;
          [~, pvals_taskhull2, ~, stats] = ttest2(rt1_all, rt2_all)
    end
end


% --- task hull plot ---
subplot(1,2,1); cla; hold on;

barData = [rate_task_taskhull, rate_nontask_taskhull];  % [nRats x 2]
nRats = size(barData, 1);
barHandle = bar(barData, 'grouped');

% Get bar X positions robustly
if numel(barHandle) == 2 && isprop(barHandle(1), 'XEndPoints')
    x1 = barHandle(1).XEndPoints;
    x2 = barHandle(2).XEndPoints;
else
    % barHandle is scalar, access children patches
    % Order of Children: rightmost series first, so flip
    barChildren = flip(barHandle.Children);  % 2 elements
    x1 = mean(barChildren(1).Vertices(:,1));  % first group (column 1)
    x2 = mean(barChildren(2).Vertices(:,1));  % second group (column 2)
end

% Plot error bars directly at bar centers
errorbar(x1, barData(:,1), sem_task_taskhull, 'k', ...
    'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 8);
errorbar(x2, barData(:,2), sem_nontask_taskhull, 'k', ...
    'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 8);

% Labels and title
xticks(1:nRats); xticklabels(ratNames);
ylabel('Firing Rate (Hz)');
title('Hull from task spikes');
legend({'In Hull', 'Out of Hull'}, 'Location', 'northwest');
set(gca, 'Box', 'off');

% Significance star
for r = 1:nRats
    yMax = max(barData(r,:)) + max([sem_task_taskhull(r), sem_nontask_taskhull(r)]);
    xStar = mean([x1(r), x2(r)]);
    if pvals_taskhull(r) < 0.05
        text(xStar, yMax * 1.05, '***', ...
            'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end

end


% --- Pos hull plot ---
subplot(1,2,2); cla; hold on;

barData2 = [rate_nontask_poshull, rate_task_poshull];  % [nRats x 2]
nRats = size(barData2, 1);
barHandle2 = bar(barData2, 'grouped');

% Set bar colors
if numel(barHandle2) >= 2
    barHandle2(1).FaceColor = [0.2 0.4 0.8];  % In Hull
    barHandle2(2).FaceColor = [0.7 0.7 0.7];  % Out of Hull
end

% Get bar X positions robustly
if numel(barHandle2) == 2 && isprop(barHandle2(1), 'XEndPoints')
    x1 = barHandle2(1).XEndPoints;
    x2 = barHandle2(2).XEndPoints;
else
    barChildren = flip(barHandle2.Children);  % Rightmost = column 1
    x1 = unique(barChildren(1).XData);  % In Hull
    x2 = unique(barChildren(2).XData);  % Out of Hull
end

% Plot error bars
errorbar(x1, barData2(:,1), sem_nontask_poshull, 'k', ...
    'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 8);
errorbar(x2, barData2(:,2), sem_task_poshull, 'k', ...
    'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 8);

% Labels and title
xticks(1:nRats); xticklabels(ratNames);
ylabel('Firing Rate (Hz)');
title('Hull from non-task spikes');
legend({'In Hull', 'Out of Hull'}, 'Location', 'northwest');
set(gca, 'Box', 'off');

% Significance stars
for r = 1:nRats
    yMax = max(barData2(r,:)) + max([sem_task_poshull(r), sem_nontask_poshull(r)], [], 'omitnan');
    xStar = mean([x1(r), x2(r)]);
    if pvals_poshull(r) < 0.05
        text(xStar, yMax * 1.05, '***', ...
            'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment','center');
    end
end
