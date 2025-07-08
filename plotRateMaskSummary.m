function plotRateMaskSummary(N)
% Summary plot using rate masks defined as N std above mean
% For An-2, An-1, and An across all rats

if nargin < 1, N = 1; end  % std threshold

ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};
nRats = numel(ratNames);

% Preallocate stats
rate_task_mask = nan(nRats,1);
rate_nontask_mask = nan(nRats,1);
rate_task_mask_non = nan(nRats,1);
rate_nontask_mask_non = nan(nRats,1);
sem_task_mask = nan(nRats,1);
sem_nontask_mask = nan(nRats,1);
sem_task_mask_non = nan(nRats,1);
sem_nontask_mask_non = nan(nRats,1);
pvals_taskmask = nan(nRats,1);
pvals_nonmask = nan(nRats,1);

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
        nNeurons = size(spikeMat,1);

        for ni = 1:nNeurons
            spikes = spikeMat(ni,:); spikes = spikes(~isnan(spikes));
            if numel(spikes) < 3, continue; end

            [rt1, rt2, rp1, rp2] = RateMaskVsTask_summary(rat, ni, dateStr, N);
            if any(isnan([rt1, rt2, rp1, rp2])), continue; end

            rt1_all(end+1) = rt1;
            rt2_all(end+1) = rt2;
            rp1_all(end+1) = rp1;
            rp2_all(end+1) = rp2;
        end
    end

    % Compute means and SEMs
    rate_task_mask(r) = mean(rt1_all);  sem_task_mask(r) = std(rt1_all)/sqrt(numel(rt1_all));
    rate_nontask_mask(r) = mean(rt2_all);  sem_nontask_mask(r) = std(rt2_all)/sqrt(numel(rt2_all));
    rate_task_mask_non(r) = mean(rp1_all); sem_task_mask_non(r) = std(rp1_all)/sqrt(numel(rp1_all));
    rate_nontask_mask_non(r) = mean(rp2_all); sem_nontask_mask_non(r) = std(rp2_all)/sqrt(numel(rp2_all));

    % Paired t-tests
    if numel(rt1_all) >= 3 && numel(rt2_all) >= 3
        fprintf('task mask')
        [~, p, ~, stats] = ttest(rt1_all, rt2_all)
        pvals_taskmask(r) = p;
        [~, p2, ~, stats2] = ttest2(rt1_all, rt2_all)
    end
    if numel(rp1_all) >= 3 && numel(rp2_all) >= 3
      fprintf('non-task mask')
        [~, p, ~, stats] = ttest(rp1_all, rp2_all)
        pvals_nonmask(r) = p;
        [~, p2, ~, stats2] = ttest(rp1_all, rp2_all)
    end
end

% --- Plot Task-based mask ---
subplot(1,2,1); cla; hold on;
barData = [rate_task_mask, rate_nontask_mask];
barHandle = bar(barData, 'grouped');
x1 = barHandle(1).XEndPoints;
x2 = barHandle(2).XEndPoints;
errorbar(x1, barData(:,1), sem_task_mask, 'k', 'LineStyle', 'none');
errorbar(x2, barData(:,2), sem_nontask_mask, 'k', 'LineStyle', 'none');
xticks(1:nRats); xticklabels(ratNames);
ylabel('Firing Rate (Hz)');
title('Task vs Non-task in task-defined mask');
legend({'In Mask','Out of Mask'}, 'Location', 'northwest');
set(gca,'Box','off');
barHandle2(1).FaceColor = [0.2 0.4 0.8];  % Task
barHandle2(2).FaceColor = [0.7 0.7 0.7];  % Non-task
for r = 1:nRats
    if pvals_taskmask(r) < 0.001
        text(mean([x1(r), x2(r)]), max(barData(r,:))*1.05, '***', 'FontSize', 16, 'HorizontalAlignment','center');
    elseif pvals_taskmask(r) < 0.01
      text(mean([x1(r), x2(r)]), max(barData(r,:))*1.05, '**', 'FontSize', 16, 'HorizontalAlignment','center');
    elseif pvals_taskmask(r) < 0.05
      text(mean([x1(r), x2(r)]), max(barData(r,:))*1.05, '*', 'FontSize', 16, 'HorizontalAlignment','center');
    end
end



% --- Plot Non-task-based mask ---
subplot(1,2,2); cla; hold on;
barData2 = [rate_task_mask_non, rate_nontask_mask_non];
barHandle2 = bar(barData2, 'grouped');
x1 = barHandle2(1).XEndPoints;
x2 = barHandle2(2).XEndPoints;
errorbar(x1, barData2(:,1), sem_task_mask_non, 'k', 'LineStyle', 'none');
errorbar(x2, barData2(:,2), sem_nontask_mask_non, 'k', 'LineStyle', 'none');
xticks(1:nRats); xticklabels(ratNames);
ylabel('Firing Rate (Hz)');
title('Task vs Non-task in non-task-defined mask');
legend({'Task', 'Non-task'}, 'Location', 'northwest');
set(gca,'Box','off');
barHandle2(1).FaceColor = [0.2 0.4 0.8];  % Task
barHandle2(2).FaceColor = [0.7 0.7 0.7];  % Non-task

for r = 1:nRats
    if pvals_nonmask(r) < 0.001
        text(mean([x1(r), x2(r)]), max(barData2(r,:))*1.05, '***', 'FontSize', 16, 'HorizontalAlignment','center');
    elseif pvals_nonmask(r) < 0.01
      text(mean([x1(r), x2(r)]), max(barData2(r,:))*1.05, '**', 'FontSize', 16, 'HorizontalAlignment','center');
    elseif pvals_nonmask(r) < 0.05
      text(mean([x1(r), x2(r)]), max(barData2(r,:))*1.05, '*', 'FontSize', 16, 'HorizontalAlignment','center');

    end
end


end



function [rate1, rate2, rate3, rate4] = RateMaskVsTask_summary(animal, neuronIdx, dateStr, N)

rate1 = NaN; rate2 = NaN; rate3 = NaN; rate4 = NaN;

pos = animal.pos.(['pos_' dateStr]);
spikes = animal.Ca_peaks.(['CA_peaks_' dateStr])(neuronIdx, :);
cs_times = animal.CS_times.(['CS_' dateStr]);
ts = pos(:,1);
pos_xy = pos(:,2:3);
spikes = spikes(~isnan(spikes));

% Skip if no spikes
if isempty(spikes), return; end

% Task mask
taskMask = false(size(ts));
for i = 1:numel(cs_times)
    taskMask = taskMask | (ts >= cs_times(i) & ts <= cs_times(i) + 2);
end
taskTimes = ts(taskMask);
nonTaskTimes = ts(~taskMask);

% Split spikes
taskSpikes = spikes(ismembertol(spikes, taskTimes, 1e-3));
nonTaskSpikes = setdiff(spikes, taskSpikes, 'stable');

% Skip if no spikes in one or both groups
if isempty(taskSpikes) || isempty(nonTaskSpikes)
    return;
end

% --- Task mask ---
rate_task = CA_normalizePosData(taskSpikes(:), pos, 2.5, 1);
rate_task = imgaussfilt(rate_task, .75);
rate_task(rate_task == 0 | isnan(rate_task)) = NaN;
mu_task = nanmean(rate_task(:));
sigma_task = nanstd(rate_task(:));
mask_task = rate_task > mu_task + N * sigma_task;

% --- Non-task mask ---
rate_non = CA_normalizePosData(nonTaskSpikes(:), pos, 2.5, 1);
rate_non = imgaussfilt(rate_non, .75);
rate_non(rate_non == 0 | isnan(rate_non)) = NaN;
mu_non = nanmean(rate_non(:));
sigma_non = nanstd(rate_non(:));
mask_non = rate_non > mu_non + N * sigma_non;

% --- Bin mapping
Xedges = linspace(min(pos_xy(:,1)), max(pos_xy(:,1)), size(rate_task,2)+1);
Yedges = linspace(min(pos_xy(:,2)), max(pos_xy(:,2)), size(rate_task,1)+1);
getBinIdx = @(xy) deal(discretize(xy(:,1), Xedges), discretize(xy(:,2), Yedges));

[xT_task, yT_task] = getBinIdx(pos_xy(ismembertol(ts, taskSpikes, 1e-3),:));
[xN_task, yN_task] = getBinIdx(pos_xy(ismembertol(ts, nonTaskSpikes, 1e-3),:));

% --- Task-defined mask
validT = xT_task > 0 & yT_task > 0 & xT_task <= size(rate_task,2) & yT_task <= size(rate_task,1);
validN = xN_task > 0 & yN_task > 0 & xN_task <= size(rate_task,2) & yN_task <= size(rate_task,1);
if ~any(validT) || ~any(validN), return; end

in_task_T = mask_task(sub2ind(size(rate_task), yT_task(validT), xT_task(validT)));
in_task_N = mask_task(sub2ind(size(rate_task), yN_task(validN), xN_task(validN)));
rate1 = sum(in_task_T) / sum(validT);  % Task in mask
rate2 = sum(in_task_N) / sum(validN);  % Non-task in task-defined mask

% --- Non-task-defined mask
validT = xT_task > 0 & yT_task > 0 & xT_task <= size(rate_non,2) & yT_task <= size(rate_non,1);
validN = xN_task > 0 & yN_task > 0 & xN_task <= size(rate_non,2) & yN_task <= size(rate_non,1);
if ~any(validT) || ~any(validN), return; end

in_non_T = mask_non(sub2ind(size(rate_non), yT_task(validT), xT_task(validT)));
in_non_N = mask_non(sub2ind(size(rate_non), yN_task(validN), xN_task(validN)));
rate3 = sum(in_non_T) / sum(validT);  % Task in non-task-defined mask
rate4 = sum(in_non_N) / sum(validN);  % Non-task in non-task-defined mask
end
