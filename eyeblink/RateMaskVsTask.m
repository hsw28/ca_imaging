function RateMaskVsTask(animal, neuronVec, dateStr, N, dim)
% Compare spike rates within high-rate masks defined by task or non-task spikes

if nargin < 4, N = 1; end %STD above mean to use
if nargin < 5, dim = 2.5; end

for ni = 1:numel(neuronVec)
    neuronIdx = neuronVec(ni);
    pos = animal.pos.(['pos_' dateStr]);
    spikes = animal.Ca_peaks.(['CA_peaks_' dateStr])(neuronIdx,:);
    cs_times = animal.CS_times.(['CS_' dateStr]);
    ts = pos(:,1);
    pos_xy = pos(:,2:3);
    spikes = spikes(~isnan(spikes));

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

    % -------------------
    % Create Task Rate Map
    % -------------------
    % For task-based mask
    rate_task = CA_normalizePosData(taskSpikes, pos, dim, 1.0);
    rate_task(isnan(rate_task)) = 0;
    rate_task(rate_task == 0) = NaN;
    mu_task = nanmean(rate_task(:));
    sigma_task = nanstd(rate_task(:));
    thresh_task = mu_task + N * sigma_task;  % e.g., N = 1.5
    mask_task = rate_task > thresh_task;



    % -------------------
    % Create Non-task Rate Map
    % -------------------
    rate_non = CA_normalizePosData(nonTaskSpikes, pos, dim, 1.0);
    rate(isnan(rate_non)) = 0;
    rate_non(rate_non == 0) = NaN;
    mu_non = nanmean(rate_non(:));
    sigma_non = nanstd(rate_non(:));
    thresh_non = mu_non + N * sigma_non;
    mask_non = rate_non > thresh_non;



    % Bin grid
    Xedges = linspace(min(pos_xy(:,1)), max(pos_xy(:,1)), size(rate_task,2)+1);
    Yedges = linspace(min(pos_xy(:,2)), max(pos_xy(:,2)), size(rate_task,1)+1);

    % --- Map spikes to bin positions
    getBinIdx = @(posList) deal( ...
        discretize(posList(:,1), Xedges), ...
        discretize(posList(:,2), Yedges));


    [xT_task, yT_task] = getBinIdx(pos_xy(ismembertol(ts, taskSpikes, 1e-3),:));
    [xN_task, yN_task] = getBinIdx(pos_xy(ismembertol(ts, nonTaskSpikes, 1e-3),:));


    % --- Copy bin indices for non-task mask reuse
    xT_non = xT_task;
    yT_non = yT_task;
    xN_non = xN_task;
    yN_non = yN_task;


    % --- Compute rates for task-defined mask
    valid_T = xT_task > 0 & yT_task > 0 & xT_task <= size(rate_task,2) & yT_task <= size(rate_task,1);
    valid_N = xN_task > 0 & yN_task > 0 & xN_task <= size(rate_task,2) & yN_task <= size(rate_task,1);
    in_task_T = mask_task(sub2ind(size(rate_task), yT_task(valid_T), xT_task(valid_T)));
    in_task_N = mask_task(sub2ind(size(rate_task), yN_task(valid_N), xN_task(valid_N)));

    rate1 = sum(in_task_T) / sum(valid_T);
    rate2 = sum(in_task_N) / sum(valid_N);

    % --- Compute rates for non-task-defined mask
    valid_T = xT_non > 0 & yT_non > 0 & xT_non <= size(rate_non,2) & yT_non <= size(rate_non,1);
    valid_N = xN_non > 0 & yN_non > 0 & xN_non <= size(rate_non,2) & yN_non <= size(rate_non,1);
    in_non_T = mask_non(sub2ind(size(rate_non), yT_non(valid_T), xT_non(valid_T)));
    in_non_N = mask_non(sub2ind(size(rate_non), yN_non(valid_N), xN_non(valid_N)));

    rate3 = sum(in_non_T) / sum(valid_T);
    rate4 = sum(in_non_N) / sum(valid_N);

    fprintf('Neuron %d (%s)\n', neuronIdx, dateStr);
    fprintf('  Task-defined mask: Task = %.2f%%, Non-task = %.2f%%\n', 100*rate1, 100*rate2);
    fprintf('  Non-task-defined mask: Task = %.2f%%, Non-task = %.2f%%\n', 100*rate3, 100*rate4);

    % --- Visualization
    figure('Name', sprintf('Neuron %d – %s', neuronIdx, strrep(dateStr,'_','-')), ...
           'Color','w','Position',[300 300 1100 400]);
           hold on

    subplot(1,3,1); imagesc(rate_task); axis image off;
    colormap(gca, 'hot'); title('Task-only rate');
    colorbar;

    subplot(1,3,2); imagesc(mask_task); axis image off;
    colormap(gca, 'gray'); title(sprintf('Task Mask (top %d%% std above mean)', N));
    colorbar;


    subplot(1,3,3); imagesc(mask_non); axis image off;
    colormap(gca, 'gray'); title(sprintf('Non-task Mask (%d%% std above mean)', N));
    colorbar;

    % --- Firing rates per bin (Hz)
  dt = median(diff(ts));  % Time resolution in sec
  fr_in_task_T = in_task_T / dt;
  fr_in_task_N = in_task_N / dt;
  fr_in_non_T = in_non_T / dt;
  fr_in_non_N = in_non_N / dt;

  % --- Mean + SEM
  m1 = mean(fr_in_task_T);   s1 = std(fr_in_task_T) / sqrt(numel(fr_in_task_T));
  m2 = mean(fr_in_task_N);   s2 = std(fr_in_task_N) / sqrt(numel(fr_in_task_N));
  m3 = mean(fr_in_non_T);    s3 = std(fr_in_non_T) / sqrt(numel(fr_in_non_T));
  m4 = mean(fr_in_non_N);    s4 = std(fr_in_non_N) / sqrt(numel(fr_in_non_N));

  % Group: [InMask:Task, InMask:NonTask] × [MaskType]
  barData = [m1, m2;  % Task-defined mask
             m3, m4]; % Non-task-defined mask
  semData = [s1, s2;
             s3, s4];

  % --- Plot
  figure('Color','w','Position',[500 300 550 400]);
  bh = bar(barData, 'grouped'); hold on;

  % Colors
  bh(1).FaceColor = [0.2 0.4 0.8];  % Task
  bh(2).FaceColor = [0.7 0.6 0.6];  % Non-task

  % Error bars
  [ngroups, nbars] = size(barData);
  groupwidth = min(0.8, nbars/(nbars + 1.5));
  for i = 1:nbars
      x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
      errorbar(x, barData(:,i), semData(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
  end

  % Axes/labels
  xticks(1:2);
  xticklabels({'Task-defined mask', 'Non-task-defined mask'});
  ylabel('Firing Rate (Hz)');
  legend({'During task', 'During non-task'}, 'Location', 'northwest');
  title(sprintf('Neuron %d – %s', neuronIdx, strrep(dateStr,'_','-')));
  box off;

  % Stats: paired t-tests
  [p1, ~] = ttest(fr_in_task_T, fr_in_task_N);
  [p2, ~] = ttest(fr_in_non_T, fr_in_non_N);

  % Significance stars
  %the star attached to indicate P < 0.05 and two stars to indicate P < 0.01. Occasionally three stars were used to indicate P < 0.001.
  yMax = max(barData(:) + semData(:)) * 1.1;
  if p1 < 0.001
      text(1, yMax, '***', 'HorizontalAlignment','center', 'FontSize', 16);
  elseif p1 < 0.01
      text(1, yMax, '**', 'HorizontalAlignment','center', 'FontSize', 16);
  elseif p1< 0.05
    text(1, yMax, '**', 'HorizontalAlignment','center', 'FontSize', 16);
  end

  if p2 < 0.001
      text(2, yMax, '***', 'HorizontalAlignment','center', 'FontSize', 16);
  elseif p2 < 0.01
      text(2, yMax, '**', 'HorizontalAlignment','center', 'FontSize', 16);
  elseif p2< 0.05
    text(2, yMax, '**', 'HorizontalAlignment','center', 'FontSize', 16);
  end
end
end
