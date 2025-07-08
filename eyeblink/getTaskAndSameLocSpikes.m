function [rate_task_taskhull, rate_nontask_taskhull, ...
          rate_task_poshull, rate_nontask_poshull] = ...
          getTaskAndSameLocSpikes(animal, neuronIdx, dateStr, win, plot_me)

% Initialize outputs
rate_task_taskhull = NaN;
rate_nontask_taskhull = NaN;
rate_task_poshull = NaN;
rate_nontask_poshull = NaN;

% Determine plotting flag
if nargin < 5
    plot_me = numel(neuronIdx) < 5;
else
    plot_me = logical(plot_me);  % convert 'false'/'true' string or 0/1 to logical
end


% Extract data
pos = animal.pos.(['pos_' dateStr]);
spikeMat = animal.Ca_peaks.(['CA_peaks_' dateStr]);
cs_times = animal.CS_times.(['CS_' dateStr]);

ts = pos(:,1);
pos_xy = pos(:,2:3);

% Create task epoch mask
taskMask = false(size(ts));
for i = 1:length(cs_times)
    taskMask = taskMask | (ts >= cs_times(i) & ts <= cs_times(i) + win(2));
end

if plot_me
    figure('Color','w','Position',[100 300 1000 400]);
end

for k = 1:2
    if k == 1
        refLabel = 'Task-based hull';
        refMask = taskMask;
        targetMask = taskMask;
        matchMask = ~taskMask;
        subplotIndex = 1;
    else
        refLabel = 'Non-task-based hull';
        refMask = ~taskMask;
        targetMask = ~taskMask;
        matchMask = taskMask;
        subplotIndex = 2;
    end

    % Get spikes
    spikes = spikeMat(neuronIdx, :);
    spikes = spikes(~isnan(spikes));

    relevantTimes = ts(refMask);
    refSpikes = spikes(ismembertol(spikes, relevantTimes, 1e-3));
    [~, posIdxs] = min(abs(refSpikes' - ts'), [], 2);

    if numel(posIdxs) < 3
        warning('Not enough spike positions to compute convex hull. Skipping...');
        continue
    end

    hullPos = pos_xy(posIdxs, :);
    k_hull = convhull(hullPos(:,1), hullPos(:,2));
    inPoly = inpolygon(pos_xy(:,1), pos_xy(:,2), hullPos(k_hull,1), hullPos(k_hull,2));
    matchMask = matchMask & inPoly;

    taskTimes = ts(taskMask & inPoly);
    nonTaskTimes = ts(~taskMask & inPoly);

    dt = median(diff(ts));  % time step
    taskDur = numel(taskTimes) * dt;
    nonTaskDur = numel(nonTaskTimes) * dt;

    spikes_task = sum(ismembertol(spikes, taskTimes, 1e-3));
    spikes_nontask = sum(ismembertol(spikes, nonTaskTimes, 1e-3));



    rate_task = spikes_task / taskDur;
    rate_nontask = spikes_nontask / nonTaskDur;

    if k == 1
        rate_task_taskhull = rate_task;
        rate_nontask_taskhull = rate_nontask;
    else
        rate_task_poshull = rate_task;
        rate_nontask_poshull = rate_nontask;
    end

    % Optional plotting
    if plot_me
        subplot(1,2,subplotIndex); hold on;
        bar(1, rate_task, 0.5, 'FaceColor', [0.5 0.2 0.8]);
        bar(2, rate_nontask, 0.5, 'FaceColor', [0.2 0.6 0.7]);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Task', 'Non-task'});
        ylabel('Firing Rate (Hz)');
        title(sprintf('Neuron %d – %s\n(%s)', neuronIdx, strrep(dateStr, '_','-'), refLabel));
    end
end
