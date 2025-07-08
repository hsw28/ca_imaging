function plotExampleHullOverlay(animal, neuronVec, dateStr, win)
% Plots position with convex hulls from task and non-task spikes
% for one or more neurons on a given session.

% Load session data
pos = animal.pos.(['pos_' dateStr]);
ts = pos(:,1);
xy = pos(:,2:3);
cs_times = animal.CS_times.(['CS_' dateStr]);

% Build task mask (0–win(2) after CS)
taskMask = false(size(ts));
for i = 1:numel(cs_times)
    taskMask = taskMask | (ts >= cs_times(i) & ts <= cs_times(i) + win(2));
end

taskTimes = ts(taskMask);
nonTaskTimes = ts(~taskMask);

for ni = 1:numel(neuronVec)
    neuronIdx = neuronVec(ni);
    spikes = animal.Ca_peaks.(['CA_peaks_' dateStr])(neuronIdx, :);
    spikes = spikes(~isnan(spikes));

    % Task spikes
    taskSpikes = spikes(ismembertol(spikes, taskTimes, 1e-3));

    % Non-task spikes: remove task first by tol, then check against non-task
    nonTaskCandidates = spikes(~ismembertol(spikes, taskSpikes, 1e-3));
    nonTaskSpikes = nonTaskCandidates(ismembertol(nonTaskCandidates, nonTaskTimes, 1e-3));

    % Position lookup
    [~, taskIdx] = min(abs(taskSpikes' - ts'), [], 2);
    taskPos = xy(taskIdx, :);
    [~, nonTaskIdx] = min(abs(nonTaskSpikes' - ts'), [], 2);
    nonTaskPos = xy(nonTaskIdx, :);

    % === Plot ===
    figure('Color','w','Position',[300 300 600 500]); hold on;

    % Trajectory
    plot(xy(:,1), xy(:,2), 'Color', [0.8 0.8 0.8]);

    % Convex hull from task spikes
    if size(taskPos,1) >= 3
        k1 = boundary(taskPos(:,1), taskPos(:,2));
        fill(taskPos(k1,1), taskPos(k1,2), [0.5 0.2 0.8], ...
             'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end

    % Convex hull from non-task spikes
    if size(nonTaskPos,1) >= 3
        k2 = boundary(nonTaskPos(:,1), nonTaskPos(:,2));
        fill(nonTaskPos(k2,1), nonTaskPos(k2,2), [0.2 0.6 0.7], ...
             'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end

    % Plot spikes (task on top to avoid hiding)
    scatter(nonTaskPos(:,1), nonTaskPos(:,2), 20, ...
        'MarkerFaceColor', [0.2 0.6 0.7], 'MarkerEdgeColor', 'none');

    scatter(taskPos(:,1), taskPos(:,2), 20, ...
        'MarkerFaceColor', [0.5 0.2 0.8], 'MarkerEdgeColor', 'none');

    % Labeling
    legend({'Trajectory', 'Task Hull', 'Non-task Hull', 'Non-task Spikes', 'Task Spikes'}, ...
           'Location', 'northeastoutside');
    title(sprintf('Neuron %d — %s', neuronIdx, strrep(dateStr, '_','-')));
    xlabel('X position'); ylabel('Y position');
    axis equal; box off;

    % Report overlap check
    n_overlap = sum(ismembertol(taskSpikes, nonTaskSpikes, 1e-3));
    if n_overlap > 0
        fprintf('[Neuron %d] WARNING: %d spikes assigned to both task and non-task!\n', ...
                neuronIdx, n_overlap);
    end
end
end
