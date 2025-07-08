function plotTaskVsNonTaskSpikes(animal, neuronIdx, dateStr)
% Visualizes spike positions during vs. outside of task+2s for a neuron
% Task-period spikes shown as colored Xs (red → blue by phase), others as gray dots

    % Extract data
    pos = animal.pos.(['pos_' dateStr]);               % [n x 3] time, x, y
    spikeTimes = animal.Ca_peaks.(['CA_peaks_' dateStr])(neuronIdx,:); % spike timestamps
    cs_times = animal.CS_times.(['CS_' dateStr]);       % CS onsets

    ts = pos(:,1); x = pos(:,2); y = pos(:,3);          % trajectory data

    % Get task-period spikes (0–2s after each CS)
    task_spikes = [];
    trial_phase = [];
    for i = 1:length(cs_times)
        t0 = cs_times(i);
        t1 = t0 + 2;
        in_trial = spikeTimes >= t0 & spikeTimes <= t1;
        task_spikes = [task_spikes, spikeTimes(in_trial)];
        trial_phase = [trial_phase, (spikeTimes(in_trial) - t0)/2];  % normalize 0–1
    end

    % Non-task spikes = not within any task window
    is_task = false(size(spikeTimes));
    for i = 1:length(cs_times)
        is_task = is_task | (spikeTimes >= cs_times(i) & spikeTimes <= cs_times(i) + 2);
    end
    non_task_spikes = spikeTimes(~is_task);

    % Interpolate positions
    task_pos = interp1(ts, [x y], task_spikes, 'linear', 'extrap');
    non_task_pos = interp1(ts, [x y], non_task_spikes, 'linear', 'extrap');

    % === Plot ===
    figure('Color','w'); hold on;
    scatter(non_task_pos(:,1), non_task_pos(:,2), 14, [0.6 0.6 0.6], 'filled');  % gray dots
    scatter(task_pos(:,1), task_pos(:,2), 45, trial_phase, 'x', 'LineWidth', 2.5);

    % Reverse colormap: red (early) → blue (late)
    colormap(flipud(parula));
    cb = colorbar;
    cb.Label.String = 'Trial Phase (0 = CS, 1 = +2s)';
    caxis([0 1]);

    title(sprintf('Neuron %d: Task vs Non-Task Spiking (%s)', neuronIdx, strrep(dateStr,'_','-')));
    xlabel('X'); ylabel('Y');
    legend({'Non-task spikes', 'Task+2s spikes'}, 'Location', 'best');
    axis equal tight;
end
