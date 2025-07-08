function f = velEMG_examples(t_interp, vel_trials, emg_trials)
% After function:
% [t_interp, mean_emg, spread_emg, mean_vel, spread_vel, vel_trials, emg_trials] = alignEMGandVelocity(ratName, spreadType)

example_idx = [1 3 24 32 39 46];
nTrials = length(example_idx);

f = figure('Color', 'w', 'Position', [300 300 700 100*nTrials]);

for i = 1:nTrials
    trial = example_idx(i);

    % === Velocity subplot ===
    if i == 1
      ax_vel = subplot(nTrials, 2, 3); hold on;
    elseif i == 2
      ax_vel = subplot(nTrials, 2, 4); hold on;
    elseif i == 3
      ax_vel = subplot(nTrials, 2, 7); hold on;
    elseif i == 4
      ax_vel = subplot(nTrials, 2, 8); hold on;
    elseif i == 5
      ax_vel = subplot(nTrials, 2, 11); hold on;
    elseif i == 6
      ax_vel = subplot(nTrials, 2, 12); hold on;
    end

    plot(t_interp, vel_trials(trial,:), 'Color', [0.1 0.6 0.9], 'LineWidth', 1.5);
    yL = ylim;
    fill([0 0.75 0.75 0], [yL(1) yL(1) yL(2) yL(2)], ...
         [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
    ylabel('Speed (cm/s)');
    if i < 6, set(gca, 'XTickLabel', []); else, xlabel('Time from CS onset (s)'); end

      scale = 1000;  % adjust to compress more/less

    % === EMG (top for this trial) ===
    if i == 1
      ax_emg = subplot(nTrials, 2, 1); hold on;
    elseif i == 2
      ax_emg = subplot(nTrials, 2, 2); hold on;
    elseif i == 3
      ax_emg = subplot(nTrials, 2, 5); hold on;
    elseif i == 4
      ax_emg = subplot(nTrials, 2, 6); hold on;
    elseif i == 5
      ax_emg = subplot(nTrials, 2, 9); hold on;
    elseif i == 6
      ax_emg = subplot(nTrials, 2, 10); hold on;
    end

    % Transform EMG using asinh scale
    y_scaled = asinh(emg_trials(trial,:) / scale);
    plot(t_interp, y_scaled, 'Color', [0.6 0.1 0.8], 'LineWidth', 1.5);

    % Shade trace period
    yL = ylim;
    fill([0 0.75 0.75 0], [yL(1) yL(1) yL(2) yL(2)], ...
         [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
    ylabel('asinh-scaled EMG');

    % Custom y-axis ticks
    ytickvals = [-10000 -5000 0 5000 10000];
    set(gca, 'YTick', asinh(ytickvals / scale), ...
             'YTickLabel', arrayfun(@num2str, ytickvals, 'UniformOutput', false));

    title(sprintf('Trial %d', trial));
    if i < 6
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time from CS onset (s)');
    end
    set(gca, 'Box', 'off', 'FontSize', 10);

end
