function plotEMG_withDoubleYBreak(t, y, midRange, highClip, lowClip)
% t: time vector (1 x T)
% y: EMG trace (1 x T)
% midRange: [low high] for main axis (e.g. [-5000 5000])
% highClip: upper y-limit for top axis (e.g. 12000)
% lowClip: lower y-limit for bottom axis (e.g. -12000)

% Create 3 stacked axes with shared x-axis
fig = figure('Color', 'w', 'Position', [300 300 600 500]);
tiledlayout(3,1, 'TileSpacing', 'none', 'Padding', 'compact');

% === Top axis (high spike range)
nexttile; hold on;
plot(t, y, 'Color', [0.6 0.1 0.8], 'LineWidth', 1);
ylim([midRange(2), highClip]);
set(gca, 'XTick', [], 'Box', 'off');
ylabel('High');

% Trace period shading
yl = ylim;
fill([0 0.75 0.75 0], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% === Middle axis (main view)
nexttile; hold on;
plot(t, y, 'Color', [0.6 0.1 0.8], 'LineWidth', 1);
ylim(midRange);
set(gca, 'XTick', [], 'Box', 'off');
ylabel('EMG (a.u.)');

% Trace period shading
yl = ylim;
fill([0 0.75 0.75 0], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% === Bottom axis (low spike range)
nexttile; hold on;
plot(t, y, 'Color', [0.6 0.1 0.8], 'LineWidth', 1);
ylim([lowClip, midRange(1)]);
xlabel('Time from CS onset (s)');
ylabel('Low');
set(gca, 'Box', 'off');

% Trace period shading
yl = ylim;
fill([0 0.75 0.75 0], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

sgtitle('EMG with Dual Y-Axis Break');
end
