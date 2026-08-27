function f = plotEMGVel()


ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};  %


spreadType = 'std';  % 'std' or 'sem' to pass to align + plot

velMat = {};
emgMat = {};
tMat   = {};

% --- layout: nRats + 1 panel for summary ---
nRats = numel(ratNames);
nPanels = nRats + 1;

figure('Color','w');
tl = tiledlayout(1,6,'TileSpacing','compact','Padding','compact');  % adjust grid as needed

ax = gobjects(1,nPanels);
for k = 1:nPanels
    ax(k) = nexttile(tl,k);
end

% --- per-rat plots + collect traces ---
for i = 1:nRats
    [t_interp, mean_emg, spread_emg, mean_vel, spread_vel] = alignEMGandVelocity(ratNames{i}, spreadType);

    plotEMGandVelocity_PerAnimal(t_interp, mean_emg, spread_emg, mean_vel, spread_vel, ratNames{i}, spreadType, 'Axes', ax(i), 'MakeFigure', false);

    % store for summary (row vectors)
    tMat{i}   = t_interp(:).';
    emgMat{i} = mean_emg(:).';
    velMat{i} = mean_vel(:).';
end

% --- build pooled summary on a common time grid ---
[t_common, emg_all, emg_spread, vel_all, vel_spread] = poolTracesToCommonGrid(tMat, emgMat, velMat, spreadType);

% --- plot summary into last axis ---
plotEMGandVelocity_PerAnimal(t_common, emg_all, emg_spread, vel_all, vel_spread, ...
    'All rats', spreadType, 'Axes', ax(end), 'MakeFigure', false);

title(tl, sprintf('EMG + speed aligned to CS (%s)', upper(spreadType)), 'Interpreter','none');
