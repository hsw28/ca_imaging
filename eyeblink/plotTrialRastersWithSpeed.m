function plotRastersWithSpeedCompactRasterLines(spikeTimesAll, trialCS, win, pos, neuronIDs)
% plotRastersWithSpeedCompactRasterLines
% Overlay traditional raster lines on a semi-transparent speed heatmap.
%
% INPUTS:
%   spikeTimesAll : cell {nNeurons×1} of spike-time vectors or numeric [nNeurons×?]
%   trialCS       : [1×nTrials] CS onset times (absolute)
%   win           : [1 2] window around CS, e.g. [-1 2]
%   pos           : whatever your ca_velocity() expects
%   neuronIDs     : vector of neuron indices to plot

% 1) get continuous speed
velData   = ca_velocity(pos);    % 2×T: [speed; absolute_time]
speed     = velData(1,:);
speed_ts  = velData(2,:);

% 2) define a uniform time grid
dt         = median(diff(speed_ts));
aligned_ts = win(1):dt:win(2);
nT         = numel(aligned_ts);
nTrials    = numel(trialCS);

for k = 1:numel(neuronIDs)
    neuronIdx = neuronIDs(k);
    % extract spikes for this neuron
    if iscell(spikeTimesAll)
        spkTimes = spikeTimesAll{neuronIdx};
    else
        spkTimes = spikeTimesAll(neuronIdx,:);
    end

    % build speed matrix (nTrials×nT)
    speedMat = nan(nTrials, nT);
    for t = 1:nTrials
        rel_ts = speed_ts - trialCS(t);
        sel    = rel_ts>=win(1) & rel_ts<=win(2);
        speedMat(t,:) = interp1(rel_ts(sel), speed(sel), aligned_ts, 'linear', NaN);
    end

    % —— plot heatmap ——
    figure('Color','w','Position',[200 200 800 600]);
    hImg = imagesc(aligned_ts, 1:nTrials, speedMat);
    set(gca, 'YDir','normal');
    colormap(parula);
    hImg.AlphaData = 0.5;            % 40% transparent
    cb = colorbar('EastOutside');
    ylabel(cb, 'Speed');

    hold on;
    % —— overlay raster lines ——
    for t = 1:nTrials
        relSpk = spkTimes - trialCS(t);
        sel    = relSpk>=win(1) & relSpk<=win(2);
        xSpk   = relSpk(sel);
        y0     = t - 0.4;
        y1     = t + 0.4;
        % draw each spike as a thick vertical line
        for xi = xSpk
            line([xi xi], [y0 y1], 'Color', 'k', 'LineWidth', 2.5);
        end
    end

    % zero‐line
    line([0 0], [0.5 nTrials+0.5], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);

    xlim(win);
    ylim([0.5 nTrials+0.5]);
    xlabel('Time from CS (s)');
    ylabel('Trial');
    title(sprintf('Neuron %d — raster over speed', neuronIdx));
    set(gca, 'Box', 'off');
end
end
