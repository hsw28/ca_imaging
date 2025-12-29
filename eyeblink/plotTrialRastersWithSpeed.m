function plotRastersWithSpeedCompactRasterLines(spikeTimesAll, trialCS, win, pos, neuronIDs)
% Overlay raster on semi-transparent speed heatmap.
% If a single neuron is requested, also plot:
%   (2) trial-aligned mean speed (±SEM)
%   (3) PSTH (mean ± SEM across trials)

% ---- speed ----
velData   = ca_velocity(pos);    % 2×T: [speed; absolute_time]
speed     = velData(1,:);
speed_ts  = velData(2,:);

% ---- common trial-aligned grid ----
dt         = median(diff(speed_ts));
aligned_ts = win(1):dt:win(2);          % centers for display
nT         = numel(aligned_ts);
nTrials    = numel(trialCS);

% ---- PSTH binning (NEW: independent of dt) ----
psthBin    = 0.05;                      % 20 ms bins (change as you like)
edgesP     = win(1):psthBin:win(2);
centersP   = edgesP(1:end-1) + psthBin/2;
nP         = numel(centersP);

singleCell = (numel(neuronIDs) == 1);

% PSTH binning (use dt as default bin)
binSize = dt;
edges   = [aligned_ts - binSize/2, aligned_ts(end) + binSize/2]; % length nT+1

for k = 1:numel(neuronIDs)
    neuronIdx = neuronIDs(k);

    % ---- spikes ----
    if iscell(spikeTimesAll)
        spkTimes = spikeTimesAll{neuronIdx};
    else
        spkTimes = spikeTimesAll(neuronIdx,:);
        spkTimes = spkTimes(~isnan(spkTimes));
    end

    % ---- speedMat: nTrials × nT ----
    speedMat = nan(nTrials, nT);
    for t = 1:nTrials
        rel_ts = speed_ts - trialCS(t);
        sel    = rel_ts >= win(1) & rel_ts <= win(2);
        if nnz(sel) >= 2
            speedMat(t,:) = interp1(rel_ts(sel), speed(sel), aligned_ts, 'linear', NaN);
        end
    end

    meanSpeed = mean(speedMat, 1, 'omitnan');
    nEffS     = sum(~isnan(speedMat), 1);
    semSpeed  = std(speedMat, 0, 1, 'omitnan') ./ sqrt(max(nEffS,1));

    % ---- PSTH: per-trial binned spikes -> Hz (using psthBin) ----
    spkRateMat = zeros(nTrials, nP);  % Hz

    for t = 1:nTrials
        relSpk = spkTimes - trialCS(t);
        relSpk = relSpk(relSpk >= win(1) & relSpk <= win(2));
        counts = histcounts(relSpk, edgesP);     % 1×nP
        spkRateMat(t,:) = counts ./ psthBin;     % Hz
    end

%    if nansum(spkRateMat(:)) < 1000
%      continue
%    end

    meanPSTH = mean(spkRateMat, 1, 'omitnan');
    semPSTH  = std(spkRateMat, 0, 1, 'omitnan') ./ sqrt(max(nTrials,1));

    % ---- figure layout ----
    if singleCell
        fig = figure('Color','w','Position',[200 150 950 900]);
        tl = tiledlayout(fig, 3, 1, 'TileSpacing','compact', 'Padding','compact');
        ax1 = nexttile(tl, 1);
    else
        figure('Color','w','Position',[200 200 800 600]);
        ax1 = gca;
    end

    % ========= Row 1: speed heatmap + raster =========
    axes(ax1); %#ok<LAXES>
    hImg = imagesc(aligned_ts, 1:nTrials, speedMat);
    set(gca, 'YDir','normal');
    colormap(parula);
    hImg.AlphaData = 0.5;
    cb = colorbar('EastOutside');
    ylabel(cb, 'Speed');

    hold on;
    for t = 1:nTrials
        relSpk = spkTimes - trialCS(t);
        sel    = relSpk >= win(1) & relSpk <= win(2);
        xSpk   = relSpk(sel);
        y0     = t - 0.4;
        y1     = t + 0.4;
        for xi = xSpk(:)'
            line([xi xi], [y0 y1], 'Color','k', 'LineWidth',2.5);
        end
    end
    line([0 0], [0.5 nTrials+0.5], 'Color','r', 'LineStyle','--', 'LineWidth',1);

    xlim(win);
    ylim([0.5 nTrials+0.5]);
    xlabel('Time from CS (s)');
    ylabel('Trial');
    title(sprintf('Neuron %d — raster over speed', neuronIdx));
    set(gca, 'Box','off');

    if singleCell
        % ========= Row 2: mean speed (±SEM) =========
        ax2 = nexttile(tl, 2);
        axes(ax2); %#ok<LAXES>
        hold on;

        x  = aligned_ts;
        lo = meanSpeed - semSpeed;
        hi = meanSpeed + semSpeed;

        % SEM shading (draw first, opaque enough to see)
        fill([x fliplr(x)], [lo fliplr(hi)], ...
             [0 0 0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        % Mean speed trace
        plot(x, meanSpeed, 'k', 'LineWidth', 2);

        % CS line
        yl = [min(lo) max(hi)];
        line([0 0], yl, 'Color','r', 'LineStyle','--', 'LineWidth',1);

        xlim(win);
        ylim(yl);              % <-- forces SEM to be visible
        xlabel('Time from CS (s)');
        ylabel('Speed (mean ± SEM)');
        title('Trial-aligned average speed');
        set(gca, 'Box','off');


        % ========= Row 3: PSTH as bars (small bins) =========
        ax3 = nexttile(tl, 3);
        axes(ax3); %#ok<LAXES>
        hold on;

        bar(centersP, meanPSTH, 1.0);   % full PSTH-bin width (since bins are now small)

        % optional SEM error bars
    %    errorbar(centersP, meanPSTH, semPSTH, 'k', 'LineStyle','none', ...
    %             'LineWidth',1, 'CapSize',0);

        yl = ylim;
        line([0 0], yl, 'Color','r', 'LineStyle','--', 'LineWidth',1);

        xlim(win);
        xlabel('Time from CS (s)');
        ylabel(sprintf('Rate (Hz), bin=%.3gs', psthBin));
        title('PSTH across trials');
        set(gca, 'Box','off');
    end
end
end
