function plotTraceRasterAndCalcium(spikeTimesAll, trialCS, win, calciumAll, calcium_ts, neuronIDs, part_of_diff_fig)
% Plots raster + calcium trace for each neuron in neuronIDs
% Inputs:
%   spikeTimesAll: matrix [nNeurons x nTimepoints], binary or inferred spikes
%   trialCS: vector of CS onset times
%   win: 1x2 time window around CS (e.g., [-1 2])
%   calciumAll: [nNeurons x nTrials x nTimepoints] calcium activity
%   calcium_ts: time vector for calcium trace (in seconds)
%   neuronIDs: vector of neuron indices to plot

nNeurons = numel(neuronIDs);
nTrials = numel(trialCS);

% Ensure calcium_ts is column vector in seconds
if size(calcium_ts,2) > 1
    calcium_ts = calcium_ts(:,2) ./ 1000;
    calcium_ts = calcium_ts(2:2:end);
end


traceStart = 0;
traceEnd = 0.75;

plotsPerFig = 3;
figCounter = 0;

pi = 0;
for ni = 1:nNeurons
    neuronIdx = neuronIDs(ni);
    spikeTimes = spikeTimesAll(neuronIdx,:);
    % Preallocate


    % Determine number of timepoints in the calcium trace per trial
    Fs = 1 / median(diff(calcium_ts));  % estimate sampling rate
    nTimepoints = round((win(2) - win(1)) * Fs);
    aligned_ts = linspace(win(1), win(2), nTimepoints);  % relative to CS

    % Preallocate
    calciumAligned = nan(nTrials, nTimepoints);
    caTrace = calciumAll(neuronIdx, :);  % [1 x time]


    for i = 1:nTrials

        cs = trialCS(i);
        % Find indices where calcium_ts falls within [cs + win(1), cs + win(2)]
        mask = calcium_ts >= cs + win(1) & calcium_ts <= cs + win(2);
        segment = caTrace(mask);


        % Handle mismatch in length (e.g., at end of recording)
        if length(segment) == nTimepoints
            calciumAligned(i,:) = segment;
        elseif length(segment) < nTimepoints && length(segment) > 0
            calciumAligned(i,1:length(segment)) = segment;
        end
    end


    % Align spikes to CS
    alignedSpikes = cell(nTrials,1);
    for i = 1:nTrials
        cs = trialCS(i);
        relSpikes = spikeTimes(spikeTimes >= cs + win(1) & spikeTimes <= cs + win(2)) - cs;
        alignedSpikes{i} = relSpikes;
    end

    % Count total spikes
    numspikes = sum(cellfun(@numel, alignedSpikes));
    if numspikes < 5
        continue
    else
      fprintf('Cell %d has %d spikes\n', neuronIdx, numspikes);

    end

    % increment plot index
    pi = pi + 1;

    % Start new figure every 6 neurons

    if nargin < 7
      if mod(pi-1, plotsPerFig) == 0
          figure('Color','w','Position',[300 100 1000 1000]);
      end
    end
    
    subplotRow = mod(pi-1, plotsPerFig) + 1;



    % --- Raster ---
    subplot(plotsPerFig, 3, (subplotRow-1)*3 + 1); hold on;

    % 💡 Draw trace-period shading first!
    fill([traceStart traceEnd traceEnd traceStart], ...
         [0 0 nTrials+1 nTrials+1], [0.9 0.9 0.9], ...
         'EdgeColor', 'none');

    % Now plot raster on top
    for i = 1:nTrials
        spks = alignedSpikes{i};
        if isempty(spks), continue; end
        y = i * ones(size(spks));
        line([spks; spks], [y - 0.4; y + 0.4], 'Color', 'k', 'LineWidth', 2);

    end

    line([0 0], ylim, 'Color', 'r', 'LineStyle', '--');
    xlim(win);
    ylim([0.5 nTrials + 0.5]);
    ylabel('Trial');
    title(sprintf('Neuron %d Raster', neuronIdx));
    set(gca, 'Box', 'off');

    % --- PSTH ---

    % Bin spike times across all trials into a histogram
    binSize = 0.05;  % 50 ms
    edges = win(1):binSize:win(2);

    allSpikes = [];
    for i = 1:nTrials
        spks = alignedSpikes{i};
        if ~isempty(spks)
            allSpikes = [allSpikes; spks(:)];  % force column, concat
        end
    end
    spikeCounts = histcounts(allSpikes, edges);

    psth = spikeCounts / (nTrials * binSize);  % spikes/sec
    binCenters = edges(1:end-1) + binSize/2;

    subplot(plotsPerFig, 3, (subplotRow-1)*3 + 2); hold on;
    yL = [0, max(psth) * 1.1 + 1e-6];  % safe default

    fill([traceStart traceEnd traceEnd traceStart], ...
         [yL(1) yL(1) yL(2) yL(2)], [0.9 0.9 0.9], 'EdgeColor', 'none');

    bar(binCenters, psth, 1, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');

    if any(isnan(yL)) || diff(yL) <= 0
        yL = [0 1];  % fallback range
    end

    line([0 0], yL, 'Color', 'r', 'LineStyle', '--');
    xlim(win); ylim(yL);
    ylabel('Firing Rate (Hz)');
    title(sprintf('Neuron %d PSTH', neuronIdx));
    set(gca, 'Box', 'off');


    % --- Calcium ---

    subplot(plotsPerFig, 3, (subplotRow-1)*3 + 3); hold on;

    meanCal = nanmean(calciumAligned, 1);
    semCal = nanstd(calciumAligned, [], 1) ./ sqrt(sum(~isnan(calciumAligned(:,1))));

    %  FIX: Use aligned_ts, not calcium_ts
    if length(aligned_ts) ~= length(meanCal)
        warning('aligned_ts (%d) ≠ meanCal (%d)', length(aligned_ts), length(meanCal));
        continue
    end

    aligned_ts = aligned_ts(:)';   % force row vector
    meanCal = meanCal(:)';         % force row vector
    semCal  = semCal(:)';          % force row vector

    % Plot trace-period shading first (so it stays in background)
    yL = [min(meanCal - semCal), max(meanCal + semCal)];
    fill([traceStart traceEnd traceEnd traceStart], ...
         [yL(1) yL(1) yL(2) yL(2)], ...
         [0.9 0.9 0.9], 'EdgeColor', 'none');

    % Plot calcium with error band
    fill([aligned_ts fliplr(aligned_ts)], ...
         [meanCal+semCal fliplr(meanCal-semCal)], ...
         [0.6 0.1 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(aligned_ts, meanCal, 'Color', [0.6 0.1 0.8], 'LineWidth', 1.5);


    line([0 0], yL, 'Color', 'r', 'LineStyle', '--');
    xlim([aligned_ts(1), aligned_ts(end)]);

    ylim(yL);
    xlabel('Time from CS (s)');
    ylabel('Calcium (a.u.)');
    title(sprintf('Neuron %d Calcium', neuronIdx));
    set(gca, 'Box', 'off');





end
