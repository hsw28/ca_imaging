function f = plot_place_vs_CSUS(animal, neuronVec, dateStr, heatmapLOW, heatmapHIGH)



nNeurons = numel(neuronVec);
plotsPerFig = 3;
plotCounter = 0;

if nNeurons>20
CSUSspikes_lim = 35;
highspeedspikes_lim = 10;
else
  CSUSspikes_lim = 0;
  highspeedspikes_lim = 0;
end

for ni = 1:nNeurons
    neuronIdx = neuronVec(ni);
    pos = animal.pos.(['pos_' dateStr]);
    peaks_time = animal.Ca_peaks.(['CA_peaks_' dateStr])(neuronIdx,:);
    calcium_ts = animal.Ca_ts.(['CA_time_' dateStr]);
    spikeTimes = peaks_time;
    cs_times = animal.CS_times.(['CS_' dateStr]);
    CSUS_id = animal.CSUS_id.(['CSUS_id_' dateStr]);
    nTrials = numel(cs_times);
    postimes = pos(:,1);

    if size(calcium_ts,2) > 1
        calcium_ts = calcium_ts(:,2) ./ 1000;
        calcium_ts = calcium_ts(2:2:end);
    end

    win = [0,2];
    Fs = 1 / median(diff(calcium_ts));
    nTimepoints = round((win(2) - win(1)) * Fs);
    aligned_ts = linspace(win(1), win(2), nTimepoints);

    alignedSpikes = [];
    csusPosIdx = [];
    for i = 1:nTrials
        cs = cs_times(i);
        relSpikes = spikeTimes(spikeTimes >= cs + win(1) & spikeTimes <= cs + win(2));
        csusIdx = find(postimes >= cs + win(1) & postimes <= cs + win(2));
        alignedSpikes = [alignedSpikes, relSpikes];
        csusPosIdx = [csusPosIdx, csusIdx'];
    end

    CSUSspikes = alignedSpikes;

    velthreshold = 4;
    dim = 2.5;
    vel = ca_velocity(pos);
    goodvel = find(vel(1,:)>=velthreshold);
    goodtime = pos(goodvel, 1);
    goodpos = pos(goodvel,:);
    goodvel = setdiff(goodvel, csusPosIdx);

    highspeedspikes = [];
    for ii = 1:length(peaks_time)
        [minValue_vel,closestIndex] = min(abs(peaks_time(ii)-goodtime));
        if minValue_vel <= 1/15 && ~isnan(peaks_time(ii))
            highspeedspikes(end+1) = peaks_time(ii);
        end
    end

    highspeedspikes = setdiff(highspeedspikes, CSUSspikes);


    if length(CSUSspikes) <= CSUSspikes_lim || length(highspeedspikes) <= highspeedspikes_lim
        continue
    end

        figure('Color','w','Position',[300 300 1000 800]);
        subplotRow = mod(plotCounter, plotsPerFig) + 1;


    % --- Plot high-speed spikes ---
    rate = CA_normalizePosData(highspeedspikes, pos, dim, 1.000);
    rate(isnan(rate)) = 0;
    rate = imgaussfilt(rate, .75);

    subplot(plotsPerFig, 2, (subplotRow-1)*2 + 1);


    numrate = sort(rate(~isnan(rate)),'descend');

    maxratefive = min(numrate(1:ceil(length(numrate)*0.01)));  % this tempers extreme high outliers

    idx = find(numrate == 0, 1, 'first');
    if length(idx)<1
      idx = length(numrate)+1;
    end
    numrate = numrate(1:idx-1);
    minratefive = max(numrate(floor(length(numrate)*0.18):end));
    if nargin<4
      imagesc(rate, [minratefive, maxratefive*1]);
    else
      imagesc(rate, [heatmapLOW, heatmapHIGH]);
    end

    colormap('parula'); colorbar;
    title(sprintf('Neuron %d (high-speed)', neuronIdx));

    % --- Plot CSUS spikes ---
    rate = CA_normalizePosData(CSUSspikes, pos, dim, 1.000);
    rate(isnan(rate)) = 0;


    subplot(plotsPerFig, 2, (subplotRow-1)*2 + 2);


    numrate = sort(rate(~isnan(rate)),'descend');
    maxratefive = min(numrate(1:ceil(length(numrate)*0.01)));  % this tempers extreme high outliers

  %  idx = find(numrate == 0, 1, 'first')
  %  numrate = numrate(1:idx-1);
  %  minratefive = max(numrate(floor(length(numrate)*0.9):end));
    imagesc(rate, [mean(numrate), maxratefive]);
    colormap('parula'); colorbar;
    title(sprintf('Neuron %d (CSUS)', neuronIdx));

    plotCounter = plotCounter + 1;
end
end
