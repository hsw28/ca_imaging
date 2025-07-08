function [shuffMeans, shuffSEMs] = plotRateMapSimilaritySummary(nShuffle)
  %shuffles spike label assignments
if nargin < 1, nShuffle = 500; end

  %metrics = {'Pearson', 'Spearman', 'Cosine'};
ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};
nRats = numel(ratNames);
metrics = {'Pearson', 'Spearman', 'Cosine'};

nMetrics = numel(metrics);

actualMeans = nan(nRats, nMetrics);
shuffMeans = nan(nRats, nMetrics);
actualSEMs = nan(nRats, nMetrics);
shuffSEMs = nan(nRats, nMetrics);
pvals = nan(nRats, nMetrics);

for r = 1:nRats
    rat = evalin('base', ratNames{r});
    dateList = autoDateList(rat);
    idx = find(strcmp(dateList, rat.An));
    theseDays = dateList(idx-2:idx);

    actualCorrs = [];
    shuffCorrs = [];

    for d = 1:3
        dateStr = theseDays{d};
        spikesMat = rat.Ca_peaks.(['CA_peaks_' dateStr]);
        nNeurons = size(spikesMat,1);
        pos = rat.pos.(['pos_' dateStr]);
        cs = rat.CS_times.(['CS_' dateStr]);
        ts = pos(:,1);

        taskMask = false(size(ts));
        for i = 1:numel(cs)
            taskMask = taskMask | (ts >= cs(i) & ts <= cs(i) + 2);
        end
        taskTS = ts(taskMask);
        nonTaskTS = ts(~taskMask);

        for ni = 1:nNeurons
            spikes = spikesMat(ni,:); spikes = spikes(~isnan(spikes));
            if numel(spikes) < 5, continue; end

            taskSpikes = spikes(ismembertol(spikes, taskTS, 1e-3));
            nonTaskSpikes = setdiff(spikes, taskSpikes, 'stable');
            if numel(taskSpikes) < 3 || numel(nonTaskSpikes) < 3, continue; end

            % Actual
            R1 = CA_normalizePosData(taskSpikes(:), pos, 2.5, 1);
            R2 = CA_normalizePosData(nonTaskSpikes(:), pos, 2.5, 1);
            R1 = imgaussfilt(R1, 0.75); R2 = imgaussfilt(R2, 0.75);
            R1(R1==0)=NaN; R2(R2==0)=NaN;
            v1 = R1(:); v2 = R2(:);
            valid = ~isnan(v1) & ~isnan(v2);
            if nnz(valid) < 10, continue; end
            v1 = v1(valid); v2 = v2(valid);

            actual = [corr(v1, v2, 'type','Pearson'), ...
                      corr(v1, v2, 'type','Spearman'), ...
                      dot(v1,v2)/(norm(v1)*norm(v2))];
            actualCorrs(end+1,:) = actual;

            % Shuffled
            allSpikes = [taskSpikes(:); nonTaskSpikes(:)];
            labels = [ones(numel(taskSpikes),1); zeros(numel(nonTaskSpikes),1)];

            shuffThis = nan(nShuffle, nMetrics);
            parfor s = 1:nShuffle
                perm = labels(randperm(numel(labels)));
                shuff1 = allSpikes(perm==1);
                shuff2 = allSpikes(perm==0);
                if numel(shuff1) < 3 || numel(shuff2) < 3, continue; end

                R1s = CA_normalizePosData(shuff1(:), pos, 2.5, 1);
                R2s = CA_normalizePosData(shuff2(:), pos, 2.5, 1);
                R1s = imgaussfilt(R1s, 0.75); R2s = imgaussfilt(R2s, 0.75);
                R1s(R1s==0)=NaN; R2s(R2s==0)=NaN;
                vs1 = R1s(:); vs2 = R2s(:);
                validS = ~isnan(vs1) & ~isnan(vs2);
                if nnz(validS) < 10, continue; end
                vs1 = vs1(validS); vs2 = vs2(validS);
                shuffThis(s,:) = [corr(vs1, vs2, 'type','Pearson'), ...
                                  corr(vs1, vs2, 'type','Spearman'), ...
                                  dot(vs1,vs2)/(norm(vs1)*norm(vs2))];
            end
            shuffCorrs(end+1,:) = nanmean(shuffThis,1);
        end
    end

    if isempty(actualCorrs), continue; end
    actualMeans(r,:) = mean(actualCorrs,1);
    shuffMeans(r,:) = mean(shuffCorrs,1);
    actualSEMs(r,:) = std(actualCorrs,0,1)/sqrt(size(actualCorrs,1));
    shuffSEMs(r,:) = std(shuffCorrs,0,1)/sqrt(size(shuffCorrs,1));

    for m = 1:nMetrics
        if size(actualCorrs,1) >= 3
          pval = mean(shuffCorrs(:,m) >= actualMeans(m), 'omitnan');
          pvals(r,m) = pval;
        end
    end
end

actualMeans;
actualSEMs;

% === Plotting ===
figure('Color','w','Position',[200 300 1200 400]);
for m = 1:nMetrics
    subplot(1,3,m); cla; hold on;
    bData = [actualMeans(:,m), shuffMeans(:,m)];
    barHandle = bar(bData, 'grouped');
    barHandle(1).FaceColor = [0.3 0.5 0.9];  % Actual
    barHandle(2).FaceColor = [0.7 0.7 0.7];  % Shuffled

    x1 = barHandle(1).XEndPoints; x2 = barHandle(2).XEndPoints;
    errorbar(x1, bData(:,1), actualSEMs(:,m), 'k', 'LineStyle','none');
    errorbar(x2, bData(:,2), shuffSEMs(:,m), 'k', 'LineStyle','none');
    xticks(1:nRats); xticklabels(ratNames);
    ylabel('Spatial similarity');
    title(sprintf('%s Corr (Task vs Non-task)', metrics{m}));
    legend({'Actual','Shuffled'}, 'Location','southwest');
    set(gca, 'Box','off', 'YLim', [0 1]);

    for r = 1:nRats
        if pvals(r,m) < 0.001
            text(mean([x1(r), x2(r)]), max(bData(r,:))*1.05, '***', 'HorizontalAlignment','center');
        elseif pvals(r,m) < 0.01
            text(mean([x1(r), x2(r)]), max(bData(r,:))*1.05, '**', 'HorizontalAlignment','center');
        elseif pvals(r,m) < 0.05
            text(mean([x1(r), x2(r)]), max(bData(r,:))*1.05, '*', 'HorizontalAlignment','center');
        end
    end
end
