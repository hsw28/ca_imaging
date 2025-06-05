function UMAP_trajectory_variability_evolution(animalName, win, binSize)
% =======================================
% UMAP Trajectory Variability Evolution
% (within-class variability across days)
% =======================================

%% Setup
Fs = 7.5;
animal = evalin('base', animalName);
dateList = autoDateList(animal);
G = animal.alignmentALL;
nDays = numel(dateList);
nAligned = size(G, 2);
if nAligned < nDays
    fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
    dateList = dateList(1:nAligned);
    nDays = nAligned;
end

minDayFrac = 0.00001;
neuronDayCounts = sum(G > 0, 2);
sharedNeurons = find(neuronDayCounts >= round(minDayFrac * nAligned));
fprintf('Using %d shared neurons for trajectory embedding.\n', numel(sharedNeurons));

nBins = round(diff(win) / binSize);

%% Assemble trial-by-bin matrix
trialVecs = {};
labels = [];
trialID = [];
trialBinIdx = [];
trialDayIdx = [];
tid = 1;

for d = 1:nAligned
    dateStr = dateList{d};
    [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
    if isempty(X), continue; end

    sharedIdx = G(sharedNeurons, d);
    validShared = sharedIdx > 0 & sharedIdx <= size(X,1);
    sharedIdx = sharedIdx(validShared);

    for t = 1:size(X,3)
        trialMat = nan(numel(sharedNeurons), nBins);
        for n = 1:numel(sharedIdx)
            trialMat(n,:) = X(sharedIdx(n), :, t);
        end
        trialMat(isnan(trialMat)) = 0;
        for b = 1:nBins
            trialVecs{end+1,1} = trialMat(:,b);
            labels(end+1,1) = y(t);
            trialID(end+1,1) = tid;
            trialBinIdx(end+1,1) = b;
            trialDayIdx(end+1,1) = d;
        end
        tid = tid + 1;
    end
end

%% UMAP embedding (global)
X_all = cell2mat(cellfun(@(x) x(:)', trialVecs, 'UniformOutput', false));
X_all = zscore(X_all, 0, 2);

%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 50, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 50, 'min_dist', 0.3, 'metric', 'euclidean', 'randomize', true, 'verbose', false);


[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.3, 'metric', 'euclidean', 'randomize', true, 'verbose', false);

%% Compute within-class variability
nBins = max(trialBinIdx);
nClasses = 2;

variability = nan(nDays, nBins, nClasses);

for d = 1:nAligned
    dayIdx = trialDayIdx == d;

    for b = 1:nBins
        for cls = 0:1
            classIdx = dayIdx & trialBinIdx==b & labels==cls;
            if sum(classIdx) < 2
                variability(d,b,cls+1) = NaN;
                continue;
            end
            mu = mean(embedding(classIdx,:), 1);
            dists = sqrt(sum((embedding(classIdx,:) - mu).^2, 2));
            variability(d,b,cls+1) = mean(dists); % you can also use std(dists)
        end
    end
end

%% Plot: variability vs time bin (averaged across days)
meanVar = squeeze(nanmean(variability,1)); % [bins x class]

figure;
hold on;
plot(1:nBins, meanVar(:,1), '-o', 'LineWidth', 2, 'Color', [0.7 0 0]);
plot(1:nBins, meanVar(:,2), '-o', 'LineWidth', 2, 'Color', [0 0 0.6]);
legend('Incorrect', 'Correct');
xlabel('Time Bin'); ylabel('Mean Within-Class Variability');
title('Within-Class Variability vs Time Bin (Avg Across Days)');
grid on;

%% Plot: variability across days
figure;
tiledlayout(ceil(nDays/5),5);

vmax = max(variability(:), [], 'omitnan');

for d = 1:nAligned
    nexttile;
    plot(1:nBins, squeeze(variability(d,:,1)), '-o', 'Color', [0.7 0 0], 'LineWidth', 2);
    hold on;
    plot(1:nBins, squeeze(variability(d,:,2)), '-o', 'Color', [0 0 0.6], 'LineWidth', 2);
    ylim([0 vmax*1.05]);
    title(sprintf('Day %s', dateList{d}));
    xlabel('Time Bin'); ylabel('Variability');
end
sgtitle('Within-Class Variability Across Days (Global UMAP)');

%% Plot: variability trend over days
meanVarDay = squeeze(nanmean(variability,2)); % [day x class]

figure;
hold on;
plot(1:nDays, meanVarDay(:,1), '-o', 'LineWidth', 2, 'Color', [0.7 0 0]);
plot(1:nDays, meanVarDay(:,2), '-o', 'LineWidth', 2, 'Color', [0 0 0.6]);
set(gca, 'XTick', 1:nDays, 'XTickLabel', dateList);
xlabel('Day'); ylabel('Mean Variability Across Bins');
title('Within-Class Variability Evolution Across Days');
grid on;
xtickangle(45);
legend('Incorrect', 'Correct');

end
