% ========================
% UMAP Trajectories by Outcome and Significance
% ========================

function UMAP_trajectory_significance(animalName, win, binSize)
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

minDayFrac = 0.5;
neuronDayCounts = sum(G > 0, 2);
sharedNeurons = find(neuronDayCounts >= round(minDayFrac * nAligned));
maxNeurons = numel(sharedNeurons);
fprintf('Using %d shared neurons for trajectory embedding.\n', maxNeurons);

nBins = round(diff(win) / binSize);

trialVecs = {};
labels = [];
trialDayIdx = [];
trialBinIdx = [];

for d = 1:nAligned
    dateStr = dateList{d};
    [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
    if isempty(X), continue; end

    sharedIdx = G(sharedNeurons, d);
    validShared = sharedIdx > 0 & sharedIdx <= size(X,1);
    sharedIdx = sharedIdx(validShared);
    neuronsUsed = sharedNeurons(validShared);

    for t = 1:size(X,3)
        trialMat = nan(maxNeurons, nBins);
        for n = 1:numel(sharedIdx)
            trialData = X(sharedIdx(n), :, t);
            if ~isempty(trialData)
                trialMat(n,:) = trialData;
            end
        end
        trialMat(isnan(trialMat)) = 0;
        for b = 1:nBins
            trialVecs{end+1,1} = trialMat(:,b);
            labels(end+1,1) = y(t);
            trialDayIdx(end+1,1) = d;
            trialBinIdx(end+1,1) = b;
        end
    end
end

if isempty(trialVecs)
    error('No valid trial data to embed.');
end

X_all = cell2mat(cellfun(@(x) x(:)', trialVecs, 'UniformOutput', false));
X_all = zscore(X_all, 0, 2);

[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.3, ...
    'metric', 'cosine', 'randomize', true, 'verbose', false);

if isempty(embedding)
    error('UMAP embedding failed.');
end

nDims = size(embedding, 2);
nTrials = length(labels) / nBins;

embed0 = reshape(embedding(labels==0,:), nBins, [], nDims);
embed1 = reshape(embedding(labels==1,:), nBins, [], nDims);

if size(embed0,2) < 2 || size(embed1,2) < 2
    warning('Insufficient trials for one or both classes.');
    return;
end

mu0 = squeeze(mean(embed0, 2));  % [nBins x nDims]
mu1 = squeeze(mean(embed1, 2));

dists = vecnorm(mu1 - mu0, 2, 2);  % Euclidean distance per bin

% Permutation test
nPerms = 1000;
permDists = nan(nBins, nPerms);
for p = 1:nPerms
    yperm = labels(randperm(length(labels)));
    e0 = reshape(embedding(yperm==0,:), nBins, [], nDims);
    e1 = reshape(embedding(yperm==1,:), nBins, [], nDims);
    if size(e0,2) < 2 || size(e1,2) < 2
        continue;
    end
    pm0 = squeeze(mean(e0, 2));
    pm1 = squeeze(mean(e1, 2));
    permDists(:,p) = vecnorm(pm1 - pm0, 2, 2);
end

pvals = mean(permDists >= dists, 2, 'omitnan');
sigBins = find(pvals < 0.05);

% Plot separation and p-values
figure;
subplot(2,1,1);
plot(1:nBins, dists, 'k-', 'LineWidth', 2); hold on;
if ~isempty(sigBins)
    sigBins = sigBins(sigBins <= numel(dists));
    scatter(sigBins, dists(sigBins), 60, 'r', 'filled');
end
xlabel('Time Bin'); ylabel('Distance'); grid on;
title('Separation between correct and incorrect trajectories');

subplot(2,1,2);
plot(1:nBins, pvals, 'b-', 'LineWidth', 2); hold on;
yline(0.05, '--r');
xlabel('Time Bin'); ylabel('p-value');
title('p-values per time bin'); grid on;

% Average trajectory plot with arrows
if nDims >= 3
    figure; hold on;
    for b = 1:(nBins-1)
        quiver3(mu0(b,1), mu0(b,2), mu0(b,3), mu0(b+1,1)-mu0(b,1), mu0(b+1,2)-mu0(b,2), mu0(b+1,3)-mu0(b,3), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 4);
        quiver3(mu1(b,1), mu1(b,2), mu1(b,3), mu1(b+1,1)-mu1(b,1), mu1(b+1,2)-mu1(b,2), mu1(b+1,3)-mu1(b,3), 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 4);
    end
    title('Average Trajectories in UMAP 3D Space');
    legend('Incorrect', 'Correct');
    xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3'); grid on;
elseif nDims >= 2
    figure; hold on;
    for b = 1:(nBins-1)
        quiver(mu0(b,1), mu0(b,2), mu0(b+1,1)-mu0(b,1), mu0(b+1,2)-mu0(b,2), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 4);
        quiver(mu1(b,1), mu1(b,2), mu1(b+1,1)-mu1(b,1), mu1(b+1,2)-mu1(b,2), 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 4);
    end
    title('Average Trajectories in UMAP 2D Space');
    legend('Incorrect', 'Correct');
    xlabel('UMAP 1'); ylabel('UMAP 2'); grid on;
else
    warning('Not enough dimensions to plot average trajectories.');
end

end
