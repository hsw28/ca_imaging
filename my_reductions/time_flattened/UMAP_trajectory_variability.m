% ===========================
% UMAP Trajectory Variability by Outcome
% ===========================
%trial flattened
%visualize trial-by-trial trajectories with mean overlays and variability shading for correct and incorrect trials.

function UMAP_trajectory_variability(animalName, win, binSize)
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

minDayFrac = 0.000001;
neuronDayCounts = sum(G > 0, 2);
sharedNeurons = find(neuronDayCounts >= round(minDayFrac * nAligned));
maxNeurons = numel(sharedNeurons);
fprintf('Using %d shared neurons for trajectory embedding.\n', maxNeurons);

nBins = round(diff(win) / binSize);

trialVecs = {};
labels = [];
trialID = [];
trialBinIdx = [];
tid = 1;

for d = 1:nAligned
    dateStr = dateList{d};
    [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
    if isempty(X), continue; end

    sharedIdx = G(sharedNeurons, d);
    validShared = sharedIdx > 0 & sharedIdx <= size(X,1);
    sharedIdx = sharedIdx(validShared);

    for t = 1:size(X,3)
        trialMat = nan(maxNeurons, nBins);
        for n = 1:numel(sharedIdx)
            trialMat(n,:) = X(sharedIdx(n), :, t);
        end
        trialMat(isnan(trialMat)) = 0;
        for b = 1:nBins
            trialVecs{end+1,1} = trialMat(:,b);
            labels(end+1,1) = y(t);
            trialID(end+1,1) = tid;
            trialBinIdx(end+1,1) = b;
        end
        tid = tid + 1;
    end
end

X_all = cell2mat(cellfun(@(x) x(:)', trialVecs, 'UniformOutput', false));
X_all = zscore(X_all, 0, 2);

nComp = 3;
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.2, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.1, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 2, 'n_neighbors', 40, 'min_dist', 0.6, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 5, 'n_neighbors', 40, 'min_dist', 0.1, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.4, 'metric', 'cosine', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 30, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.6, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.4, 'metric', 'cosine', 'randomize', true, 'verbose', false);
%
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 50, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 50, 'min_dist', 0.3, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.3, 'metric', 'euclidean', 'randomize', true, 'verbose', false);




nBins = max(trialBinIdx);
nTrials = max(trialID);

colors = lines(2);

% Figure with STD and SEM side-by-side
figure;
tiledlayout(1,2);

% STD subplot
nexttile;
title('Trajectory Variability (STD)');
hold on;
view(3);
for t = 1:nTrials
    idx = trialID == t;
    emb = embedding(idx,:);
    if all(labels(idx)==0)
        plot3(emb(:,1), emb(:,2), emb(:,3), '-', 'Color', colors(1,:), 'LineWidth', 0.5);
    elseif all(labels(idx)==1)
        plot3(emb(:,1), emb(:,2), emb(:,3), '-', 'Color', colors(2,:), 'LineWidth', 0.5);
    end
end

mu0 = nan(nBins, nComp);
mu1 = nan(nBins, nComp);
std0 = nan(nBins, nComp);
std1 = nan(nBins, nComp);
sem0 = nan(nBins, nComp);
sem1 = nan(nBins, nComp);
for b = 1:nBins
    bin0 = embedding(trialBinIdx==b & labels==0,:);
    bin1 = embedding(trialBinIdx==b & labels==1,:);
    t = mean(bin0, 1);
    mu0(b,:) = t(1:3);
    t = mean(bin1, 1);
    mu1(b,:) = t(1:3);
    t = std(bin0, 0, 1);
    std0(b,:) = t(1:3);
    t = std(bin1, 0, 1);
    std1(b,:) = t(1:3);
    t = std0(b,:) / sqrt(size(bin0,1));
    sem0(b,:) = t(1:3);
    t = std1(b,:) / sqrt(size(bin1,1));
    sem1(b,:) = t(1:3);
end

plot3(mu0(:,1), mu0(:,2), mu0(:,3), '-', 'Color', colors(1,:), 'LineWidth', 3);
plot3(mu1(:,1), mu1(:,2), mu1(:,3), '-', 'Color', colors(2,:), 'LineWidth', 3);
for b = 1:nBins
    [sx, sy, sz] = sphere;
    surf(std0(b,1)*sx + mu0(b,1), std0(b,2)*sy + mu0(b,2), std0(b,3)*sz + mu0(b,3), 'FaceColor', colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    surf(std1(b,1)*sx + mu1(b,1), std1(b,2)*sy + mu1(b,2), std1(b,3)*sz + mu1(b,3), 'FaceColor', colors(2,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3');
grid on;

% SEM subplot
nexttile;
title('Trajectory Variability (SEM)');
hold on;
view(3);
for t = 1:nTrials
    idx = trialID == t;
    emb = embedding(idx,:);
    if all(labels(idx)==0)
        plot3(emb(:,1), emb(:,2), emb(:,3), '-', 'Color', colors(1,:), 'LineWidth', 0.5);
    elseif all(labels(idx)==1)
        plot3(emb(:,1), emb(:,2), emb(:,3), '-', 'Color', colors(2,:), 'LineWidth', 0.5);
    end
end
plot3(mu0(:,1), mu0(:,2), mu0(:,3), '-', 'Color', colors(1,:), 'LineWidth', 3);
plot3(mu1(:,1), mu1(:,2), mu1(:,3), '-', 'Color', colors(2,:), 'LineWidth', 3);
for b = 1:nBins
    [sx, sy, sz] = sphere;
    surf(sem0(b,1)*sx + mu0(b,1), sem0(b,2)*sy + mu0(b,2), sem0(b,3)*sz + mu0(b,3), 'FaceColor', colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    surf(sem1(b,1)*sx + mu1(b,1), sem1(b,2)*sy + mu1(b,2), sem1(b,3)*sz + mu1(b,3), 'FaceColor', colors(2,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3');
grid on;

legend('Incorrect Trials', 'Correct Trials');
end
