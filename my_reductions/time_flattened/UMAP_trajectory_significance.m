% ========================
% UMAP Trajectories by Outcome and Significance
% ========================
% Usage examples:
%   UMAP_trajectory_significance('rat0314', [0 .75], 1/7.5);
%   UMAP_trajectory_significance('rat0314', [0 .75], 1/7.5, 'perDayUMAP', true);

function UMAP_trajectory_significance(animalName, win, binSize, varargin)

%% Args
p = inputParser;
addParameter(p, 'perDayUMAP', false);
parse(p, varargin{:});
perDayUMAP = p.Results.perDayUMAP;

%% Basic setup
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

minDayFrac = 0.001;
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

%% Run UMAP on ALL days together
X_all = cell2mat(cellfun(@(x) x(:)', trialVecs, 'UniformOutput', false));
X_all = zscore(X_all, 0, 2);

%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.2, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.1, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 2, 'n_neighbors', 40, 'min_dist', 0.6, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 5, 'n_neighbors', 40, 'min_dist', 0.1, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%%%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.4, 'metric', 'cosine', 'randomize', true, 'verbose', false);
%%%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 30, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.6, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.4, 'metric', 'cosine', 'randomize', true, 'verbose', false);

%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 50, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 50, 'min_dist', 0.3, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.3, 'metric', 'euclidean', 'randomize', true, 'verbose', false);

%% Compute mean traj and separation (GLOBAL UMAP)
nBins = max(trialBinIdx);
colors = lines(2);

mu0 = nan(nBins, 3);
mu1 = nan(nBins, 3);
for b = 1:nBins
    mu0(b,:) = mean(embedding(trialBinIdx==b & labels==0,:), 1);
    mu1(b,:) = mean(embedding(trialBinIdx==b & labels==1,:), 1);
end

dists = sqrt(sum((mu0 - mu1).^2, 2));

%% Permutation test (GLOBAL UMAP)
nPerm = 1000;
permDists = nan(nPerm, nBins);
for p = 1:nPerm
    permLabel = labels(randperm(length(labels)));
    for b = 1:nBins
        permMu0 = mean(embedding(trialBinIdx==b & permLabel==0,:), 1);
        permMu1 = mean(embedding(trialBinIdx==b & permLabel==1,:), 1);
        permDists(p,b) = sqrt(sum((permMu0 - permMu1).^2));
    end
end
pVals = mean(permDists >= dists', 1);

%% Plot GLOBAL separation + p-values
figure;
subplot(2,1,1);
plot(1:nBins, dists, '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('Distance'); xlabel('Time Bin');
title('Separation between correct and incorrect trajectories (GLOBAL UMAP)');
subplot(2,1,2);
plot(1:nBins, pVals, '-ob', 'LineWidth', 2);
hold on; yline(0.05, 'r--');
xlabel('Time Bin'); ylabel('p-value');
title('p-values per time bin (GLOBAL UMAP)');

%% Plot GLOBAL mean trajectory
figure;
hold on;
plot3(mu0(:,1), mu0(:,2), mu0(:,3), '-o', 'Color', colors(1,:), 'LineWidth', 3);
plot3(mu1(:,1), mu1(:,2), mu1(:,3), '-o', 'Color', colors(2,:), 'LineWidth', 3);
legend('Incorrect', 'Correct');
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3');
grid on;
title('Average Trajectories in GLOBAL UMAP 3D Space');

%% Plot how GLOBAL separation changes across days
dayMeans = nan(nAligned,1);
for d = 1:nAligned
    % Subset for this day:
    dayIdx = (trialDayIdx == d);
    Y_day = labels(dayIdx);
    binIdx_day = trialBinIdx(dayIdx);
    embed_day = embedding(dayIdx,:);

    mu0_d = nan(nBins, size(embedding,2));
    mu1_d = nan(nBins, size(embedding,2));
    for b = 1:nBins
        mu0_d(b,:) = mean(embed_day(binIdx_day==b & Y_day==0,:), 1);
        mu1_d(b,:) = mean(embed_day(binIdx_day==b & Y_day==1,:), 1);
    end
    distsD = sqrt(sum((mu0_d - mu1_d).^2, 2));
    dayMeans(d) = mean(distsD, 'omitnan');
end

figure;
plot(1:nAligned, dayMeans, '-o', 'LineWidth', 2);
xticks(1:nAligned);
xticklabels(strrep(dateList,'_','\_'));
xtickangle(45);
ylabel('Mean separation (GLOBAL UMAP)');
xlabel('Day');
title('Mean separation across days (GLOBAL UMAP)');
grid on;


%% ======== Optional: Per-Day UMAP =========
if perDayUMAP
    fprintf('=== Running per-day UMAPs ===\n');
    allDayDists = nan(nAligned, nBins);
    allDayPvals = nan(nAligned, nBins);
    dayMeanSep = nan(nAligned,1);

    for d = 1:nAligned
        X_day = X_all(trialDayIdx==d,:);
        Y_day = labels(trialDayIdx==d);
        binIdx_day = trialBinIdx(trialDayIdx==d);

        if numel(unique(Y_day)) < 2
            fprintf('Skipping %s: only one class.\n', dateList{d});
            continue;
        end



        %[embedDay, ~] = run_umap(X_day, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.2, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
        %[embedDay, ~] = run_umap(X_day, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.1, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
        %[embedDay, ~] = run_umap(X_day, 'n_components', 2, 'n_neighbors', 40, 'min_dist', 0.6, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
        %[embedDay, ~] = run_umap(X_day, 'n_components', 5, 'n_neighbors', 40, 'min_dist', 0.1, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
        [embedDay, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.4, 'metric', 'cosine', 'randomize', true, 'verbose', false);
        %[embedDay, ~] = run_umap(X_day, 'n_components', 3, 'n_neighbors', 30, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
        %[embedDay, ~] = run_umap(X_day, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.6, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
        %[embedDay, ~] = run_umap(X_day, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
        %[embedDay, ~] = run_umap(X_day, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.4, 'metric', 'cosine', 'randomize', true, 'verbose', false);

        mu0 = nan(nBins, 3);
        mu1 = nan(nBins, 3);
        for b = 1:nBins
            mu0(b,:) = mean(embedDay(binIdx_day==b & Y_day==0,:), 1);
            mu1(b,:) = mean(embedDay(binIdx_day==b & Y_day==1,:), 1);
        end
        distsD = sqrt(sum((mu0 - mu1).^2, 2));
        dayMeanSep(d) = mean(distsD);

        % Perm test
        permD = nan(nPerm, nBins);
        for p = 1:nPerm
            permL = Y_day(randperm(length(Y_day)));
            for b = 1:nBins
                pmu0 = mean(embedDay(binIdx_day==b & permL==0,:), 1);
                pmu1 = mean(embedDay(binIdx_day==b & permL==1,:), 1);
                permD(p,b) = sqrt(sum((pmu0 - pmu1).^2));
            end
        end
        pV = mean(permD >= distsD', 1);

        % Save
        allDayDists(d,:) = distsD;
        allDayPvals(d,:) = pV;
    end

    % Plot separation per day
    figure;
    tiledlayout(ceil(nAligned/5), 5);
    for d = 1:nAligned
        nexttile;
        plot(1:nBins, allDayDists(d,:), '-or', 'LineWidth', 2);
        ylim([0 max(allDayDists(:))*1.1]);
        xlabel('Time Bin'); ylabel('Dist');
        title(sprintf('Day %s', dateList{d}));
    end
    sgtitle('Separation (per-day UMAP; axes NOT comparable)');

    % Plot p-values per day
    figure;
    tiledlayout(ceil(nAligned/5), 5);
    for d = 1:nAligned
        nexttile;
        plot(1:nBins, allDayPvals(d,:), '-ob', 'LineWidth', 2);
        hold on; yline(0.05, 'r--');
        ylim([0 1]);
        xlabel('Time Bin'); ylabel('p');
        title(sprintf('Day %s', dateList{d}));
    end
    sgtitle('p-values (per-day UMAP; axes NOT comparable)');

    % NEW: Plot how separation changes across days (per-day UMAP)
    figure;
    plot(1:nAligned, dayMeanSep, '-o', 'LineWidth', 2);
    xticks(1:nAligned);
    xticklabels(strrep(dateList,'_','\_'));
    xtickangle(45);
    ylabel('Mean separation (PER-DAY UMAP)');
    xlabel('Day');
    title('Mean separation across days (PER-DAY UMAP)');
    grid on;

end

%% ======== New: Per-Day Trajectories in Global UMAP ========
fprintf('=== Plotting per-day trajectories in global UMAP ===\n');

nCols = 5;
nRows = ceil(nAligned / nCols);

figure;
for d = 1:nAligned
    subplot(nRows, nCols, d);

    dayIdx = trialDayIdx == d;
    if sum(dayIdx) == 0
        title(sprintf('Day %s (no data)', dateList{d}));
        continue;
    end

    for b = 1:nBins
        mu0 = mean(embedding(dayIdx & trialBinIdx==b & labels==0,:), 1);
        mu1 = mean(embedding(dayIdx & trialBinIdx==b & labels==1,:), 1);

        % Plot as points and connecting lines
        if b == 1
            prev_mu0 = mu0;
            prev_mu1 = mu1;
        else
            % Plot trajectory line segment
            plot3([prev_mu0(1), mu0(1)], [prev_mu0(2), mu0(2)], [prev_mu0(3), mu0(3)], '-', 'Color', colors(1,:), 'LineWidth', 2); hold on;
            plot3([prev_mu1(1), mu1(1)], [prev_mu1(2), mu1(2)], [prev_mu1(3), mu1(3)], '-', 'Color', colors(2,:), 'LineWidth', 2);
            % Update
            prev_mu0 = mu0;
            prev_mu1 = mu1;
        end
    end

    xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3');
    title(sprintf('Day %s', dateList{d}));
    grid on;
end
sgtitle('Per-day mean trajectories (Global UMAP)');

%% ======== New: Per-Day Separation vs Time Bin (Global UMAP) ========
fprintf('=== Computing per-day separation vs time bin (Global UMAP) ===\n');

% Initialize storage for global axis standardization
allDayDists_store = nan(nAligned, nBins);
allDayPvals_store = nan(nAligned, nBins);

nCols = 5;
nRows = ceil(nAligned / nCols);


for d = 1:nAligned
    subplot(nRows, nCols, d);

    dayIdx = trialDayIdx == d;
    if sum(dayIdx) == 0
        title(sprintf('Day %s (no data)', dateList{d}));
        continue;
    end

    mu0_day = nan(nBins, size(embedding,2));
    mu1_day = nan(nBins, size(embedding,2));
    dist_day = nan(1, nBins);
    pvals_day = nan(1, nBins);

    % Compute per-bin means and distances
    for b = 1:nBins
        mu0_day(b,:) = mean(embedding(dayIdx & trialBinIdx==b & labels==0,:), 1);
        mu1_day(b,:) = mean(embedding(dayIdx & trialBinIdx==b & labels==1,:), 1);
        dist_day(b) = sqrt(sum((mu0_day(b,:) - mu1_day(b,:)).^2));
    end

    % Permutation test for day
    permDist_day = nan(nPerm, nBins);
    idx_in_day = find(dayIdx);
    permL_orig = labels(idx_in_day);

    for p = 1:nPerm
        permL = permL_orig(randperm(length(permL_orig)));

        for b = 1:nBins
            bin_mask_in_day = trialBinIdx(idx_in_day)==b;

            pmu0 = mean(embedding(idx_in_day(bin_mask_in_day & permL==0),:), 1);
            pmu1 = mean(embedding(idx_in_day(bin_mask_in_day & permL==1),:), 1);

            permDist_day(p,b) = sqrt(sum((pmu0 - pmu1).^2));
        end
    end

    pvals_day = mean(permDist_day >= dist_day, 1);
    allDayDists_store(d,:) = dist_day;
    allDayPvals_store(d,:) = pvals_day;


end


%% ======= Plot standardized per-day separation =======
fprintf('=== Plotting per-day separation with GLOBAL Y-Axis ===\n');

globalYmax = max(allDayDists_store(:)) * 1.1;

figure;
tiledlayout(ceil(nAligned/5), 5);
for d = 1:nAligned
    nexttile;
    plot(1:nBins, allDayDists_store(d,:), '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    ylim([0 globalYmax]);
    xlabel('Time Bin'); ylabel('Distance');
    title(sprintf('Day %s', strrep(dateList{d},'_','\_')));
end
sgtitle('Per-day Separation (Global UMAP, Standardized Y-Axis)');

%% ======= Plot standardized per-day p-values =======
fprintf('=== Plotting per-day p-values with GLOBAL Y-Axis ===\n');
globalPmax = min(1, max(allDayPvals_store(:)) * 1.1);  % cap at 1

figure;
tiledlayout(ceil(nAligned/5), 5);
for d = 1:nAligned
    nexttile;
    plot(1:nBins, allDayPvals_store(d,:), '-ob', 'LineWidth', 2);
    hold on; yline(0.05, 'r--');
    ylim([0 globalPmax]);
    xlabel('Time Bin'); ylabel('p-value');
    title(sprintf('Day %s', strrep(dateList{d},'_','\_')));
end
sgtitle('Per-day p-values (Global UMAP, Standardized Y-Axis)');


end
