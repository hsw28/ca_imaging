function driftStats = computeCellProjectionDrift(Z_cells_corr, alignment_data, validDays, AUROC)
% can get output directly from UMAP_cellsByOutcome
% Computes per-cell projection drift across days after UMAP alignment
% Inputs:
%   - Z_cells_corr: [1 x nDays] cell array of [latentDim x numCells]
%   - alignment_data: [nGlobalCells x nDays] matrix with index into Z_cells_corr{d}
%   - validDays: indices of days with usable data
% Output:
%   - driftStats: struct with fields:
%       - driftPerCell: [nTrackedCells x 1] average drift per global cell
%       - driftMatrix:  [nTrackedCells x nValidDays] Euclidean drift per day
%       - stableIdx:    indices of cells with low drift

G = alignment_data;

latentDim = size(Z_cells_corr{validDays(1)}, 1);
nDays = length(validDays);
nGlobal = size(G, 1);
driftMatrix = nan(nGlobal, nDays);

% Choose reference day as last
refDay = validDays(end);
refEmbeddings = Z_cells_corr{refDay};

for i = 1:nGlobal
    if G(i,refDay) == 0, continue; end
    refIdx = G(i,refDay);
    refVec = refEmbeddings(:, refIdx);

    for d = validDays
        if d == refDay || G(i,d) == 0, continue; end
        idx = G(i,d);
        if idx > size(Z_cells_corr{d}, 2), continue; end
        vec = Z_cells_corr{d}(:, idx);
        drift = norm(vec - refVec);
        driftMatrix(i, d) = drift;
    end
end

driftPerCell = mean(driftMatrix, 2, 'omitnan');
stableIdx = find(driftPerCell < prctile(driftPerCell, 25));

% Plot average drift per cell
figure;
subplot(3,2,2)
plot(driftPerCell, '.');
ylabel('Average Drift Across Days'); xlabel('Global Cell Index');
title('UMAP Projection Drift Across Days'); grid on;

% Histogram of drift
subplot(3,2,1)
histogram(driftPerCell, 36);
xlabel('Average Drift'); ylabel('Cell Count');
title('Distribution of Projection Drift');
grid on;

% Mean drift per day (aggregate across cells)
subplot(3,2,3)
meanDriftPerDay = nanmean(driftMatrix, 1);
plot(meanDriftPerDay, '-o');
xlabel('Day'); ylabel('Mean drift'); title('Drift over learning');

% For selected stable/unstable cells
subplot(3,2,4)
[~, sortIdx] = sort(nanmean(driftMatrix, 2)); % sort cells by stability
imagesc(driftMatrix(sortIdx,:));
xlabel('Day'); ylabel('Cell'); title('Drift per cell across learning');
colorbar;

subplot(3,2,5:6)
scatter(AUROC, meanDriftPerDay, 60, 'filled');
xlabel('Decoder AUROC');
ylabel('Mean cell drift');
title('Drift vs. decoding performance');

% Optionally fit a regression
lsline;
[r, p] = corr(meanDriftPerDay(:), AUROC(:), 'rows','complete');
text(0.1, 0.9, sprintf('r = %.2f, p = %.3f', r, p), 'Units','normalized');

driftStats.driftPerCell = driftPerCell;
driftStats.driftMatrix = driftMatrix;
driftStats.stableIdx = stableIdx;
