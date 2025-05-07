function results = UMAP_gridSearch(animalName, refSpec, latentDim, neighborVals, minDistVals)
% Grid search over UMAP parameters to optimize decoding AUROC
%
% Inputs:
%   animalName: string name of animal struct in base workspace (eg 'rat0314')
%   refSpec: 'last', 'end', or index for reference day
%   latentDim: UMAP latent dimensionality (e.g., 3)
%   neighborVals: array of values for 'n_neighbors' eg [10, 20]
%   minDistVals: array of values for 'min_dist' eg [.2, .5]
%
% Output:
%   results: struct with AUROC matrix, corresponding parameter values, etc.

if nargin < 3 || isempty(latentDim)
    latentDim = 3;
end
if nargin < 4 || isempty(neighborVals)
    neighborVals = 40;
end
if nargin < 5 || isempty(minDistVals)
    minDistVals = 0.6;
end

results = struct();
results.AUROCs = nan(numel(neighborVals), numel(minDistVals));
results.params = struct('n_neighbors', neighborVals, 'min_dist', minDistVals);
results.allAUROC = cell(numel(neighborVals), numel(minDistVals));

fprintf('Running grid search over %d neighbor values and %d min_dist values...\n', numel(neighborVals), numel(minDistVals));

for i = 1:numel(neighborVals)
    for j = 1:numel(minDistVals)
        n_neighbors = neighborVals(i);
        min_dist = minDistVals(j);
        fprintf('\n[Grid %d,%d] n_neighbors=%d, min_dist=%.2f\n', i, j, n_neighbors, min_dist);

        try
            [AUROC, ~, ~] = runUMAP_single(animalName, refSpec, latentDim, n_neighbors, min_dist);
            results.AUROCs(i,j) = mean(AUROC,'omitnan');
            results.allAUROC{i,j} = AUROC;
        catch ME
            warning('Grid search failed at i=%d, j=%d: %s', i, j, ME.message);
        end
    end
end


fprintf('\nGrid search complete.\n');

% Plot results
figure;
imagesc(results.AUROCs);
colorbar;
xlabel('min\_dist');
ylabel('n\_neighbors');
title('Mean AUROC (last 3 days) across UMAP grid search');
xticks(1:numel(minDistVals));
xticklabels(string(minDistVals));
yticks(1:numel(neighborVals));
yticklabels(string(neighborVals));
end



function [AUROC, Z_trials, labels] = runUMAP_single(animalName, refSpec, latentDim, n_neighbors, min_dist)
    % Core logic extracted from runUMAP_fromStruct, parametrized by n_neighbors and min_dist

    Fs  = 7.5;
    win = [0, 1.3];

    animal = evalin('base', animalName);
    G = animal.alignmentALL;
    dateList = autoDateList(animal);
    nDays = numel(dateList);

    nAligned = size(G, 2);
    if nAligned < nDays
        dateList = dateList(1:nAligned);
        nDays = nAligned;
    end

    nBins = round(diff(win) * Fs);
    Z_trials = cell(1, nDays);
    labels = cell(1, nDays);

    for d = 1:nDays
        dateStr = dateList{d};
        if isBadDay(animal, dateStr), continue; end

        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
        if isempty(X), continue; end

        X2d_trials = squeeze(nanmean(X, 2));
        validTrials = all(isfinite(X2d_trials), 1);
        Xclean_trials = zscore(X2d_trials(:, validTrials)', 0, 2);

        embedding = run_umap(double(Xclean_trials), 'n_components', latentDim, ...
            'n_neighbors', n_neighbors, 'min_dist', min_dist, 'verbose', false, ...
            'cluster_output', 'none', 'metric', 'cosine');

        if istable(embedding)
            embedding = table2array(embedding);
        end

        Z_trials{d} = embedding(:, 1:min(latentDim, size(embedding,2)))';
        labels{d} = y(validTrials);
    end

    validDays = find(~cellfun(@isempty, Z_trials));
    if strcmpi(refSpec, 'last') || strcmpi(refSpec, 'end')
        refDay = validDays(end);
    else
        refDay = validDays(refSpec);
    end

    AUROC = nan(1, nDays);
    refZ = Z_trials{refDay};
    if istable(refZ), refZ = table2array(refZ); end

    mdl = fitcsvm(refZ', labels{refDay}, 'KernelFunction','linear');

    for d = 1:nDays
        if isempty(Z_trials{d}) || isempty(labels{d}), continue; end

        valid = ~isnan(labels{d});
        yt = labels{d}(valid);
        Ztest = Z_trials{d}(:, valid);

        if istable(Ztest), Ztest = table2array(Ztest); end
        if numel(unique(yt)) < 2, continue; end

        [~, score] = predict(mdl, Ztest');
        [~,~,~,AUROC(d)] = perfcurve(yt, score(:,2), 1);
    end
end
