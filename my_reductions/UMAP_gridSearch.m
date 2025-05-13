function results = UMAP_gridSearch(animalName, refSpec, latentDim, neighborVals, minDistVals, metricVals, kernelVals, kernelScaleVals)
  % Grid search over UMAP + SVM parameters to optimize decoding AUROC
  %
  % Inputs:
  %   animalName: string name of animal struct in base workspace (eg 'rat0314')
  %   refSpec: 'last', 'end', or index for reference day
  %   latentDim: UMAP latent dimensionality (e.g., 3)
  %   neighborVals: array of values for 'n_neighbors' eg [10, 20]
  %   minDistVals: array of values for 'min_dist' eg [.2, .5]
  %   metricVals: cell array of distance metrics {'cosine','euclidean'}
  %   kernelVals: cell array of SVM kernels {'linear','rbf'}
  %   kernelScaleVals: cell array of KernelScale values {'auto', 1, 0.5, etc.}
  %
  % ex: rat0314.UMAPgrid = UMAP_gridSearch('rat0314', 'last', 3, [10, 20], [0.2, 0.6], {'cosine'}, {'linear'}, kernelScales);
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
  if nargin < 6 || isempty(metricVals)
      metricVals = {'cosine'};
  end
  if nargin < 7 || isempty(kernelVals)
      kernelVals = {'linear'};
  end
  if nargin < 8 || isempty(kernelScaleVals)
      kernelScaleVals = {'auto'};
  end

  results = struct();
  results.AUROCs = nan(numel(neighborVals), numel(minDistVals), numel(metricVals), numel(kernelVals), numel(kernelScaleVals));
  results.params = struct('n_neighbors', neighborVals, 'min_dist', minDistVals, 'metric', {metricVals}, 'kernel', {kernelVals}, 'kernelScale', {kernelScaleVals});
  results.allAUROC = cell(size(results.AUROCs));

  labelsOut = strings(numel(neighborVals)*numel(minDistVals)*numel(metricVals)*numel(kernelVals)*numel(kernelScaleVals),1);
  aucsOut = nan(size(labelsOut));
  counter = 1;

  fprintf('Running grid search over neighbors x min_dist x metric x kernel x kernelScale (%d x %d x %d x %d x %d)...\n', numel(neighborVals), numel(minDistVals), numel(metricVals), numel(kernelVals), numel(kernelScaleVals));

  for i = 1:numel(neighborVals)
      for j = 1:numel(minDistVals)
          for k = 1:numel(metricVals)
              for l = 1:numel(kernelVals)
                  for m = 1:numel(kernelScaleVals)
                      n_neighbors = neighborVals(i);
                      min_dist = minDistVals(j);
                      metric = metricVals{k};
                      kernel = kernelVals{l};
                      kernelScale = kernelScaleVals{m};

                      fprintf('\n[Grid %d,%d,%d,%d,%d] n_neighbors=%d, min_dist=%.2f, metric=%s, kernel=%s, scale=%s\n', ...
                          i, j, k, l, m, n_neighbors, min_dist, metric, kernel, mat2str(kernelScale));

                      try
                          [AUROC, ~, ~] = runUMAP_single(animalName, refSpec, latentDim, n_neighbors, min_dist, metric, kernel, kernelScale);
                          last3 = AUROC(end-2:end);
                          meanAUROC = mean(last3,'omitnan');
                          results.AUROCs(i,j,k,l,m) = meanAUROC;
                          results.allAUROC{i,j,k,l,m} = AUROC;
                          AUROC
                          meanAUROC
                          perc = [NaN	0.410000000000000	0.620000000000000	0.510000000000000	0.500000000000000	0.570000000000000	0.890000000000000	0.760000000000000	0.440000000000000	0.690000000000000	0.950000000000000	0.630000000000000	0.860000000000000	0.860000000000000	0.880000000000000];
                          meanDif = nanmean(abs(AUROC-perc))
                          stdDif = nanstd(abs(AUROC-perc))

                          labelsOut(counter) = sprintf('k=%s | m=%s | n=%d | d=%.2f | s=%s', kernel, metric, n_neighbors, min_dist, mat2str(kernelScale));
                          aucsOut(counter) = meanAUROC;
                          counter = counter + 1;
                      catch ME
                          warning('Grid search failed at i=%d, j=%d, k=%d, l=%d, m=%d: %s', i, j, k, l, m, ME.message);
                      end
                  end
              end
          end
      end
  end

  fprintf('\nGrid search complete.\n');

  % Trim unused entries before plotting
  validIdx = ~isnan(aucsOut);
  labelsOut = labelsOut(validIdx);
  aucsOut = aucsOut(validIdx);

  % Bar plot of average AUROC over last 3 days
  figure;
  bar(aucsOut);
  ylabel('Mean AUROC (last 3 days)');
  title('UMAP Grid Search Performance');
  xticks(1:numel(aucsOut));
  xticklabels(labelsOut);
  xlabel('Param Combo (kernel | metric | neighbors | min\_dist | scale)');
  xtickangle(45);
  grid on;
end

function [AUROC, Z_trials, labels] = runUMAP_single(animalName, refSpec, latentDim, n_neighbors, min_dist, metric, kernel, kernelScale)
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
            'cluster_output', 'none', 'metric', metric);

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

    mdl = fitcsvm(refZ', labels{refDay}, 'KernelFunction', kernel, 'KernelScale', kernelScale);

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
