function results = UMAPcells_gridSearch(animalName, latentDim, neighborVals, minDistVals, metricVals)

  % Grid search over UMAP parameters for decoding trial outcome
  % Returns centroid distance, SVM AUROC, and within-class variance
      %Centroid distance between correct and incorrect trial embeddings (want high)
      %SVM AUROC to quantify separability (want high)
      %Within-class variance (want low)



  if nargin < 2, latentDim = 3; end
  if nargin < 3, neighborVals = [10, 15, 30]; end
  if nargin < 4, minDistVals = [0.1, 0.5, 0.8]; end
  if nargin < 5, metricVals = {'euclidean'}; end

  Fs  = 7.5;
  win = [0, 1.3];

  animal = evalin('base', animalName);
  dateList = autoDateList(animal);
  nDays = numel(dateList);

  G = animal.alignmentALL;
  nAligned = size(G, 2);
  if nAligned < nDays
      dateList = dateList(1:nAligned);
      nDays = nAligned;
  end

  nN = length(neighborVals);
  nD = length(minDistVals);
  nM = length(metricVals);

  centroidDist = nan(nN, nD, nM);
  meanAUROC    = nan(nN, nD, nM);
  withinVar    = nan(nN, nD, nM);

  results = struct();
  results.metrics = struct();

  for i = 1:nN
      for j = 1:nD
          for k = 1:nM
              n_neighbors = neighborVals(i);
              min_dist = minDistVals(j);
              metric = metricVals{k};
              fprintf("Running UMAP with n_neighbors=%d, min_dist=%.2f, metric=%s\n", ...
                  n_neighbors, min_dist, metric);

              centroid_dists = [];
              svm_aurocs = [];
              var_corr = [];
              var_incorr = [];

              %for d = 1:nDays
              for d = nDays-2:nDays

                  dateStr = dateList{d};
                  if isBadDay(animal, dateStr), continue; end
                  [X, y] = getDayMatrixFromStruct(animal, dateStr, win, round(diff(win)*Fs), Fs);
                  if isempty(X), continue; end

                  correct_trials = y == 1;
                  incorrect_trials = y == 0;
                  if sum(correct_trials) < 2 || sum(incorrect_trials) < 2, continue; end

                  X_corr = squeeze(nanmean(X(:,:,correct_trials), 3));
                  X_incorr = squeeze(nanmean(X(:,:,incorrect_trials), 3));
                  X_corr = zscore(X_corr, 0, 2);
                  X_incorr = zscore(X_incorr, 0, 2);

                  valid_corr = all(isfinite(X_corr), 2);
                  valid_incorr = all(isfinite(X_incorr), 2);
                  X_corr_valid = X_corr(valid_corr, :);
                  X_incorr_valid = X_incorr(valid_incorr, :);

                  Xcat = double([X_corr_valid; X_incorr_valid]);
                  labels = [ones(size(X_corr_valid,1),1); zeros(size(X_incorr_valid,1),1)];

                  [embedding, ~, ~] = run_umap(Xcat, 'n_components', latentDim, ...
                      'n_neighbors', n_neighbors, 'min_dist', min_dist, 'metric', metric, ...
                      'verbose', false, 'cluster_output', 'none');

                  corr_pts = embedding(1:size(X_corr_valid,1), :);
                  incorr_pts = embedding(size(X_corr_valid,1)+1:end, :);

                  dist = norm(mean(corr_pts) - mean(incorr_pts));
                  centroid_dists(end+1) = dist;

                  mdl = fitcsvm(embedding, labels, 'KernelFunction','linear');
                  [~, score] = predict(mdl, embedding);
                  [~,~,~,auroc] = perfcurve(labels, score(:,2), 1);
                  svm_aurocs(end+1) = auroc;

                  var_corr(end+1) = mean(var(corr_pts));
                  var_incorr(end+1) = mean(var(incorr_pts));
              end

              centroidDist(i,j,k) = mean(centroid_dists);
              meanAUROC(i,j,k)    = mean(svm_aurocs);
              withinVar(i,j,k)    = mean([var_corr var_incorr]);

              results.metrics(i,j,k).n_neighbors = n_neighbors;
              results.metrics(i,j,k).min_dist = min_dist;
              results.metrics(i,j,k).metric = metric;
              results.metrics(i,j,k).centroid_dist = centroidDist(i,j,k);
              results.metrics(i,j,k).svm_auroc = meanAUROC(i,j,k);
              results.metrics(i,j,k).within_var = withinVar(i,j,k);
          end
      end
  end

  %% Plotting per metric
  for k = 1:nM
      figure;
      imagesc(squeeze(centroidDist(:,:,k))); colorbar;
      title(sprintf('Centroid Distance (%s)', metricVals{k}));
      xlabel('min\_dist'); ylabel('n\_neighbors');
      xticks(1:nD); xticklabels(string(minDistVals));
      yticks(1:nN); yticklabels(string(neighborVals));

      figure;
      imagesc(squeeze(meanAUROC(:,:,k))); colorbar;
      title(sprintf('SVM AUROC (%s)', metricVals{k}));
      xlabel('min\_dist'); ylabel('n\_neighbors');
      xticks(1:nD); xticklabels(string(minDistVals));
      yticks(1:nN); yticklabels(string(neighborVals));

      figure;
      imagesc(squeeze(withinVar(:,:,k))); colorbar;
      title(sprintf('Within-Class Variance (%s)', metricVals{k}));
      xlabel('min\_dist'); ylabel('n\_neighbors');
      xticks(1:nD); xticklabels(string(minDistVals));
      yticks(1:nN); yticklabels(string(neighborVals));
  end
  end
