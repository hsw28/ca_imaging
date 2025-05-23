function UMAP_trajectory_byOutcome(animalName, win, binSize)
  %time flattened, looks at trajectories of individual cells throughout days
  %compresses all time-binned population vectors across all days and trials into a single manifold,
  %where each line = one trial and each point = one time bin within that trial.
  %These are trial-level trajectories, not cell-level.
  % ========================
  % UMAP Trajectories by Outcome and by Day
  % ========================

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

  % First pass: aggregate for all-day embedding
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

  X_all = cell2mat(cellfun(@(x) x(:)', trialVecs, 'UniformOutput', false));
  X_all = zscore(X_all, 0, 2);

  % Global UMAP embedding
  %[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.2, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.1, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 2, 'n_neighbors', 40, 'min_dist', 0.6, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 5, 'n_neighbors', 40, 'min_dist', 0.1, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 40, 'min_dist', 0.4, 'metric', 'cosine', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 30, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.6, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
%[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.4, 'metric', 'euclidean', 'randomize', true, 'verbose', false);
[embedding, ~] = run_umap(X_all, 'n_components', 3, 'n_neighbors', 10, 'min_dist', 0.4, 'metric', 'cosine', 'randomize', true, 'verbose', false);

  nTrials = length(labels) / nBins;
  cmap = parula(nBins);

  % ==============================
  % (1) Trajectories by Time, Per Day
  % ==============================
  figure;
  nCols = 5;
  nRows = ceil(nAligned / nCols);
  for d = 1:nAligned
      subplot(nRows, nCols, d); hold on;
      idx = find(trialDayIdx == d);
      if isempty(idx), continue; end
      for i = 1:sum(trialDayIdx == d)/nBins
          ix = idx((i-1)*nBins + (1:nBins));
          traj = embedding(ix,:);
          for j = 1:nBins-1
              plot(traj(j:j+1,1), traj(j:j+1,2), '-', 'Color', cmap(j,:), 'LineWidth', 1);
          end
      end
      title(dateList{d}, 'Interpreter','none');
      axis tight; grid on;
  end

  % ==============================
  % (2) Trajectories by Time, All Days
  % ==============================
  figure; hold on;
  for i = 1:nTrials
      idxs = (i-1)*nBins + (1:nBins);
      traj = embedding(idxs, :);
      for j = 1:nBins-1
          plot(traj(j:j+1,1), traj(j:j+1,2), '-', 'Color', cmap(j,:), 'LineWidth', 1.2);
      end
  end
  title('UMAP Trajectories: All Days (Colored by Time)');
  xlabel('UMAP 1'); ylabel('UMAP 2'); grid on;

  % ==============================
  % (3) Separation Between CR and no-CR (All Days)
  % ==============================
  embed0 = reshape(embedding(labels==0,:), nBins, []);
  embed1 = reshape(embedding(labels==1,:), nBins, []);
  mu0 = mean(embed0, 2);
  mu1 = mean(embed1, 2);
  separation = vecnorm(mu1 - mu0, 2, 2);

  figure;
  plot(1:nBins, separation, 'k-', 'LineWidth', 2);
  title('Separation Between Correct and Incorrect Trajectories (All Days)');
  xlabel('Time Bin'); ylabel('Euclidean Distance'); grid on;

  % ==============================
  % (4) Separation by Day (Subplots)
  % ==============================
  figure;
  for d = 1:nAligned
      idx = find(trialDayIdx == d);
      if isempty(idx), continue; end
      embed_d = embedding(idx,:);
      labels_d = labels(idx);

      embed0 = embed_d(labels_d==0,:);
      embed1 = embed_d(labels_d==1,:);
      if isempty(embed0) || isempty(embed1), continue; end

      embed0 = reshape(embed0, nBins, []);
      embed1 = reshape(embed1, nBins, []);
      mu0 = mean(embed0, 2);
      mu1 = mean(embed1, 2);
      sep = vecnorm(mu1 - mu0, 2, 2);

      subplot(nRows, nCols, d);
      plot(1:nBins, sep, 'k-', 'LineWidth', 1.5);
      title(dateList{d}, 'Interpreter','none');
      xlabel('Time'); ylabel('Dist');
      axis([1 nBins 0 max(sep)*1.1]);
  end

  % ==============================
  % (5) Outcome Trajectories by Day (Subplots)
  % ==============================
  figure;
  for d = 1:nAligned
      subplot(nRows, nCols, d);
      hold on;
      idx = find(trialDayIdx == d);
      if isempty(idx), continue; end
      uT = sum(trialDayIdx == d) / nBins;
      for i = 1:uT
          ix = idx((i-1)*nBins + (1:nBins));
          traj = embedding(ix,:);
          lbl = labels(ix(1));
          if lbl == 1
              col = [0 0 1];  % blue for correct
          else
              col = [1 0 0];  % red for incorrect
          end
          for j = 1:nBins-1
              plot(traj(j:j+1,1), traj(j:j+1,2), '-', 'Color', col, 'LineWidth', 1);
          end
      end
      title(dateList{d}, 'Interpreter','none');
      axis tight; grid on;
  end

  % ==============================
  % (6) Outcome Trajectories All Days
  % ==============================
  figure; hold on;
  for i = 1:nTrials
      idxs = (i-1)*nBins + (1:nBins);
      traj = embedding(idxs, :);
      label = labels(idxs(1));
      if label == 1
          col = [0 0 1];  % blue
      else
          col = [1 0 0];  % red
      end
      for j = 1:nBins-1
          plot(traj(j:j+1,1), traj(j:j+1,2), '-', 'Color', col, 'LineWidth', 1.2);
      end
  end
  title('UMAP Trajectories by Trial Outcome (All Days)');
  xlabel('UMAP 1'); ylabel('UMAP 2'); grid on;

  end
