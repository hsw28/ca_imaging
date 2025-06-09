function [Z_cells_corr alignment_data validDays svm_aurocs] = UMAP_cellsByOutcome(animalName, latentDim, neighborVals, minDistVals, metricVals)
    % Performs UMAP embedding of cells based on their activity during
    % correct and incorrect trials separately, per day. Optionally aligns and overlays in common space.
    % ex    runUMAP_cellsByTrialOutcome('rat0314', 3, 190, .2, {'cosine'})


%It performs neuron-level embedding, not trial-level.
%For each day:
%It averages over all correct trials, and separately over all incorrect trials.
%So each condition (CR or no-CR) produces a single vector per neuron (its mean activity over time).
%These neuron vectors are used as input to UMAP.
%🧾 UMAP Input Matrix:
%Shape = [2 × neurons] × timepoints
%Then z-scored across time and passed to UMAP.
%It’s neuron-flattened.


    if nargin < 2, latentDim = 3; end
    if nargin < 3, neighborVals = 15; end
    if nargin < 4, minDistVals = .5; end
    if nargin < 5, metricVals = {'euclidean'}; end
      metric = metricVals{1};

    Fs  = 7.5;
    win = [0, 1.3];

    animal = evalin('base', animalName);
    dateList = autoDateList(animal);
    nDays = numel(dateList);

    G = animal.alignmentALL;
    alignment_data = G;
    nAligned = size(G, 2);
    if nAligned < nDays
        fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
        dateList = dateList(1:nAligned);
        nDays = nAligned;
    end

    Z_cells_corr = cell(1, nDays);
    Z_cells_incorr = cell(1, nDays);
    validDays = [];
    svm_aurocs = [];

    fprintf('--- Extracting cell embeddings per trial outcome (%d days) ---\n', nDays);

    for d = 1:nDays
        dateStr = dateList{d};
        if isBadDay(animal, dateStr)
            fprintf('Skipping %s: no trials or invalid CS times.\n', dateStr);
            continue;
        end

        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, round(diff(win)*Fs), Fs);
        if isempty(X), continue; end

        correct_trials = y == 1;
        incorrect_trials = y == 0;

        if sum(correct_trials) < 2 || sum(incorrect_trials) < 2
            fprintf('Skipping %s: not enough trials per condition.\n', dateStr);
            continue;
        end

        % average over trials on a given day
        X_corr = squeeze(nanmean(X(:,:,correct_trials), 3));
        X_incorr = squeeze(nanmean(X(:,:,incorrect_trials), 3));



        X_corr = zscore(X_corr, 0, 2);
        X_incorr = zscore(X_incorr, 0, 2);

        valid_corr = all(isfinite(X_corr), 2);
        valid_incorr = all(isfinite(X_incorr), 2);
        X_corr_valid = X_corr(valid_corr, :);
        X_incorr_valid = X_incorr(valid_incorr, :);


          Xcat = double([X_corr_valid; X_incorr_valid]);
          labels_cat = [ones(size(X_corr_valid,1),1); zeros(size(X_incorr_valid,1),1)];

          [embedding, ~, ~] = run_umap(Xcat, 'n_components', latentDim, ...
          'n_neighbors', neighborVals, 'min_dist', minDistVals, 'metric', metric, ...
          'verbose', false, 'cluster_output', 'none');

          embed_corr = embedding(1:size(X_corr_valid,1), :);
          embed_incorr = embedding(size(X_corr_valid,1)+1:end, :);

          mdl = fitcsvm(embedding, labels_cat, 'KernelFunction','linear');
          [~, score] = predict(mdl, embedding);
          [~,~,~,auroc] = perfcurve(labels_cat, score(:,2), 1);
          svm_aurocs(end+1) = auroc;

        Z_cells_corr{d} = embed_corr';
        Z_cells_incorr{d} = embed_incorr';
        validDays = [validDays, d];
    end



    %% Align correct trials across days
       Z_corr_aligned = Z_cells_corr;
       refDay = validDays(end);

       for d = validDays
           if d == refDay || isempty(Z_cells_corr{d}), continue; end
           sharedGlobal = find(G(:,refDay) > 0 & G(:,d) > 0);
           if numel(sharedGlobal) < latentDim, continue; end

           idxRef = arrayfun(@(g) find(G(g,refDay)==G(:,refDay)), sharedGlobal);
           idxNew = arrayfun(@(g) find(G(g,d)==G(:,d)), sharedGlobal);
           if any(idxRef > size(Z_cells_corr{refDay},2)) || any(idxNew > size(Z_cells_corr{d},2)), continue; end

           A = Z_cells_corr{d}(:, idxNew);
           B = Z_cells_corr{refDay}(:, idxRef);
           [U,~,V] = svd(A*B','econ');
           R = U*V';
           Z_corr_aligned{d} = R * Z_cells_corr{d};
       end

       %% Align incorrect trials across days
       Z_incorr_aligned = Z_cells_incorr;

       for d = validDays
           if d == refDay || isempty(Z_cells_incorr{d}), continue; end
           sharedGlobal = find(G(:,refDay) > 0 & G(:,d) > 0);
           if numel(sharedGlobal) < latentDim, continue; end

           idxRef = arrayfun(@(g) find(G(g,refDay)==G(:,refDay)), sharedGlobal);
           idxNew = arrayfun(@(g) find(G(g,d)==G(:,d)), sharedGlobal);
           if any(idxRef > size(Z_cells_incorr{refDay},2)) || any(idxNew > size(Z_cells_incorr{d},2)), continue; end

           A = Z_cells_incorr{d}(:, idxNew);
           B = Z_cells_incorr{refDay}(:, idxRef);
           [U,~,V] = svd(A*B','econ');
           R = U*V';
           Z_incorr_aligned{d} = R * Z_cells_incorr{d};
       end



       figure
    %% Plot aligned correct and incorrect embeddings together per day
    for d = validDays
        if isempty(Z_corr_aligned{d}) || isempty(Z_incorr_aligned{d}), continue; end
          subplot_tight(3,round((nAligned-1)/3),d-1)

        dataCat = [Z_corr_aligned{d}'; Z_incorr_aligned{d}'];
        labelsCat = [ones(size(Z_corr_aligned{d},2),1); zeros(size(Z_incorr_aligned{d},2),1)];

        colors = [1 0 0; 0 0 1];  % red for 0, blue for 1
        pointColors = colors(labelsCat + 1, :);  % convert labels (0/1) to row indices (1/2)

        if latentDim >= 3
            scatter3(dataCat(:,1), dataCat(:,2), dataCat(:,3), 25, pointColors, 'filled');
            xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3');
        else
            gscatter(dataCat(:,1), dataCat(:,2), labelsCat, 'br', 'ox');
            xlabel('UMAP 1'); ylabel('UMAP 2');
        end
        title(sprintf('Aligned cell embeddings for (day %s)', dateList{d}));
        legend('Correct', 'Incorrect');
        grid on;
    end

%{
        %% NEW: Plot combined trajectory across days with rainbow gradient
        figure; hold on;
        cmap = jet(numel(validDays));
        colors = jet(numel(validDays));  % Gradient by day
        for i = 1:numel(validDays)
            d = validDays(i);
            Zc = Z_corr_aligned{d}';
            Zi = Z_incorr_aligned{d}';

            scatter3(Zc(:,1), Zc(:,2), Zc(:,3), 10, colors(i,:), 'o', 'filled', 'MarkerFaceAlpha', 0.3);
            scatter3(Zi(:,1), Zi(:,2), Zi(:,3), 10, colors(i,:), '^', 'filled', 'MarkerFaceAlpha', 0.3);
        end

        xlabel('UMAP 1'); ylabel('UMAP 2');
        if latentDim == 3
            zlabel('UMAP 3');
        end
        title('Trajectory of aligned cell embeddings across days');
        grid on; view(3);
        scatter3(nan, nan, nan, 36, 'b', 'filled');
%}

quantifyUMAPManifolds(Z_corr_aligned, Z_incorr_aligned, validDays)
