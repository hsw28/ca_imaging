function [AUROC, Z_trials, labels] = runUMAP_fromStruct(animalName, refSpec, latentDim)
    % Performs UMAP-based manifold extraction, Procrustes alignment, and decoding
    % of trial correctness from calcium imaging data stored in a structured format.

    if nargin < 3
        latentDim = 15;
    end

    Fs  = 7.5;
    win = [0, 1.3];

    animal = evalin('base', animalName);
    G = animal.alignmentALL;

    dateList = autoDateList(animal);
    nDays = numel(dateList);

    nAligned = size(G, 2);
    if nAligned < nDays
        fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
        dateList = dateList(1:nAligned);
        nDays = nAligned;
    end

    nBins = round(diff(win) * Fs);
    Z_cells = cell(1, nDays);
    Z_trials = cell(1, nDays);
    labels = cell(1, nDays);

    fprintf('--- Extracting latent spaces (%d days) ---\n', nDays);
    for d = 1:nDays
        dateStr = dateList{d};
        if isBadDay(animal, dateStr)
            fprintf('Skipping %s: no trials or invalid CS times.\n', dateStr);
            continue;
        end

        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
        if isempty(X), continue; end

        X2d = squeeze(nanmean(X, 3));
        validCells = sum(isfinite(X2d), 2) >= round(0.8 * size(X2d,2));
        X2d = X2d(validCells, :);
        Xclean_cells = zscore(X2d, 0, 2);


        [embedding_cells, ~, ~] = run_umap(double(Xclean_cells), 'n_components', latentDim, ...
            'n_neighbors', 40, 'min_dist', 0.6, ...
            'verbose', false, 'cluster_output', 'none', 'metric', 'cosine');



            if istable(embedding_cells)
                embedding_trials = table2array(embedding_cells);
            elseif isstruct(embedding_cells)
                error('UMAP returned a struct instead of an array. Please update run_umap or check arguments.');
            end
            Z_cells{d} = embedding_cells(:, 1:min(latentDim, size(embedding_cells,2)))';

        X2d_trials = squeeze(nanmean(X, 2));
        validTrials = all(isfinite(X2d_trials), 1);
        Xclean_trials = zscore(X2d_trials(:, validTrials)', 0, 2);

        [embedding_trials, ~, ~] = run_umap(double(Xclean_trials), 'n_components', latentDim, ...
            'n_neighbors', 40, 'min_dist', 0.6, ...
            'verbose', false, 'cluster_output', 'none', 'metric', 'cosine');
            if istable(embedding_trials)
                embedding_trials = table2array(embedding_trials);
            elseif isstruct(embedding_trials)
                error('UMAP returned a struct instead of an array. Please update run_umap or check arguments.');
            end
            Z_trials{d} = embedding_trials(:, 1:min(latentDim, size(embedding_trials,2)))';


        labels{d} = y(validTrials);
    end

    validDays = find(~cellfun(@isempty, Z_cells));
    if nargin < 2 || strcmpi(refSpec, 'last') || strcmpi(refSpec, 'end')
        refDay = validDays(end);
    else
        refDay = validDays(refSpec);
    end
    if isempty(refDay), error('No valid reference day.'); end
    fprintf('Using day %d (%s) as alignment reference.\n', refDay, dateList{refDay});

    Z_aligned = Z_cells;
    for d = 1:nDays
        if d == refDay || isempty(Z_cells{d}), continue; end
        sharedGlobal = find(G(:,refDay) > 0 & G(:,d) > 0);
        if numel(sharedGlobal) < latentDim, continue; end

        idxRef = arrayfun(@(g) find(G(g,refDay)==G(:,refDay)), sharedGlobal);
        idxNew = arrayfun(@(g) find(G(g,d)==G(:,d)), sharedGlobal);
        if any(idxRef > size(Z_cells{refDay},2)) || any(idxNew > size(Z_cells{d},2)), continue; end

        A = Z_cells{d}(:, idxNew);
        B = Z_cells{refDay}(:, idxRef);
        [U,~,V] = svd(A*B','econ');
        R = U*V';
        Z_aligned{d} = R * Z_cells{d};
    end

    AUROC = nan(1,nDays);

    fprintf('\n--- Decoding correctness from reference day %d (%s) ---\n', refDay, dateList{refDay});

    refZ = Z_trials{refDay};

    disp(class(Z_trials{refDay}))
    whos Z_trials{refDay}

    if istable(refZ)
        refZ = table2array(refZ);
    elseif ~isnumeric(refZ)
        error('Z_trials{%d} is not numeric and cannot be used with SVM.', refDay);
    end
    mdl = fitcsvm(refZ', labels{refDay}, 'KernelFunction','linear');

    for d = 1:nDays
        if isempty(Z_trials{d}) || isempty(labels{d})
            fprintf('Day %d: skipped (empty data)\n', d);
            continue;
        end

        valid = ~isnan(labels{d});
        yt = labels{d}(valid);
        Ztest = Z_trials{d}(:, valid);

        if istable(Ztest)
            Ztest = table2array(Ztest);
        elseif ~isnumeric(Ztest)
            warning('Skipping day %d: Z_trials{%d} is not numeric.', d, d);
            continue
        end

        if numel(unique(yt)) < 2
            fprintf('Day %d: not enough label diversity — skipping AUROC.\n', d);
            continue;
        end

        [~, score] = predict(mdl, Ztest');
        [~,~,~,AUROC(d)] = perfcurve(yt, score(:,2), 1);
        fprintf('Day %d: AUROC = %.2f\n', d, AUROC(d));
    end


    figure;
    plot(1:nDays, AUROC, '-o','LineWidth',1.5);
    xlabel('Day index (in dateList)'); ylabel('AUROC');
    title('Cross-day correctness decoder (UMAP-based)');
    ylim([0, 1]); grid on;
    yline(0.5, '--k');
end
