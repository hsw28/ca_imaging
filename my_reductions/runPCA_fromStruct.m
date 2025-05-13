
function AUROC = runPCA_fromStruct(animalName, refSpec, latentDim)
    % animalName  : string (e.g., 'rat0314') — must match workspace variable
    % refSpec     : integer index or 'last'/'end' (optional) — alignment reference day

    %ex   runPCA_fromStruct('rat0314', 1)
    %     runPCA_fromStruct('rat0314', 1)
    %     runPCA_fromStruct('rat0314', 'last')

    %% ----------------------------------------------------------------------
    %  runPCA_fromStruct.m
    %  Performs PCA-based manifold extraction, Procrustes alignment, and decoding
    %  of trial correctness from calcium imaging data stored in a structured format.
    %  ----------------------------------------------------------------------



    %% ---------- PARAMETERS -------------------------------------------------
    if nargin < 3;
      latentDim = 15;
    end

    Fs       = 7.5;          % frame rate (Hz)
    win      = [0, 1.3];    % time window around CS (in seconds)
    %latentDim = 25;          % number of principal components to keep

    % Access the structured data from the base workspace
    animal = evalin('base', animalName);
    G = animal.alignmentALL;    % cell-to-day global index matrix

    % Automatically get list of date strings
    dateList = autoDateList(animal);
    nDays = numel(dateList);

    % Only analyze aligned sessions
    nAligned = size(G, 2);
    if nAligned < nDays
        fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
        dateList = dateList(1:nAligned);
        nDays = nAligned;
    end

    nBins = round(diff(win) * Fs);  % number of time bins per trial
    Z_cells = cell(1, nDays);       % PC projections across cells (for alignment)
    Z_trials = cell(1, nDays);      % PC projections across trials (for decoding)
    labels = cell(1, nDays);        % trial correctness labels

    %% ---------- EXTRACT LATENT SPACES --------------------------------------
    fprintf('--- Extracting latent spaces (%d days) ---\n', nDays);
    figure
    for d = 1:nDays
        dateStr = dateList{d};

        if isBadDay(animal, dateStr)
            fprintf('Skipping %s: no trials or invalid CS times.\n', dateStr);
            continue
        end

        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
        if isempty(X), continue; end

        %across cells (average across trials/time bins)
        %What it does: Reduces the dimensionality of each cell’s temporal activity (e.g., average trace over time per condition).
        %Focus: How cells vary in temporal dynamics during a trial.
        %Use: Alignment across days using Procrustes or similar transforms to assess representational stability.
        %Interpretation: Reveals whether individual cells maintain similar functional roles (e.g., ramping, phasic responses) across days.
        %Helps answer:
        %Are the same cells encoding similar dynamics over time across different days or conditions?


        X2d = squeeze(nanmean(X, 3));
        validCells = sum(isfinite(X2d), 2) >= round(0.8 * size(X2d,2));
        X2d = X2d(validCells, :);
        Xclean_cells = zscore(X2d, 0, 2);
        [~, score_cells] = pca(Xclean_cells);
        Z_cells{d} = score_cells(:, 1:min(latentDim, end))';
        numPCs_cells = size(score_cells, 2);  % how many were returned
        fprintf('Day %d | PCA (cells): %d components returned\n', d, numPCs_cells);



        % PCA across cells (average across time for each trial)
        %What it does: Reduces dimensionality of the population vector per trial (spatial pattern across cells).
        %Focus: How trials differ in overall cell activation patterns.
        %Use: Trial-by-trial decoding (e.g., predict trial type or correctness).
        %Interpretation: Shows the structure of trial identity (e.g., correct vs incorrect) in population space.
        %Helps answer:
        %Can I predict trial outcome or type from the population activity pattern?

        % Collapse time to get a [nCells × nTrials] matrix
        X2d_trials = squeeze(nanmean(X, 2));

        % Select trials with full data
        validTrials = all(isfinite(X2d_trials), 1);

        % Normalize per trial across cells
        Xclean_trials = zscore(X2d_trials(:, validTrials)', 0, 2);  % [nTrials × nCells]

        % Run PCA and store top components
        [coeff, score_trials] = pca(Xclean_trials);
        Z_trials{d} = score_trials(:, 1:min(latentDim, size(score_trials,2)))';
        labels{d} = y(validTrials);

        fprintf('Day %d | PCA (trials): %d components returned\n', d, size(score_trials,2));

        [coeff, ~, ~] = pca(Xclean_cells);  % coeff is [nBins × nPCs]

          subplot(4, ceil((nDays-1)./4), d-1)
          length(coeff(:,1))
          plot(coeff(:,1:3), 'LineWidth', 2);  % Plot first 3 PCs
          xlabel('Time bin'); ylabel('PC amplitude');
          legend({'PC1','PC2','PC3'});
          if d-1>1
            legend('off')
          end
          title('Top Principal Components Across Cells (Temporal Modes)');


    end

    %% ---------- SELECT REFERENCE DAY --------------------------------------
    validDays = find(~cellfun(@isempty, Z_cells));
    if nargin < 2 || strcmpi(refSpec, 'last') || strcmpi(refSpec, 'end')
        refDay = validDays(end);
    else
        refDay = validDays(refSpec);
    end
    if isempty(refDay), error('No valid reference day.'); end
    fprintf('Using day %d (%s) as alignment reference.\n', refDay, dateList{refDay});

    %% ---------- ALIGN Z_cells VIA PROCRUSTES ------------------------------
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

    %% ---------- DECODE CORRECTNESS ACROSS DAYS ---------------------------
    mdl = fitcsvm(Z_trials{refDay}', labels{refDay}, 'KernelFunction','linear');
    AUROC = nan(1,nDays);

    fprintf('\n--- Decoding correctness from reference day %d (%s) ---\n', refDay, dateList{refDay});
    for d = 1:nDays
        if isempty(Z_trials{d}) || isempty(labels{d})
            fprintf('Day %d: skipped (empty data)\n', d);
            continue;
        end

        valid = ~isnan(labels{d});
        yt = labels{d}(valid);
        Ztest = Z_trials{d}(:, valid)';

        %fprintf('Day %d: %d valid trials, %d unique labels\n', d, sum(valid), numel(unique(yt)));
        %disp(tabulate(yt));

        if numel(unique(yt)) < 2
            fprintf('Day %d: not enough label diversity — skipping AUROC.\n', d);
            continue;
        end

        [~, score] = predict(mdl, Ztest);
        [~,~,~,AUROC(d)] = perfcurve(yt, score(:,2), 1);
    end

    % Plot AUROC across days


    figure;
    plot(1:nDays, AUROC, '-o','LineWidth',1.5);
    xlabel('Day index (in dateList)'); ylabel('AUROC');
    title('Cross-day correctness decoder (trial-level PCA)');
    ylim([0, 1]); grid on;
    hline(.5)
end
