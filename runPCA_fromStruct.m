
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

        % PCA across cells (average across trials) for decoding
        %How trials differ from each other based on the spatial pattern of cell activity.
        %The main dimensions along which trials vary — often related to trial type, performance (CR vs no CR), or learning stage.
        %Decoding trial identity — by training a classifier in trial space.
        %Helps answer: Can I tell what kind of trial this is, based on cell activation pattern?
        X2d = squeeze(nanmean(X, 3));
        validCells = sum(isfinite(X2d), 2) >= round(0.8 * size(X2d,2));
        X2d = X2d(validCells, :);
        Xclean_cells = zscore(X2d, 0, 2);
        [~, score_cells] = pca(Xclean_cells);
        Z_cells{d} = score_cells(:, 1:min(latentDim, end))';
        numPCs_cells = size(score_cells, 2);  % how many were returned
        fprintf('Day %d | PCA (cells): %d components returned\n', d, numPCs_cells);

        % PCA across trials (average across time bins) for alignment
        %How cells vary in their temporal activity patterns during a trial.
        %Finds components that capture shared dynamics across time (e.g., response to CS).
        %Alignment: Procrustes or other transformations align neural manifolds day to day.
        %Helps answer: Are the same cells encoding similar temporal dynamics across days?

        % This step extracts how whole-trial responses vary across cells.
        % Each trial is compressed into a vector (mean across time bins),
        % and we use PCA to find patterns of trial-to-trial variability
        % (e.g., related to behavioral outcome).


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

%% ==================== Helper functions ===============================
function [X, y, nC, nT, CSon, t] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs)
    csVar = ['CS_' dateStr];
    caVar = ['CA_traces_' dateStr];
    ftVar = ['CA_time_' dateStr];
    ttVar = ['EMGts_' dateStr];

    if ~isfield(animal.CS_times, csVar) || ~isfield(animal.Ca_traces, caVar) || ~isfield(animal.Ca_ts, ftVar) || ~isfield(animal.CRs, ttVar)
        X = []; y = [];
        return
    end

    CSon = animal.CS_times.(csVar);
    F = animal.Ca_traces.(caVar);
    t = animal.Ca_ts.(ftVar); t = t(:,2)/1000; t = t(1:2:end); t = t(1:size(F,2));
    y = animal.CRs.(ttVar);

    %fprintf('Raw y labels for %s: %d zeros, %d ones, %d NaNs\n', dateStr, sum(y==0), sum(y==1), sum(isnan(y)));

    if isempty(CSon) || numel(CSon) < 2 || all(isnan(CSon))
        X = []; y = [];
        return
    end

    nC = size(F,1);
    nT = numel(CSon);
    X = nan(nC, nBins, nT, 'single');
    usedTrialIdx = false(nT,1);

    for k = 1:nT
        t0 = CSon(k) + win(1);
        t1 = CSon(k) + win(2);
        if t0 < t(1) || t1 > t(end), continue; end

        idx = find(t >= t0 & t < t1);
        if numel(idx) < 2, continue; end

        trace = F(:, idx);
        if size(trace,2) < nBins
            trace(:, end+1:nBins) = NaN;
        elseif size(trace,2) > nBins
            trace = trace(:,1:nBins);
        end

        X(:,:,k) = trace;
        usedTrialIdx(k) = true;
    end

    X = X(:,:,usedTrialIdx);
    y = y(usedTrialIdx);

    if ~isempty(X)
        reshaped = reshape(X, nC, []);
        mu = mean(reshaped, 2, 'omitnan');
        X = bsxfun(@minus, X, mu);
    end
end

function dateList = autoDateList(animal)
    getDate = @(s) regexp(s,'\d{4}_\d{2}_\d{2}','match','once');
    fNames  = fieldnames(animal.Ca_traces);
    dateList = unique(cellfun(getDate, fNames, 'uni',0));
    dateList = sort(dateList);
end

function skip = isBadDay(animal, dateStr)
    csField = ['CS_' dateStr];
    skip = ~isfield(animal.CS_times, csField) || isempty(animal.CS_times.(csField)) || numel(animal.CS_times.(csField)) < 2 || all(isnan(animal.CS_times.(csField)));
end
