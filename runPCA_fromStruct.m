%% ----------------------------------------------------------------------
%  runPCA_fromStruct.m
%  Performs PCA-based manifold extraction + alignment AND trial decoding
%  from structured animal data with G index matrix.
%  ----------------------------------------------------------------------

%% ---------- USER INPUT -------------------------------------------------
animal   = rat0313;                         % master struct
Fs       = 7.5;                             % frame rate (Hz)
G        = rat0313.alignmentALL;           % [nGlobalCells × nDays] index matrix

% Automatically get list of date strings
dateList = autoDateList(animal);
nDays    = numel(dateList);

% Only analyze aligned sessions
nAligned = size(G,2);
if nAligned < nDays
    fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
    dateList = dateList(1:nAligned);
    nDays = nAligned;
end

%% ---------- PARAMETERS -------------------------------------------------
win        = [-1 , 1.5];              % sec around CS onset
nBins      = round(diff(win)*Fs);     % time bins
latentDim  = 10;                      % # of PCs

Z_cells   = cell(1,nDays);            % latentDim × nCells per day (for alignment)
Z_trials  = cell(1,nDays);            % latentDim × nTrials per day (for decoding)
labels    = cell(1,nDays);            % trial correctness vectors

fprintf('--- Extracting latent spaces (%d days) ---\n', nDays);
for d = 1:nDays
    dateStr = dateList{d};

    if isBadDay(animal, dateStr)
        fprintf('Skipping %s: no trials or invalid CS times.\n', dateStr);
        continue
    end

    [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
    if isempty(X)
        continue
    end

    % --- PCA across cells (for alignment)
    X2d = squeeze(nanmean(X, 3));              % [nCells × nBins]
    validCells = sum(isfinite(X2d), 2) >= round(0.8 * size(X2d,2));
    X2d = X2d(validCells, :);                  % filter cells
    Xclean_cells = zscore(X2d, 0, 2);          % normalize per cell across time
    [~, score_cells] = pca(Xclean_cells);
    score_cells = score_cells(:, 1:min(latentDim, size(score_cells,2)));
    Z_cells{d} = score_cells';                 % [latentDim × nCells]

    % --- PCA across trials (for decoding)
    X2d_trials = squeeze(nanmean(X, 2));       % [nCells × nTrials]
    validTrials = all(isfinite(X2d_trials), 1);

    fprintf('Day %d raw: %d trials, valid: %d trials after NaN filter\n', ...
        d, size(X,3), sum(validTrials));

    Xclean_trials = zscore(X2d_trials(:, validTrials)', 0, 2);  % [nTrials × nCells]
    [~, score_trials] = pca(Xclean_trials);
    score_trials = score_trials(:, 1:min(latentDim, size(score_trials, 2)));
    Z_trials{d} = score_trials';              % [latentDim × nTrials]
    labels{d} = y(validTrials);

    fprintf('Day %s  |  %d cells  %d trials\n', dateStr, size(X,1), numel(y));
end

% Pick first valid ref day for alignment
refDay = find(~cellfun(@isempty, Z_cells), 1, 'last');
if isempty(refDay), error('No valid reference day.'); end
fprintf('Using day %d (%s) as alignment reference.\n', refDay, dateList{refDay});

%% ---------- Align Z_cells via Procrustes ---------------------
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

%% ---------- Decode correctness across days -------------------
mdl = fitcsvm(Z_trials{refDay}', labels{refDay}, 'KernelFunction','linear');
AUROC = nan(1,nDays);

fprintf('\n--- Decoding correctness from reference day %d (%s) ---\n', ...
        refDay, dateList{refDay});

for d = 1:nDays
    if isempty(Z_trials{d}) || isempty(labels{d})
        fprintf('Day %d: skipped (empty data)\n', d);
        continue;
    end

    valid = ~isnan(labels{d});
    yt = labels{d}(valid);
    Ztest = Z_trials{d}(:, valid)';

    fprintf('Day %d: %d valid trials, %d unique labels\n', ...
            d, sum(valid), numel(unique(yt)));
    disp(tabulate(yt));

    if numel(unique(yt)) < 2
        fprintf('Day %d: not enough label diversity — skipping AUROC.\n', d);
        continue;
    end

    [~, score] = predict(mdl, Ztest);
    [~,~,~,AUROC(d)] = perfcurve(yt, score(:,2), 1);
end

figure;
plot(1:nDays, AUROC, '-o','LineWidth',1.5);
xlabel('Day index (in dateList)'); ylabel('AUROC');
title('Cross-day correctness decoder (trial-level PCA)');
ylim([0.5, 1]); grid on;

%% ==================== Helper functions ===============================
function [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs)
    csVar = ['CS_' dateStr];
    caVar = ['CA_traces_' dateStr];
    ftVar = ['CA_time_' dateStr];
    ttVar = ['EMGts_' dateStr];

    if ~isfield(animal.CS_times, csVar) || ...
       ~isfield(animal.Ca_traces, caVar) || ...
       ~isfield(animal.Ca_ts, ftVar) || ...
       ~isfield(animal.CRs, ttVar)
        X = []; y = [];
        return
    end

    CSon = animal.CS_times.(csVar);
    F    = animal.Ca_traces.(caVar);
    t    = animal.Ca_ts.(ftVar); t = t(:,2)/1000; t = t(1:2:end); t = t(1:size(F,2));
    y    = animal.CRs.(ttVar);

    fprintf('Raw y labels for %s: %d zeros, %d ones, %d NaNs\n', ...
        dateStr, sum(y==0), sum(y==1), sum(isnan(y)));

    if isempty(CSon) || numel(CSon) < 2 || all(isnan(CSon))
        X = []; y = [];
        return
    end

    nC = size(F,1); nT = numel(CSon);
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
    skip = ~isfield(animal.CS_times, csField) || ...
           isempty(animal.CS_times.(csField)) || ...
           numel(animal.CS_times.(csField)) < 2 || ...
           all(isnan(animal.CS_times.(csField)));
end
