function neuralWithinSession(csusNum, varargin)
% neuralWithinSession(csusNum, 'K', 5, 'MakePlots', true, 'RatNames', {'rat0222'})
% - csusNum: 15/30/45/60/90 (selects rat.csusXX)
% - K: number of early and late trials to average (default 5)
% - MakePlots: show per-day curves in addition to the per-rat summary (default true)
% - RatNames: which rats to run (default {'rat0222'})

% ---------------- args ----------------
p = inputParser;
p.addParameter('K', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('MakePlots', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('RatNames', {'rat0222','rat0307','rat0313','rat0314','rat0816'}, @(c)iscell(c)||isstring(c));
p.parse(varargin{:});
K          = p.Results.K;
MakePlots  = p.Results.MakePlots;
ratNames   = cellstr(p.Results.RatNames);

SR = 7.5; %#ok<NASGU>  % available if needed downstream

for ii = 1:numel(ratNames)
    ratVar = ratNames{ii};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratVar))
        warning('Variable %s not found in base workspace. Skipping.', ratVar);
        continue;
    end
    rat = evalin('base', ratVar);

    fprintf('\n=== %s : within-session analysis (neural) ===\n', ratVar);

    % ---------------- 1) day selection (CS days up to An; last 3) ----------------
    dateListRaw = getAllDayKeys(rat);
    [csKeys, csDates] = keepCSdays(dateListRaw);
    anDate = datenum(strrep(rat.An,'_','-'));
    useMask = csDates <= anDate;
    csKeysUpToAn = csKeys(useMask);
    csDatesUpToAn = csDates(useMask);
    if numel(csKeysUpToAn) < 3
        warning('%s: Not enough CS days up to An (%s). Skipping.', ratVar, rat.An);
        continue;
    end
    [~, ord] = sort(csDatesUpToAn, 'ascend');
    daysToUse = csKeysUpToAn(ord(end-2:end));  % last 3 CS days up to An

    % storage for per-day early/late summaries and TT trend
    earlyLate.TT = nan(numel(daysToUse),2);
    earlyLate.TN = nan(numel(daysToUse),2);
    earlyLate.NN = nan(numel(daysToUse),2);
    rhoTT        = nan(numel(daysToUse),1);

    % ---------------- 2) iterate days ----------------
    for d = 1:numel(daysToUse)
        dayKey = daysToUse{d};              % e.g., 'CS_2023_05_08'
        usKey  = regexprep(dayKey,'^CS_','US_');
        dateTok = regexp(dayKey,'(\d{4}_\d{2}_\d{2})','tokens','once');  % 'YYYY_MM_DD'
        if isempty(dateTok)
            fprintf('  [%s] bad key — skipping.\n', dayKey);
            continue;
        end
        dateTok = dateTok{1};

        % --- CS/US onsets ---
        if ~isfield(rat.CS_times, dayKey) || ~isfield(rat.US_times, usKey)
            fprintf('  [%s] missing CS/US times — skipping.\n', dayKey);
            continue;
        end
        cs_on = valueVec(rat.CS_times.(dayKey));
        us_on = valueVec(rat.US_times.(usKey));
        nT    = min(numel(cs_on), numel(us_on));
        if nT < 1
            fprintf('  [%s] no trials — skipping.\n', dayKey);
            continue;
        end
        cs_on = cs_on(1:nT);  us_on = us_on(1:nT);
        good  = (us_on - cs_on) > 0.1;      % positive TEBC length
        cs_on = cs_on(good);  us_on = us_on(good);
        nT    = numel(cs_on);
        if nT < 6
            fprintf('  [%s] only %d valid trials — skipping.\n', dayKey, nT);
            continue;
        end

        % --- Ca_peaks for this date (field could be CA_peaks_YYYY_MM_DD etc.) ---
        peaksFN = findFieldWithDate(rat.Ca_peaks, dateTok);
        if isempty(peaksFN)
            fprintf('  [%s] no Ca_peaks for %s — skipping.\n', dayKey, dateTok);
            continue;
        end
        S_pk = rat.Ca_peaks.(peaksFN);      % [cells x events], NaN/0 padded times (sec)

        % ---------------- windows (relative to CS/US) ----------------
        csWin   = [0,    0.50];   % CS: 0–0.50 s after CS
        trWin   = [0.50, NaN];    % Trace: 0.50 s after CS up to US
        usWin   = [0.00, 0.15];   % US: 0–0.15 s after US
        baseWin = [-1.0, 0.0];    % baseline: 1 s pre-CS

        % use first up-to-10 trials for the “learning-ish” look
        maxK = min(10, nT);
        trIdx   = 1:maxK;
        cs_on_K = cs_on(trIdx);
        us_on_K = us_on(trIdx);

        % Per-neuron, per-trial rates
        rate = perNeuronRatesInWindows(S_pk, cs_on_K, us_on_K, csWin, trWin, usWin, baseWin);

        % ---------------- population-vector similarity, binned across trials ----------------
        nBins = 5;  % show ~early→late evolution within first ~10 trials
        [TT, TN, NN, trialBins] = popVecSim_acrossTrials(rate, cs_on_K, us_on_K, baseWin, nBins);

        % quick TT monotonic trend over the first ~10 trials
        tPairs = (1:numel(TT))';
        if numel(TT) >= 3
            rhoTT(d,1) = corr(tPairs, TT, 'type','Spearman','rows','pairwise');
        end

        % ---------------- Early vs Late (uses user K; auto-clipped) ----------------
        Keff_TN = min(K, numel(TN));
        Keff_TT = min(K, numel(TT));
        Keff_NN = min(K, numel(NN));

        % TT/NN are defined on pairs (length nT-1), TN on trials (length nT)
        idxTT_E = 1:max(1, Keff_TT-1);
        idxTT_L = max(1, numel(TT)-Keff_TT+1) : numel(TT);

        idxTN_E = 1:Keff_TN;
        idxTN_L = max(1, numel(TN)-Keff_TN+1) : numel(TN);

        idxNN_E = 1:max(1, Keff_NN-1);
        idxNN_L = max(1, numel(NN)-Keff_NN+1) : numel(NN);

        earlyLate.TT(d,:) = [mean(TT(idxTT_E),'omitnan'), mean(TT(idxTT_L),'omitnan')];
        earlyLate.TN(d,:) = [mean(TN(idxTN_E),'omitnan'), mean(TN(idxTN_L),'omitnan')];
        earlyLate.NN(d,:) = [mean(NN(idxNN_E),'omitnan'), mean(NN(idxNN_L),'omitnan')];

        % ---------------- optional per-day figure ----------------
        if MakePlots
            % For plotting curves, ensure x and y lengths match per series
            figure('Name',sprintf('%s | %s | pop vec sim',ratVar,dayKey), ...
                   'Color','w','Position',[100 100 640 420]); hold on
            plot(trialBins(1:numel(TT)), TT, '-o','LineWidth',2,'DisplayName','Trial–Trial');
            plot(trialBins(1:numel(TN)), TN, '-s','LineWidth',2,'DisplayName','Trial–Nontrial');
            plot(trialBins(1:numel(NN)), NN, '-^','LineWidth',2,'DisplayName','Nontrial–Nontrial');
            xlabel('Trial bin'); ylabel('Cosine similarity');
            title(sprintf('%s — %s', ratVar, dayKey));
            legend('Location','best'); grid on
        end
    end

    % ---------------- ONE FIGURE PER RAT: Early vs Late summary ----------------
    if any(isfinite(earlyLate.TT(:)) | isfinite(earlyLate.TN(:)) | isfinite(earlyLate.NN(:)))
        Kshow = K;  % just for title
        figure('Name',sprintf('%s | pop vec sim: early vs late',ratVar), ...
               'Color','w','Position',[80 80 980 380]);
        tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        labels = {'Early','Late'};

        % Helper for bars+dots
        panelPlot = @(sub, dat, ttl) ...
            plotBarWithDots(sub, dat, labels, ttl);

        nexttile; panelPlot(1, earlyLate.TT, 'TT (trial–trial)');
        nexttile; panelPlot(2, earlyLate.TN, 'TN (trial–nontrial)');
        nexttile; panelPlot(3, earlyLate.NN, 'NN (nontrial–nontrial)');

        title(tl, sprintf('%s — early (%d) vs late (%d) trials | TT trend \\rho (per day): %s', ...
            ratVar, Kshow, Kshow, strtrim(sprintf('%.2f ', rhoTT))))
    end
end
end

% ================= helpers =================

function out = valueVec(x)
% Turn whatever is stored in CS_times/US_times (vector or struct) into a numeric row vector
    if isnumeric(x), out = x(:)'; return; end
    if isstruct(x)
        f = fieldnames(x);
        for k = 1:numel(f)
            if isnumeric(x.(f{k})), out = x.(f{k})(: )'; return; end
        end
    end
    error('CS/US times are not numeric.');
end

function fn = findFieldWithDate(S, dateTok)
% Return the FIRST field in struct S whose name contains the date token 'YYYY_MM_DD'
    fn = '';
    if ~isstruct(S), return; end
    f = fieldnames(S);
    pat = strrep(dateTok,'_','[_-]?'); % tolerate '_' or '-'
    for i = 1:numel(f)
        if ~isempty(regexp(f{i}, pat, 'once'))
            fn = f{i}; return;
        end
    end
end

function out = getAllDayKeys(rat)
% Union of keys across day-structured fields (we'll filter to CS days later)
    buckets = {'CS_times','US_times','Ca_peaks','Ca_ts','pos','csus15','csus30','csus45','csus60','csus90'};
    out = {};
    for b = 1:numel(buckets)
        if isfield(rat, buckets{b}) && isstruct(rat.(buckets{b}))
            out = [out; fieldnames(rat.(buckets{b}))]; %#ok<AGROW>
        end
    end
    out = unique(out);
end

function [keysCS, datesCS] = keepCSdays(keysIn)
% Keep only 'CS_YYYY_MM_DD' (ignore extinction etc.)
    keysCS = {};
    datesCS = [];
    for i = 1:numel(keysIn)
        k = keysIn{i};
        tok = regexp(k, '^CS_(\d{4})_(\d{2})_(\d{2})$', 'tokens', 'once');
        if ~isempty(tok)
            dnum = datenum(sprintf('%s-%s-%s', tok{1}, tok{2}, tok{3}));
            keysCS{end+1,1} = k; %#ok<AGROW>
            datesCS(end+1,1) = dnum; %#ok<AGROW>
        end
    end
end

function rate = perNeuronRatesInWindows(S_pk, cs_on, us_on, csWin, trWin, usWin, baseWin)
% S_pk: [cells x nEventsTimes] (NaN/0 padded event times in seconds)
% Returns fields with size [cells x nTrials] of rates (events/s)
    nC = size(S_pk,1); nT = numel(cs_on);
    rate.cs    = nan(nC,nT);
    rate.trace = nan(nC,nT);
    rate.us    = nan(nC,nT);
    rate.base  = nan(nC,nT);
    rate.csus_evtMask = false(nC,nT); % any event in [CS,US)?

    for t=1:nT
        tCS = cs_on(t); tUS = us_on(t);
        wCS  = tCS + csWin;                       % [CS + a, CS + b]
        wTR  = [tCS + trWin(1), tUS];             % [CS+0.50, US)
        wUS  = [tUS + usWin(1), tUS + usWin(2)];  % [US+0, US+0.15]
        wBL  = tCS + baseWin;                     % [CS-1, CS)

        for c=1:nC
            x = S_pk(c,:); x = x(~isnan(x) & x~=0);
            if isempty(x), continue; end
            nCS = sum(x>=wCS(1) & x<wCS(2));
            nTR = sum(x>=wTR(1) & x<wTR(2));
            nUS = sum(x>=wUS(1) & x<wUS(2));
            nBL = sum(x>=wBL(1) & x<wBL(2));

            rate.cs(c,t)    = nCS / max(1e-6, diff(wCS));
            rate.trace(c,t) = nTR / max(1e-6, diff(wTR));
            rate.us(c,t)    = nUS / max(1e-6, diff(wUS));
            rate.base(c,t)  = nBL / max(1e-6, diff(wBL));

            rate.csus_evtMask(c,t) = any(x>=tCS & x<tUS);
        end
    end
end

function [TT, TN, NN, trialBins] = popVecSim_acrossTrials(rate, ~, ~, ~, nBins)
% Compute population vector similarities across a day and bin by trial index.
% Returns TT, TN, NN as vectors of length nBins (or unbinned if nBins==1).
% - TT: trial–trial cosine (consecutive trials within each bin)
% - TN: trial–nontrial cosine (trials within each bin vs their baseline)
% - NN: nontrial–nontrial cosine (consecutive baselines within each bin)

if nargin < 5, nBins = 5; end

nT = size(rate.cs,2);

% ----- build population vectors (rows = trials) -----
V_trials    = rate.cs';
V_trials    = V_trials ./ max(vecnorm(V_trials,2,2), eps);

V_nontrials = rate.base';
V_nontrials = V_nontrials ./ max(vecnorm(V_nontrials,2,2), eps);

% ----- raw similarities -----
TT_raw = nan(nT-1,1);
for i = 1:nT-1
    TT_raw(i) = dot(V_trials(i,:), V_trials(i+1,:));
end

TN_raw = nan(nT,1);
for i = 1:nT
    TN_raw(i) = dot(V_trials(i,:), V_nontrials(i,:));
end

NN_raw = nan(nT-1,1);
for i = 1:nT-1
    NN_raw(i) = dot(V_nontrials(i,:), V_nontrials(i+1,:));
end

% ----- binning -----
if nBins > 1 && nT >= nBins
    edges = round(linspace(1, nT+1, nBins+1));   % [1 ... nT+1]
    TT = nan(nBins,1); TN = nan(nBins,1); NN = nan(nBins,1);
    for b = 1:nBins
        tr1 = edges(b); tr2 = edges(b+1)-1;
        if tr1 > tr2, continue; end

        % TT/NN are on pairs starting within the bin
        if tr2-1 >= tr1
            pairIdx = tr1 : tr2-1;
            TT(b) = mean(TT_raw(pairIdx), 'omitnan');
            NN(b) = mean(NN_raw(pairIdx), 'omitnan');
        end
        TN(b) = mean(TN_raw(tr1:tr2), 'omitnan');
    end
    trialBins = (1:nBins)';                      % for plotting
else
    TT = TT_raw; TN = TN_raw; NN = NN_raw;
    trialBins = (1:numel(TN))';
end
end

function plotBarWithDots(~, dat, labels, ttl)
% dat: [nDays x 2] (col1 early, col2 late)
    m = [mean(dat(:,1),'omitnan'), mean(dat(:,2),'omitnan')];
    b = bar(m); %#ok<NASGU>
    hold on
    xE = 0.85; xL = 2.15;
    plot(xE*ones(size(dat,1),1), dat(:,1), 'ko','MarkerFaceColor','k');
    plot(xL*ones(size(dat,1),1), dat(:,2), 'ko','MarkerFaceColor','k');
    set(gca,'XTick',1:2,'XTickLabel',labels); ylabel('cosine');
    title(ttl); grid on
end
