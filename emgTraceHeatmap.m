function emgTraceHeatmap(csusNum, varargin)
% PLOTS EMG POWER

%emgTraceHeatmap(csusNum, 'RatNames', {'rat0222'}, 'Pre',0.5, 'Post',2, ...
%                 'K',5, 'ZscorePerDay',true, 'MakePlots',true, ...
%                 'EmgVarSuffix','emg')
%
% Mirrors neuralTraceHeatmap but plots EMG instead of calcium.
% Assumes:
%   - CS/US times live in base var 'ratXXXX' (e.g., rat0222)
%   - EMG lives in base var 'ratXXXXemg' (or configurable suffix)
%       EMG data:   rat0222emg.EMG.EMG_YYYY_MM_DD          -> vector [1 x nTime] (or [nTime x 1])
%       EMG time:   rat0222emg.EMG_ts.EMGts_YYYY_MM_DD     -> vector time axis (s)
%
% Panels per RAT:
%   (1) Trial x Time heatmap of EMG (hot colormap), averaged across last 3 CS days
%   (2) “PV” similarity across trials (cosine between per‑trial EMG time‑segments)
%       TT: trial–trial (consecutive EMG segments in [CS,US))
%       TN: trial–nontrial (EMG segment [CS,US) vs baseline [−1,0))
%       NN: nontrial–nontrial (consecutive baselines)
%   (3) Early (first K) vs Late (last K) mean EMG traces (±SEM), with individual early traces
%
% Notes:
%   - Handles arbitrary EMG sampling; segments are interpolated onto a common
%     relative time axis per day for heatmaps and per‑trial vectors.
%   - If 'SkipFirst2022' is true, drops the first CS/US pair on 2022 dates.

% -------- inputs --------
p = inputParser;
p.addParameter('RatNames', {'rat0222'}, @(c)iscell(c)||isstring(c));
p.addParameter('Pre', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Post', 2.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('K', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('ZscorePerDay', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('MakePlots', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('EmgVarSuffix', 'emg', @(s)ischar(s)||isstring(s));
p.addParameter('SkipFirst2022', true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});

ratNames     = cellstr(p.Results.RatNames);
preWin       = p.Results.Pre;
postWin      = p.Results.Post;
K            = p.Results.K;
doZ          = p.Results.ZscorePerDay;
makePlots    = p.Results.MakePlots;
emgSuffix    = char(p.Results.EmgVarSuffix);
skip2022     = p.Results.SkipFirst2022;

for ii = 1:numel(ratNames)
    ratVar   = ratNames{ii};
    emgVar   = [ratVar emgSuffix];      % e.g., 'rat0222emg'
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratVar))
        warning('Variable %s not found in base workspace. Skipping.', ratVar);
        continue;
    end
    if ~evalin('base', sprintf('exist(''%s'',''var'')', emgVar))
        warning('Variable %s not found in base workspace. Skipping.', emgVar);
        continue;
    end
    rat    = evalin('base', ratVar);
    ratEMG = evalin('base', emgVar);

    fprintf('\n=== %s : within-session analysis (EMG) ===\n', ratVar);

    % ---- choose last 3 CS days up to An ----
    [csKeys, csDates] = keepCSdays(getAllDayKeys(rat));
    anDate = datenum(strrep(rat.An,'_','-'));
    useMask = csDates <= anDate;
    csKeysUpToAn   = csKeys(useMask);
    csDatesUpToAn  = csDates(useMask);
    if numel(csKeysUpToAn) < 3
        warning('%s: Not enough CS days up to An (%s). Skipping.', ratVar, rat.An);
        continue;
    end
    [~, ord] = sort(csDatesUpToAn, 'ascend');
    daysToUse = csKeysUpToAn(ord(end-2:end));   % last 3 CS days

    % sanity on csusXX existence (just like neural)
    csusField = sprintf('csus%d', csusNum);
    if ~isfield(rat, csusField)
        warning('%s missing field %s. Skipping...', ratVar, csusField);
        continue;
    end

    % ----- collectors -----
    H_days = {}; tRel_days = {}; usLag_days = [];           % Panel 1
    TT_days = {}; TN_days = {}; NN_days = {};               % Panel 2
    earlyTrace_days = {}; lateTrace_days = {};              % Panel 3 (per-day mean traces)

    for d = 1:numel(daysToUse)
        dayKey = daysToUse{d};               % 'CS_YYYY_MM_DD'
        usKey  = regexprep(dayKey,'^CS_','US_');
        tok    = regexp(dayKey,'(\d{4}_\d{2}_\d{2})','tokens','once');
        if isempty(tok), fprintf('  [%s] bad key — skipping.\n', dayKey); continue; end
        dateTok = tok{1};

        % ---- CS/US times (from rat) ----
        if ~isfield(rat.CS_times, dayKey) || ~isfield(rat.US_times, usKey)
            fprintf('  [%s] missing CS/US times — skipping.\n', dayKey);
            continue;
        end
        cs_on = valueVec(rat.CS_times.(dayKey));
        us_on = valueVec(rat.US_times.(usKey));

        % optional: drop first trial for 2022 dates
        if skip2022 && contains(dayKey, '2022')
            if numel(cs_on)>1 && numel(us_on)>1
                cs_on = cs_on(2:end);
                us_on = us_on(2:end);
            end
        end

        nT = min(numel(cs_on), numel(us_on));
        if nT < 6, fprintf('  [%s] only %d trials — skipping.\n', dayKey, nT); continue; end
        cs_on = cs_on(1:nT); us_on = us_on(1:nT);
        good  = (us_on - cs_on) > 0.1;
        cs_on = cs_on(good);  us_on = us_on(good);
        nT    = numel(cs_on);
        if nT < 6, fprintf('  [%s] <6 valid trials — skipping.\n', dayKey); continue; end

        % ---- EMG data & time (from ratEMG) ----

        [emgVec, tEmg] = getDayEMG(ratEMG, dateTok);




        if isempty(emgVec) || isempty(tEmg) || numel(tEmg) < 2
            fprintf('  [%s] no EMG/ts for %s — skipping.\n', dayKey, dateTok);
            continue;
        end

        emgVec = emgVec(:)';      % row
        tEmg   = tEmg(:);         % column
        % optional per-day z-score
        if doZ
            mu = mean(emgVec,'omitnan'); sd = std(emgVec,[],'omitnan'); if sd==0, sd=1; end
            emgVec = (emgVec - mu) / sd;
        end

        % ===== Panel 1: heatmap (trials x time) of EMG =====
        [H_day, tRel, usLag] = makeEMGTrialHeatmap(emgVec, tEmg, cs_on, us_on, preWin, postWin);
        H_days{end+1}    = H_day;       %#ok<AGROW>
        tRel_days{end+1} = tRel;        %#ok<AGROW>
        usLag_days(end+1)= usLag;       %#ok<AGROW>

        % per-day early/late mean traces for Panel 3
        K_use = min(K, size(H_day,1));
        earlyTrace_days{end+1} = mean(H_day(1:K_use ,:), 1, 'omitnan'); %#ok<AGROW>
        lateTrace_days{end+1}  = mean(H_day(end-K_use+1:end,:), 1, 'omitnan'); %#ok<AGROW>

        % ===== Panel 2: EMG “PV” similarity across trials =====
        % define vectors as EMG samples in [CS,US) (trial) and baseline [-1,0) (nontrial),
        % interpolated to a fixed relative axis per day for consistent cosine.
        [TT_raw, TN_raw, NN_raw] = emgPvSimilarity(emgVec, tEmg, cs_on, us_on);
        if any(isfinite(TT_raw)), TT_days{end+1} = TT_raw; end %#ok<AGROW>
        if any(isfinite(TN_raw)), TN_days{end+1} = TN_raw; end %#ok<AGROW>
        if any(isfinite(NN_raw)), NN_days{end+1} = NN_raw; end %#ok<AGROW>
    end

    if ~makePlots, continue; end
    if isempty(H_days)
        fprintf('  No valid days to plot for %s.\n', ratVar);
        continue;
    end

    % ===== FIGURE per RAT =====
    figure('Name', sprintf('%s | EMG heatmap + similarity + early/late', ratVar), ...
           'Color','w', 'Position', [80 80 1200 420]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % -------- Panel 1: Heatmap averaged across days --------
    nexttile; hold on;
    [H_mean, tRel_mean] = meanHeatmapAcrossDays(H_days, tRel_days);

    imagesc(tRel_mean, 1:size(H_mean,1), H_mean);
    set(gca,'YDir','normal');
    colormap(gca, hot);

    colorbar;

    clim([nanmin(H_mean(:)) nanmean(H_mean(:))+(7*nanstd(H_mean(:)))]);

    xline(0, '--k', 'CS', 'LabelVerticalAlignment','bottom');
    if ~isempty(usLag_days)
        xline(nanmean(usLag_days), ':k', 'US (mean)', 'LabelVerticalAlignment','bottom');
    end
    xlabel('Time from CS (s)'); ylabel('Trial #');
    title(sprintf('%s — EMG POWER heatmap (z=%d) ', ratVar, doZ));

    % -------- Panel 2: EMG similarity vs trial (mean ± SEM across days) --------

    %nexttile; hold on; grid on;
    %[TTm, TTsem, xTT] = meanCurveAcrossDays(TT_days);
    %[TNm, TNsem, xTN] = meanCurveAcrossDays(TN_days);
    %[NNm, NNsem, xNN] = meanCurveAcrossDays(NN_days);
    %co = lines(3);
    %ebfill(xTT, TTm, TTsem, co(1,:)); plot(xTT, TTm, 'o-', 'LineWidth',2, 'Color',co(1,:), 'DisplayName','TT');
    %ebfill(xTN, TNm, TNsem, co(2,:)); plot(xTN, TNm, 's-', 'LineWidth',2, 'Color',co(2,:), 'DisplayName','TN');
    %ebfill(xNN, NNm, NNsem, co(3,:)); plot(xNN, NNm, '^-', 'LineWidth',2, 'Color',co(3,:), 'DisplayName','NN');
    %xlabel('Trial (or pair index)'); ylabel('Cosine similarity');
    %title('EMG similarity across trials (mean ± SEM across days)');
    %legend('Location','best');

    % -------- Panel 3: Early vs Late mean EMG (±SEM) + individual early trials --------
    nexttile; hold on; grid on;
    [Emean, Esem] = meanTracesAcrossDays(earlyTrace_days, tRel_days, tRel_mean);
    [Lmean, Lsem] = meanTracesAcrossDays(lateTrace_days,  tRel_days, tRel_mean);

    % Also draw the individual early trial means across days for the first up to 5 trials

    Kplot = min(5, min(numel(earlyTrace_days), K));  % same K cap as before
    if Kplot>0
        cEarly = [0 0.45 0.74];
        % reconstruct per-day trial matrices from H_days to pick individual early trials:
        % (We can reuse H_days directly since each H_day is trials x time on its own tRel)
        for et = 1:Kplot
            % collect trial et across days (skip if day has < et trials)
            D = nan(numel(H_days), numel(tRel_mean));
            have = false(numel(H_days),1);
            for d=1:numel(H_days)
                H = H_days{d}; tr = tRel_days{d};
                if size(H,1) >= et
                    yi = interp1(tr, H(et,:), tRel_mean, 'linear','extrap');
                    D(d,:) = yi; have(d)=true;
                end
            end
            if any(have)
                mi = mean(D(have,:), 1, 'omitnan');
      %          plot(tRel_mean, mi-nanmean(mi), '-', 'Color', cEarly, 'LineWidth', 1, 'HandleVisibility','off');
            end
        end
    end

    c = [0 0.45 0.74];
    ebfill(tRel_mean, Emean, Esem, lighten(c,0.6));
    plot(tRel_mean, Emean, '-', 'LineWidth',2, 'Color',c, 'DisplayName', sprintf('Early mean (first %d)',K));

    c = [0.45 0.74 0];
    ebfill(tRel_mean, Lmean, Lsem, lighten(c,0.4));
    plot(tRel_mean, Lmean, '-',  'LineWidth',2, 'Color',c, 'DisplayName', sprintf('Late mean (last %d)',K));

    xline(0, '--k');
    if ~isempty(usLag_days), xline(nanmean(usLag_days), ':k'); end
    xlabel('Time from CS (s)'); ylabel('EMG (a.u., z if enabled)');
    title('EMG: early vs late (±SEM) with individual early trials');
    legend('Location','best');

    title(tl, sprintf('%s — last 3 CS days up to %s', ratVar, rat.An));
end
end

% ===================== helpers =====================

function [emgVec, tEmg] = getDayEMG(ratEMG, dateTok)
    emgVec = []; tEmg = [];
    emgFN = findFieldWithDate(ratEMG.EMG,   dateTok);   % 'EMG_YYYY_MM_DD'
    tsFN  = findFieldWithDate(ratEMG.EMG_ts,dateTok);   % 'EMGts_YYYY_MM_DD'
    if isempty(emgFN) || isempty(tsFN), return; end

    v = double(ratEMG.EMG.(emgFN)(:).');  % row
    t = double(ratEMG.EMG_ts.(tsFN)(:));  % col (sec)

    % keep time strictly increasing and keep matching samples in v
    [t, uniqIdx] = unique(t, 'stable');
    v = v(uniqIdx).';

    % drop trailing NaNs safely (time or signal)
    nanIdx = find(isnan(t) | isnan(v), 1, 'first');
    if ~isempty(nanIdx)
        t = t(1:nanIdx-1);
        v = v(1:nanIdx-1);
    end
    if numel(t) < 2 || numel(v) ~= numel(t), return; end

    % --- power transform (rectified-squared) ---
    emgVec = v.^2;             % <-- square the SIGNAL, not the time
    % Optional: smooth power (e.g., 50 ms boxcar)
%     dt = median(diff(t), 'omitnan');
%     if isfinite(dt) && dt > 0
%         w = max(1, round(0.005/dt));
%         emgVec = movmean(emgVec, w, 'omitnan');
%     end

    tEmg = t;
end


function [H, tRel, usLag] = makeEMGTrialHeatmap(emgVec, tEmg, cs_on, us_on, pre, post)
% Build trials x time EMG matrix aligned to CS, spanning [-pre .. (US+post)].
    usLag = nanmean(us_on - cs_on);
    T_afterCS = max(us_on - cs_on);
    T_end = T_afterCS + post;
    % Choose a relative axis resolution close to native sampling (median dt)
    dt = median(diff(tEmg),'omitnan'); if ~isfinite(dt) || dt<=0, dt = 1/500; end
    tRel = (-pre : dt : T_end)';
    H = nan(numel(cs_on), numel(tRel));
    % Interpolate EMG onto each trial window
    for k=1:numel(cs_on)
        tSeg = cs_on(k) + tRel;
        H(k,:) = interp1(tEmg, emgVec, tSeg, 'linear', 'extrap');
    end
end

function [TT_raw, TN_raw, NN_raw] = emgPvSimilarity(emgVec, tEmg, cs_on, us_on)
% Cosine similarity using time‑sample vectors:
%   Vt(i,:) = EMG samples in [CS_i, US_i)
%   Vn(i,:) = EMG samples in [CS_i-1, CS_i)
% Interpolate each trial to a fixed number of samples so dot products are valid.
    nT = numel(cs_on);
    if nT<2, TT_raw=[]; TN_raw=[]; NN_raw=[]; return; end

    % choose per‑trial relative axes
    % set a fixed number of samples for vectors for numerical stability
    ns = 200;  % adjust if you want finer/coarser
    Vt = nan(nT, ns);
    Vn = nan(nT, ns);

    for i=1:nT
        tCS = cs_on(i); tUS = us_on(i);
        if ~(isfinite(tCS)&&isfinite(tUS)&& tUS>tCS), continue; end

        % trial segment [CS,US)
        t_rel_t = linspace(0, tUS-tCS, ns);
        Vt(i,:) = interp1(tEmg, emgVec, tCS + t_rel_t, 'linear','extrap');

        % baseline segment [CS-1, CS)
        t_rel_n = linspace(-1, 0, ns);
        Vn(i,:) = interp1(tEmg, emgVec, tCS + t_rel_n, 'linear','extrap');
    end

    % L2 normalize rows
    Vt = Vt ./ max(vecnorm(Vt,2,2), eps);
    Vn = Vn ./ max(vecnorm(Vn,2,2), eps);

    % similarities
    TT_raw = nan(nT-1,1);
    for i=1:nT-1, TT_raw(i) = dot(Vt(i,:), Vt(i+1,:)); end

    TN_raw = nan(nT,1);
    for i=1:nT,   TN_raw(i) = dot(Vt(i,:), Vn(i,:));   end

    NN_raw = nan(nT-1,1);
    for i=1:nT-1, NN_raw(i) = dot(Vn(i,:), Vn(i+1,:)); end
end

% ===== shared helpers (borrowed from your neural file) =====

function [Hmean, tRel_mean] = meanHeatmapAcrossDays(H_days, tRel_days)
    [~,ix] = max(cellfun(@numel, tRel_days));
    tRel_mean = tRel_days{ix};
    maxT = max(cellfun(@(H) size(H,1), H_days));
    Hstack = nan(maxT, numel(tRel_mean), numel(H_days));
    for d=1:numel(H_days)
        H = H_days{d};
        tRel = tRel_days{d};
        Hrg = nan(size(H,1), numel(tRel_mean));
        for tr=1:size(H,1)
            %Hrg(tr,:) = (interp1(tRel, H(tr,:), tRel_mean, 'linear','extrap')).^2;
            Hrg(tr,:) = (interp1(tRel, H(tr,:), tRel_mean, 'linear','extrap'));

        end
        Hstack(1:size(H,1),:,d) = Hrg;
    end
    Hmean = mean(Hstack, 3, 'omitnan');
end

function [m, s, x] = meanCurveAcrossDays(CURVES)
    if isempty(CURVES), m=[]; s=[]; x=[]; return; end
    Lmax = max(cellfun(@numel, CURVES));
    M = nan(numel(CURVES), Lmax);
    for d=1:numel(CURVES)
        v = CURVES{d}(:)'; M(d,1:numel(v)) = v;
    end
    m = mean(M,1,'omitnan')'; s = sem(M); x = (1:Lmax)';
end

function [m, s] = meanTracesAcrossDays(TRACES, tRel_days, tRef)
    if isempty(TRACES), m = []; s = []; return; end
    nD = numel(TRACES);
    A = nan(nD, numel(tRef));
    for d=1:nD
        y = TRACES{d};
        tr = tRel_days{d};
        yi = interp1(tr, y, tRef, 'linear', 'extrap');
        A(d,:) = yi;
    end
    m = mean(A,1,'omitnan'); s = sem(A);
end

function s = sem(X)
    if isvector(X)
        x = X(:); n = sum(isfinite(x)); s = std(x,'omitnan')/max(1,sqrt(n));
    else
        n = sum(isfinite(X),1);
        s = std(X,0,1,'omitnan') ./ max(1,sqrt(n)); s = s(:);
    end
end

function ebfill(x, m, se, col)
    if isempty(x) || isempty(m) || isempty(se), return; end
    xv = [x(:); flipud(x(:))];
    yv = [m(:)-se(:); flipud(m(:)+se(:))];
    patch('XData',xv,'YData',yv,'FaceColor',col,'FaceAlpha',0.15,'EdgeColor','none');
end

function c2 = lighten(c, frac)
    c2 = c + (1-c)*min(max(frac,0),1);
end

function out = getAllDayKeys(rat)
    buckets = {'CS_times','US_times','Ca_peaks','Ca_ts','Ca_traces','pos','csus15','csus30','csus45','csus60','csus90'};
    out = {};
    for b = 1:numel(buckets)
        if isfield(rat, buckets{b}) && isstruct(rat.(buckets{b}))
            out = [out; fieldnames(rat.(buckets{b}))]; %#ok<AGROW>
        end
    end
    out = unique(out);
end

function [keysCS, datesCS] = keepCSdays(keysIn)
    keysCS = {}; datesCS = [];
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

function fn = findFieldWithDate(S, dateTok)
    fn = ''; if ~isstruct(S), return; end
    f = fieldnames(S); pat = strrep(dateTok,'_','[_-]?');
    for i = 1:numel(f)
        if ~isempty(regexp(f{i}, pat, 'once')), fn = f{i}; return; end
    end
end

function out = valueVec(x)
    if isnumeric(x), out = x(:)'; return; end
    if isstruct(x)
        f = fieldnames(x);
        for k = 1:numel(f)
            if isnumeric(x.(f{k})), out = x.(f{k})(: )'; return; end
        end
    end
    error('CS/US times are not numeric.');
end
