function R = speedDependence_taskVsRun_population(ratNames, varargin)
% speedDependence_taskVsRun_population
% Compare POPULATION speed dependence in TRACE (task) vs NON-TASK running.
% Aggregates across cells (mean or sum), matches by speed *for stats only*,
% and now stores BOTH pre-match ("all") and post-match ("matched") samples.
%
% Usage:
%   R = speedDependence_taskVsRun_population({'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%        'SpeedBinWidth',2,'MinDurPerBin',1,'MinBins',5,'binSize',1/7.5, ...
%        'SpeedThresh',4,'PopAgg','mean','PopNorm','none','DoPlots',true, ...
%        'MinBinsPolicy','lenient','ScatterUse','all');
%
% Name-Value options (key):
%   'Win'            [0 2]
%   'binSize'        1
%   'SpeedThresh'    4
%   'SpeedBinWidth'  1
%   'MinDurPerBin'   5
%   'MinBins'        1
%   'MinBinsPolicy'  'lenient'   % 'day' | 'lenient'
%   'PopAgg'         'mean'      % 'mean' | 'sum'
%   'PopNorm'        'none'      % 'none' | 'mean' | 'demean'
%   'CapUseSpeedGatedTask' false % if true, count task spikes only when v>=SpeedThresh
%   'DoPlots'        true
%   'ScatterUse'     'all'       % 'all' | 'matched'  (what to show in scatter figures)
%
% Outputs (key fields):
%   R.r_task_combined, R.p_task_combined, R.r_run_combined, R.p_run_combined
%   R.pooled.r_task, R.pooled.r_run
%   R.stats.p_task_vs_run (per-rat), R.stats.p_task_vs_run_pooled (pooled)
%   R.raw.task/run{r}{d}.speed/pop          (matched)
%   R.raw.task/run{r}{d}.speed_all/pop_all  (pre-match)

% ---------- args ----------
p = inputParser;
addParameter(p,'Win',[0 2]);
addParameter(p,'binSize',1);
addParameter(p,'SpeedThresh',0);
addParameter(p,'SpeedBinWidth',3);
addParameter(p,'MinDurPerBin',1/7.5);
addParameter(p,'MinBins',1);
addParameter(p,'MinBinsPolicy','day',@(s) any(validatestring(lower(s),{'day','lenient'})));
addParameter(p,'PopAgg','mean',@(s) any(validatestring(lower(s),{'mean','sum'})));
addParameter(p,'PopNorm','none',@(s) any(validatestring(lower(s),{'none','mean','demean'})));
addParameter(p,'CapUseSpeedGatedTask',false);
addParameter(p,'DoPlots',true);
addParameter(p,'ScatterUse','all',@(s) any(validatestring(lower(s),{'all','matched'})));
parse(p,varargin{:});

Win           = p.Results.Win;
binSize       = p.Results.binSize;
vThresh       = p.Results.SpeedThresh;
SBW           = p.Results.SpeedBinWidth;
minDurPerBin  = p.Results.MinDurPerBin;
minBins       = p.Results.MinBins;
minBinsPolicy = lower(p.Results.MinBinsPolicy);
popAgg        = lower(p.Results.PopAgg);
popNorm       = lower(p.Results.PopNorm);
capTaskGate   = p.Results.CapUseSpeedGatedTask;
doplot        = p.Results.DoPlots;
scatterUse    = lower(p.Results.ScatterUse);

nR = numel(ratNames);

% Per-rat containers (each cell = vector over days)
R = struct( ...
  'r_task',{cell(1,nR)}, 'p_task',{cell(1,nR)}, 'n_task',{cell(1,nR)}, ...
  'r_run' ,{cell(1,nR)}, 'p_run' ,{cell(1,nR)}, 'n_run' ,{cell(1,nR)}, ...
  'r_task_combined',nan(nR,1), 'p_task_combined',nan(nR,1), ...
  'r_run_combined' ,nan(nR,1), 'p_run_combined' ,nan(nR,1), ...
  'rateRatio_match',{cell(1,nR)}, ...
  'stats',struct(), 'pooled',struct(), 'meta',[]);

% Raw x–y capture per rat/day (for plotting)
R.raw = struct('task',{cell(1,nR)}, 'run',{cell(1,nR)});

% Pooled collectors (across rats/days)
pool_task_r = []; pool_task_n = [];
pool_run_r  = []; pool_run_n  = [];
pool_rratio = [];

for r = 1:nR
    rat = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    iAn   = find(strcmp(dates,rat.An),1);
    if isempty(iAn) || iAn < 3
        warning('%s: not enough days before An; skipping.', ratNames{r});
        continue;
    end
    days = dates(iAn-2:iAn);  % last 3 days up to An

    r_task_d = nan(numel(days),1);
    p_task_d = nan(numel(days),1);
    n_task_d = nan(numel(days),1);

    r_run_d  = nan(numel(days),1);
    p_run_d  = nan(numel(days),1);
    n_run_d  = nan(numel(days),1);

    rratio_d = nan(numel(days),1);  % population speed-matched rate ratio

    if isempty(R.raw.task{r}), R.raw.task{r} = cell(numel(days),1); end
    if isempty(R.raw.run{r}),  R.raw.run{r}  = cell(numel(days),1); end

    for d = 1:numel(days)
        D = days{d};

        % ---------- data ----------
        Sraw = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));   % n×? cell or matrix
        S    = normalizeCApeaks(Sraw);                    % {n,1} spike vectors

        pos = smoothpos(rat.pos.(sprintf('pos_%s',D)));   % match speedBinMatched
        [vt, vv] = velocityFromPos(pos);

        CS  = rat.CS_times.(sprintf('CS_%s',D));          % CS onset times

        % ---------- task side (trial windows split into binSize chunks) ----------
        if capTaskGate
            [trialSpeed, trialRatePerCell, trialEffDur] = taskBinsSpeedAndRates_withDur( ...
                S, CS, pos(:,1), vt, vv, Win, binSize, vThresh);
        else
            [trialSpeed, trialRatePerCell, trialEffDur] = taskBinsSpeedAndRates_withDur( ...
                S, CS, pos(:,1), vt, vv, Win, binSize, []);   % no speed gate
        end
        trialRatePerCell_use = normalizePerCell(trialRatePerCell, popNorm);
        trialPop = aggregateAcrossCells(trialRatePerCell_use, popAgg);

        % ---------- non-task side: bins of length binSize, v>=vThresh ----------
        [runSpeed, runRatePerCell, runEffDur] = runBinsSpeedAndRates_withDur( ...
            S, CS, pos(:,1), vt, vv, Win, binSize, vThresh);
        runRatePerCell_use = normalizePerCell(runRatePerCell, popNorm);
        runPop = aggregateAcrossCells(runRatePerCell_use, popAgg);

        if isempty(trialSpeed) || isempty(runSpeed)
            continue;
        end

        % ---- capture ALL (pre-matching) samples for plotting ----
        t_ok_all = isfinite(trialSpeed) & isfinite(trialPop);
        r_ok_all = isfinite(runSpeed)  & isfinite(runPop);

        t_speed_all = trialSpeed(t_ok_all);
        t_pop_all   = trialPop(t_ok_all);

        r_speed_all = runSpeed(r_ok_all);
        r_pop_all   = runPop(r_ok_all);

        % ---------- SPEED-BIN matched selection (for stats) ----------
        [taskMask, runMaskSel, qualBins] = selectBySpeedBins_secondsMatched( ...
            trialSpeed, trialEffDur, runSpeed, runEffDur, SBW, minDurPerBin, minBins);

        % Strict day-level enforcement if requested
        if strcmp(minBinsPolicy,'day') && qualBins < minBins
            % still store the ALL samples for plotting context
            R.raw.task{r}{d} = struct('speed',[], 'pop',[], ...
                                      'speed_all',t_speed_all,'pop_all',t_pop_all);
            R.raw.run{r}{d}  = struct('speed',[], 'pop',[], ...
                                      'speed_all',r_speed_all,'pop_all',r_pop_all);
            continue;
        end

        % Apply masks to get MATCHED samples
        t_ok = taskMask & isfinite(trialSpeed) & isfinite(trialPop);
        r_ok = runMaskSel & isfinite(runSpeed) & isfinite(runPop);

        if sum(t_ok) < 3 || sum(r_ok) < 3
            % store ALL samples, but skip stats this day
            R.raw.task{r}{d} = struct('speed',[], 'pop',[], ...
                                      'speed_all',t_speed_all,'pop_all',t_pop_all);
            R.raw.run{r}{d}  = struct('speed',[], 'pop',[], ...
                                      'speed_all',r_speed_all,'pop_all',r_pop_all);
            continue;
        end

        t_speed = trialSpeed(t_ok);
        t_pop   = trialPop(t_ok);
        r_speed = runSpeed(r_ok);
        r_pop   = runPop(r_ok);

        % --------- capture BOTH matched (speed/pop) AND pre-match (_all) ---------
        R.raw.task{r}{d} = struct('speed',t_speed,'pop',t_pop, ...
                                  'speed_all',t_speed_all,'pop_all',t_pop_all);
        R.raw.run{r}{d}  = struct('speed',r_speed,'pop',r_pop, ...
                                  'speed_all',r_speed_all,'pop_all',r_pop_all);

        % ---------- correlations (population; on MATCHED samples) ----------
        [rt, pt, nt] = corr_with_counts(t_speed, t_pop); % task
        [rr, pr, nr] = corr_with_counts(r_speed, r_pop); % run

        r_task_d(d) = rt;  p_task_d(d) = pt;  n_task_d(d) = nt;
        r_run_d(d)  = rr;  p_run_d(d)  = pr;  n_run_d(d)  = nr;

        % ---------- population speed-matched rate ratio ----------
        rratio_d(d) = populationSpeedMatchedRateRatio_nearest( ...
            t_speed, t_pop, r_speed, r_pop);
    end

    % store per rat (per-day vectors)
    R.r_task{r} = r_task_d;  R.p_task{r} = p_task_d;  R.n_task{r} = n_task_d;
    R.r_run{r}  = r_run_d;   R.p_run{r}  = p_run_d;   R.n_run{r}  = n_run_d;
    R.rateRatio_match{r} = rratio_d;

    % combine per rat across its days (fixed-effects Fisher-z, weighted by n-3)
    [R.r_task_combined(r), R.p_task_combined(r)] = fisherZ_fixed_effects(r_task_d, n_task_d);
    [R.r_run_combined(r) , R.p_run_combined(r) ] = fisherZ_fixed_effects(r_run_d , n_run_d );

    % pooled collectors for pooled FEZ
    okT = isfinite(r_task_d) & (n_task_d>=5);
    okR = isfinite(r_run_d)  & (n_run_d >=5);
    pool_task_r = [pool_task_r; r_task_d(okT)];
    pool_task_n = [pool_task_n; n_task_d(okT)];
    pool_run_r  = [pool_run_r ; r_run_d(okR)];
    pool_run_n  = [pool_run_n ; n_run_d(okR)];
    pool_rratio = [pool_rratio; rratio_d(isfinite(rratio_d))];
end

% ----- Task vs Run paired tests (per rat) and pooled, using MATCHED r -----
R.stats.p_task_vs_run = nan(nR,1);
R.stats.n_pairs       = nan(nR,1);

all_rt = [];   % collect day-level Task r across rats
all_rr = [];   % collect day-level Run  r across rats

for r = 1:nR
    if isempty(R.r_task{r}), continue; end
    rt = R.r_task{r}(:);
    rr = R.r_run{r}(:);
    ok = isfinite(rt) & isfinite(rr);           % paired valid days
    if sum(ok) >= 2
        [~, p] = ttest(rt(ok), rr(ok));         % two-sided paired
        R.stats.p_task_vs_run(r) = p;
        R.stats.n_pairs(r)       = sum(ok);

        all_rt = [all_rt; rt(ok)];
        all_rr = [all_rr; rr(ok)];
    end
end

% Pooled paired test across all day-level pairs from all rats
if numel(all_rt) >= 2
    [~, pPool] = ttest(all_rt, all_rr);         % paired across concatenated day-pairs
else
    pPool = NaN;
end
R.stats.p_task_vs_run_pooled = pPool;

% -------- pooled summaries across rats (FEZ on day-level r) --------
[R.pooled.r_task, R.pooled.p_task] = fisherZ_fixed_effects(pool_task_r, pool_task_n);
[R.pooled.r_run , R.pooled.p_run ] = fisherZ_fixed_effects(pool_run_r , pool_run_n);
R.pooled.rateRatio_match_geoMean   = exp(mean(log(pool_rratio(pool_rratio>0)), 'omitnan'));

% -------- meta --------
R.meta = struct('Win',Win,'binSize',binSize,'SpeedThresh',vThresh, ...
                'SpeedBinWidth',SBW,'MinDurPerBin',minDurPerBin,'MinBins',minBins, ...
                'MinBinsPolicy',minBinsPolicy, ...
                'ratNames',{ratNames}, 'PopAgg',popAgg,'PopNorm',popNorm, ...
                'CapUseSpeedGatedTask',capTaskGate, 'ScatterUse',scatterUse);

% -------- plots (scatter uses pre-match by default) --------
if doplot
    plotPopulationSummaryBars(R);
    plot_scatter_per_rat_and_pooled(R, 'task', 'ScatterUse', scatterUse);
    plot_all_rats_one_graph(R, 'task', 'ScatterUse', scatterUse);
    plot_scatter_per_rat_and_pooled(R, 'run',  'ScatterUse', scatterUse);
    plot_all_rats_one_graph(R, 'run',  'ScatterUse', scatterUse);
    plot_task_vs_run_per_rat_and_pooled(R, 'ScatterUse', scatterUse);
    paired_connected_plot(R, ...
       'Title','Population speed–activity correlation (per rat)', ...
       'YLabel','population r(speed, activity)', ...
       'XTickLabels',{'Run (Non-trial)','Task (Trial)'});
end
end

% ================== HELPERS ==================

function [vt, vv] = velocityFromPos(pos)
t = pos(:,1); xy = pos(:,2:3);
dt = diff(t); dt(dt==0)=NaN;
v  = [0; sqrt(sum(diff(xy).^2,2))./dt]; % cm/s
v(~isfinite(v)) = 0;
vt = t; vv = v;
end

function Scell = normalizeCApeaks(Sraw)
if iscell(Sraw)
    if isvector(Sraw), Scell = Sraw(:); return; end
    K = size(Sraw,2);
    scores = zeros(1,K);
    for j = 1:K
        col = Sraw(:,j);
        scores(j) = mean(cellfun(@(x) isnumeric(x) && isvector(x), col));
    end
    [bestScore, jbest] = max(scores);
    if bestScore == 0, error('normalizeCApeaks:NoNumericColumn','No numeric spike-time column.'); end
    Scell = Sraw(:, jbest);
    Scell = cellfun(@(v) v(:), Scell, 'uni', false);
    return;
end
if isnumeric(Sraw)
    n = size(Sraw,1);
    Scell = cell(n,1);
    for i = 1:n
        v = Sraw(i,:).';
        v = v(isfinite(v) & v>0);
        Scell{i,1} = v;
    end
    return;
end
error('normalizeCApeaks:BadType','Unsupported CA_peaks type: %s.', class(Sraw));
end

function [binSpeed, rateMat, effDur] = runBinsSpeedAndRates_withDur(S, CS, tPos, vt, vv, Win, binSec, vThr)
% Build non-task bins of length binSec with v>=vThr; count spikes only during masked time.
nCells = size(S,1);
% mask out trial windows
csMask = false(size(tPos));
for k = 1:numel(CS)
    csMask = csMask | ((tPos >= CS(k)+Win(1)) & (tPos < CS(k)+Win(2)));
end
dtPos = median(diff(tPos),'omitnan');
vOnPos = interp1(vt, vv, tPos, 'linear','extrap');
runMask = ~csMask & (vOnPos >= vThr);

tmin = tPos(find(runMask,1,'first'));
tmax = tPos(find(runMask,1,'last'));
if isempty(tmin) || isempty(tmax)
    binSpeed = []; rateMat = []; effDur = []; return
end

edges = tmin:binSec:tmax; if edges(end) < tmax, edges = [edges tmax]; end
nB = numel(edges)-1;

binSpeed = nan(nB,1);
effDur   = zeros(nB,1);
rateMat  = nan(nB,nCells);

for b = 1:nB
    iv = (tPos>=edges(b)) & (tPos<edges(b+1)) & runMask;   % masked samples in bin
    effDur(b) = sum(iv) * dtPos;
    if effDur(b) <= 0, continue; end
    binSpeed(b) = mean(vOnPos(iv),'omitnan');

    % count spikes only during masked time in this bin
    intervals = maskToIntervals(tPos, iv, dtPos);
    for c = 1:nCells
        sp = S{c,1};
        if isempty(sp), rateMat(b,c) = 0; continue; end
        nSp = countInIntervals(sp, intervals);
        rateMat(b,c) = nSp / effDur(b);  % Hz
    end
end

ok = isfinite(binSpeed) & (effDur>0) & any(isfinite(rateMat),2);
binSpeed = binSpeed(ok);
rateMat  = rateMat(ok,:);
effDur   = effDur(ok);
end

function intervals = maskToIntervals(t, mask, dt)
if nargin<3 || isempty(dt), dt = median(diff(t),'omitnan'); end
mask = mask(:); t = t(:);
if numel(t)~=numel(mask), error('maskToIntervals: size mismatch'); end
dm = diff([false; mask; false]);
starts = find(dm==1);
ends   = find(dm==-1)-1;
intervals = [t(starts)  t(ends)+dt];   % include sample; extend ~dt
end

function n = countInIntervals(spikes, intervals)
if isempty(spikes) || isempty(intervals), n = 0; return; end
spikes = spikes(:);
inside = false(size(spikes));
for k = 1:size(intervals,1)
    inside = inside | (spikes>=intervals(k,1) & spikes<intervals(k,2));
end
n = sum(inside);
end

function [binSpeed, rateMat, effDur] = taskBinsSpeedAndRates_withDur( ...
        S, CS, tPos, vt, vv, Win, binSec, vThr)
% Split each trial window [CS+Win(1), CS+Win(2)) into binSec chunks.
% If vThr is provided, count time/spikes only where speed>=vThr inside each chunk.

nCells  = numel(S);
dtPos   = median(diff(tPos),'omitnan');
vOnPos  = interp1(vt, vv, tPos, 'linear','extrap');

binSpeed = []; effDur = []; rateMat  = [];

for k = 1:numel(CS)
    t0 = CS(k)+Win(1);  t1 = CS(k)+Win(2);
    if ~isfinite(t0) || ~isfinite(t1) || t1<=t0, continue; end

    edges = t0:binSec:t1;
    if edges(end) < t1, edges = [edges t1]; end

    for b = 1:numel(edges)-1
        b0 = edges(b); b1 = edges(b+1);
        inBin = (tPos>=b0) & (tPos<b1);

        if isempty(vThr)
            mask = inBin;                    % full subwindow
        else
            mask = inBin & (vOnPos>=vThr);   % gated samples only
        end

        dur = sum(mask) * dtPos;
        if dur <= 0, continue; end

        sbin = mean(vOnPos(mask), 'omitnan');

        if isempty(vThr)
            intervals = [b0 b1];
        else
            intervals = maskToIntervals(tPos, mask, dtPos);
        end

        row = NaN(1,nCells);
        for c = 1:nCells
            sp = S{c};
            if isempty(sp), row(c) = 0; continue; end
            if isempty(vThr)
                nSp = sum(sp>=b0 & sp<b1);
            else
                nSp = countInIntervals(sp, intervals);
            end
            row(c) = nSp / dur;           % Hz
        end

        binSpeed = [binSpeed; sbin];
        effDur   = [effDur;   dur];
        rateMat  = [rateMat;  row];
    end
end

ok = isfinite(binSpeed) & (effDur>0) & any(isfinite(rateMat),2);
binSpeed = binSpeed(ok);
effDur   = effDur(ok);
rateMat  = rateMat(ok,:);
end

function Y = normalizePerCell(X, mode)
switch lower(mode)
    case 'none'
        Y = X;
    case 'mean'     % scale by per-cell mean
        mu = mean(X,1,'omitnan');
        mu(mu==0 | ~isfinite(mu)) = NaN;
        Y = X ./ mu;
    case 'demean'   % subtract per-cell mean
        mu = mean(X,1,'omitnan');
        Y = X - mu;
    otherwise
        error('PopNorm must be ''none'', ''mean'', or ''demean''.');
end
end

function y = aggregateAcrossCells(rateMat, modeStr)
switch lower(modeStr)
    case 'mean'
        y = mean(rateMat, 2, 'omitnan');
    case 'sum'
        y = sum(rateMat, 2, 'omitnan');
    otherwise
        error('PopAgg must be ''mean'' or ''sum''.');
end
end

function [r, p, n] = corr_with_counts(x, y)
ok = isfinite(x) & isfinite(y);
n  = sum(ok);
if n >= 5
    [r, p] = corr(x(ok), y(ok), 'type','Pearson');
else
    r = NaN; p = NaN;
end
end

function geo = populationSpeedMatchedRateRatio_nearest(t_speed, t_pop, r_speed, r_pop)
% For each task sample, match nearest RUN sample by speed; geo-mean ratio.
if isempty(t_speed) || isempty(r_speed)
    geo = NaN; return;
end
ratios = nan(numel(t_speed),1);
for k = 1:numel(t_speed)
    s = t_speed(k);
    if ~isfinite(s), continue; end
    [~, idx] = min(abs(r_speed - s));
    denom = r_pop(idx);
    nume  = t_pop(k);
    if isfinite(denom) && denom>0 && isfinite(nume)
        ratios(k) = nume/denom;
    end
end
geo = exp(mean(log(ratios(ratios>0 & isfinite(ratios))), 'omitnan'));
end

function [rComb, pComb] = fisherZ_fixed_effects(rVec, nVec)
ok = isfinite(rVec) & isfinite(nVec) & (nVec>=5);
rVec = rVec(ok); nVec = nVec(ok);
if isempty(rVec)
    rComb = NaN; pComb = NaN; return;
end
w   = max(nVec - 3, 1);
z   = atanh(rVec);
zbar= sum(w .* z) / sum(w);
rComb = tanh(zbar);
se = 1/sqrt(sum(w));
pComb = 2*normcdf(-abs(zbar)/se,0,1);
end

% ================== PLOTTING (scatter uses ALL by default) ==================

function plot_scatter_per_rat_and_pooled(R, side, varargin)
% One figure:
%   • one scatter per rat (all days concatenated) + ONE best-fit line per rat
%   • one final subplot with ALL rats pooled + ONE best-fit line
% side: 'task' | 'run'
q = inputParser; addParameter(q,'ScatterUse','all'); parse(q,varargin{:});
use = lower(q.Results.ScatterUse);

if nargin<2, side='task'; end
side = validatestring(lower(side), {'task','run'});
rats = R.meta.ratNames; nR = numel(rats);

[X, Y] = collect_xy(R, side, use);    % <-- uses 'all' or 'matched'

% layout: grid + pooled
nCols = 3; nRows = ceil((nR+1)/nCols);
f = figure('Color','w','Position',[120 120 1100 680]); %#ok<NASGU>
tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');

col = lines(nR);
for r=1:nR
    nexttile; hold on;
    x = X{r}; y = Y{r};
    if isempty(x), title(sprintf('%s (no data)',rats{r})); axis off; continue; end
    scatter(x, y, 10, col(r,:), 'filled', 'MarkerFaceAlpha',0.35, 'MarkerEdgeAlpha',0.10);
  %%%  set(gca,'YScale','log')
    p  = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 200);
    plot(xx, polyval(p, xx), '-', 'Color', col(r,:), 'LineWidth', 2);
    rr = corr(x, y, 'type','Pearson', 'rows','complete');
    title(sprintf('%s  r=%.3f (n=%d)', rats{r}, rr, numel(x)));
    xlabel('speed (cm/s)'); ylabel('population activity'); grid on; box off;
end

% pooled
nexttile; hold on;
XP = vertcat(X{:}); YP = vertcat(Y{:});
scatter(XP, YP, 8, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.08);
%%% set(gca,'YScale','log')
pp  = polyfit(XP, YP, 1);
xx  = linspace(min(XP), max(XP), 250);
plot(xx, polyval(pp, xx), 'k-', 'LineWidth', 2.4);
rP  = corr(XP, YP, 'type','Pearson', 'rows','complete');
title(sprintf('All rats  r=%.3f (n=%d) — %s', rP, numel(XP), use));
xlabel('speed (cm/s)'); ylabel('population activity'); grid on; box off;
sgtitle(sprintf('Speed vs population activity — %s (scatter uses %s)', upper(side), use));
end

function plot_all_rats_one_graph(R, side, varargin)
% ALL rats on ONE axes, one best-fit line per rat (+ mean line).
q = inputParser; addParameter(q,'ScatterUse','all'); parse(q,varargin{:});
use = lower(q.Results.ScatterUse);

if nargin<2, side='task'; end
side = validatestring(lower(side), {'task','run'});
rats = R.meta.ratNames; nR = numel(rats);

[X, Y] = collect_xy(R, side, use);

% global x-range for consistent lines
allX = vertcat(X{:}); if isempty(allX), return; end
xmin = min(allX); xmax = max(allX);
xx   = linspace(xmin, xmax, 300);

col = lines(nR);
P = nan(nR,2);   % [slope intercept] per rat

figure('Color','w','Position',[140 140 800 520]); hold on;
for r=1:nR
    x = X{r}; y = Y{r};
    if isempty(x), continue; end
    scatter(x, y, 8, col(r,:), 'filled', 'MarkerFaceAlpha',0.12, 'MarkerEdgeAlpha',0.05, 'HandleVisibility','off');
%%%    set(gca,'YScale','log')
    p  = polyfit(x, y, 1);
    P(r,:) = p;
    plot(xx, polyval(p, xx), '-', 'Color', col(r,:), 'LineWidth', 2, 'DisplayName', rats{r});
end

P = P(all(isfinite(P),2),:);
if ~isempty(P)
    pMean = mean(P,1);
    plot(xx, polyval(pMean, xx), 'r--', 'LineWidth', 2.5, 'DisplayName', 'Mean across rats');
%%%    set(gca,'YScale','log')
end

xlabel('speed (cm/s)'); ylabel('population activity');
title(sprintf('All rats — %s (scatter uses %s)', upper(side), use));
grid on; box off; legend('Location','best'); legend boxoff;
end

function plot_task_vs_run_per_rat_and_pooled(R, varargin)
% For each rat: plot TASK and RUN on the same axes (+ pooled).
q = inputParser; addParameter(q,'ScatterUse','all'); parse(q,varargin{:});
use = lower(q.Results.ScatterUse);

rats = R.meta.ratNames; nR = numel(rats);
[XT, YT] = collect_xy(R, 'task', use);
[XR, YR] = collect_xy(R, 'run',  use);

nCols = 3; nRows = ceil((nR+1)/nCols);
figure('Color','w','Position',[100 100 1200 700]);
tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');

cols = [0.2 0.45 0.85; 0.9 0.4 0.2]; % task=blue, run=red
for r=1:nR
    nexttile; hold on;
    xt = XT{r}; yt = YT{r};
    xr = XR{r}; yr = YR{r};
    if isempty(xt) || isempty(xr)
        title(sprintf('%s (no data)', rats{r})); axis off; continue;
    end
    scatter(xt, yt, 10, cols(1,:), 'filled', 'MarkerFaceAlpha',0.3);
    scatter(xr, yr, 10, cols(2,:), 'filled', 'MarkerFaceAlpha',0.3);
%%%    set(gca,'YScale','log')

    pt = polyfit(xt, yt, 1);
    pr = polyfit(xr, yr, 1);
    xx = linspace(min([xt;xr]), max([xt;xr]), 200);
    plot(xx, polyval(pt,xx), '-', 'Color', cols(1,:), 'LineWidth',2);
    plot(xx, polyval(pr,xx), '-', 'Color', cols(2,:), 'LineWidth',2);

    [rt, ptval] = corr(xt, yt, 'type','Pearson', 'rows','complete');
    [rr, prval] = corr(xr, yr, 'type','Pearson', 'rows','complete');

    title(sprintf('%s — %s\nTask: r=%.3f, p=%.3g\nRun:  r=%.3f, p=%.3g', ...
        rats{r}, use, rt, ptval, rr, prval));
    xlabel('speed (cm/s)'); ylabel('population activity');
    grid on; box off;
end

% pooled
nexttile; hold on;
XP = vertcat(XT{:}); YP = vertcat(YT{:});
XRp= vertcat(XR{:}); YRp= vertcat(YR{:});
scatter(XP, YP, 8, cols(1,:), 'filled', 'MarkerFaceAlpha',0.25);
scatter(XRp,YRp, 8, cols(2,:), 'filled', 'MarkerFaceAlpha',0.25);
%%% set(gca,'YScale','log')

pt = polyfit(XP, YP, 1); pr = polyfit(XRp, YRp, 1);
xx = linspace(min([XP;XRp]), max([XP;XRp]), 250);
plot(xx, polyval(pt,xx), '-', 'Color', cols(1,:), 'LineWidth',2.5);
plot(xx, polyval(pr,xx), '-', 'Color', cols(2,:), 'LineWidth',2.5);

[rt, ptval] = corr(XP, YP, 'type','Pearson','rows','complete');
[rr, prval] = corr(XRp, YRp, 'type','Pearson','rows','complete');

title(sprintf('Pooled — %s\nTask: r=%.3f, p=%.3g\nRun:  r=%.3f, p=%.3g', use, rt, ptval, rr, prval));
xlabel('speed (cm/s)'); ylabel('population activity'); grid on; box off;
sgtitle(sprintf('Task vs Run: speed–population correlations (scatter uses %s)', use));
end

function [X, Y] = collect_xy(R, side, use)
% Helper to pull x,y per rat for 'all' or 'matched' samples.
rats = R.meta.ratNames; nR = numel(rats);
X = cell(nR,1); Y = cell(nR,1);
for r=1:nR
    X{r} = []; Y{r} = [];
    if r>numel(R.raw.(side)) || isempty(R.raw.(side){r}), continue; end
    for d=1:numel(R.raw.(side){r})
        S = R.raw.(side){r}{d};
        if isempty(S), continue; end
        switch use
            case 'all'
                if ~isfield(S,'speed_all'), continue; end
                x = S.speed_all; y = S.pop_all;
            otherwise % 'matched'
                x = S.speed; y = S.pop;
        end
        ok = isfinite(x) & isfinite(y);
        X{r} = [X{r}; x(ok)];
        Y{r} = [Y{r}; y(ok)];
    end
end
end

% ===== other plotting/stat helpers (unchanged except where called) =====

function s = p2stars_local(p)
if ~isfinite(p), s = ''; return; end
if p < 1e-3, s = '***';
elseif p < 1e-2, s = '**';
elseif p < 5e-2, s = '*';
else, s = '';
end
end

function plotPopulationSummaryBars(R)
xlab = [R.meta.ratNames, {'pooled'}];
task = [R.r_task_combined(:); R.pooled.r_task];
run  = [R.r_run_combined(:) ; R.pooled.r_run];

pvec = nan(numel(xlab),1);
if isfield(R,'stats') && isfield(R.stats,'p_task_vs_run')
    pvec(1:numel(R.stats.p_task_vs_run)) = R.stats.p_task_vs_run(:);
end
if isfield(R,'stats') && isfield(R.stats,'p_task_vs_run_pooled')
    pvec(end) = R.stats.p_task_vs_run_pooled;
end

fBars = figure('Color','w','Position',[100 100 860 500]); %#ok<NASGU>
hold on;
b = bar([task run],'grouped'); % two bars per group
b(1).FaceColor = [0.20 0.45 0.85]; % Task
b(2).FaceColor = [0.90 0.40 0.20]; % Run

ax = gca; ax.Box='off'; ax.YGrid='on'; ax.FontSize=12;
ax.XTick = 1:numel(xlab); ax.XTickLabel = xlab; ax.XTickLabelRotation = 25;
ylabel('population r(speed, activity)');
title('Population speed–activity correlation');
legend({'Task','Run'},'Location','northwest');

if isprop(b(1),'XEndPoints')
    x1 = b(1).XEndPoints;  y1 = b(1).YEndPoints;
    x2 = b(2).XEndPoints;  y2 = b(2).YEndPoints;
else
    gx = 1:numel(xlab); off = 0.15; %#ok<NASGU>
    error('This MATLAB version needs XEndPoints for annotation.');
end

yl = ylim; pad = 0.04 * range(yl);
i = numel(x1);                         % last group = pooled
yh = max([y1(i) y2(i)]) + pad;
plot([x1(i) x1(i) x2(i) x2(i)], [yh yh+pad/2 yh+pad/2 yh], 'k','LineWidth',1.2);

% pooled p from paired t-test on day pairs (already computed)
pp = R.stats.p_task_vs_run_pooled;
if     pp < 1e-3, s='***';
elseif pp < 1e-2, s='**';
elseif pp < 5e-2, s='*';
else,  s='n.s.';
end
text(mean([x1(i) x2(i)]), yh+pad*0.65, s, 'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', 'FontWeight','bold', 'FontSize',12);

end

function paired_connected_plot(R, varargin)
p = inputParser;
addParameter(p,'X',[]);                 % non-trial
addParameter(p,'Y',[]);                 % trial
addParameter(p,'RatNames',[]);
addParameter(p,'Title','');
addParameter(p,'YLabel','');
addParameter(p,'DoLogY',false,@islogical);
addParameter(p,'ShowMeanLine',true,@islogical);
addParameter(p,'DoPairedT',true,@islogical);
addParameter(p,'XTickLabels',{'Non-trial','Trial'});
parse(p,varargin{:});
X = p.Results.X; Y = p.Results.Y;

if isempty(X) || isempty(Y)
    if nargin<1 || isempty(R)
        error('Provide R or explicit X/Y.');
    end
    X = R.r_run_combined(:);            % Non-trial (Run)
    Y = R.r_task_combined(:);           % Trial (Task)
end

if isempty(p.Results.RatNames)
    if ~isempty(R) && isfield(R,'meta') && isfield(R.meta,'ratNames')
        names = R.meta.ratNames(:);
    else
        names = arrayfun(@(k) sprintf('rat%d',k), 1:numel(X), 'uni',0).';
    end
else
    names = p.Results.RatNames(:);
end

ok = isfinite(X) & isfinite(Y);
X = X(ok); Y = Y(ok); names = names(ok);
nR = numel(X);
if nR==0, warning('No finite pairs to plot.'); return; end

figure('Color','w','Position',[260 220 520 520]); hold on;
ax = gca; ax.Box='off'; ax.FontSize=12;
cNon = [0.50 0.70 0.95];
cTri = [0.95 0.60 0.50];

for i = 1:nR
    plot([1 2], [X(i) Y(i)], '-', 'Color',[0.5 0.5 0.5], 'LineWidth',2.2);
    plot(1, X(i), 'o', 'MarkerFaceColor',cNon, 'MarkerEdgeColor','k', 'MarkerSize',6);
    plot(2, Y(i), 'o', 'MarkerFaceColor',cTri, 'MarkerEdgeColor','k', 'MarkerSize',6);
    text(1-0.06, X(i), string(names{i}), 'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', 'Color',[0.1 0.1 0.1], 'FontAngle','italic');
end

if p.Results.ShowMeanLine
    muX = mean(X,'omitnan');
    muY = mean(Y,'omitnan');
    plot([1 2], [muX muY], 'r--', 'LineWidth',3);
end

statStr = '';
if p.Results.DoPairedT && nR>=2
    try
        [~,pt,~,st] = ttest(Y, X);
        statStr = sprintf('paired t: p=%.3g, t(%d)=%.2f', pt, st.df, st.tstat);
    catch
        pSR = signrank(Y, X);
        statStr = sprintf('signrank: p=%.3g', pSR);
    end
    yl = ylim;
    text(1.5, yl(2) - 0.06*(yl(2)-yl(1)), statStr, 'HorizontalAlignment','center', ...
        'FontSize',11, 'FontWeight','bold');
end

xlim([0.7 2.3]);
set(gca,'XTick',[1 2],'XTickLabel',p.Results.XTickLabels);
if p.Results.DoLogY, set(gca,'YScale','log'); end
ylabel(p.Results.YLabel);
title(p.Results.Title);
grid on; box off;
end


% ------------------ SPEED-BIN MATCH CORE ------------------
function [taskMask, runMaskSel, qualBins] = selectBySpeedBins_secondsMatched( ...
        trialSpeed, trialEffDur, runSpeed, runEffDur, binWidth, minDurPerBin, minBins)
% Returns logical masks of which TASK trials and RUN bins are included after:
% 1) Bin speeds with fixed width binWidth (cm/s),
% 2) Require >= minDurPerBin seconds on BOTH sides (task and run) in that bin,
% 3) For each qualifying bin, select RUN bins deterministically to match TASK seconds
%    (closest to bin center first; tie-break by smaller effDur),
% 4) qualBins = number of bins that passed step (2). The caller decides whether to skip day.

% Define bin edges from combined speeds
smin = min([trialSpeed(:); runSpeed(:)], [], 'omitnan');
smax = max([trialSpeed(:); runSpeed(:)], [], 'omitnan');
if ~isfinite(smin) || ~isfinite(smax)
    taskMask = false(size(trialSpeed));
    runMaskSel = false(size(runSpeed));
    qualBins = 0;
    return;
end
% align edges to binWidth
a = floor(smin/binWidth)*binWidth;
b = ceil(smax/binWidth)*binWidth;
edges = a:binWidth:b;
if numel(edges) < 2
    edges = [a a+binWidth];
end
centers = edges(1:end-1) + binWidth/2;

% Assign bins
[~,~,tBinIdx] = histcounts(trialSpeed, edges);
[~,~,rBinIdx] = histcounts(runSpeed,   edges);

taskMask   = false(size(trialSpeed));
runMaskSel = false(size(runSpeed));

qualBins = 0;
for bi = 1:numel(centers)
    t_in = (tBinIdx == bi) & isfinite(trialEffDur);
    r_in = (rBinIdx == bi);

    taskSec = sum(trialEffDur(t_in), 'omitnan');
    runSec  = sum(runEffDur(r_in),   'omitnan');

    if taskSec >= minDurPerBin && runSec >= minDurPerBin
        qualBins = qualBins + 1;
        % keep all task trials in this bin
        taskMask(t_in) = true;

        % select run bins to match taskSec (deterministic greedy)
        cand = find(r_in);
        if ~isempty(cand)
            % rank by |speed - center|, then by smaller effDur
            [~,o1] = sort(abs(runSpeed(cand) - centers(bi)), 'ascend');
            [~,o2] = sort(runEffDur(cand(o1)), 'ascend');
            ord    = o1(o2);

            acc = 0;
            for j = 1:numel(ord)
                ii = cand(ord(j));
                runMaskSel(ii) = true;
                acc = acc + runEffDur(ii);
                if acc >= taskSec, break; end
            end
        end
    end
end

% NOTE: No hard enforcement of MinBins here.
% The caller chooses to skip the day (policy=='day') or proceed (policy=='lenient').
end

% --------- small plotting helper ---------
function ymax = plotGroup_local(i, rt, rr, off, ms)
% returns max y among points to help autoscale
x1 = i - off; x2 = i + off;
scatter(repmat(x1,numel(rt),1), rt, ms, [0.20 0.45 0.85],'filled','MarkerFaceAlpha',0.5);
scatter(repmat(x2,numel(rr),1), rr, ms, [0.90 0.40 0.20],'filled','MarkerFaceAlpha',0.5);
ymax = max([rt(:); rr(:)], [], 'omitnan');
end





function [lme, anovTbl, betaTbl, summaryStr] = lme_speed_vs_activity_task_vs_run(R)
% Build a long table of the ACTUAL points used for r (after your matching)
% and test if the speed->activity slope differs between Task and Run.

rats = R.meta.ratNames;
y = []; x = []; cond = {}; rat = {}; day = {};

for r = 1:numel(rats)
    % ----- Task -----
    if r<=numel(R.raw.task) && ~isempty(R.raw.task{r})
        for d = 1:numel(R.raw.task{r})
            S = R.raw.task{r}{d}; if isempty(S), continue; end
            ok = isfinite(S.speed) & isfinite(S.pop);
            n  = sum(ok); if n<3, continue; end
            x   = [x;   S.speed(ok)];
            y   = [y;   S.pop(ok)];
            cond= [cond; repmat({'task'}, n,1)];
            rat = [rat; repmat(rats(r), n,1)];
            day = [day; repmat(string(d), n,1)];
        end
    end
    % ----- Run -----
    if r<=numel(R.raw.run) && ~isempty(R.raw.run{r})
        for d = 1:numel(R.raw.run{r})
            S = R.raw.run{r}{d}; if isempty(S), continue; end
            ok = isfinite(S.speed) & isfinite(S.pop);
            n  = sum(ok); if n<3, continue; end
            x   = [x;   S.speed(ok)];
            y   = [y;   S.pop(ok)];
            cond= [cond; repmat({'run'}, n,1)];
            rat = [rat; repmat(rats(r), n,1)];
            day = [day; repmat(string(d), n,1)];
        end
    end
end

tbl = table(y, x, categorical(cond), categorical(rat), categorical(day), ...
            'VariableNames', {'y','x','cond','rat','day'});

% Fit LME (random intercepts for rat and day nested in rat).
% If you want, try adding a random slope for speed by rat: (1 + x|rat)
lme = fitlme(tbl, 'y ~ x*cond + (1|rat) + (1|rat:day)');

anovTbl = anova(lme);              % for overall x:cond term
betaTbl = lme.Coefficients;        % estimates and p-values

% Extract the interaction row (x:cond_run if 'task' is reference)
ix = find(strcmp(betaTbl.Name, 'x:cond_run'));
if isempty(ix)                        % if 'run' became reference
    ix = find(contains(betaTbl.Name, 'x:cond_'));
end
b  = betaTbl.Estimate(ix);
p  = betaTbl.pValue(ix);

summaryStr = sprintf('Mixed-effects slope difference (Run vs Task): \\beta=%.4g, p=%.3g, N=%d pts, R=%d rats', ...
                     b, p, height(tbl), numel(unique(tbl.rat)));
end
