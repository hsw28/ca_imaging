function R = speedDependence_taskVsRun(ratNames, varargin)
% speedDependence_taskVsRun
% Compare per-cell speed dependence in TRACE (task) vs NON-TASK running.
%
% Usage:
%   R = speedDependence_taskVsRun({'rat0222','rat0307','rat0313','rat0314','rat0816'});
%
% Inputs:
%   ratNames : cellstr of base-workspace variable names for rats
%
% Name-Value options:
%   'Win'         [0 2]     % trace window relative to CS (s)
%   'BinSec'      2         % non-task running bin length (s)
%   'SpeedThresh' 4         % cm/s cutoff for running bins and (optional) task gating
%   'QFDR'        0.05      % Benjamini-Hochberg FDR q
%   'NShuffle'    1000      % (placeholder) permutations for label-shuffle nulls
%   'DoPlots'     false     % quick histograms
%
%   % NEW: cap non-task analyzed time to match task analyzed time
%   'CapNonTaskToTask'     false      % if true, cap non-task total seconds ~= task total
%   'CapMode'              'random'   % 'random' | 'velmatch'
%   'CapSeed'              []         % RNG seed for reproducibility ([] -> do not set)
%   'VelMatchNBins'        10         % # speed bins for 'velmatch' mode
%
% Data requirements per day:
%   rat.pos.pos_DATE           -> [t x y] with t in seconds
%   rat.vel.vel_DATE (optional)-> velocity time series [t v]; else derived from pos
%   rat.CS_times.CS_DATE       -> CS onset times (s)
%   rat.Ca_peaks.CA_peaks_DATE -> 1xN cell, spike times per neuron (s)
%
% Returns R with per-rat and pooled fields (same as before) plus meta flags.

% ---------- args ----------
p = inputParser;
addParameter(p,'Win',[0 2]);
addParameter(p,'BinSec',2);
addParameter(p,'SpeedThresh',4);
addParameter(p,'QFDR',0.05);
addParameter(p,'NShuffle',1000);
addParameter(p,'DoPlots',true);
addParameter(p,'CapNonTaskToTask',false);
addParameter(p,'CapMode','random');
addParameter(p,'CapUseSpeedGatedTask',false);
addParameter(p,'CapSeed',[]);
addParameter(p,'VelMatchNBins',10);
parse(p,varargin{:});
Win         = p.Results.Win;
BinSec      = p.Results.BinSec;
vThresh     = p.Results.SpeedThresh;
qFDR        = p.Results.QFDR;
Nsh         = p.Results.NShuffle; %#ok<NASGU>  % placeholder
doplot      = p.Results.DoPlots;
capOn       = p.Results.CapNonTaskToTask;
capMode     = lower(p.Results.CapMode);
capSeed     = p.Results.CapSeed;
capTaskGate = p.Results.CapUseSpeedGatedTask;
nBinsVM     = p.Results.VelMatchNBins;

nR = numel(ratNames);
R = struct('r_task',{cell(1,nR)},'r_run',{cell(1,nR)}, ...
           'p_task',{cell(1,nR)},'p_run',{cell(1,nR)}, ...
           'sig_task_FDR',{cell(1,nR)},'sig_run_FDR',{cell(1,nR)}, ...
           'delta_r',{cell(1,nR)},'p_delta',{cell(1,nR)}, ...
           'sig_delta_FDR',{cell(1,nR)}, ...
           'rateRatio_match',{cell(1,nR)}, ...
           'fracSig_task',[],'fracSig_run',[],'fracSig_delta',[],'meta',[]);

all_task = []; all_run = []; all_pt = []; all_pr = []; all_d = []; all_pd = [];
all_sig_t = []; all_sig_r = []; all_sig_d = [];
all_rat_idx = [];

for r = 1:nR
    rat = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    iAn   = find(strcmp(dates,rat.An),1);
    if isempty(iAn) || iAn < 3
        warning('%s: not enough days before An; skipping.', ratNames{r});
        continue;
    end
    days = dates(iAn-2:iAn);  % 3 days

    % ---------- collect per-day trials & bins ----------
    r_task = []; r_run = []; p_task = []; p_run = [];
    rateRatio_match = [];
    allTrialSpeed = [];  % for vel-match cap
    taskEffSec_allDays = 0;

    % Hold non-task bins across days; we'll cap AFTER pooling
    binSpeed_all = cell(numel(days),1);
    binRate_all  = cell(numel(days),1);
    effDur_all   = cell(numel(days),1);

    % ---------- collect per-day trials & bins ----------
    r_task = []; r_run = []; p_task = []; p_run = [];
    rateRatio_match = [];
    allTrialSpeed = [];  % only used for pooled summaries if desired (not needed for capping now)

    for d = 1:numel(days)
        D = days{d};
        Sraw = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));   % n×m cell (n=cells)
        S    = normalizeCApeaks(Sraw);                    % {n,1} spike vectors

        pos = rat.pos.(sprintf('pos_%s',D));              % [t x y]
        [vt, vv] = velocityFromPos(pos);

        CS  = rat.CS_times.(sprintf('CS_%s',D));          % CS onset times

        % ---------- trial side ----------
        if capTaskGate
            [trialSpeed_d, trialRatePerCell_d] = trialSpeedAndRates_speedGated( ...
                S, CS, vt, vv, Win, vThresh, pos(:,1));
        else
            [trialSpeed_d, trialRatePerCell_d] = trialSpeedAndRates_basic( ...
                S, CS, vt, vv, Win);
        end
        allTrialSpeed = [allTrialSpeed; trialSpeed_d]; %#ok<AGROW>

        % ---------- non-task side (uncapped bins first) ----------
        [binSpeed_d, binRatePerCell_d, effDur_d] = runBinsSpeedAndRates_withDur( ...
            S, CS, pos(:,1), vt, vv, Win, BinSec, vThresh);

        % ---------- compute this day's task effective seconds ----------
        dtPos  = median(diff(pos(:,1)),'omitnan');
        vOnPos = interp1(vt, vv, pos(:,1), 'linear','extrap');
        if capTaskGate
            inWin = false(size(pos,1),1);
            for k = 1:numel(CS)
                t0 = CS(k)+Win(1); t1 = CS(k)+Win(2);
                inWin = inWin | (pos(:,1)>=t0 & pos(:,1)<t1);
            end
            taskEffSec_d = sum(inWin & (vOnPos >= vThresh)) * dtPos;
        else
            taskEffSec_d = numel(CS) * diff(Win);
        end

        % ---------- optionally cap THIS DAY'S non-task seconds ----------
        if capOn & ~isempty(binSpeed_d)
            if ~isempty(capSeed), rng(capSeed); end
            switch capMode
                case 'random'
                    sel = selectBinsToMatchSeconds_random(effDur_d, taskEffSec_d);
                case 'velmatch'
                    sel = selectBinsToMatchSeconds_velmatch(binSpeed_d, effDur_d, ...
                            trialSpeed_d, taskEffSec_d, nBinsVM);
                otherwise
                    error('CapMode must be ''random'' or ''velmatch''.');
            end
            binSpeed_d       = binSpeed_d(sel);
            binRatePerCell_d = binRatePerCell_d(sel,:);
            effDur_d         = effDur_d(sel); %#ok<NASGU> % kept for clarity
        end

        % ---------- correlations for this day ----------
        [rt, pt] = perCellCorr(trialSpeed_d, trialRatePerCell_d);  % task
        [rr, pr] = perCellCorr(binSpeed_d , binRatePerCell_d );    % run

        % ---------- speed-matched trial vs run ratios (use the same selected bins) ----------
        rr_match = perCellSpeedMatchedRateRatio(trialSpeed_d, trialRatePerCell_d, ...
                                                binSpeed_d, binRatePerCell_d);

        % ---------- append across days (per rat) ----------
        r_task = [r_task; rt];
        r_run  = [r_run ; rr];
        p_task = [p_task; pt];
        p_run  = [p_run ; pr];
        rateRatio_match = [rateRatio_match; rr_match];
    end




    % FDR for task & run
    sig_task = fdr_mask(p_task, qFDR);
    sig_run  = fdr_mask(p_run , qFDR);

    % delta r and approximate per-cell p via Fisher z (conservative)
    [p_delta, sig_delta] = deltaPermuteP(r_task, r_run);

    % store per rat
    R.r_task{r} = r_task; R.r_run{r} = r_run;
    R.p_task{r} = p_task; R.p_run{r} = p_run;
    R.sig_task_FDR{r} = sig_task; R.sig_run_FDR{r} = sig_run;
    R.delta_r{r} = r_task - r_run; R.p_delta{r} = p_delta;
    R.sig_delta_FDR{r} = sig_delta;
    R.rateRatio_match{r} = rateRatio_match;

    R.fracSig_task(r,1)  = mean(sig_task,'omitnan');
    R.fracSig_run(r,1)   = mean(sig_run ,'omitnan');
    R.fracSig_delta(r,1) = mean(sig_delta,'omitnan');

    % pooled collectors
    all_task   = [all_task; r_task];
    all_run    = [all_run ; r_run ];
    all_pt     = [all_pt  ; p_task];
    all_pr     = [all_pr  ; p_run ];
    all_d      = [all_d   ; r_task - r_run];
    all_pd     = [all_pd  ; p_delta];
    all_sig_t  = [all_sig_t; sig_task];
    all_sig_r  = [all_sig_r; sig_run];
    all_sig_d  = [all_sig_d; sig_delta];
    all_rat_idx = [all_rat_idx; r*ones(numel(r_task),1)];
end

% pooled summaries
R.pooled.r_task = all_task;
R.pooled.r_run  = all_run;
R.pooled.delta  = all_d;
R.pooled.fracSig_task  = mean(all_sig_t,'omitnan');
R.pooled.fracSig_run   = mean(all_sig_r,'omitnan');
R.pooled.fracSig_delta = mean(all_sig_d,'omitnan');

R.meta = struct('Win',Win,'BinSec',BinSec,'SpeedThresh',vThresh,'QFDR',qFDR, ...
                'NShuffle',Nsh,'ratNames',{ratNames}, ...
                'CapNonTaskToTask',capOn,'CapMode',capMode, ...
                'CapUseSpeedGatedTask',capTaskGate,'VelMatchNBins',nBinsVM);

% Optional: McNemar test per rat comparing sig fractions (run vs task)
for r = 1:numel(R.r_task)
    if isempty(R.r_task{r}), continue; end
    sigTask = R.sig_task_FDR{r};
    sigRun  = R.sig_run_FDR{r};
    n12 = sum(~sigTask & sigRun);
    n21 = sum(sigTask & ~sigRun);
    nTot = n12 + n21;
    if nTot > 0
        pval = 2*binocdf(min(n12,n21), nTot, 0.5);
        fprintf('%s: run vs task sig fraction, McNemar p = %.4g\n', R.meta.ratNames{r}, pval);
    end
end

% inside speedDependence_taskVsRun, replace your current "if doplot" block with:
if doplot
    plotFracSigBars(R, true);   % set second arg to false if you want ONLY Task/Run
end


end

% ---------- helpers ----------
function [vt, vv] = velocityFromPos(pos)
% pos: [t x y], t in seconds, x,y in cm
t = pos(:,1); xy = pos(:,2:3);
dt = diff(t); dt(dt==0)=NaN;
v  = [0; sqrt(sum(diff(xy).^2,2))./dt]; % cm/s assuming cm units
v(~isfinite(v)) = 0;
vt = t; vv = v;
end

function Scell = normalizeCApeaks(Sraw)
% normalizeCApeaks  Convert CA_peaks to {n,1} cell of spike-time vectors.
if iscell(Sraw)
    if isvector(Sraw), Scell = Sraw(:); return; end
    K = size(Sraw,2);
    scores = zeros(1,K);
    for j = 1:K
        col = Sraw(:,j);
        scores(j) = mean(cellfun(@(x) isnumeric(x) && isvector(x), col));
    end
    [bestScore, jbest] = max(scores);
    if bestScore == 0
        error('normalizeCApeaks:NoNumericColumn', ...
              'No column in CA_peaks looks like numeric spike-time vectors.');
    end
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
error('normalizeCApeaks:BadType', ...
      'Unsupported CA_peaks type: %s. Expected cell or numeric matrix.', class(Sraw));
end

function [trialSpeed, rateMat] = trialSpeedAndRates_basic(S, CS, vt, vv, Win)
% Per-trial mean speed in Win and event rate (spikes/s) using full window length.
nCells = size(S,1);
nTr    = numel(CS);
trialSpeed = nan(nTr,1);
rateMat    = nan(nTr,nCells);
wlen = diff(Win);
for k = 1:nTr
    t0 = CS(k)+Win(1); t1 = CS(k)+Win(2);
    ivt = vt>=t0 & vt<t1;
    trialSpeed(k) = mean(vv(ivt),'omitnan');
    for c = 1:nCells
        sp = S{c,1};
        if isempty(sp), rateMat(k,c)=0; continue; end
        rateMat(k,c) = sum(sp>=t0 & sp<t1) / wlen;
    end
end
end

function [binSpeed, rateMat, effDur] = runBinsSpeedAndRates_withDur(S, CS, tPos, vt, vv, Win, binSec, vThr)
% Build non-task bins outside [CS+Win] with v >= vThr; compute mean speed, rates, and effective duration.
nCells = size(S,1);
% mask out trial windows
csMask = false(size(tPos));
for k = 1:numel(CS)
    csMask = csMask | ((tPos >= CS(k)+Win(1)) & (tPos < CS(k)+Win(2)));
end
vOnPos = interp1(vt, vv, tPos, 'linear','extrap');
runMask = ~csMask & (vOnPos >= vThr);

tmin = tPos(find(runMask,1,'first'));
tmax = tPos(find(runMask,1,'last'));
if isempty(tmin) || isempty(tmax)
    binSpeed = []; rateMat = []; effDur = [];
    return
end

edges = tmin:binSec:tmax;
if edges(end) < tmax, edges = [edges tmax]; end
nB = numel(edges)-1;

binSpeed = nan(nB,1);
effDur   = zeros(nB,1);
rateMat  = nan(nB,nCells);

dtPos = median(diff(tPos), 'omitnan');

for b = 1:nB
    iv = (tPos>=edges(b)) & (tPos<edges(b+1)) & runMask;
    effDur(b) = sum(iv) * dtPos;     % seconds with v>=thr
    if effDur(b) <= 0, continue; end
    binSpeed(b) = mean(vOnPos(iv), 'omitnan');

    % spike counts in the bin; divide by effective duration
    for c = 1:nCells
        sp = S{c,1};
        if isempty(sp), rateMat(b,c) = 0; continue; end
        nSp = sum(sp>=edges(b) & sp<edges(b+1));
        rateMat(b,c) = nSp / effDur(b);
    end
end

ok = isfinite(binSpeed) & (effDur>0) & any(isfinite(rateMat),2);
binSpeed = binSpeed(ok);
rateMat  = rateMat(ok,:);
effDur   = effDur(ok);
end

function [r, p] = perCellCorr(x, Y)
% x: n x 1 speeds; Y: n x cells rates
nC = size(Y,2); r = nan(nC,1); p = nan(nC,1);
for c = 1:nC
    xc = x; yc = Y(:,c);
    ok = isfinite(xc) & isfinite(yc);
    if sum(ok) >= 5
        [rc, pc] = corr(xc(ok), yc(ok), 'type','Pearson');
        r(c)=rc; p(c)=pc;
    end
end
end

function thr = fdr_bh(p, q)
% Return a single BH threshold in [0,1] for vector p (ignores NaNs).
p = p(:); p = p(isfinite(p)); n = numel(p);
if n == 0, thr = 0; return; end
p = sort(p);
line = ((1:n).' / n) * q;
k = find(p <= line, 1, 'last');
if isempty(k), thr = 0; else, thr = p(k); end
end

function mask = fdr_mask(p, q)
thr = fdr_bh(p, q);
mask = false(size(p));
ii = isfinite(p);
mask(ii) = p(ii) < thr;
end

function [p_delta, sigDelta] = deltaPermuteP(r_task, r_run)
% Approximate per-cell p for delta using Fisher z difference (two-sided normal).
z = atanh(r_task) - atanh(r_run);
p_delta = 2*normcdf(-abs(z),0,1);
sigDelta = fdr_mask(p_delta, 0.05);
end

function rr_match = perCellSpeedMatchedRateRatio(trialSpeed, trialRate, binSpeed, binRate)
% For each trial, pick nearest-speed run bin; compute (trial rate)/(matched bin rate).
% Return per-cell geometric mean of ratios across trials.
nC = size(trialRate,2);
rr_match = nan(nC,1);
if isempty(binSpeed), return; end
for c = 1:nC
    ratios = nan(numel(trialSpeed),1);
    for k = 1:numel(trialSpeed)
        s = trialSpeed(k);
        if ~isfinite(s), continue; end
        [~, idx] = min(abs(binSpeed - s));
        denom = binRate(idx,c);
        nume  = trialRate(k,c);
        if isfinite(denom) && denom>0 && isfinite(nume)
            ratios(k) = nume/denom;
        end
    end
    rr_match(c) = exp(mean(log(ratios(ratios>0 & isfinite(ratios))), 'omitnan')); % geo-mean
end
end

function sel = selectBinsToMatchSeconds_random(effDur, targetSec)
% Randomly sample bins (without replacement) until accumulated effDur ~ targetSec.
idx = find(effDur > 0);
if isempty(idx), sel = idx; return; end
ord = randperm(numel(idx));
acc = 0; keep = false(size(effDur));
for j = 1:numel(ord)
    b = idx(ord(j));
    if acc + effDur(b) > targetSec
        % choose the closer of including or stopping
        if (targetSec - acc) < (acc + effDur(b) - targetSec)
            % stop before adding b
        else
            keep(b) = true; acc = acc + effDur(b);
        end
        break;
    else
        keep(b) = true; acc = acc + effDur(b);
    end
end
% Fallback: ensure at least one bin
if ~any(keep), keep(idx(ord(1))) = true; end
sel = find(keep);
end

function sel = selectBinsToMatchSeconds_velmatch(binSpeed, effDur, trialSpeed, targetSec, nBins)
% Match non-task seconds to the distribution of trial speeds.
% 1) Build speed quantile bins from trialSpeed (finite only).
% 2) Allocate target seconds across quantile bins proportional to trial counts.
% 3) Within each bin, randomly sample non-task bins until its per-bin target is reached.

% Clean inputs
okTr = isfinite(trialSpeed);
trialSpeed = trialSpeed(okTr);
if isempty(trialSpeed) || all(~isfinite(binSpeed))
    % fallback to random if we cannot define bins
    sel = selectBinsToMatchSeconds_random(effDur, targetSec);
    return;
end

% Quantile edges from trial speeds
q = linspace(0,1,nBins+1);
edges = quantile(trialSpeed, q);
% Ensure monotonic non-decreasing edges (quantile can repeat if degenerate)
edges(1)   = -inf;
edges(end) = inf;

% Assign trials to bins & compute per-bin targets
trialCounts = histcounts(trialSpeed, edges);
if sum(trialCounts)==0
    sel = selectBinsToMatchSeconds_random(effDur, targetSec);
    return;
end
targetPerBin = targetSec * (trialCounts / sum(trialCounts));

% Assign non-task bins to the same edges based on binSpeed
[~,~,binIdx] = histcounts(binSpeed, edges);

keep = false(size(effDur));
for b = 1:nBins
    targ = targetPerBin(b);
    if targ <= 0, continue; end
    cand = find(binIdx==b & effDur>0);
    if isempty(cand), continue; end
    cand = cand(randperm(numel(cand))); % randomize within bin
    acc = 0;
    for j = 1:numel(cand)
        ii = cand(j);
        if acc + effDur(ii) > targ
            % choose closer of including or not
            if (targ - acc) < (acc + effDur(ii) - targ)
                % stop
            else
                keep(ii) = true; acc = acc + effDur(ii);
            end
            break;
        else
            keep(ii) = true; acc = acc + effDur(ii);
        end
    end
    % If still short (no single bin perfectly fills), that's okay; leftovers will remain.
end

% If we're severely short (e.g., scarce non-task data in some speed ranges), top up randomly.
accTot = sum(effDur(keep));
if accTot < 0.9*targetSec
    topup = selectBinsToMatchSeconds_random(effDur(~keep), targetSec - accTot);
    iiAll = find(~keep);
    keep(iiAll(topup)) = true;
end

% If we overshot due to rounding, trim a little from the largest-effDur bins
accTot = sum(effDur(keep));
if accTot > 1.1*targetSec
    ii = find(keep);
    [~,order] = sort(effDur(ii),'descend');
    for j = 1:numel(order)
        if accTot <= targetSec, break; end
        keep(ii(order(j))) = false;
        accTot = sum(effDur(keep));
    end
end

% Fallback: ensure at least one bin
if ~any(keep)
    cand = find(effDur>0);
    if ~isempty(cand), keep(cand(1)) = true; end
end
sel = find(keep);

end

function plotFracSigBars(R, includeDelta)
% plotFracSigBars  Grouped bars of fraction FDR-significant per rat (+ pooled)
% includeDelta: if true, add a 3rd bar for Δ (r_task - r_run)

    if nargin < 2, includeDelta = true; end

    % Assemble per-rat + pooled vectors
    frac_task = [R.fracSig_task(:);  R.pooled.fracSig_task];
    frac_run  = [R.fracSig_run(:);   R.pooled.fracSig_run];

    if includeDelta
        frac_delta = [R.fracSig_delta(:); R.pooled.fracSig_delta];
        Y = [frac_task frac_run frac_delta];
        leg = {'Task','Run','\Delta'};
    else
        Y = [frac_task frac_run];
        leg = {'Task','Run'};
    end

    % X-labels = rats + pooled
    xlab = [R.meta.ratNames, {'pooled'}];

    % Plot
    figure('Color','w','Position',[100 100 760 420]);
    b = bar(Y,'grouped'); hold on;

    % Colors (nice, distinct)
    if size(Y,2)>=1, b(1).FaceColor = [0.20 0.45 0.85]; end   % Task (blue)
    if size(Y,2)>=2, b(2).FaceColor = [0.90 0.40 0.20]; end   % Run  (orange-red)
    if size(Y,2)>=3, b(3).FaceColor = [0.95 0.85 0.20]; end   % Δ    (gold)

    ax = gca;
    ax.YLim = [0 1];
    ax.YGrid = 'on';
    ax.Box = 'off';
    ax.FontSize = 12;
    ax.XTick = 1:numel(xlab);
    ax.XTickLabel = xlab;
    ax.XTickLabelRotation = 25;

    ylabel('fraction significant (FDR q=0.05)');
    title('FDR-significant proportions');
    legend(leg,'Location','northwest');

    hold off;
end
