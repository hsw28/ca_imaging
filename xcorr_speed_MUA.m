function OUT = xcorr_speed_MUA(ratNames, varargin)
% xcorr_speed_MUA
% Cross-correlate running speed and population multi-unit rate (MUA)
% within each rat, then compare to swapped across-rat pairs.
%
% Uses the same data flow as speedDependence_taskVsRun_population:
%   CA_peaks -> normalizeCApeaks -> smoothpos -> velocityFromPos ->
%   taskBinsSpeedAndRates_withDur / runBinsSpeedAndRates_withDur ->
%   normalizePerCell -> aggregateAcrossCells
%
% INPUT:
%   ratNames : cell array of rat variable names in base workspace
%
% Name-Value options:
%   'Win'         [0 2]      % trial window for tasks; also excluded from Run
%   'binSize'     1/7.5      % bin length (s)
%   'SpeedThresh' 4          % v>=SpeedThresh for Run bins (ignored for task)
%   'PopAgg'      'mean'     % 'mean' or 'sum' across cells
%   'PopNorm'     'none'     % 'none' | 'mean' | 'demean' per cell
%   'MaxLagSec'   1          % xcorr lags ±MaxLagSec
%   'MinSamples'  50         % minimum samples per rat to keep
%   'Mode'        'run'      % 'run' (non-task) or 'task'
%   'TrialDemean' true       % if true and Mode='task', demean per trial
%   'DoPlot'      true       % plot curves
%
%   Grouping:
%   'GroupMode'   'none'     % 'none' | 'on'
%   'GroupVar'    ''         % field name in rat struct, e.g. 'rates' or 'MI_noCSUS30'
%                             (must be a struct with per-day fields)
%   'GroupNum'    1          % # groups (e.g. 3 for tertiles, 4 for labels 1–4)
%   'GroupType'   'quantile' % 'quantile' (continuous) | 'discrete' (integer labels)
%
% OUTPUT struct OUT:
%   .lagsSec             : [L x 1] lag axis (s)
%   .withinCurves        : [L x nR x nG] xcorr per rat and group
%   .withinPeakR         : [nR x nG] peak corr (by |r|)
%   .withinPeakLagSec    : [nR x nG] lag-of-peak (s)
%   .withinR0            : [nR x nG] r at lag 0
%
%   .acrossCurves        : [L x nPairs x nG] xcorr for swapped pairs
%   .acrossPairs         : [nPairs x 2] [i j] indices (speed_i, pop_j)
%   .acrossPeakR         : [nPairs x nG] peak corr (by |r|)
%   .acrossPeakLagSec    : [nPairs x nG] lag-of-peak (s)
%   .acrossR0            : [nPairs x nG] r at lag 0
%
%   .meta                : struct of settings and ratNames

% ---------- options ----------
p = inputParser;
p.addParameter('Win',[0 2]);
p.addParameter('binSize',1/7.5);
p.addParameter('SpeedThresh',0);
p.addParameter('PopAgg','mean',@(s)any(strcmpi(s,{'mean','sum'})));
p.addParameter('PopNorm','none',@(s)any(strcmpi(s,{'none','mean','demean'})));
p.addParameter('MaxLagSec',1);
p.addParameter('MinSamples',50);
p.addParameter('Mode','task',@(s)any(strcmpi(s,{'run','task'})));
p.addParameter('TrialDemean',false,@islogical); %%dont think i ever want this
p.addParameter('DoPlot',true,@islogical);

% NEW grouping params
p.addParameter('GroupMode','none',@(s)any(strcmpi(s,{'none','on'})));
p.addParameter('GroupVar','',@(s)ischar(s) || isstring(s));
p.addParameter('GroupNum',1,@(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('GroupType','quantile',@(s)any(strcmpi(s,{'quantile','discrete'})));

p.parse(varargin{:});

Win         = p.Results.Win;
binSize     = p.Results.binSize;
vThresh     = p.Results.SpeedThresh;
popAgg      = lower(p.Results.PopAgg);
popNorm     = lower(p.Results.PopNorm);
MaxLagSec   = p.Results.MaxLagSec;
MinSamples  = p.Results.MinSamples;
Mode        = lower(p.Results.Mode);
TrialDemean = p.Results.TrialDemean;
DoPlot      = p.Results.DoPlot;

GroupMode   = lower(p.Results.GroupMode);
GroupVar    = char(p.Results.GroupVar);
GroupNum    = round(p.Results.GroupNum);
GroupType   = lower(p.Results.GroupType);

doGroup = strcmp(GroupMode,'on');
if doGroup
    if isempty(GroupVar)
        error('GroupMode ''on'' requires GroupVar (e.g. ''rates'' or ''MI_noCSUS30'').');
    end
    if GroupNum < 1
        error('GroupNum must be >= 1.');
    end
else
    GroupNum = 1;   % single group
end

if TrialDemean && ~strcmp(Mode,'task')
    warning('TrialDemean is only meaningful for Mode=''task''; ignoring.');
end

nR = numel(ratNames);

% ---------- collect per-rat speed & pop traces ----------
if doGroup
    speed_all = cell(nR, GroupNum);   % per rat x group
    pop_all   = cell(nR, GroupNum);
else
    speed_all = cell(nR, 1);
    pop_all   = cell(nR, 1);
end

for r = 1:nR
    rat = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    iAn   = find(strcmp(dates, rat.An), 1);
    if isempty(iAn) || iAn < 3
        warning('%s: not enough days before An; skipping.', ratNames{r});
        continue;
    end
    days = dates(iAn-2:iAn);   % last 3 days up to An

    if doGroup
        thisSpeed = cell(GroupNum,1);
        thisPop   = cell(GroupNum,1);
    else
        thisSpeed = [];
        thisPop   = [];
    end

    for d = 1:numel(days)
        D = days{d};

        % CA peaks -> cell array of spike times
        Sraw   = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));
        S      = normalizeCApeaks_local(Sraw);
        nCells = numel(S);

        % smoothed position and velocity
        pos       = smoothpos(rat.pos.(sprintf('pos_%s',D)));
        [vt, vv]  = velocityFromPos_local(pos);
        CS        = rat.CS_times.(sprintf('CS_%s',D));

        switch Mode
            case 'run'
                [binSpeed, ratePerCell, effDur] = runBinsSpeedAndRates_withDur_local( ...
                    S, CS, pos(:,1), vt, vv, Win, binSize, vThresh);
            case 'task'
                [binSpeed, ratePerCell, effDur] = taskBinsSpeedAndRates_withDur_local( ...
                    S, CS, pos(:,1), vt, vv, Win, binSize, []);  % no speed gate
        end

        if isempty(binSpeed)
            continue;
        end

        % per-cell normalization
        ratePerCell_use = normalizePerCell_local(ratePerCell, popNorm);

        if doGroup
            % --------- grouping by per-cell variable (e.g. firing rate / MI / label) ---------
            if ~isfield(rat, GroupVar)
                error('Rat %s is missing field ''%s'' for grouping.', ratNames{r}, GroupVar);
            end

            GVstruct = rat.(GroupVar);   % e.g. rat.MI_noCSUS30

            if ~isstruct(GVstruct)
                error('Rat %s: %s must be a struct containing per-day fields.', ...
                      ratNames{r}, GroupVar);
            end

            % find a field whose name contains the date string D, e.g. 'MI_2023_05_28'
            fns = fieldnames(GVstruct);
            idxFld = contains(fns, D);
            if ~any(idxFld)
                error('Rat %s: no field in %s contains date %s.', ...
                    ratNames{r}, GroupVar, D);
            end
            if sum(idxFld) > 1
                warning('Rat %s: multiple fields in %s match date %s; using first.', ...
                        ratNames{r}, GroupVar, D);
            end

            fname = fns{find(idxFld,1,'first')};
            gvDay = GVstruct.(fname);
            gvDay = gvDay(:);   % force column

            if numel(gvDay) < nCells
                error('Rat %s: %s.%s has only %d entries for %d cells.', ...
                    ratNames{r}, GroupVar, fname, numel(gvDay), nCells);
            end
            gvDay = gvDay(1:nCells);

            groupID = nan(nCells,1);

            switch GroupType
                case 'quantile'
                    % ignore NaNs when computing quantiles
                    validMask = isfinite(gvDay);
                    if ~any(validMask)
                        warning('%s: all NaN in %s.%s; skipping grouping for this day.', ...
                            ratNames{r}, GroupVar, fname);
                        continue;
                    end
                    edges = prctile(gvDay(validMask), linspace(0,100,GroupNum+1));

                    for g = 1:GroupNum
                        if g < GroupNum
                            inG = gvDay >= edges(g) & gvDay <  edges(g+1);
                        else
                            inG = gvDay >= edges(g) & gvDay <= edges(g+1);
                        end
                        groupID(inG) = g;
                    end

                case 'discrete'
                    % labels already encode group numbers (e.g. 1,2,3,4)
                    validMask = isfinite(gvDay);
                    if ~any(validMask)
                        warning('%s: all NaN in %s.%s; skipping grouping for this day.', ...
                            ratNames{r}, GroupVar, fname);
                        continue;
                    end
                    groupID(validMask) = gvDay(validMask);

                    bad = groupID(validMask) < 1 | groupID(validMask) > GroupNum;
                    if any(bad)
                        warning(['%s: %s.%s has labels outside [1..%d]; ' ...
                                 'ignoring those cells for grouping.'], ...
                                ratNames{r}, GroupVar, fname, GroupNum);
                        tmp = groupID;
                        tmp(groupID < 0 | groupID > GroupNum) = NaN;
                        groupID = tmp;
                    end

                otherwise
                    error('Unknown GroupType: %s', GroupType);
            end

            % build group-specific population vectors for this day
            for g = 1:GroupNum
                cols = find(groupID==g);
                if isempty(cols)
                    % no cells in this group for this day
                    continue;
                end
                popVec_g = aggregateAcrossCells_local(ratePerCell_use(:,cols), popAgg);

                binSpeed_g = binSpeed;
                effDur_g   = effDur;

                % optional TRIAL-DEMEAN (task only), per group
                if TrialDemean && strcmp(Mode,'task')
                    [binSpeed_g, popVec_g, effDur_g] = trial_demean_bins( ...
                        binSpeed_g, popVec_g, effDur_g, CS, Win, binSize);
                end

                ok = isfinite(binSpeed_g) & isfinite(popVec_g) & (effDur_g>0);
                thisSpeed{g} = [thisSpeed{g}; binSpeed_g(ok)];
                thisPop{g}   = [thisPop{g};   popVec_g(ok)];
            end

        else
            % --------- aggregate across ALL cells (no grouping) ----------
            popVec = aggregateAcrossCells_local(ratePerCell_use, popAgg);

            % optional TRIAL-DEMEAN (task only)
            if TrialDemean && strcmp(Mode,'task')
                [binSpeed, popVec, effDur] = trial_demean_bins( ...
                    binSpeed, popVec, effDur, CS, Win, binSize);
            end

            ok = isfinite(binSpeed) & isfinite(popVec) & (effDur>0);
            thisSpeed = [thisSpeed; binSpeed(ok)];
            thisPop   = [thisPop;   popVec(ok)];
        end
    end

    % ---------- per-rat sample threshold & store ----------
    if doGroup
        for g = 1:GroupNum
            nS = numel(thisSpeed{g});
            if nS < MinSamples
                if nS > 0
                    fprintf('%s: only %d %s samples in group %d, skipping this group.\n', ...
                        ratNames{r}, nS, Mode, g);
                end
                continue;
            end
            speed_all{r,g} = thisSpeed{g}(:);
            pop_all{r,g}   = thisPop{g}(:);
            fprintf('%s: kept %d %s samples for group %d.\n', ...
                ratNames{r}, nS, Mode, g);
        end
    else
        if numel(thisSpeed) < MinSamples
            fprintf('%s: only %d %s samples, skipping this rat.\n', ...
                ratNames{r}, numel(thisSpeed), Mode);
            continue;
        end
        speed_all{r,1} = thisSpeed(:);
        pop_all{r,1}   = thisPop(:);
        fprintf('%s: kept %d %s samples.\n', ratNames{r}, numel(thisSpeed), Mode);
    end
end

% ---------- drop rats with no groups ----------
if doGroup
    keep = any(~cellfun(@isempty, speed_all), 2);
else
    keep = ~cellfun(@isempty, speed_all(:,1));
end

ratNames  = ratNames(keep);
speed_all = speed_all(keep,:);
pop_all   = pop_all(keep,:);
nR        = sum(keep);

if nR < 2
    warning('Need at least 2 rats with data for within vs across; found %d.', nR);
end

% effective number of groups among kept rats (still 1..GroupNum conceptually)
nG = GroupNum;

% ---------- xcorr params ----------
dt          = binSize;
maxLagBins  = round(MaxLagSec/dt);
lagsBins    = (-maxLagBins:maxLagBins).';
lagsSec     = lagsBins * dt;
nL          = numel(lagsBins);

% ---------- within-rat xcorr (truncate to same length) ----------
withinCurves     = nan(nL, nR, nG);
withinPeakR      = nan(nR, nG);
withinPeakLagSec = nan(nR, nG);
withinR0         = nan(nR, nG);

for r = 1:nR
    for g = 1:nG
        x = speed_all{r,g};
        y = pop_all{r,g};

        if isempty(x) || isempty(y)
            continue;
        end

        L = min(numel(x), numel(y));   % truncate both
        x = x(1:L);
        y = y(1:L);

        if L < 5
            continue;
        end

        xc = xcorr(x, y, maxLagBins, 'coeff');
        withinCurves(:,r,g) = xc(:);

        corr_direct = corr(x, y, 'rows','complete');
        if doGroup
            fprintf('%s (group %d) Pearson r = %.3f\n', ratNames{r}, g, corr_direct);
        else
            fprintf('%s Pearson r = %.3f\n', ratNames{r}, corr_direct);
        end

        [~,idx] = max(abs(xc));
        withinPeakR(r,g)      = xc(idx);
        withinPeakLagSec(r,g) = lagsSec(idx);

        i0 = find(lagsBins==0,1);
        if ~isempty(i0)
            withinR0(r,g) = xc(i0);
        end
    end
end

% ---------- across-rat (swapped) xcorr ----------
pairs = [];
for i = 1:nR
    for j = 1:nR
        if i ~= j
            pairs = [pairs; i j]; %#ok<AGROW>
        end
    end
end
nPairs = size(pairs,1);

acrossCurves     = nan(nL, nPairs, nG);
acrossPeakR      = nan(nPairs, nG);
acrossPeakLagSec = nan(nPairs, nG);
acrossR0         = nan(nPairs, nG);

for k = 1:nPairs
    i = pairs(k,1);  % speed source
    j = pairs(k,2);  % pop source

    for g = 1:nG
        x = speed_all{i,g};
        y = pop_all{j,g};

        if isempty(x) || isempty(y)
            continue;
        end

        L = min(numel(x), numel(y));   % truncate both
        x = x(1:L);
        y = y(1:L);

        if L < 5
            continue;
        end

        xc = xcorr(x, y, maxLagBins, 'coeff');
        acrossCurves(:,k,g) = xc(:);

        [~,idx] = max(abs(xc));
        acrossPeakR(k,g)      = xc(idx);
        acrossPeakLagSec(k,g) = lagsSec(idx);

        i0 = find(lagsBins==0,1);
        if ~isempty(i0)
            acrossR0(k,g) = xc(i0);
        end
    end
end

% ---------- pack output ----------
OUT = struct();
OUT.lagsSec          = lagsSec;
OUT.withinCurves     = withinCurves;
OUT.withinPeakR      = withinPeakR;
OUT.withinPeakLagSec = withinPeakLagSec;
OUT.withinR0         = withinR0;

OUT.acrossCurves     = acrossCurves;
OUT.acrossPairs      = pairs;
OUT.acrossPeakR      = acrossPeakR;
OUT.acrossPeakLagSec = acrossPeakLagSec;
OUT.acrossR0         = acrossR0;

OUT.meta = struct('Win',Win,'binSize',binSize,'SpeedThresh',vThresh, ...
                  'PopAgg',popAgg,'PopNorm',popNorm, ...
                  'MaxLagSec',MaxLagSec,'MinSamples',MinSamples, ...
                  'Mode',Mode,'TrialDemean',TrialDemean, ...
                  'GroupMode',GroupMode,'GroupVar',GroupVar, ...
                  'GroupNum',GroupNum,'GroupType',GroupType, ...
                  'ratNames',{ratNames});

if DoPlot
    plot_xcorr_curves(OUT);
end
end

% ================== LOCAL HELPERS ==================

function Scell = normalizeCApeaks_local(Sraw)
% Same idea as normalizeCApeaks in your code.
if iscell(Sraw)
    if isvector(Sraw)
        Scell = Sraw(:);
        return;
    end
    K = size(Sraw,2);
    scores = zeros(1,K);
    for j = 1:K
        col = Sraw(:,j);
        scores(j) = mean(cellfun(@(x) isnumeric(x) && isvector(x), col));
    end
    [bestScore, jbest] = max(scores);
    if bestScore == 0
        error('normalizeCApeaks_local:NoNumericColumn','No numeric spike-time column.');
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
error('normalizeCApeaks_local:BadType','Unsupported CA_peaks type: %s.', class(Sraw));
end

function [vt, vv] = velocityFromPos_local(pos)
% pos: [T x 3], columns: [t x y]
t  = pos(:,1);
xy = pos(:,2:3);
dt = diff(t);
dt(dt==0) = NaN;
v  = [0; sqrt(sum(diff(xy).^2,2)) ./ dt];  % cm/s
v(~isfinite(v)) = 0;
vt = t;
vv = v;
end

function [binSpeed, rateMat, effDur] = runBinsSpeedAndRates_withDur_local( ...
        S, CS, tPos, vt, vv, Win, binSec, vThr)
% Non-task bins of length binSec with v>=vThr; exclude [CS+Win(1), CS+Win(2)).

nCells = numel(S);

% mask out trial windows
csMask = false(size(tPos));
for k = 1:numel(CS)
    csMask = csMask | ((tPos >= CS(k)+Win(1)) & (tPos < CS(k)+Win(2)));
end
dtPos   = median(diff(tPos),'omitnan');
vOnPos  = interp1(vt, vv, tPos, 'linear','extrap');
runMask = ~csMask & (vOnPos >= vThr);

tmin = tPos(find(runMask,1,'first'));
tmax = tPos(find(runMask,1,'last'));
if isempty(tmin) || isempty(tmax)
    binSpeed = []; rateMat = []; effDur = []; return;
end

edges = tmin:binSec:tmax;
if edges(end) < tmax
    edges = [edges tmax];
end
nB = numel(edges)-1;

binSpeed = nan(nB,1);
effDur   = zeros(nB,1);
rateMat  = nan(nB,nCells);

for b = 1:nB
    iv = (tPos>=edges(b)) & (tPos<edges(b+1)) & runMask;   % masked samples in bin
    effDur(b) = sum(iv) * dtPos;
    if effDur(b) <= 0
        continue;
    end
    binSpeed(b) = mean(vOnPos(iv),'omitnan');

    intervals = maskToIntervals_local(tPos, iv, dtPos);
    for c = 1:nCells
        sp = S{c,1};
        if isempty(sp)
            rateMat(b,c) = 0;
            continue;
        end
        nSp = countInIntervals_local(sp, intervals);
        rateMat(b,c) = nSp / effDur(b);  % Hz
    end
end

ok = isfinite(binSpeed) & (effDur>0) & any(isfinite(rateMat),2);
binSpeed = binSpeed(ok);
rateMat  = rateMat(ok,:);
effDur   = effDur(ok);
end

function [binSpeed, rateMat, effDur] = taskBinsSpeedAndRates_withDur_local( ...
        S, CS, tPos, vt, vv, Win, binSec, vThr)
% Split each trial window [CS+Win(1), CS+Win(2)) into binSec chunks.

nCells  = numel(S);
dtPos   = median(diff(tPos),'omitnan');
vOnPos  = interp1(vt, vv, tPos, 'linear','extrap');

binSpeed = [];
effDur   = [];
rateMat  = [];

for k = 1:numel(CS)
    t0 = CS(k)+Win(1);
    t1 = CS(k)+Win(2);
    if ~isfinite(t0) || ~isfinite(t1) || t1<=t0
        continue;
    end

    edges = t0:binSec:t1;
    if edges(end) < t1
        edges = [edges t1];
    end

    for b = 1:numel(edges)-1
        b0 = edges(b);
        b1 = edges(b+1);
        inBin = (tPos>=b0) & (tPos<b1);

        if isempty(vThr)
            mask = inBin;                    % full subwindow
        else
            mask = inBin & (vOnPos>=vThr);   % gated samples only
        end

        dur = sum(mask) * dtPos;
        if dur <= 0
            continue;
        end

        sbin = mean(vOnPos(mask), 'omitnan');

        if isempty(vThr)
            intervals = [b0 b1];
        else
            intervals = maskToIntervals_local(tPos, mask, dtPos);
        end

        row = NaN(1,nCells);
        for c = 1:nCells
            sp = S{c};
            if isempty(sp)
                row(c) = 0;
                continue;
            end
            if isempty(vThr)
                nSp = sum(sp>=b0 & sp<b1);
            else
                nSp = countInIntervals_local(sp, intervals);
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

function [binSpeedOut, popVecOut, effDurOut] = trial_demean_bins( ...
        binSpeed, popVec, effDur, CS, Win, binSec)
% Demean speed and population activity within each trial.

nBinsPerTrial = round(diff(Win) / binSec);
if nBinsPerTrial <= 0
    binSpeedOut = binSpeed; popVecOut = popVec; effDurOut = effDur;
    return;
end

nTotal = numel(binSpeed);
nFull  = floor(nTotal / nBinsPerTrial);
if nFull == 0
    binSpeedOut = binSpeed; popVecOut = popVec; effDurOut = effDur;
    return;
end

idxKeep = 1:(nFull*nBinsPerTrial);

Smat = reshape(binSpeed(idxKeep),  nBinsPerTrial, nFull).';
Pmat = reshape(popVec(idxKeep),    nBinsPerTrial, nFull).';
Dmat = reshape(effDur(idxKeep),    nBinsPerTrial, nFull).';

% Demean per trial (row), ignoring NaNs
Smu = mean(Smat,2,'omitnan');
Pmu = mean(Pmat,2,'omitnan');

Smat = Smat - Smu;
Pmat = Pmat - Pmu;

binSpeedOut = Smat.'; binSpeedOut = binSpeedOut(:);
popVecOut   = Pmat.'; popVecOut   = popVecOut(:);
effDurOut   = Dmat.'; effDurOut   = effDurOut(:);

% tack on any leftover bins (partial final trial) unchanged
if nFull*nBinsPerTrial < nTotal
    tailIdx = (nFull*nBinsPerTrial+1):nTotal;
    binSpeedOut = [binSpeedOut; binSpeed(tailIdx)];
    popVecOut   = [popVecOut;   popVec(tailIdx)];
    effDurOut   = [effDurOut;   effDur(tailIdx)];
end
end

function intervals = maskToIntervals_local(t, mask, dt)
if nargin<3 || isempty(dt)
    dt = median(diff(t),'omitnan');
end
mask = mask(:);
t    = t(:);
if numel(t)~=numel(mask)
    error('maskToIntervals_local: size mismatch');
end
dm     = diff([false; mask; false]);
starts = find(dm==1);
ends   = find(dm==-1)-1;
intervals = [t(starts)  t(ends)+dt];   % include sample; extend ~dt
end

function n = countInIntervals_local(spikes, intervals)
if isempty(spikes) || isempty(intervals)
    n = 0;
    return;
end
spikes = spikes(:);
inside = false(size(spikes));
for k = 1:size(intervals,1)
    inside = inside | (spikes>=intervals(k,1) & spikes<intervals(k,2));
end
n = sum(inside);
end

function Y = normalizePerCell_local(X, mode)
switch lower(mode)
    case 'none'
        Y = X;
    case 'mean'     % divide by per-cell mean
        mu = mean(X,1,'omitnan');
        mu(mu==0 | ~isfinite(mu)) = NaN;
        Y = X ./ mu;
    case 'demean'   % subtract per-cell mean
        mu = mean(X,1,'omitnan');
        Y  = X - mu;
    otherwise
        error('PopNorm must be ''none'', ''mean'', or ''demean''.');
end
end

function y = aggregateAcrossCells_local(rateMat, modeStr)
switch lower(modeStr)
    case 'mean'
        y = mean(rateMat, 2, 'omitnan');
    case 'sum'
        y = sum(rateMat, 2, 'omitnan');
    otherwise
        error('PopAgg must be ''mean'' or ''sum''.');
end
end

function plot_xcorr_curves(OUT)
% Plot mean ± SEM xcorr curves for within- and across-animal results.

lags = OUT.lagsSec(:);
W    = OUT.withinCurves;   % [nL x nR x nG] or [nL x nR]
A    = OUT.acrossCurves;   % [nL x nPairs x nG] or [nL x nPairs]

hasGrouping = isfield(OUT.meta,'GroupMode') && strcmpi(OUT.meta.GroupMode,'on') ...
              && OUT.meta.GroupNum > 1;

if ~hasGrouping
    % ----- single-group behavior -----
    if ndims(W) == 3
        W = W(:,:,1);
    end
    if ndims(A) == 3
        A = A(:,:,1);
    end

    Wmean = mean(W, 2, 'omitnan');
    Wsem  = std(W, 0, 2, 'omitnan') ./ max(1, sqrt(sum(isfinite(W),2)));

    Amean = mean(A, 2, 'omitnan');
    Asem  = std(A, 0, 2, 'omitnan') ./ max(1, sqrt(sum(isfinite(A),2)));

    figure('Color','w','Position',[200 200 700 450]); hold on;

    % faint individual within curves
    if ~isempty(W)
        plot(lags, W, 'Color',[0 0.4 1 0.12], 'LineWidth',1);
    end

    % faint individual across curves
    if ~isempty(A)
        plot(lags, A, 'Color',[1 0 0 0.08], 'LineWidth',1);
    end

    % shaded SEM (within)
    if any(isfinite(Wmean))
        fill([lags; flipud(lags)], ...
             [Wmean-Wsem; flipud(Wmean+Wsem)], ...
             [0 0.4 1], 'FaceAlpha',0.25, 'EdgeColor','none');
    end

    % shaded SEM (across)
    if any(isfinite(Amean))
        fill([lags; flipud(lags)], ...
             [Amean-Asem; flipud(Amean+Asem)], ...
             [1 0 0], 'FaceAlpha',0.20, 'EdgeColor','none');
    end

    % mean curves
    plot(lags, Wmean, 'Color',[0 0.2 1], 'LineWidth',3);
    plot(lags, Amean, 'Color',[0.8 0 0], 'LineWidth',3);

    % zero-lag line
    yl = ylim;
    plot([0 0], yl, '--k');
    ylim(yl);

    xlabel('Lag (s)');
    ylabel('xcorr(speed, MUA)');
    title(sprintf('Within vs Across speed–MUA xcorr (%s, TrialDemean=%d)', ...
        OUT.meta.Mode, OUT.meta.TrialDemean));
    legend({'within (indiv)','across (indiv)', ...
            'within mean','across mean'}, ...
           'Location','best');
    box off;

else
    % ----- grouped version: skip groups with no data -----
    [nL, ~, nG] = size(W);

    hasData = false(1,nG);
    for g = 1:nG
        Wg = W(:,:,g);
        Ag = A(:,:,g);
        if any(isfinite(Wg),'all') || any(isfinite(Ag),'all')
            hasData(g) = true;
        end
    end
    gList  = find(hasData);
    nGplot = numel(gList);

    if nGplot == 0
        warning('No groups with finite data to plot.');
        return;
    end

    figure('Color','w','Position',[100 100 350*nGplot 450]);

    for ii = 1:nGplot
        g = gList(ii);
        subplot(1,nGplot,ii); hold on;

        Wg = W(:,:,g);          % [nL x nR]
        Ag = A(:,:,g);          % [nL x nPairs]

        Wmean = mean(Wg, 2, 'omitnan');
        Wsem  = std(Wg, 0, 2, 'omitnan') ./ max(1, sqrt(sum(isfinite(Wg),2)));

        Amean = mean(Ag, 2, 'omitnan');
        Asem  = std(Ag, 0, 2, 'omitnan') ./ max(1, sqrt(sum(isfinite(Ag),2)));

        % faint individual within curves
        if ~isempty(Wg)
            plot(lags, Wg, 'Color',[0 0.4 1 0.10], 'LineWidth',1);
        end

        % faint individual across curves
        if ~isempty(Ag)
            plot(lags, Ag, 'Color',[1 0 0 0.06], 'LineWidth',1);
        end

        % shaded SEM (within)
        if any(isfinite(Wmean))
            fill([lags; flipud(lags)], ...
                 [Wmean-Wsem; flipud(Wmean+Wsem)], ...
                 [0 0.4 1], 'FaceAlpha',0.25, 'EdgeColor','none');
        end

        % shaded SEM (across)
        if any(isfinite(Amean))
            fill([lags; flipud(lags)], ...
                 [Amean-Asem; flipud(Amean+Asem)], ...
                 [1 0 0], 'FaceAlpha',0.20, 'EdgeColor','none');
        end

        % mean curves
        plot(lags, Wmean, 'Color',[0 0.2 1], 'LineWidth',3);
        plot(lags, Amean, 'Color',[0.8 0 0], 'LineWidth',3);

        % zero-lag line
        yl = ylim;
        plot([0 0], yl, '--k');
        ylim(yl);

        xlabel('Lag (s)');
        if ii == 1
            ylabel('xcorr(speed, MUA)');
        end
        title(sprintf('Group %d (%s, TrialDemean=%d)', ...
            g, OUT.meta.Mode, OUT.meta.TrialDemean));

        if ii == nGplot
            legend({'within (indiv)','across (indiv)', ...
                    'within mean','across mean'}, ...
                   'Location','best');
        end
        box off;
    end
end
end
