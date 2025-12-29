function plotTaskVar_InOut_vs_Distance(ratNames, varargin)
% plotTaskVar_InOut_vs_Distance
% Core outputs (per rat + pooled):
%   1) Paired per-cell means: in-PF vs out-PF (cells that have both)
%   2) Distance between trial centroids vs variability (|Δ|/mean) (all trial pairs)
%   3) Trial firing rate vs distance to PF center (all trials)
%
% RATEMASK RULE:
%   For each day 'YYYY_MM_DD', use rat.ratemask.ratemask_YYYY_MM_DD (Ncells x 1).
%   Keep cells where mv==1. If ratemask missing or size mismatch, skip that day.

% ---------- args ----------
p = inputParser;
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty0(d));
addParameter(p,'WinSecs',[0 2], @(v) isnumeric(v)&&numel(v)==2&&v(2)>v(1));
addParameter(p,'SpeedThresh',4, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'NSD',1, @(x) isnumeric(x)&&isscalar(x)&&isfinite(x));
addParameter(p,'MinTrials',10, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'MinSpikes',0, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'CentroidMethod','mode', ...
    @(s) any(validatestring(lower(s),{'mode','medoid','mean','speed'})));
addParameter(p,'CentroidBinCm',2.5, @(x) isnumeric(x)&&isscalar(x)&&x>0);
parse(p,varargin{:});
daysArg        = p.Results.Days;
winSecs        = p.Results.WinSecs;
vMin           = p.Results.SpeedThresh;
Nsd            = p.Results.NSD;
minT           = p.Results.MinTrials;
minSpk         = p.Results.MinSpikes;
centroidMethod = lower(p.Results.CentroidMethod);
centroidBinCm  = p.Results.CentroidBinCm;

if ischar(ratNames) || isstring(ratNames), ratNames = cellstr(ratNames); end
nRats = numel(ratNames);

% ---------- compute R per rat ----------
ratResults = cell(1, nRats);
for r = 1:nRats
    ratName = ratNames{r};
    ratResults{r} = getRatR(ratName, daysArg, winSecs, vMin, Nsd, minT, minSpk, centroidMethod, centroidBinCm);
end

% ---------- pooled containers ----------
pooled_in_means  = [];
pooled_out_means = [];

pooled_dist      = [];
pooled_cv        = [];

pooled_dist_pf   = [];
pooled_rates_all = [];

% ---------- FIG 1: paired in/out per cell ----------
nRows = nRats + 1;
figure('Color','w','Position',[100 50 900 250 + 180*nRows]);
t1 = tiledlayout(nRows,1,'TileSpacing','compact','Padding','compact');

for r = 1:nRats
    R = ratResults{r};
    ratName = ratNames{r};

    [inVec, outVec] = collectInOutPerCell(R, ratName);

    nexttile; hold on; box on
    if isempty0(inVec)
        text(0.5,0.5,'No cells with both in/out-PF (after ratemask)','HorizontalAlignment','center');
        axis off
    else
        pairedDotPlot(inVec, outVec);
        [~,p_pair,~,stats] = ttest(inVec, outVec);
        title(sprintf('%s | paired in/out  (p=%.3g, t=%.2f, n=%d)', ratName, p_pair, stats.tstat, numel(inVec)));
    end
    ylabel(ratName);
    if r==nRats, xlabel('Condition'); end

    pooled_in_means  = [pooled_in_means;  inVec];   %#ok<AGROW>
    pooled_out_means = [pooled_out_means; outVec];  %#ok<AGROW>
end

% pooled row for Fig 1
nexttile; hold on; box on
if isempty0(pooled_in_means)
    text(0.5,0.5,'No pooled cells with both in/out-PF (after ratemask)','HorizontalAlignment','center');
    axis off
else
    pairedDotPlot(pooled_in_means, pooled_out_means);
    [~,p_pair,~,stats] = ttest(pooled_in_means, pooled_out_means);
    title(sprintf('POOLED | paired in/out  (p=%.3g, t=%.2f, n=%d)', p_pair, stats.tstat, numel(pooled_in_means)));
end
ylabel('POOLED');
title(t1,'Paired per-cell means (ratemask-filtered)','FontWeight','bold');

% ---------- FIG 2: dist vs var (rows) ----------
figure('Color','w','Position',[100 50 900 250 + 180*nRows]);
t2 = tiledlayout(nRows,1,'TileSpacing','compact','Padding','compact');

for r = 1:nRats
    R = ratResults{r};
    ratName = ratNames{r};

    [distAll, cvAll] = collectDistVsVarAllTrials(R, ratName);

    nexttile; hold on; box on
    if isempty0(distAll)
        axis off
        text(0.5,0.5,'No trial pairs (after ratemask)','HorizontalAlignment','center');
        title(sprintf('%s | dist vs var', ratName));
    else
        tallBinScatterFit(gca, distAll, cvAll, 60, 'Nshow', 5e4);
        xlabel('Centroid distance between trials (cm)');
        ylabel('|Δ rate| / mean');
        [m,b,rho,pv] = fitLineAndCorr(distAll, cvAll);
        xx = linspace(min(distAll), max(distAll), 200);
        plot(xx, m*xx + b, 'k-', 'LineWidth',1.2);
        title(sprintf('%s | dist vs var  r=%.2f, p=%.3g  (n=%d pairs)', ratName, rho, pv, numel(cvAll)));
    end

    pooled_dist = [pooled_dist; distAll]; %#ok<AGROW>
    pooled_cv   = [pooled_cv;   cvAll];   %#ok<AGROW>
end

% pooled row for Fig 2
nexttile; hold on; box on
if isempty0(pooled_dist)
    axis off
    text(0.5,0.5,'No pooled trial pairs (after ratemask)','HorizontalAlignment','center');
else
    tallBinScatterFit(gca, pooled_dist, pooled_cv, 80, 'Nshow', 1e5);
    xlabel('Centroid distance between trials (cm)');
    ylabel('|Δ rate| / mean');
    [m,b,rho,pv] = fitLineAndCorr(pooled_dist, pooled_cv);
    xx = linspace(min(pooled_dist), max(pooled_dist), 300);
    plot(xx, m*xx + b, 'k-', 'LineWidth',1.2);
    title(sprintf('POOLED | dist vs var  r=%.2f, p=%.3g  (n=%d pairs)', rho, pv, numel(pooled_cv)));
end
title(t2,'Distance vs variability (ratemask-filtered)','FontWeight','bold');

% ---------- FIG 3: rate vs PF distance (rows) ----------
figure('Color','w','Position',[100 50 900 250 + 180*nRows]);
t3 = tiledlayout(nRows,1,'TileSpacing','compact','Padding','compact');

for r = 1:nRats
    R = ratResults{r};
    ratName = ratNames{r};

    [distPF, trialRates] = collectRateVsPFdistance(R, ratName);

    nexttile; hold on; box on
    if isempty0(distPF)
        axis off
        text(0.5,0.5,'No rate vs PF-distance data (after ratemask)','HorizontalAlignment','center');
        title(sprintf('%s | rate vs PF-dist', ratName));
    else
        tallBinScatterFit(gca, distPF, trialRates, 60, 'Nshow', 5e4);
        xlabel('Distance to PF center (cm)');
        ylabel('In-trial firing rate (Hz)');
        [m2,b2,r2,p2] = fitLineAndCorr(distPF, trialRates);
        xx2 = linspace(min(distPF), max(distPF), 200);
        plot(xx2, m2*xx2 + b2, 'k-', 'LineWidth',1.2);
        title(sprintf('%s | rate vs PF-dist  r=%.2f, p=%.3g  (n=%d trials)', ratName, r2, p2, numel(trialRates)));
    end

    pooled_dist_pf   = [pooled_dist_pf;   distPF(:)];      %#ok<AGROW>
    pooled_rates_all = [pooled_rates_all; trialRates(:)];  %#ok<AGROW>
end

% pooled row for Fig 3
nexttile; hold on; box on
if isempty0(pooled_dist_pf)
    axis off
    text(0.5,0.5,'No pooled rate vs PF-distance data (after ratemask)','HorizontalAlignment','center');
else
    tallBinScatterFit(gca, pooled_dist_pf, pooled_rates_all, 80, 'Nshow', 1e5);
    xlabel('Distance to PF center (cm)');
    ylabel('In-trial firing rate (Hz)');
    [m2,b2,r2,p2] = fitLineAndCorr(pooled_dist_pf, pooled_rates_all);
    xx2 = linspace(min(pooled_dist_pf), max(pooled_dist_pf), 300);
    plot(xx2, m2*xx2 + b2, 'k-', 'LineWidth',1.2);
    title(sprintf('POOLED | rate vs PF-dist  r=%.2f, p=%.3g  (n=%d trials)', r2, p2, numel(pooled_rates_all)));
end
title(t3,'Rate vs distance-to-PF (ratemask-filtered)','FontWeight','bold');

end % main


% ======================================================================
% Helpers
% ======================================================================

function R = getRatR(ratName, daysArg, winSecs, vMin, Nsd, minT, minSpk, centroidMethod, centroidBinCm)
R = taskFiringVariability(ratName, ...
        'Days',daysArg, ...
        'WinSecs',winSecs, 'SpeedThresh',vMin, 'NSD',Nsd, ...
        'MinTrials',minT, 'MinSpikes',minSpk, ...
        'CentroidMethod',centroidMethod, 'CentroidBinCm',centroidBinCm, ...
        'Plot',false);
end

function [inVec, outVec] = collectInOutPerCell(R, ratName)
inVec  = [];
outVec = [];
if ~isfield(R,'perDay') || isempty0(R.perDay), return, end

for d = 1:numel(R.perDay)
    if ~isfield(R.perDay(d),'perCell'), continue, end
    PC = R.perDay(d).perCell;
    if isempty0(PC), continue, end

    dayStr  = getDayStr(R.perDay(d));
    keepIdx = ratemask_keepIdx(ratName, dayStr, numel(PC));
    if isempty0(keepIdx), continue, end

    for c = keepIdx
        if ~isfield(PC(c),'rates') || isempty0(PC(c).rates) || ~isfield(PC(c),'inPF_trial') || isempty0(PC(c).inPF_trial)
            continue
        end
        rates = PC(c).rates(:);
        inpf  = PC(c).inPF_trial(:);
        if numel(inpf) ~= numel(rates), continue, end

        mi = mean(rates(inpf==1), 'omitnan');
        mo = mean(rates(inpf==0), 'omitnan');
        ni = nnz(inpf==1 & isfinite(rates));
        no = nnz(inpf==0 & isfinite(rates));

        if isfinite(mi) && isfinite(mo) && ni>=1 && no>=1
            inVec  = [inVec;  mi]; %#ok<AGROW>
            outVec = [outVec; mo]; %#ok<AGROW>
        end
    end
end
end

function [distAll, cvAll] = collectDistVsVarAllTrials(R, ratName)
distAll = [];
cvAll   = [];
if ~isfield(R,'perDay') || isempty0(R.perDay), return, end

for d = 1:numel(R.perDay)
    if ~isfield(R.perDay(d),'perCell'), continue, end
    PC = R.perDay(d).perCell;
    if isempty0(PC), continue, end

    dayStr  = getDayStr(R.perDay(d));
    keepIdx = ratemask_keepIdx(ratName, dayStr, numel(PC));
    if isempty0(keepIdx), continue, end

    for c = keepIdx
        if ~isfield(PC(c),'rates')     || isempty0(PC(c).rates) || ...
           ~isfield(PC(c),'centroids') || isempty0(PC(c).centroids)
            continue
        end

        rates = PC(c).rates(:);
        cents = PC(c).centroids;

        if size(cents,1) ~= numel(rates)
            n = min(numel(rates), size(cents,1));
            if n < 2, continue, end
            rates = rates(1:n);
            cents = cents(1:n,:);
        end

        keep  = isfinite(rates) & all(isfinite(cents),2);
        rates = rates(keep);
        cents = cents(keep,:);
        if numel(rates) < 2, continue, end

        D = squareform(pdist(cents,'euclidean'));
        [I,J] = find(triu(true(size(D)),1));
        dvec  = D(sub2ind(size(D), I, J));
        diffA = abs(rates(I) - rates(J));
        cv    = diffA ./ max(eps, 0.5*(rates(I)+rates(J)));

        distAll = [distAll; dvec(:)]; %#ok<AGROW>
        cvAll   = [cvAll;   cv(:)];   %#ok<AGROW>
    end
end
end

function [distPF, rates] = collectRateVsPFdistance(R, ratName)
distPF = [];
rates  = [];
if ~isfield(R,'perDay') || isempty0(R.perDay), return, end

for d = 1:numel(R.perDay)
    if ~isfield(R.perDay(d),'perCell'), continue, end
    PC = R.perDay(d).perCell;
    if isempty0(PC), continue, end

    dayStr  = getDayStr(R.perDay(d));
    keepIdx = ratemask_keepIdx(ratName, dayStr, numel(PC));
    if isempty0(keepIdx), continue, end

    for c = keepIdx
        if ~isfield(PC(c),'rates') || isempty0(PC(c).rates)
            continue
        end
        r = PC(c).rates(:);

        % Prefer stored dist_to_pf, else recompute from centroids + pf_center
        dp = [];
        have_dp     = isfield(PC(c),'dist_to_pf') && ~isempty0(PC(c).dist_to_pf);
        have_cent   = isfield(PC(c),'centroids') && ~isempty0(PC(c).centroids);
        have_center = isfield(PC(c),'pf_center') && numel(PC(c).pf_center)==2 && all(isfinite(PC(c).pf_center));

        if have_dp
            dp = PC(c).dist_to_pf(:);
        end

        if isempty0(dp) || numel(dp) ~= numel(r)
            if have_cent && have_center
                cents = PC(c).centroids;
                n = min(numel(r), size(cents,1));
                if n < 1, continue, end
                r = r(1:n);
                cents = cents(1:n,:);
                ctr = PC(c).pf_center(:).';
                dp = sqrt(sum((cents - ctr).^2, 2));
            else
                if ~isempty0(dp)
                    n = min(numel(r), numel(dp));
                    if n < 1, continue, end
                    r  = r(1:n);
                    dp = dp(1:n);
                else
                    continue
                end
            end
        end

        good = isfinite(r) & isfinite(dp);
        if any(good)
            distPF = [distPF; dp(good)]; %#ok<AGROW>
            rates  = [rates;  r(good)];  %#ok<AGROW>
        end
    end
end
end

function pairedDotPlot(inVec, outVec)
n = numel(inVec);
x1 = ones(n,1);
x2 = 2*ones(n,1);
for k = 1:n
    plot([x1(k) x2(k)], [inVec(k) outVec(k)], '-', 'Color',[0 0 0 0.2]); hold on
end
scatter(x1, inVec, 18, 'filled');
scatter(x2, outVec,18, 'filled');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'In-PF','Out-PF'});
ylabel('In-trial firing rate (Hz)');
end

function [m,b,rho,pv] = fitLineAndCorr(x,y)
x = x(:); y = y(:);
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);
if numel(x) < 3
    m = NaN; b = NaN; rho = NaN; pv = NaN; return
end
P = polyfit(x, y, 1);
m = P(1); b = P(2);
[rho,pv] = corr(x,y,'Rows','complete','Type','Pearson');
end

function tallBinScatterFit(ax, x, y, nbins, varargin)
ip = inputParser;
addParameter(ip,'Nshow',5e4,@(z)isnumeric(z)&&isscalar(z)&&z>=0);
parse(ip,varargin{:});
Nshow = ip.Results.Nshow;

hold(ax,'on'); box(ax,'on');

x = x(:); y = y(:);
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);
n = numel(x);

if Nshow>0 && n>Nshow
    idx = randperm(n, Nshow);
else
    idx = 1:n;
end
scatter(ax, x(idx), y(idx), 10, 'filled', 'MarkerFaceAlpha',0.25);

% (bin overlay intentionally omitted, like your current version)
end

% ---------- ratemask + day helpers ----------

function dayStr = getDayStr(perDayStruct)
% Pull a day string like '2023_05_09' from common perDay fields.
dayStr = '';
cands = {'day','date','dayStr','dateStr','Day','Date'};
for i = 1:numel(cands)
    f = cands{i};
    if isfield(perDayStruct,f) && ~isempty0(perDayStruct.(f))
        v = perDayStruct.(f);
        if isstring(v) || ischar(v)
            dayStr = char(v);
            return
        end
    end
end
end

function keepIdx = ratemask_keepIdx(ratName, dayStr, nCells)
% rat.ratemask.ratemask_YYYY_MM_DD is Ncells x 1, with 1=keep, 0=drop.
% Strict: missing/mismatch => [] (skip day).
keepIdx = [];

if isempty0(dayStr), return, end
if ~evalin('base', sprintf('exist(''%s'',''var'')', ratName)), return, end
rat = evalin('base', ratName);

maskField = sprintf('ratemask_%s', dayStr);
if ~isfield(rat,'ratemask') || ~isfield(rat.ratemask, maskField), return, end

mv = rat.ratemask.(maskField);
mv = mv(:);

if numel(mv) ~= nCells
    return
end

keepIdx = find(mv == 1);
end


function tf = isempty0(varargin)
% Robust wrapper around builtin isempty to avoid path shadowing.
% Accepts any #args so it never throws "Too many input arguments".
if nargin < 1
    tf = true;
else
    tf = builtin('isempty', varargin{1});
end
end
