function R = tracecellactivity(ratNames, varargin)
% tracecellactivity  Plot pooled CS-aligned activity for trace-modulated cells.
%
% Uses rat.traceneurons.tn_<day> vectors created by epochModulation_venn_dep:
%   1   = trace modulated
%   0   = analyzed but not trace modulated
%   NaN = excluded by rate filter / not analyzed
%
% Usage:
%   R = tracecellactivity;
%   R = tracecellactivity({'rat0222','rat0307','rat0313','rat0314','rat0816'});
%   R = tracecellactivity(ratNames, 'Normalize', 'zscore');

if nargin < 1 || isempty(ratNames)
    ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end

p = inputParser;
p.addParameter('Window',[-1 2],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Bin',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Normalize','max',@(s)any(strcmpi(s,{'max','zscore','none'})));
p.addParameter('SmoothBins',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinTrials',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('MakePlot',true,@islogical);
p.parse(varargin{:});
opt = p.Results;

ratNames = cellstr(ratNames);
win = opt.Window;
bin = opt.Bin;
edges = win(1):bin:win(2);
centers = edges(1:end-1) + bin/2;

epoch.CS    = [0 0.25];
epoch.Trace = [0.25 0.75];
epoch.US    = [0.75 0.85];

P = [];
meta = struct('rat',{},'day',{},'cell',{},'peakTime',{});

for r = 1:numel(ratNames)
    ratName = ratNames{r};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratName))
        warning('%s not found in base workspace. Skipping.', ratName);
        continue
    end

    rat = evalin('base', ratName);
    days = lastThreeAnalysisDays(rat);

    for d = 1:numel(days)
        dayStr = days{d};
        peakFld = sprintf('CA_peaks_%s', dayStr);
        csFld = sprintf('CS_%s', dayStr);
        tnFld = sprintf('tn_%s', dayStr);

        if ~isfield(rat,'Ca_peaks') || ~isfield(rat.Ca_peaks, peakFld) || ...
           ~isfield(rat,'CS_times') || ~isfield(rat.CS_times, csFld)
            warning('%s %s missing Ca_peaks or CS_times. Skipping.', ratName, dayStr);
            continue
        end
        if ~isfield(rat,'traceneurons') || ~isfield(rat.traceneurons, tnFld)
            warning('%s missing traceneurons.%s. Run epochModulation_venn_dep first. Skipping %s.', ...
                ratName, tnFld, dayStr);
            continue
        end

        spk = rat.Ca_peaks.(peakFld);
        csTimes = rat.CS_times.(csFld);
        csTimes = csTimes(:);
        tn = rat.traceneurons.(tnFld);
        tn = tn(:);

        nCells = min(size(spk,1), numel(tn));
        if nCells < size(spk,1) || nCells < numel(tn)
            warning('%s %s cell count mismatch between CA_peaks and traceneurons; using first %d cells.', ...
                ratName, dayStr, nCells);
        end
        if numel(csTimes) < opt.MinTrials
            warning('%s %s has only %d trials. Skipping.', ratName, dayStr, numel(csTimes));
            continue
        end

        traceCells = find(tn(1:nCells) == 1);
        for ii = 1:numel(traceCells)
            c = traceCells(ii);
            psth = cellPSTH(spk(c,:), csTimes, edges, bin);
            if all(~isfinite(psth)) || all(psth == 0)
                continue
            end

            psth = smoothVector(psth, opt.SmoothBins);
            psthNorm = normalizePSTH(psth, centers, opt.Normalize);
            peakTime = tracePeakTime(psthNorm, centers, epoch.Trace);

            P = [P; psthNorm]; %#ok<AGROW>
            meta(end+1).rat = ratName; %#ok<AGROW>
            meta(end).day = dayStr;
            meta(end).cell = c;
            meta(end).peakTime = peakTime;
        end
    end
end

peakTimes = [meta.peakTime]';
if isempty(P)
    warning('No trace-modulated cell PSTHs were found.');
    R.opts = opt;
    R.ratNames = ratNames;
    R.time = centers;
    R.psth = P;
    R.psthSorted = P;
    R.peakTimes = peakTimes;
    R.peakTimesSorted = peakTimes;
    R.meta = meta;
    R.metaSorted = meta;
    R.epochs = epoch;
    return
end

peakForSort = peakTimes;
peakForSort(~isfinite(peakForSort)) = inf;
[~, sortIdx] = sort(peakForSort, 'ascend');
Psort = P(sortIdx,:);
metaSort = meta(sortIdx);
peakTimesSort = peakTimes(sortIdx);

R.opts = opt;
R.ratNames = ratNames;
R.time = centers;
R.psth = P;
R.psthSorted = Psort;
R.peakTimes = peakTimes;
R.peakTimesSorted = peakTimesSort;
R.meta = meta;
R.metaSorted = metaSort;
R.epochs = epoch;

if opt.MakePlot
    plotTraceCellActivity(R);
end
end

function days = lastThreeAnalysisDays(rat)
dates = autoDateList(rat);
idx = find(strcmp(dates, rat.An), 1);
if isempty(idx)
    idx = numel(dates);
end
days = dates(max(1,idx-2):idx);
end

function psth = cellPSTH(st, csTimes, edges, bin)
st = st(~isnan(st) & st>0);
nTrials = numel(csTimes);
counts = zeros(1, numel(edges)-1);
for t = 1:nTrials
    rel = st(st >= csTimes(t)+edges(1) & st < csTimes(t)+edges(end)) - csTimes(t);
    if ~isempty(rel)
        counts = counts + histcounts(rel, edges);
    end
end
psth = counts ./ (nTrials * bin);
end

function y = smoothVector(x, smoothBins)
if smoothBins <= 0
    y = x;
    return
end
w = ones(1, 2*round(smoothBins)+1);
w = w ./ sum(w);
y = conv(x, w, 'same');
end

function y = normalizePSTH(x, t, method)
switch lower(method)
    case 'max'
        base = min(x(isfinite(x)));
        y = x - base;
        mx = max(y(isfinite(y)));
        if mx > 0
            y = y ./ mx;
        end
    case 'zscore'
        baseMask = t < 0;
        mu = mean(x(baseMask), 'omitnan');
        sd = std(x(baseMask), 0, 'omitnan');
        y = (x - mu) ./ (sd + eps);
    otherwise
        y = x;
end
end

function pt = tracePeakTime(psth, t, traceWin)
traceMask = t >= traceWin(1) & t < traceWin(2);
if ~any(traceMask)
    pt = NaN;
    return
end
v = psth(traceMask);
tt = t(traceMask);
if all(~isfinite(v))
    pt = NaN;
    return
end
[~, idx] = max(v);
pt = tt(idx);
end

function plotTraceCellActivity(R)
t = R.time;
P = R.psthSorted;
epoch = R.epochs;

figure('Color','w','Position',[100 100 950 430]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on
imagesc(t, 1:size(P,1), P);
set(gca,'YDir','normal');
shadeEpochs(gca, epoch, [0.02 0.02 0.02], 0.08);
xline(epoch.CS(2),'k:','Trace onset');
xline(epoch.Trace(2),'k:','US onset');
xlim(R.opts.Window);
xlabel('Time from CS (s)');
ylabel('Trace-modulated cells');
title(sprintf('Normalized PSTH heatmap (n = %d)', size(P,1)));
colormap(gca, parula);
colorbar;
box off

nexttile; hold on
pk = R.peakTimes(isfinite(R.peakTimes));
histogram(pk, epoch.Trace(1):R.opts.Bin:epoch.Trace(2), ...
    'FaceColor',[0.15 0.55 0.45], 'EdgeColor','none');
xline(epoch.Trace(1),'k:','Trace onset');
xline(epoch.Trace(2),'k:','US onset');
xlim(epoch.Trace);
xlabel('Peak time within trace (s)');
ylabel('# cells');
title('Trace peak-time distribution');
box off

title(tl, 'Trace-cell CS-aligned activity');
end

function shadeEpochs(ax, epoch, color, alphaVal)
yl = ylim(ax);
patch(ax, [epoch.CS(1) epoch.CS(2) epoch.CS(2) epoch.CS(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], color, 'FaceAlpha',alphaVal, 'EdgeColor','none');
patch(ax, [epoch.Trace(1) epoch.Trace(2) epoch.Trace(2) epoch.Trace(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], color, 'FaceAlpha',alphaVal*1.35, 'EdgeColor','none');
patch(ax, [epoch.US(1) epoch.US(2) epoch.US(2) epoch.US(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], color, 'FaceAlpha',alphaVal*1.8, 'EdgeColor','none');
img = findobj(ax,'Type','image');
if ~isempty(img)
    uistack(img,'bottom');
else
    patches = findobj(ax,'Type','patch');
    uistack(patches,'bottom');
end
end
