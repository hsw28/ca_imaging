function SpikesPerBin(ratNames, varargin)
% SpikesPerBin(ratNames, ...)
% Build a (nBins × nCells) event-rate matrix per day and store in rat.spikesperbin.
% Optionally demean each cell by its pre-CS baseline mean (Hz).
%
% Example:
%   SpikesPerBin({'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%       'Demean',true,'BaselineWindow',[-0.5 0],'StoreBoth',true,'PerDay',true,'GrandPlot',true)
%
% Name-Value options:
%   'Demean'          true
%   'BaselineWindow'  [-0.5 0]
%   'StoreBoth'       true
%   'PerDay'          false
%   'DoPlot'          true
%   'Axes'            []      (see notes below)
%   'GrandPlot'       true
%   'EndTime'         4
%   'BinWidth'        (1/10)
%
% Axes usage:
%   PerDay=false:
%       - Axes = [] -> creates a figure
%       - Axes = vector length nRats   -> plots per-rat only
%       - Axes = vector length nRats+1 -> last axis is grand
%   PerDay=true:
%       - Axes = [] -> creates per-rat figure(s) + optional grand figure
%       - Axes can be 2x3 (top row rat day1-3, bottom row grand day1-3)
%         or vector length 6 in order [ratD1 ratD2 ratD3 grandD1 grandD2 grandD3]

% ---------- options ----------
p = inputParser;
p.addParameter('Demean', false, @(x)islogical(x) && isscalar(x));
p.addParameter('BaselineWindow', [-0.5 0], @(v)isnumeric(v) && numel(v)==2 && v(1)<v(2));
p.addParameter('StoreBoth', true, @(x)islogical(x) && isscalar(x));
p.addParameter('PerDay', false, @(x)islogical(x) && isscalar(x));
p.addParameter('DoPlot', true, @(x)islogical(x) && isscalar(x));
p.addParameter('Axes', [], @(x) isempty(x) || all(ishghandle(x)));
p.addParameter('GrandPlot', true, @(x)islogical(x) && isscalar(x));
p.addParameter('EndTime', 4, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('BinWidth', (1/15), @(x)isnumeric(x) && isscalar(x) && x>0);
p.parse(varargin{:});

Demean          = p.Results.Demean;
BaselineWindow  = p.Results.BaselineWindow;
StoreBoth       = p.Results.StoreBoth;
PerDay          = p.Results.PerDay;
DoPlot          = p.Results.DoPlot;
AxesIn          = p.Results.Axes;
GrandPlot       = p.Results.GrandPlot;
endtime         = p.Results.EndTime;
binWidth        = p.Results.BinWidth;

% --- bar border style (NEW) ---
barEdgeColor = [0 0 0];
barLineWidth = 0.75;

% ---------- analysis window & binning ----------
winFull  = [-0.5 endtime];
binEdges = winFull(1):binWidth:winFull(2);
nBins    = numel(binEdges)-1;
tAxis    = binEdges(1:end-1) + binWidth/2;   % bin centers

% baseline mask
baseMask = (tAxis >= BaselineWindow(1)) & (tAxis < BaselineWindow(2));

% ---------- main loop: compute + store per day ----------
for r = 1:numel(ratNames)
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    if isempty(idx) || idx < 3
        warning('%s has fewer than 3 learned sessions (An=%s). Skipping.', ratNames{r}, rat.An);
        continue;
    end
    days  = dates(idx-2:idx);

    for d = 1:3
        dayStr   = days{d};
        spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',dayStr));
        csTimes  = rat.CS_times.(sprintf('CS_%s',dayStr));
        ratemask = (rat.ratemask.(sprintf('ratemask_%s',dayStr)) == 1);

        nTotal   = size(spk,1);
        inclMask = ratemask(:) == 1;
        inclIdx  = find(inclMask);

        counts = zeros(nBins, nTotal);
        nTrialsUse = numel(csTimes);

        for t = 1:nTrialsUse
            t0 = csTimes(t);
            for jj = 1:numel(inclIdx)
                c  = inclIdx(jj);
                st = spk(c,:);
                st = st(~isnan(st) & st>0);
                rel = st - t0;
                rel = rel(rel >= binEdges(1) & rel < binEdges(end));
                if ~isempty(rel)
                    counts(:,c) = counts(:,c) + histcounts(rel, binEdges).';
                end
            end
        end

        ratesHz = counts / (binWidth * max(1, nTrialsUse));
        ratesHz(:, ~inclMask) = NaN;

        if Demean
            if ~any(baseMask)
                warning('BaselineWindow [%.3f %.3f) does not overlap bins; skipping demeaning for %s %s.', ...
                    BaselineWindow(1), BaselineWindow(2), ratNames{r}, dayStr);
                ratesDM = ratesHz;
            else
                baseMeans = nanmean(ratesHz(baseMask, :), 1);
                ratesDM   = ratesHz - baseMeans;
            end
        end

        if ~isfield(rat,'spikesperbin'), rat.spikesperbin = struct(); end

        if StoreBoth
            rat.spikesperbin.(sprintf('cpb_%s',dayStr)) = ratesHz.';        % nCells × nBins
            if Demean
                rat.spikesperbin.(sprintf('cpbDM_%s',dayStr)) = ratesDM.';  % nCells × nBins
            end
        else
            if Demean
                rat.spikesperbin.(sprintf('cpb_%s',dayStr)) = ratesDM.';
            else
                rat.spikesperbin.(sprintf('cpb_%s',dayStr)) = ratesHz.';
            end
        end
    end

    assignin('base', ratNames{r}, rat);
    fprintf('Stored spikesperbin for %s (days: %s, %s, %s)\n', ratNames{r}, days{:});
end

% ---------- plotting ----------
if ~DoPlot
    return;
end

nRats = numel(ratNames);
nBinsLocal = nBins; %#ok<NASGU>

% Precompute per-rat day means (for grand panels)
ratMu = cell(nRats,3);
ratDayLabels = cell(nRats,3);

for rr = 1:nRats
    rat   = evalin('base', ratNames{rr});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    if isempty(idx) || idx < 3
        for d=1:3, ratMu{rr,d} = nan(1,nBins); ratDayLabels{rr,d} = ''; end
        continue;
    end
    days  = dates(idx-2:idx);

    for d = 1:3
        ratDayLabels{rr,d} = days{d};
        if Demean && isfield(rat.spikesperbin, sprintf('cpbDM_%s',days{d}))
            X = rat.spikesperbin.(sprintf('cpbDM_%s',days{d}));
        else
            X = rat.spikesperbin.(sprintf('cpb_%s',days{d}));
        end
        ratMu{rr,d} = nanmean(X,1);
    end
end

if Demean
    ylab = '\Delta rate (Hz from baseline)';
else
    ylab = 'Event rate (Hz)';
end

% =========================
% PerDay = false (pooled)
% =========================
if ~PerDay
    % axes allocation
    if isempty(AxesIn)
        figure('Color','w');
        nAx = nRats + double(GrandPlot);
        axList = gobjects(1,nAx);
        for k = 1:nAx
            axList(k) = subplot(1,nAx,k);
        end
    else
        axList = AxesIn(:).';
        if GrandPlot && numel(axList) < nRats+1
            error('Axes must be length nRats+1 when GrandPlot=true (or omit Axes).');
        elseif ~GrandPlot && numel(axList) < nRats
            error('Axes must be length nRats when GrandPlot=false (or omit Axes).');
        end
    end

    grandMat = nan(nRats, nBins);

    for rr = 1:nRats
        ax = axList(rr);
        axes(ax); cla(ax); hold(ax,'on');

        rat   = evalin('base', ratNames{rr});
        dates = autoDateList(rat);
        idx   = find(strcmp(dates, rat.An),1);
        if isempty(idx) || idx < 3
            title(ax, sprintf('%s (insufficient days)', ratNames{rr}), 'Interpreter','none');
            box(ax,'off');
            continue;
        end
        days  = dates(idx-2:idx);

        % concatenate cells across 3 days for this rat
        M = [];
        for d = 1:3
            if Demean && isfield(rat.spikesperbin, sprintf('cpbDM_%s',days{d}))
                X = rat.spikesperbin.(sprintf('cpbDM_%s',days{d}));
            else
                X = rat.spikesperbin.(sprintf('cpb_%s',days{d}));
            end
            M = [M; X]; %#ok<AGROW>
        end

        mu  = nanmean(M,1);
        sem = nanstd(M,0,1) ./ sqrt(sum(~isnan(M),1));
        grandMat(rr,:) = mu;

        yTop = max(0, max(mu+sem)*1.1);
        fill(ax, [0 2 2 0], [0 0 yTop yTop], [.9 .9 .9], 'EdgeColor','none');

        bar(ax, tAxis, mu, 1, ...
            'FaceColor',[.5 .5 .8], ...
            'EdgeColor',barEdgeColor, ...
            'LineWidth',barLineWidth);
        errorbar(ax, tAxis, mu, sem, '.k', 'CapSize',4);

        yl = ylim(ax);
        line(ax, [0 0], yl, 'Color','r','LineStyle','--');
        line(ax, [.75 .75], yl, 'Color','b','LineStyle','--');

        xlabel(ax, 'Time from CS (s)');
        ylabel(ax, ylab);
        title(ax, ratNames{rr}, 'Interpreter','none');
        box(ax,'off');
    end

    if GrandPlot
        axG = axList(nRats+1);
        axes(axG); cla(axG); hold(axG,'on');

        grandMu  = nanmean(grandMat,1);
        grandSem = nanstd (grandMat,0,1) ./ sqrt(sum(~isnan(grandMat),1));

        yTop = max(0, max(grandMu+grandSem)*1.1);
        fill(axG, [0 2 2 0], [0 0 yTop yTop], [.9 .9 .9], 'EdgeColor','none');

        bar(axG, tAxis, grandMu, 1, ...
            'FaceColor',[.8 .5 .5], ...
            'EdgeColor',barEdgeColor, ...
            'LineWidth',barLineWidth);
        errorbar(axG, tAxis, grandMu, grandSem, '.k', 'CapSize',4);

        yl = ylim(axG);
        line(axG, [0 0], yl, 'Color','r','LineStyle','--');
        line(axG, [.75 .75], yl, 'Color','b','LineStyle','--');

        xlabel(axG, 'Time from CS (s)');
        ylabel(axG, ylab);
        title(axG, 'All rats', 'Interpreter','none');
        box(axG,'off');
    end

    return;
end

% =========================
% PerDay = true
% =========================
for rr = 1:nRats
    rat   = evalin('base', ratNames{rr});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    if isempty(idx) || idx < 3
        warning('%s: insufficient days for PerDay plot. Skipping plot.', ratNames{rr});
        continue;
    end
    days  = dates(idx-2:idx);

    if isempty(AxesIn)
        figure('Color','w');
        axRat = gobjects(1,3);
        for d = 1:3, axRat(d) = subplot(1,3,d); end
        axGrand = [];
        if GrandPlot
            figure('Color','w');
            axGrand = gobjects(1,3);
            for d = 1:3, axGrand(d) = subplot(1,3,d); end
        end
    else
        % AxesIn can be 3 (rat only) or 6 / 2x3 (rat + grand)
        if ismatrix(AxesIn) && size(AxesIn,1)==2 && size(AxesIn,2)==3
            axRat   = AxesIn(1,:);
            axGrand = AxesIn(2,:);
        else
            a = AxesIn(:).';
            if numel(a) < 3
                error('Axes must have at least 3 axes for PerDay=true.');
            end
            axRat = a(1:3);
            if GrandPlot
                if numel(a) < 6
                    error('Axes must have 6 axes (rat 3 + grand 3) when GrandPlot=true.');
                end
                axGrand = a(4:6);
            else
                axGrand = [];
            end
        end
    end

    % --- rat day-by-day ---
    for d = 1:3
        ax = axRat(d);
        axes(ax); cla(ax); hold(ax,'on');

        if Demean && isfield(rat.spikesperbin, sprintf('cpbDM_%s',days{d}))
            X = rat.spikesperbin.(sprintf('cpbDM_%s',days{d}));
        else
            X = rat.spikesperbin.(sprintf('cpb_%s',days{d}));
        end

        mu  = nanmean(X,1);
        sem = nanstd(X,0,1) ./ sqrt(sum(~isnan(X),1));

        yTop = max(0, max(mu+sem)*1.1);
        fill(ax, [0 2 2 0], [0 0 yTop yTop], [.9 .9 .9], 'EdgeColor','none');

        bar(ax, tAxis, mu, 1, ...
            'FaceColor',[.5 .5 .8], ...
            'EdgeColor',barEdgeColor, ...
            'LineWidth',barLineWidth);
        errorbar(ax, tAxis, mu, sem, '.k', 'CapSize',4);

        yl = ylim(ax);
        line(ax, [0 0], yl, 'Color','r','LineStyle','--');
        line(ax, [.75 .75], yl, 'Color','b','LineStyle','--');

        xlabel(ax, 'Time from CS (s)');
        ylabel(ax, ylab);
        title(ax, sprintf('%s: %s', ratNames{rr}, days{d}), 'Interpreter','none');
        box(ax,'off');
    end

    % --- grand day-by-day across rats ---
    if GrandPlot && ~isempty(axGrand)
        for d = 1:3
            ax = axGrand(d);
            axes(ax); cla(ax); hold(ax,'on');

            thisDayMat = nan(nRats, nBins);
            for r2 = 1:nRats
                thisDayMat(r2,:) = ratMu{r2,d};
            end

            grandMu  = nanmean(thisDayMat,1);
            grandSem = nanstd (thisDayMat,0,1) ./ sqrt(sum(~isnan(thisDayMat),1));

            yTop = max(0, max(grandMu+grandSem)*1.1);
            fill(ax, [0 2 2 0], [0 0 yTop yTop], [.9 .9 .9], 'EdgeColor','none');

            bar(ax, tAxis, grandMu, 1, ...
                'FaceColor',[.8 .5 .5], ...
                'EdgeColor',barEdgeColor, ...
                'LineWidth',barLineWidth);
            errorbar(ax, tAxis, grandMu, grandSem, '.k', 'CapSize',4);

            yl = ylim(ax);
            line(ax, [0 0], yl, 'Color','r','LineStyle','--');
            line(ax, [.75 .75], yl, 'Color','b','LineStyle','--');

            xlabel(ax, 'Time from CS (s)');
            ylabel(ax, ylab);
            title(ax, sprintf('All rats: day %d', d), 'Interpreter','none');
            box(ax,'off');
        end
    end
end
end
