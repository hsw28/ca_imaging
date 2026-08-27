function SpikesPerBin_extinction(ratNames, varargin)
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
%   'StartTime'       -1
%   'EndTime'         3
%   'BinWidth'        (1/10)
%   'DropFirstTrials' 0       (use all trials by default)
%   'ExtinctionDays'  3       (1 = final day, 2 = final two days, etc.)
%
% Axes usage:
%   PerDay=false:
%       - Axes = [] -> creates a figure
%       - Axes = vector length nRats   -> plots per-rat only
%       - Axes = vector length nRats+1 -> last axis is grand
%   PerDay=true:
%       - Axes = [] -> creates per-rat figure(s) + optional grand figure
%       - Axes can be 2 x ExtinctionDays (top row rat days, bottom row grand days)
%         or vector length 2*ExtinctionDays in order [rat days, grand days]

% ---------- options ----------
p = inputParser;
p.addParameter('Demean', false, @(x)islogical(x) && isscalar(x));
p.addParameter('BaselineWindow', [-0.5 0], @(v)isnumeric(v) && numel(v)==2 && v(1)<v(2));
p.addParameter('StoreBoth', true, @(x)islogical(x) && isscalar(x));
p.addParameter('PerDay', false, @(x)islogical(x) && isscalar(x));
p.addParameter('DoPlot', true, @(x)islogical(x) && isscalar(x));
p.addParameter('Axes', [], @(x) isempty(x) || all(ishghandle(x)));
p.addParameter('GrandPlot', true, @(x)islogical(x) && isscalar(x));
p.addParameter('StartTime', -1, @(x)isnumeric(x) && isscalar(x));
p.addParameter('EndTime', 3, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('BinWidth', (1/15), @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('DropFirstTrials', 0, @(x)isnumeric(x) && isscalar(x) && x>=0 && mod(x,1)==0);
p.addParameter('ExtinctionDays', 2, @(x)isnumeric(x) && isscalar(x) && x>=1 && mod(x,1)==0);
p.parse(varargin{:});

Demean          = p.Results.Demean;
BaselineWindow  = p.Results.BaselineWindow;
StoreBoth       = p.Results.StoreBoth;
PerDay          = p.Results.PerDay;
DoPlot          = p.Results.DoPlot;
AxesIn          = p.Results.Axes;
GrandPlot       = p.Results.GrandPlot;
starttime       = p.Results.StartTime;
endtime         = p.Results.EndTime;
binWidth        = p.Results.BinWidth;
dropFirstN      = p.Results.DropFirstTrials;
nExtDays        = p.Results.ExtinctionDays;

if endtime <= starttime
    error('EndTime must be greater than StartTime.');
end

% --- bar border style (NEW) ---
barEdgeColor = [0 0 0];
barLineWidth = 0.75;

% ---------- analysis window & binning ----------
winFull  = [starttime endtime];
binEdges = winFull(1):binWidth:winFull(2);
nBins    = numel(binEdges)-1;
tAxis    = binEdges(1:end-1) + binWidth/2;   % bin centers

% baseline mask
baseMask = (tAxis >= BaselineWindow(1)) & (tAxis < BaselineWindow(2));

% ---------- main loop: compute + store per day ----------
for r = 1:numel(ratNames)
    rat   = evalin('base', ratNames{r});
    days  = extinctionDateList(rat, nExtDays);
    if numel(days) < nExtDays
        warning('%s has fewer than %d extinction sessions. Skipping.', ratNames{r}, nExtDays);
        continue;
    end

    for d = 1:nExtDays
        dayStr   = days{d};
        spkField = findFieldForDate(rat.Ca_peaks, dayStr, {'CA_peaks_', 'Ca_peaks_'});
        csField  = findFieldForDate(rat.CS_times, dayStr, ...
            {'CS_exinction_', 'CS_extinction_', 'CS_'});
        spk      = rat.Ca_peaks.(spkField);
        csTimes  = rat.CS_times.(csField);
        csTimes  = dropFirstTrials(csTimes, dropFirstN);

        nTotal   = size(spk,1);
        inclIdx  = [1:size(spk,1)];

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
    fprintf('Stored spikesperbin for %s (days: %s)\n', ratNames{r}, strjoin(days, ', '));
end

% ---------- plotting ----------
if ~DoPlot
    return;
end

nRats = numel(ratNames);
nBinsLocal = nBins; %#ok<NASGU>

% Precompute per-rat day means (for grand panels)
ratMu = cell(nRats,nExtDays);
ratDayLabels = cell(nRats,nExtDays);

for rr = 1:nRats
    rat   = evalin('base', ratNames{rr});
    days  = extinctionDateList(rat, nExtDays);
    if numel(days) < nExtDays
        for d=1:nExtDays, ratMu{rr,d} = nan(1,nBins); ratDayLabels{rr,d} = ''; end
        continue;
    end

    for d = 1:nExtDays
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
        days  = extinctionDateList(rat, nExtDays);
        if numel(days) < nExtDays
            title(ax, sprintf('%s (insufficient days)', ratNames{rr}), 'Interpreter','none');
            box(ax,'off');
            continue;
        end

        % concatenate cells across selected extinction days for this rat
        M = [];
        for d = 1:nExtDays
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
    days  = extinctionDateList(rat, nExtDays);
    if numel(days) < nExtDays
        warning('%s: insufficient days for PerDay plot. Skipping plot.', ratNames{rr});
        continue;
    end

    if isempty(AxesIn)
        figure('Color','w');
        axRat = gobjects(1,nExtDays);
        for d = 1:nExtDays, axRat(d) = subplot(1,nExtDays,d); end
        axGrand = [];
        if GrandPlot
            figure('Color','w');
            axGrand = gobjects(1,nExtDays);
            for d = 1:nExtDays, axGrand(d) = subplot(1,nExtDays,d); end
        end
    else
        % AxesIn can be nExtDays (rat only) or 2*nExtDays (rat + grand)
        if ismatrix(AxesIn) && size(AxesIn,1)==2 && size(AxesIn,2)==nExtDays
            axRat   = AxesIn(1,:);
            axGrand = AxesIn(2,:);
        else
            a = AxesIn(:).';
            if numel(a) < nExtDays
                error('Axes must have at least ExtinctionDays axes for PerDay=true.');
            end
            axRat = a(1:nExtDays);
            if GrandPlot
                if numel(a) < 2*nExtDays
                    error('Axes must have 2*ExtinctionDays axes (rat + grand) when GrandPlot=true.');
                end
                axGrand = a(nExtDays+1:2*nExtDays);
            else
                axGrand = [];
            end
        end
    end

    % --- rat day-by-day ---
    for d = 1:nExtDays
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
        for d = 1:nExtDays
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

function days = extinctionDateList(rat, nDays)
% Use the actual extinction CS fields; autoDateList only looks at Ca_traces.
days = {};
if ~isfield(rat, 'CS_times') || isempty(rat.CS_times)
    return;
end

fields = fieldnames(rat.CS_times);
isExtinction = ~cellfun('isempty', regexpi(fields, '^CS_(exinction|extinction)_'));
dateTokens = cellfun(@extractDateToken, fields(isExtinction), 'UniformOutput', false);
dateTokens = dateTokens(~cellfun('isempty', dateTokens));

if isempty(dateTokens)
    dateTokens = cellfun(@extractDateToken, fields, 'UniformOutput', false);
    dateTokens = dateTokens(~cellfun('isempty', dateTokens));
end

days = unique(dateTokens, 'stable');
if isempty(days)
    return;
end

dayNums = datenum(strrep(days, '_', '-'), 'yyyy-mm-dd');
[~, order] = sort(dayNums);
days = days(order);
days = days(max(1, numel(days)-nDays+1):end);
end

function token = extractDateToken(fieldName)
token = regexp(fieldName, '\d{4}[_-]\d{2}[_-]\d{2}', 'match', 'once');
if ~isempty(token)
    token = strrep(token, '-', '_');
end
end

function trials = dropFirstTrials(trials, nDrop)
if numel(trials) <= nDrop
    trials = trials([]);
else
    trials = trials(nDrop+1:end);
end
end

function fieldName = findFieldForDate(S, dayStr, prefixes)
fields = fieldnames(S);
for i = 1:numel(prefixes)
    candidate = [prefixes{i} dayStr];
    if isfield(S, candidate)
        fieldName = candidate;
        return;
    end
end

datePattern = strrep(dayStr, '_', '[_-]');
matches = fields(~cellfun('isempty', regexp(fields, datePattern, 'once')));
if isempty(matches)
    error('Could not find a field for date %s.', dayStr);
elseif numel(matches) > 1
    prefixMatches = false(size(matches));
    for i = 1:numel(prefixes)
        prefixMatches = prefixMatches | startsWith(matches, prefixes{i});
    end
    matches = matches(prefixMatches);
end

if isempty(matches)
    error('Could not find a field for date %s.', dayStr);
elseif numel(matches) > 1
    error('Found multiple fields for date %s: %s', dayStr, strjoin(matches, ', '));
end

fieldName = matches{1};
end
