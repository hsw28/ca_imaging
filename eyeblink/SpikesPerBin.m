function SpikesPerBin(ratNames, varargin)
% SpikesPerBin(ratNames, 'Demean',false, 'BaselineWindow',[-0.5 0], 'StoreBoth',true)
% Build a (nBins × nCells) event-rate matrix per day and store in rat.spikesperbin.
% Optionally demean each cell by its pre-CS baseline mean (Hz).
%
% Example:
%   SpikesPerBin({'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%                'Demean',true, 'BaselineWindow',[-0.5 0])

% ---------- options ----------
p = inputParser;
p.addParameter('Demean', true, @(x)islogical(x) && isscalar(x));
p.addParameter('BaselineWindow', [-0.5 0], @(v)isnumeric(v) && numel(v)==2 && v(1)<v(2));
p.addParameter('StoreBoth', true, @(x)islogical(x) && isscalar(x));   % if true, save both raw (cpb_) and de-meaned (cpbDM_)
p.parse(varargin{:});
Demean          = p.Results.Demean;
BaselineWindow  = p.Results.BaselineWindow;
StoreBoth       = p.Results.StoreBoth;

% ---------- analysis window & binning ----------
endtime = 4;                   % you can change this; plot shades still mark 0–2 s
winFull  = [-.5 endtime];      % analysis window (s)
binWidth = (1/15);             % 0.06666.. s (15 Hz sampling-aligned bins)
binEdges = winFull(1):binWidth:winFull(2);
nBins    = numel(binEdges)-1;
tAxis    = binEdges(1:end-1) + binWidth/2;   % bin centers

% indices used for baseline (for demeaning)
baseMask = (tAxis >= BaselineWindow(1)) & (tAxis < BaselineWindow(2));

% ---------- main loop ----------
for r = 1:numel(ratNames)
    rat   = evalin('base', ratNames{r});          % pull from workspace
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);        % last analyzed day
    days  = dates(idx-2:idx);                     % last 3 analyzed days

    for d = 1:3
        dayStr   = days{d};
        spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',dayStr));   % nCells × T (timestamps of events per cell)
        csTimes  = rat.CS_times.(sprintf('CS_%s',dayStr));         % 1 × nTrials (CS onsets)
        ratemask = (rat.ratemask.(sprintf('ratemask_%s',dayStr)) == 1);  % 1 × nCells keep-mask

        nTotal   = size(spk,1);                       % total cells this day
        inclMask = ratemask(:) == 1;                  % logical (nCells × 1)
        inclIdx  = find(inclMask);

        % accumulate counts across trials, then convert to Hz
        counts = zeros(nBins, nTotal);                % start at 0
        for t = 1:numel(csTimes)
            t0 = csTimes(t);
            for jj = 1:numel(inclIdx)
                c  = inclIdx(jj);
                st = spk(c,:);
                st = st(~isnan(st) & st>0);          % valid timestamps
                rel = st - t0;                        % relative to CS
                rel = rel(rel >= binEdges(1) & rel < binEdges(end));
                counts(:,c) = counts(:,c) + histcounts(rel, binEdges).';
            end
        end
        ratesHz = counts / (binWidth * max(1, numel(csTimes)));   % per-bin event rate (Hz)
        ratesHz(:, ~inclMask) = NaN;                               % blank excluded cells

        % ----- optional demeaning (per cell) by pre-CS baseline -----
        if Demean
            if ~any(baseMask)
                warning('BaselineWindow [%.3f %.3f) does not overlap bins; skipping demeaning for %s %s.', ...
                        BaselineWindow(1), BaselineWindow(2), ratNames{r}, dayStr);
                ratesDM = ratesHz;  % no-op
            else
                baseMeans = nanmean(ratesHz(baseMask, :), 1);      % 1 × nCells
                ratesDM   = ratesHz - baseMeans;                   % broadcast
            end
        end

        % ----- store back into the rat struct -----
        if ~isfield(rat,'spikesperbin'), rat.spikesperbin = struct(); end

        if StoreBoth
            rat.spikesperbin.(sprintf('cpb_%s',dayStr)) = ratesHz.';     % (nCells × nBins)
            if Demean
                rat.spikesperbin.(sprintf('cpbDM_%s',dayStr)) = ratesDM.'; % demeaned
            end
        else
            % overwrite single product depending on Demean
            if Demean
                rat.spikesperbin.(sprintf('cpb_%s',dayStr)) = ratesDM.';  % store de-meaned as cpb_
            else
                rat.spikesperbin.(sprintf('cpb_%s',dayStr)) = ratesHz.';  % raw
            end
        end
    end

    % push updated struct back to base workspace
    assignin('base', ratNames{r}, rat);
    fprintf('Stored spikesperbin for %s (days: %s, %s, %s)\n', ratNames{r}, days{:});
end

% ---------- figure: mean ± SEM per rat + grand panel ----------
nRats = numel(ratNames);
figure('Color','w','Position',[200 200 1200 300]);

grandMat = [];  % collect per-rat means for final panel

for rr = 1:nRats
    rat   = evalin('base', ratNames{rr});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);

    % pull matrices (prefer demeaned if requested and stored)
    M = [];
    for d = 1:3
        if Demean && isfield(rat.spikesperbin, sprintf('cpbDM_%s',days{d}))
            X = rat.spikesperbin.(sprintf('cpbDM_%s',days{d}));  % nCells × nBins
        else
            X = rat.spikesperbin.(sprintf('cpb_%s',days{d}));
        end
        M = [M; X]; %#ok<AGROW>
    end

    mu  = nanmean(M,1);                           % 1 × nBins
    sem = nanstd(M,0,1) ./ sqrt(sum(~isnan(M),1));

    subplot(1,nRats+1,rr); hold on;

    yTop = max(0, max(mu+sem)*1.1);
    fill([0 2 2 0], [0 0 yTop yTop], [.9 .9 .9], 'EdgeColor','none');  % trace shading

    bar(tAxis, mu, 1, 'FaceColor',[.5 .5 .8], 'EdgeColor','none');
    errorbar(tAxis, mu, sem, '.k', 'CapSize',4);

    yl = ylim; line([0 0], yl, 'Color','r','LineStyle','--');         % CS
    line([.75 .75], yl, 'Color','b','LineStyle','--');                % US (adjust if needed)

    xlabel('Time from CS (s)');
    if Demean
        ylabel('\Delta rate (Hz from baseline)');
    else
        ylabel('Event rate (Hz)');
    end
    title(ratNames{rr}); box off;

    grandMat = [grandMat; mu]; %#ok<AGROW>
end

% ----------- combined panel ------------------------------------------
grandMu  = nanmean(grandMat,1);
grandSem = nanstd (grandMat,0,1) ./ sqrt(sum(~isnan(grandMat),1));

subplot(1,nRats+1,nRats+1); hold on;
yTop = max(0, max(grandMu+grandSem)*1.1);
fill([0 2 2 0], [0 0 yTop yTop], [.9 .9 .9], 'EdgeColor','none');
bar(tAxis, grandMu, 1, 'FaceColor',[.8 .5 .5], 'EdgeColor','none');
errorbar(tAxis, grandMu, grandSem, '.k', 'CapSize',4);
yl = ylim; line([0 0], yl, 'Color','r','LineStyle','--');
line([.75 .75], yl, 'Color','b','LineStyle','--');  % US marker

xlabel('Time from CS (s)');
if Demean
    ylabel('\Delta rate (Hz from baseline)');
else
    ylabel('Event rate (Hz)');
end
title('All rats'); box off;

sgtitle(sprintf('Mean event rate per %.3f-s bin  (%.2f → %.2f s, CS at 0)', ...
        binWidth, winFull(1), winFull(2)));
end
