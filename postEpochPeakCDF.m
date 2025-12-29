function OUT = postEpochPeakCDF(ratNames, varargin)
% postEpochPeakCDF
%   Compute and plot CDF of PSTH peak times for post-epoch–modulated cells.
%
% For each rat (last up-to-3 pre-An days):
%   1) Load Ca peaks, position, CS times, ratemask, epoch labels.
%   2) Restrict to kept cells (ratemask==1).
%   3) Select post-modulated cells via epoch label (e.g., epoch==4).
%   4) Build continuous FR(t) for each cell from Ca peaks.
%   5) For each post cell, compute trial-averaged PSTH in postWin and
%      take the time of the PSTH maximum as the cell's post peak time
%      (relative to the start of the post window).
%   6) Pool peak times across rats/days and plot the CDF.
%
% Example:
%   OUT = postEpochPeakCDF({'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%                          'postWin',[0.85 2.0], 'PostLabel',4);

% ---------- options ----------
p = inputParser;
p.addParameter('postWin',[0.85 2.0]);     % post window relative to CS (s)
p.addParameter('EpochVar','epoch',@(s)ischar(s)||isstring(s)); % field with epoch labels
p.addParameter('PostLabel',4,@(x)isnumeric(x)&&isscalar(x));   % label for post epoch
p.addParameter('SmoothBins',3,@(x)isnumeric(x)&&isscalar(x));  % Gaussian smoothing width (bins)
p.addParameter('DoPlot',true,@islogical);

p.parse(varargin{:});
opt = p.Results;

if nargin < 1 || isempty(ratNames)
    error('Specify ratNames cell array.');
end

nR = numel(ratNames);

allPeakTimes = [];          % pooled across rats/days
perRat       = struct([]);  % per-rat details
ratPeakTimes = cell(nR,1);  % per-rat pooled peak times

for r = 1:nR
    ratName = ratNames{r};
    rat     = evalin('base', ratName);

    dates = autoDateList(rat);
    iAn   = find(strcmp(dates, rat.An), 1);
    if isempty(iAn) || iAn < 1
        warning('%s: no An date or too few days; skipping.', ratName);
        continue;
    end

    % last up-to-3 days before/including An
    dayIdx = max(1, iAn-2):iAn;
    days   = dates(dayIdx);

    perRat(r).ratName   = ratName;
    perRat(r).days      = days;
    perRat(r).peakTimes = cell(numel(days),1);

    pt_r = [];  % collector for this rat across days

    for d = 1:numel(days)
        D = days{d};

        % --- spikes & ratemask ---
        spkAll = rat.Ca_peaks.(sprintf('CA_peaks_%s', D));        % [nCells_all x nSpikes?]
        mask   = rat.ratemask.(sprintf('ratemask_%s', D))==1;     % logical [nCells_all x 1]

        if isempty(spkAll) || ~any(mask)
            warning('%s %s: no spikes or no kept cells; skipping day.', ratName, D);
            continue;
        end

        % --- epoch labels (to pick post-modulated cells) ---
        if ~isfield(rat, opt.EpochVar)
            error('%s: missing field %s for epoch labels.', ratName, opt.EpochVar);
        end
        GV = rat.(opt.EpochVar);
        if ~isstruct(GV)
            error('%s.%s must be a struct with per-day fields.', ratName, opt.EpochVar);
        end

        fns   = fieldnames(GV);
        idxFn = contains(fns, D);
        if ~any(idxFn)
            error('%s: no epoch field matches date %s.', ratName, D);
        end
        if sum(idxFn) > 1
            warning('%s: multiple epoch fields match %s; using first.', ratName, D);
        end
        fname     = fns{find(idxFn,1,'first')};
        epochFull = GV.(fname);
        epochFull = epochFull(:);

        nCellsAll = size(spkAll,1);
        if numel(epochFull) < nCellsAll
            error('%s %s: epoch labels (%d) < cells (%d).', ratName, D, numel(epochFull), nCellsAll);
        end

        % restrict to kept cells
        spk   = spkAll(mask,:);
        epoch = epochFull(mask);

        % logical selection of post-modulated cells
        postCell = (epoch == opt.PostLabel);
        if ~any(postCell)
            fprintf('%s %s: no post-modulated cells.\n', ratName, D);
            continue;
        end

        % --- build FR matrix over full session for kept cells ---
        pos = rat.pos.(sprintf('pos_%s', D));
        t   = pos(:,1);
        FR  = spikesToRate_local(spk, t);   % [time x nCells_kept]

        csTimes = rat.CS_times.(sprintf('CS_%s', D));
        if isempty(csTimes)
            warning('%s %s: no CS times; skipping day.', ratName, D);
            continue;
        end

        % subset to post cells
        FR_post = FR(:, postCell);

        % --- compute PSTH-based peak times for post cells (0..postDur) ---
        peakTimes_day = computePostPeakTimes_local(FR_post, t, csTimes, opt.postWin, opt.SmoothBins);

        perRat(r).peakTimes{d} = peakTimes_day(:);
        pt_r        = [pt_r; peakTimes_day(:)]; %#ok<AGROW>
        allPeakTimes = [allPeakTimes; peakTimes_day(:)]; %#ok<AGROW>

        fprintf('%s %s: %d post cells, %d peaks collected.\n', ...
            ratName, D, nnz(postCell), numel(peakTimes_day));
    end

    ratPeakTimes{r} = pt_r;
end

OUT = struct();
OUT.ratNames      = ratNames;
OUT.perRat        = perRat;
OUT.ratPeakTimes  = ratPeakTimes;
OUT.peakTimes_all = allPeakTimes(:);
OUT.postWin       = opt.postWin;
OUT.options       = opt;

% ---------- CDF plot: one line per rat + mean line ----------
if opt.DoPlot && ~isempty(allPeakTimes)
    dur = opt.postWin(2) - opt.postWin(1);     % length of post window (s)

    % rats that actually contributed peaks
    validRats = find(~cellfun(@isempty, ratPeakTimes));
    nValid    = numel(validRats);

    figure('Color','w','Position',[300 300 650 450]); hold on;

    % common x-grid for mean CDF
    xGrid = linspace(0, dur, 200);
    Fgrid = nan(nValid, numel(xGrid));

    cols = lines(nValid);
    leg  = cell(nValid+1,1);
    row  = 0;

    for i = 1:nValid
        rIdx = validRats(i);
        pt   = ratPeakTimes{rIdx};
        pt   = pt(isfinite(pt));   % safety

        if isempty(pt)
            continue;
        end

        row = row + 1;

        [f,x] = ecdf(pt);

        % make sample points unique for interp1
        [xu, ia] = unique(x);
        fu       = f(ia);

        % per-rat line
        plot(xu, fu, 'LineWidth', 1.5, 'Color', cols(row,:));
        leg{row} = ratNames{rIdx};

        % interpolate onto common grid for mean CDF
        if numel(xu) == 1
            % step from that point onward
            Fgrid(row,:) = double(xGrid >= xu);
        else
            Fgrid(row,:) = interp1(xu, fu, xGrid, 'previous', 'extrap');
        end
    end

    % trim unused legend slots if some rats had no peaks
    leg = leg(1:row);

    % mean CDF across rats
    meanF = mean(Fgrid(1:row,:), 1, 'omitnan');
    plot(xGrid, meanF, 'k-', 'LineWidth', 2.5);
    leg{end+1} = 'Mean across rats';

    xlabel('Peak time within post window (s)');
    ylabel('Cumulative fraction of cells');
    title(sprintf('Post-epoch PSTH peak times (n = %d cells total)', numel(allPeakTimes)));
    xlim([0 dur]);
    ylim([0 1]);
    box on;
    legend(leg, 'Location', 'southeast');
end

end % main


% ================= helper functions ===================

function FR = spikesToRate_local(spk, t)
% spk: [nCells x nSpikes] numeric, timestamps >0, NaN for padding
% t:   [nTime x 1] time vector (e.g. from pos)
dt = median(diff(t),'omitnan');
nT = numel(t);
nC = size(spk,1);
FR = zeros(nT,nC);
for c = 1:nC
    st = spk(c,:);
    st = st(~isnan(st) & st>0);
    if isempty(st), continue; end
    idx = discretize(st, t);
    idx = idx(isfinite(idx) & idx>=1 & idx<=nT);
    FR(idx,c) = FR(idx,c) + 1;
end
FR = FR ./ dt;   % Hz
end

function peakTimes = computePostPeakTimes_local(FR, t, csTimes, postWin, smoothBins)
% FR: [time x nCells] for a single day (kept/post cells subset)
% t:  [time x 1]
% csTimes: CS timestamps
% postWin: [a b] relative to CS (e.g. [0.85 2.0])
% smoothBins: Gaussian smoothing width (bins) for PSTH

nCells    = size(FR,2);
peakTimes = nan(nCells,1);

dt  = median(diff(t),'omitnan');
dur = postWin(2) - postWin(1);          % length of post window in s

% define RELATIVE bin edges (0 → dur)
edgesRel   = 0:dt:dur;
if edgesRel(end) < dur
    edgesRel = [edgesRel dur];
end
binCenters = edgesRel(1:end-1) + diff(edgesRel)/2;

for c = 1:nCells
    psth_sum = zeros(1, numel(binCenters));
    ntr = 0;

    for k = 1:numel(csTimes)
        % absolute window on this trial
        t0 = csTimes(k) + postWin(1);
        t1 = csTimes(k) + postWin(2);

        mask = (t >= t0 & t < t1);
        if ~any(mask), continue; end

        fr_seg = FR(mask, c);
        % time RELATIVE to t0 (0 at start of post window)
        t_rel  = t(mask) - t0;          % 0 → dur

        % bin FR into RELATIVE edges
        bIdx = discretize(t_rel, edgesRel);

        for b = 1:numel(binCenters)
            inBin = (bIdx == b);
            if any(inBin)
                psth_sum(b) = psth_sum(b) + mean(fr_seg(inBin),'omitnan');
            end
        end
        ntr = ntr + 1;
    end

    if ntr == 0
        continue;
    end

    psth = psth_sum / ntr;

    if smoothBins > 1
        psth = smoothdata(psth,'gaussian',smoothBins);
    end

    [~,imax] = max(psth);

    % quadratic interpolation around max to de-quantize peak time
    if imax > 1 && imax < numel(psth)
        x3 = binCenters(imax-1:imax+1);
        y3 = psth(imax-1:imax+1);
        p  = polyfit(x3, y3, 2);           % parabola
        t_hat = -p(2)/(2*p(1));            % vertex
        t_hat = max(0, min(dur, t_hat));   % clamp
        peakTimes(c) = t_hat;
    else
        peakTimes(c) = binCenters(imax);
    end
end

peakTimes = peakTimes(~isnan(peakTimes));
end
