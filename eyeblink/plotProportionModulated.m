function change = plotProportionModulated(varargin)
% plotProportionModulated  Fraction of neurons modulated in trace vs null,
% with observed-vs-shuffled histograms for any metric.
%
% Panels: one per rat + pooled.
%
% Options:
%   'Metric' : 'logfold' | 'fold'  (default)  | 'symdiff'
%              'logfold'  -> mean(FRt)/FRr (log x-axis)
%              'fold'     -> mean(FRt)/FRr (linear x-axis)
%              'symdiff'  -> (mean(FRt)-FRr)/(mean(FRt)+FRr)
%   'Eps'    : numeric epsilon added to both rates (default 0)
%   'Alpha'  : shuffle alpha (still computed; not used in plots) (default 0.05)
%   'NPerm'  : # shuffles per cell (default 500)
%   'Mode'   : 'allNonTrial' (default) | 'preTrial'
%              'allNonTrial' -> trial vs all non-trial baseline,
%                               shuffle: random windows vs session baseline
%              'preTrial'    -> trial vs 2 s pre-CS baseline,
%                               shuffle: post vs pre around pseudo-CS
%   'SaveMod': true/false (default false)
%              If true, save per-rat/day p-values as:
%                ratXXXX.mod.mod_YYYY_MM_DD (vector, len = nCells)
%              Value = p-value per cell; NaN = not analyzed.

% ---------- options ----------
p = inputParser;
addParameter(p,'Metric','fold',@(s) any(strcmpi(s,{'logfold','fold','symdiff'})));
addParameter(p,'Eps',0,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Alpha',0.05,@(x) isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'NPerm',500,@(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'Mode','allNonTrial',@(s) any(strcmpi(s,{'allNonTrial','preTrial'})));
addParameter(p,'SaveMod',true,@(x) (islogical(x) || isnumeric(x)));
%addParameter(p,'Tail','two-sided',@(s) any(strcmpi(s,{'two-sided','right','left'})));
addParameter(p,'Tail','right',@(s) any(strcmpi(s,{'two-sided','right','left'})));


parse(p,varargin{:});
Metric  = lower(p.Results.Metric);
Eps     = p.Results.Eps;
alpha   = p.Results.Alpha;
nPerm   = p.Results.NPerm;
Mode    = lower(p.Results.Mode);
SaveMod = logical(p.Results.SaveMod);
Tail = lower(p.Results.Tail);

% ---------- config ----------
ratNames   = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
win        = [0 2];   % trace window relative to CS (s)
minSpikes  = 0;       % optional per-cell inclusion threshold (currently unused)

nRats  = numel(ratNames);
change = NaN(nRats,3000);   % store metric per cell per rat (size generous)
fracRat = nan(nRats,1);

% Directional stats (still computed; not used in plots)
incFracRat = nan(nRats,1);
decFracRat = nan(nRats,1);
nValidRat  = zeros(nRats,1);

% ---------- figure ----------
fig = figure('Color','w','Position',[100 100 1200 600]);
tiledlayout(fig,2,3,'Padding','compact','TileSpacing','compact');

colObs   = [0.20 0.60 1.00];  % observed
colShuf  = [0.70 0.70 0.70];  % shuffled (null)

% ---------- helper: robust day TS fetch ----------
    function [T0,T1,hadTS] = fetch_day_ts(rat, dayStr, spkDay)
        hadTS = false; T0 = NaN; T1 = NaN;
        keys = {sprintf('Ca_ts_%s',dayStr), sprintf('CA_ts_%s',dayStr), ...
                sprintf('ts_%s',dayStr),      sprintf('TS_%s',dayStr)};
        if isfield(rat,'Ca_ts') && isstruct(rat.Ca_ts)
            for k = 1:numel(keys)
                if isfield(rat.Ca_ts, keys{k})
                    ts = rat.Ca_ts.(keys{k});
                    if ~isempty(ts) && numel(ts)>=2
                        T0 = ts(1); T1 = ts(end);
                        hadTS = true; return;
                    end
                end
            end
        end
        for k = 1:numel(keys)
            if isfield(rat, keys{k})
                ts = rat.(keys{k});
                if ~isempty(ts) && numel(ts)>=2
                    T0 = ts(1); T1 = ts(end);
                    hadTS = true; return;
                end
            end
        end
        % fallback to min/max spike times for the day
        if ~isempty(spkDay)
            allTimes = spkDay(:);
            allTimes = allTimes(~isnan(allTimes)&allTimes>0);
            if ~isempty(allTimes)
                T0 = min(allTimes); T1 = max(allTimes);
                hadTS = false;
            end
        end
    end

% ---------- helper to gather cells for one rat or pooled ----------
    function [cells, meta] = gather_cells(ratIdx)
        cells = struct('st',{},'cs',{},'t0',{},'t1',{});
        meta  = struct('ratIdx',{},'dayStr',{},'cellIdx',{});

        if ratIdx <= nRats
            ratsToDo = ratIdx;
        else
            ratsToDo = 1:nRats; % pooled
        end

        for rr = ratsToDo
            rat   = evalin('base', ratNames{rr});
            dates = autoDateList(rat);
            idx   = find(strcmp(dates, rat.An),1);
            if isempty(idx) || idx < 3
                warning('%s: could not center on An or not enough days; using last 3 days.', ratNames{rr});
                idx = numel(dates);
            end
            days  = dates(max(1,idx-2):idx);

            for d = 1:numel(days)
                dayStr   = days{d};
                spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',dayStr));
                csTimes  = rat.CS_times.(sprintf('CS_%s',dayStr));
                ratemask = rat.ratemask.(sprintf('ratemask_%s',dayStr));

                [T0_day, T1_day, ~] = fetch_day_ts(rat, dayStr, spk);

                [nCells,~] = size(spk);
                for c = 1:nCells
                    if ratemask(c)==0, continue, end
                    times = spk(c,:);
                    times = times(~isnan(times) & times>0);
                    if isempty(times), continue, end

                    cells(end+1).st = times(:);
                    cells(end).cs   = csTimes(:);
                    cells(end).t0   = T0_day;
                    cells(end).t1   = T1_day;

                    meta(end+1).ratIdx = rr;
                    meta(end).dayStr   = dayStr;
                    meta(end).cellIdx  = c;
                end
            end
        end
    end

% ---------- loop panels: 1..nRats then pooled ----------
for panelIdx = 1:(nRats+1)
    [cells, meta] = gather_cells(panelIdx);

    % collect all shuffles for this panel
    nullCat = [];

    % ---- compute metric + shuffle per cell ----
    nCells = numel(cells);
    obsM   = nan(nCells,1);
    hCell  = false(nCells,1);   % kept for summary printing
    incCell = false(nCells,1);  % used only for symdiff summaries
    decCell = false(nCells,1);

    % NEW: store p-values per cell in this panel (for SaveMod)
    pVals = nan(nCells,1);

    for i = 1:nCells
        st = cells(i).st;
        cs = cells(i).cs;
        t0_rec = cells(i).t0;
        t1_rec = cells(i).t1;

        if isempty(st) || isempty(cs) || ~isfinite(t0_rec) || ~isfinite(t1_rec)
            continue
        end

        % CS with full window inside recording span
        dur = diff(win);
        validCS = cs(cs+win(1) >= t0_rec & cs+win(2) <= t1_rec);
        nT = numel(validCS);
        if nT<1, continue, end

        % trial FRs (post-CS window defined by win)
        FRt = nan(nT,1);
        for t = 1:nT
            a = validCS(t) + win(1);
            b = a + dur;
            FRt(t) = sum(st>=a & st<b)/dur;
        end
        FRt_mean = mean(FRt,'omitnan');

        % ---------- baseline FR (observed) ----------
        switch Mode
            case 'allnontrial'
                % trial vs global non-trial baseline
                inWin = false(size(st));
                for t = 1:nT
                    a = validCS(t) + win(1);
                    b = a + dur;
                    inWin = inWin | (st>=a & st<b);
                end
                nSpk_non = sum(~inWin);
                totalTime = (t1_rec - t0_rec);
                totalNon  = totalTime - nT*dur;
                if totalNon <= 0
                    FRr = NaN;
                else
                    FRr = nSpk_non / totalNon;
                end

            case 'pretrial'
                % trial vs 2-sec pre-CS baseline
                preDur = 2;   % seconds
                FRr_vec = nan(nT,1);

                for t = 1:nT
                    a = validCS(t) - preDur;   % window start
                    b = validCS(t);            % window end

                    % must be inside recording
                    if a < t0_rec || b > t1_rec
                        continue
                    end
                    FRr_vec(t) = sum(st>=a & st<b) / preDur;
                end

                % mean baseline across all valid pre-CS windows
                FRr = mean(FRr_vec,'omitnan');
        end

        if ~isfinite(FRr)
            continue
        end

        % observed metric for this cell
        obsM(i) = compute_metric(FRt_mean, FRr, Metric, Eps);
        if panelIdx<=nRats
            change(panelIdx,i) = obsM(i);
        end

        % ---------- null via shuffled starts ----------
        nullM = nan(nPerm,1);

        switch Mode
            % ---------------------- Mode: allNonTrial ----------------------
            case 'allnontrial'
                % session-baseline null: random windows vs session complement
                if (t1_rec - t0_rec) <= dur
                    % too short; leave nullM as NaN
                else
                    tStarts = linspace(t0_rec, t1_rec - dur, 500);
                    for ip = 1:nPerm
                        samp = randsample(tStarts, nT, true);

                        % shuffled trial FR
                        FRs = arrayfun(@(s) sum(st >= s & st < s+dur) / dur, samp);
                        t_shuf = mean(FRs,'omitnan');

                        % shuffled non-trial FR = complement of shuffled windows
                        inWin_shuf = false(size(st));
                        for tt = 1:nT
                            a = samp(tt); b = a + dur;
                            inWin_shuf = inWin_shuf | (st >= a & st < b);
                        end
                        nSpk_non_shuf = sum(st >= t0_rec & st < t1_rec & ~inWin_shuf);
                        totalNon_shuf = (t1_rec - t0_rec) - nT*dur;
                        if totalNon_shuf <= 0
                            b_shuf = NaN;
                        else
                            b_shuf = nSpk_non_shuf / max(totalNon_shuf, eps);
                        end

                        nullM(ip) = compute_metric(t_shuf, b_shuf, Metric, Eps);
                    end
                end

            % ---------------------- Mode: preTrial ------------------------
            case 'pretrial'
                % prePostMatched null: post vs pre around pseudo-CS times
                preDur = 2;   % seconds (match observed)
                if (t1_rec - t0_rec) <= (preDur + dur)
                    % not enough span to place both pre and post windows
                else
                    % choose anchors so [cs_s-preDur, cs_s+dur] is inside recording
                    tStarts = linspace(t0_rec + preDur, t1_rec - dur, 500);
                    for ip = 1:nPerm
                        samp = randsample(tStarts, nT, true);

                        FRt_shuf = nan(nT,1);
                        FRr_shuf = nan(nT,1);

                        for tt = 1:nT
                            cs_s   = samp(tt);

                            % post (trial) window, defined by win
                            a_post = cs_s + win(1);
                            b_post = a_post + dur;

                            % pre window (2 s before pseudo-CS)
                            a_pre  = cs_s - preDur;
                            b_pre  = cs_s;

                            % these should already be inside [t0_rec, t1_rec]
                            if a_pre < t0_rec || b_post > t1_rec
                                continue
                            end

                            FRt_shuf(tt) = sum(st >= a_post & st < b_post) / dur;
                            FRr_shuf(tt) = sum(st >= a_pre  & st < b_pre ) / preDur;
                        end

                        t_bar = mean(FRt_shuf,'omitnan');
                        b_bar = mean(FRr_shuf,'omitnan');

                        nullM(ip) = compute_metric(t_bar, b_bar, Metric, Eps);
                    end
                end
        end

        % append this cell's null to the panel's pool
        nullCat = [nullCat; nullM(:)];

        % significance (kept for printed summaries only)
        nullUse = nullM(isfinite(nullM));
obsVal  = obsM(i);

if isempty(nullUse) || ~isfinite(obsVal)
    pVal = NaN;
    p_right = NaN;
    p_left  = NaN;
else

    switch Metric
        case 'symdiff'
            % centered at 0
            p_right = mean(nullUse >= obsVal);
            p_left  = mean(nullUse <= obsVal);
            p_two   = mean(abs(nullUse) >= abs(obsVal));

        otherwise
            % fold metrics (centered at 1, use log symmetry)
            if obsVal <= 0
                p_two = NaN; p_right = NaN; p_left = NaN;
            else
                obsL  = log(obsVal);
                nullL = log(nullUse);

                p_right = mean(nullUse >= obsVal);
                p_left  = mean(nullUse <= obsVal);
                p_two   = mean(abs(nullL) >= abs(obsL));
            end
    end

    % choose which tail to use
    switch Tail
        case 'two-sided'
            pVal = p_two;
        case 'right'
            pVal = p_right;
        case 'left'
            pVal = p_left;
    end
end

hCell(i)   = (pVal < alpha);
incCell(i) = (p_right < alpha);
decCell(i) = (p_left  < alpha);
pVals(i)   = pVal;
    end

    % ---- per-rat roll-up for summaries ----
    if panelIdx <= nRats
        validMask = isfinite(obsM);
        nValidRat(panelIdx) = nnz(validMask);
        if nValidRat(panelIdx) > 0
            fracRat(panelIdx)   = mean(hCell(validMask),'omitnan');
            incFracRat(panelIdx)= mean(incCell(validMask),'omitnan');
            decFracRat(panelIdx)= mean(decCell(validMask),'omitnan');
        end
    end

    % ---- SAVE per-rat/day p-values into rat.mod if requested ----
  if SaveMod && panelIdx <= nRats && nCells>0
      rat = evalin('base', ratNames{panelIdx});
      if ~isfield(rat,'mod') || ~isstruct(rat.mod)
          rat.mod = struct();
      end

      % group meta entries by day
      dayList = unique({meta.dayStr});
      for dIdx = 1:numel(dayList)
          dayStr = dayList{dIdx};
          maskDay = strcmp({meta.dayStr}, dayStr);
          if ~any(maskDay), continue, end

          % true number of cells that day from CA_peaks
          spkDay    = rat.Ca_peaks.(sprintf('CA_peaks_%s', dayStr));
          nCellsDay = size(spkDay,1);

          % default NaN for all cells (not examined)
          vec = nan(nCellsDay,1);

          % fill in p-values for examined cells
          inds = find(maskDay);
          for kIdx = 1:numel(inds)
              ii   = inds(kIdx);          % index into cells/pVals/meta
              cIdx = meta(ii).cellIdx;    % cell index in that day
              pVal = pVals(ii);
              if isfinite(pVal)
                  vec(cIdx) = pVal;       % p-value per cell
              end
          end

          fld = sprintf('mod_%s', dayStr);
          rat.mod.(fld) = vec;
      end

      assignin('base', ratNames{panelIdx}, rat);
  end
    % ---- plot into this panel: Observed vs Shuffled ----
    nexttile; hold on;

    % Use exactly what was saved for rats; pooled uses current obsM
    if panelIdx <= nRats
        dataAll = change(panelIdx,:).';
    else
        dataAll = obsM(:);
    end

    switch Metric
        case 'logfold'  % ratio; log x-axis; require positive values
            obsUse  = dataAll(isfinite(dataAll) & dataAll>0);
            shfUse  = nullCat(isfinite(nullCat) & nullCat>0);
            if isempty(obsUse) || isempty(shfUse)
                text(0.5,0.5,'No positive values to plot','HorizontalAlignment','center');
            else
                minF = max(min([obsUse(:); shfUse(:)]), realmin);
                maxF = max([obsUse(:); shfUse(:)]);
                nBins = 24;
                edges = logspace(log10(minF), log10(maxF), nBins+1);

                histogram(shfUse, 'BinEdges', edges, 'Normalization','probability', ...
                          'FaceColor',colShuf); hold on
                histogram(obsUse, 'BinEdges', edges, 'Normalization','probability', ...
                          'FaceColor',colObs, 'FaceAlpha',0.75);

                set(gca,'XScale','log');
                xline(1,'k-');
                xlabel('Fold-change (trial / nontrial)');
                ylabel('Probability');
                legend('shuffled','observed','Location','best');
            end

        case 'fold'     % ratio; linear x-axis; include zeros
            obsUse  = dataAll(isfinite(dataAll) & dataAll>=0);
            shfUse  = nullCat(isfinite(nullCat) & nullCat>=0);
            if isempty(obsUse) || isempty(shfUse)
                text(0.5,0.5,'No finite nonnegative values','HorizontalAlignment','center');
            else

                hi = max(3, prctile([obsUse(:); shfUse(:)], 99.9));
                lo = 0;
                nBins = round(hi.*5);
                edges = linspace(lo, hi, nBins+1);
                % exact edge at 1 so <1 is guaranteed to have its own bins
                if lo < 1 && 1 < hi
                    [~,k] = min(abs(edges - 1));
                    edges(k) = 1;
                end

                histogram(shfUse, 'BinEdges', edges, 'Normalization','probability', ...
                          'FaceColor',colShuf); hold on
                histogram(obsUse, 'BinEdges', edges, 'Normalization','probability', ...
                          'FaceColor',colObs, 'FaceAlpha',0.75);

                xline(1,'k-');
                xlim([lo hi]);
                xlabel('Fold-change (trial / nontrial)');
                ylabel('Probability');
                legend('shuffled','observed','Location','best');
            end

        case 'symdiff'  % symmetric difference in [-1,1]
            obsUse  = dataAll(isfinite(dataAll));
            shfUse  = nullCat(isfinite(nullCat));
            edges = linspace(-1,2.5,41);

            histogram(shfUse, 'BinEdges', edges, 'Normalization','probability', ...
                      'FaceColor',colShuf); hold on
            histogram(obsUse, 'BinEdges', edges, 'Normalization','probability', ...
                      'FaceColor',colObs, 'FaceAlpha',0.75);

            xline(0,'k-');
            xlim([-1 1]);
            xlabel('Symmetric difference (t-nt)/(t+nt)');
            ylabel('Probability');
            legend('shuffled','observed','Location','best');
    end

end

% ---- print summaries ----
fprintf('\n=== Fraction modulated by rat (%s, Mode=%s) ===\n', Metric, Mode);
for r = 1:nRats
    if nValidRat(r) > 0
        fprintf('%s: n=%d | mod=%.1f%%\n', ratNames{r}, nValidRat(r), 100*fracRat(r));
    else
        fprintf('%s: n=0 | mod=NA\n', ratNames{r});
    end
end
fracAll = mean(fracRat,'omitnan');
fprintf('All rats combined (mean of rats): %.1f%% modulated\n', 100*fracAll);

% Directional split for symdiff (printed only)
if strcmp(Metric,'symdiff')
    fprintf('\n=== Percent modulated by direction (symdiff; per rat) ===\n');
    for r = 1:nRats
        if nValidRat(r) > 0
            fprintf('%s: n=%d | %%↑=%.1f%% | %%↓=%.1f%% | %%any=%.1f%%\n', ...
                ratNames{r}, nValidRat(r), 100*incFracRat(r), 100*decFracRat(r), 100*fracRat(r));
        else
            fprintf('%s: n=0\n', ratNames{r});
        end
    end
    % pooled, weighted by # valid cells per rat
    w = nValidRat; w(~isfinite(incFracRat)|~isfinite(decFracRat)) = 0;
    incPooled = nansum(w .* incFracRat) / max(nansum(w),1);
    decPooled = nansum(w .* decFracRat) / max(nansum(w),1);
    anyPooled = nansum(w .* fracRat   ) / max(nansum(w),1);
    fprintf('ALL (weighted by n): %%↑=%.1f%% | %%↓=%.1f%% | %%any=%.1f%% | n_total=%d\n', ...
        100*incPooled, 100*decPooled, 100*anyPooled, nansum(nValidRat));
end

end

% ================= helper =================
function m = compute_metric(FRt_mean, FRr, Metric, Eps)
% compute_metric  Apply chosen metric with optional epsilon for stability.
t = FRt_mean + Eps;
b = FRr      + Eps;
switch Metric
    case {'logfold','fold'}
        m = t ./ b;
    case 'symdiff'
        denom = t + b;
        if denom==0
            m = NaN;
        else
            m = (t - b) ./ denom;
        end
    otherwise
        error('Unknown Metric: %s', Metric);
end
end
