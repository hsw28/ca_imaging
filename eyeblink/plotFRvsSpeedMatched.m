function Res = plotFRvsSpeedMatched(compareMode, varargin)
% Compare trial-window FR to speed-matched non-trial FR (paired per cell).
% Returns a results struct; can be called in grid mode with plotting off.
%
% compareMode: 'seconds' | 'bins'
%
% Name/value (common):
%   'ratNames'           cellstr   default {'rat0222','rat0307','rat0313','rat0314','rat0816'}
%   'win'                [1x2]     default [0 2] s
%   'binSize'            scalar    default 1/7.5 s (≈0.133 s)  [used in 'bins' mode]
%   'alpha'              scalar    default 0.05    (FDR across cells)
%   'test'               char      'ttest'|'signrank' (per-cell paired)
%   'MinPairsFrac'       scalar    default 0.2     (min frac trials a cell must have)
%   'AllowReuseControl'  logical   default true    (can re-use non-trial control time)
%   'doPlots'            logical   default true
%
% Grid helpers:
%   'MinDayCells'        scalar    default 0       (early quit if paired N < this)
%   'EarlyQuitOnLowN'    logical   default false
%
% OUTPUT (Res):
%   Res.pctPerDayIncluded  [nRats x 3] %sig among tested cells (FDR)
%   Res.pctPerRatIncluded  [nRats x 1]
%   Res.meanPctIncluded    scalar
%   Res.dayTestsByRat      [nRats x 3] struct: p, tstat, df, nCells, deltaMeanHz
%   Res.perRatPooled       [nRats x 1] struct: p, tstat, df, nCells, deltaMeanHz
%   Res.allRatsPooled      struct:     p, tstat, df, nCells, deltaMeanHz
%   Res.meta               parameters used
%   Res.earlyQuit          true/false when EarlyQuitOnLowN triggers (with Res.lowN)
%
% DEPENDS ON: trialVsSpeedMatched (your per-cell paired comparator)

if nargin<1 || isempty(compareMode), compareMode = 'bins'; end
compareMode = lower(string(compareMode));
assert(ismember(compareMode,["seconds","bins"]), 'compareMode must be ''seconds'' or ''bins''.');

ip = inputParser();
ip.addParameter('ratNames', {'rat0222','rat0307','rat0313','rat0314','rat0816'});
ip.addParameter('win',       [0 2]);
ip.addParameter('binSize',    1/7.5);
ip.addParameter('alpha',      0.05);
ip.addParameter('test',      'ttest');
ip.addParameter('MinPairsFrac',      0.2);
ip.addParameter('AllowReuseControl', true, @islogical);
ip.addParameter('doPlots',    true, @islogical);

% NEW (grid helpers)
ip.addParameter('MinDayCells',       0);
ip.addParameter('EarlyQuitOnLowN',   false, @islogical);

ip.parse(varargin{:});
ratNames  = ip.Results.ratNames;
win       = ip.Results.win;
binSize   = ip.Results.binSize;
alpha     = ip.Results.alpha;
testType  = ip.Results.test;
mpf       = ip.Results.MinPairsFrac;
reuseOK   = ip.Results.AllowReuseControl;
doPlots   = ip.Results.doPlots;
minDayCells = ip.Results.MinDayCells;
earlyQuit   = ip.Results.EarlyQuitOnLowN;

nRats = numel(ratNames);
dayLabels = {'day-2','day-1','day 0'};

% prealloc
pctPerDayIncluded = nan(nRats,3);
pctPerRatIncluded = nan(nRats,1);

% per-(rat,day) paired tests
template = struct('p',NaN,'tstat',NaN,'df',NaN,'nCells',0,'deltaMeanHz',NaN);
dayTestsByRat = repmat(template, nRats, 3);

% per-rat pooling collectors
perRat_trialRates = cell(nRats,1);
perRat_matchRates = cell(nRats,1);

for r = 1:nRats
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An), 1);
    if isempty(idx) || idx<3
        warning('%s missing 3 days ending at An; skipping.', ratNames{r});
        continue;
    end
    days  = dates(idx-2:idx);  % d-2, d-1, d0

    for d=1:3
        day      = days{d};
        spk      = rat.Ca_peaks.(['CA_peaks_' day]);     % cell array
        posFull  = smoothpos(rat.pos.(['pos_' day]));
        ts       = posFull(:,1);
        xy       = posFull(:,2:3); %#ok<NASGU>
        cs       = rat.CS_times.(['CS_' day]);
        mask     = rat.ratemask.(['ratemask_' day]);

        S = trialVsSpeedMatched( spk, ts, xy, cs, mask, ...
              'win',win,'binSize',binSize,'alpha',alpha,'test',testType, ...
              'CompareMode',char(compareMode), ...
              'MinPairsFrac',mpf, ...
              'AllowReuseControl',reuseOK );

        pctPerDayIncluded(r,d) = S.pctSigDiffFDR;

        % build paired vectors for this (rat,day)
        tested = isfinite(S.meanTrialRate) & isfinite(S.meanMatchedRate);
        x = double(S.meanTrialRate(tested));
        y = double(S.meanMatchedRate(tested));
        K = min(numel(x), numel(y));
        if K>0
            x = x(1:K); y = y(1:K);
        end

        % EARLY QUIT?
        if earlyQuit && minDayCells>0 && K < minDayCells
            Res = struct('earlyQuit', true, ...
                         'lowN', struct('rat', ratNames{r}, ...
                                        'dayIdx', d, 'dayLabel', dayLabels{d}, ...
                                        'nCells', K));
            return;
        end

        % console line about tested/sig (optional)
        fprintf('[%s %s] tested=%d, sig=%d (%.1f%%), median p=%.3g\n', ...
            ratNames{r}, day, sum(isfinite(S.pVal)), sum(S.sigDiffFDR), ...
            S.pctSigDiffFDR, median(S.pVal(isfinite(S.pVal)),'omitnan'));

        % per-(rat,day) paired test across cells
        if K >= 2
            [~,p,~,st] = ttest(x, y);
            dayTestsByRat(r,d).p           = p;
            dayTestsByRat(r,d).tstat       = st.tstat;
            dayTestsByRat(r,d).df          = st.df;
            dayTestsByRat(r,d).nCells      = K;
            dayTestsByRat(r,d).deltaMeanHz = mean(x - y,'omitnan');
            fprintf('   [%s %s] n=%d  Δ=%.4f Hz  p=%.3g\n', ...
                ratNames{r}, dayLabels{d}, K, dayTestsByRat(r,d).deltaMeanHz, p);
        else
            dayTestsByRat(r,d).nCells = K;
        end

        % stash for pooled-per-rat later
        perRat_trialRates{r} = [perRat_trialRates{r}; x];
        perRat_matchRates{r} = [perRat_matchRates{r}; y];
    end

    pctPerRatIncluded(r) = mean(pctPerDayIncluded(r,:),'omitnan');
end

% ---------- pooled per rat ----------
template2 = struct('p',NaN,'tstat',NaN,'df',NaN,'nCells',0,'deltaMeanHz',NaN);
perRatPooled = repmat(template2, nRats, 1);
for r=1:nRats
    xr = perRat_trialRates{r}; yr = perRat_matchRates{r};
    K  = min(numel(xr), numel(yr));
    if K >= 2
        xr = xr(1:K); yr = yr(1:K);
        mask = isfinite(xr) & isfinite(yr);
        xr = xr(mask); yr = yr(mask);
        if numel(xr) >= 2
            [~,p,~,st] = ttest(xr, yr);
            perRatPooled(r).p           = p;
            perRatPooled(r).tstat       = st.tstat;
            perRatPooled(r).df          = st.df;
            perRatPooled(r).nCells      = numel(xr);
            perRatPooled(r).deltaMeanHz = mean(xr - yr,'omitnan');
        end
    end
end

% ---------- pooled across all rats & days ----------
x_all = vertcat(perRat_trialRates{:});
y_all = vertcat(perRat_matchRates{:});
Kall  = min(numel(x_all), numel(y_all));
allRatsPooled = template2;
if Kall >= 2
    x_all = x_all(1:Kall); y_all = y_all(1:Kall);
    mask  = isfinite(x_all) & isfinite(y_all);
    x_all = x_all(mask);    y_all = y_all(mask);
    if numel(x_all) >= 2
        [~,p,~,st] = ttest(x_all, y_all);
        allRatsPooled.p           = p;
        allRatsPooled.tstat       = st.tstat;
        allRatsPooled.df          = st.df;
        allRatsPooled.nCells      = numel(x_all);
        allRatsPooled.deltaMeanHz = mean(x_all - y_all,'omitnan');
    end
end

% ---------- pack result ----------
Res = struct();
Res.pctPerDayIncluded = pctPerDayIncluded;
Res.pctPerRatIncluded = pctPerRatIncluded;
Res.meanPctIncluded   = mean(pctPerDayIncluded(:),'omitnan');
Res.dayTestsByRat     = dayTestsByRat;
Res.perRatPooled      = perRatPooled;
Res.allRatsPooled     = allRatsPooled;
Res.meta = struct('compareMode',char(compareMode), 'win',win, 'binSize',binSize, ...
                  'alpha',alpha, 'test',testType, 'MinPairsFrac',mpf, ...
                  'AllowReuseControl',reuseOK, 'ratNames',{ratNames}, ...
                  'MinDayCells',minDayCells, 'EarlyQuitOnLowN',earlyQuit);

% (optional) plots — skipped in grid mode
if doPlots
    % your plotting code (unchanged) can go here if you still want it
end
end


% ===== plotting helper (compact) =====
function local_make_plots(ratNames, pctPerDayIncluded, pooledByRat, dayTestsByRat, compareMode)
nRats = numel(ratNames); DL = {'day-2','day-1','day 0'};
figure('Color','w','Position',[120 200 1280 520]);

subplot(1,2,1); hold on;
bar(pctPerDayIncluded,'grouped'); ylim([0 100]);
xticks(1:nRats); xticklabels(ratNames); ylabel('% sig (FDR) among tested');
legend(DL,'Location','northwest');
title(sprintf('Included cells per day (%s)', upper(char(compareMode)))); box on;

subplot(1,2,2); hold on;
vals = arrayfun(@(s) s.deltaMeanHz, pooledByRat);
ns   = arrayfun(@(s) s.nCells,      pooledByRat);
bar(vals,'FaceColor',[.3 .6 .85]); yline(0,'k-'); box on;
xticks(1:nRats); xticklabels(ratNames); ylabel('\Delta rate (Hz)  Trial - Matched');
title('Pooled per rat (paired across cells)');
for r=1:nRats
    text(r, vals(r)+(vals(r)>=0)*0.02-(vals(r)<0)*0.02, sprintf('n=%d',ns(r)), ...
        'HorizontalAlignment','center','Color',[.2 .2 .2]);
end
end

%%===========================================================================%%
function stats = trialVsSpeedMatched(spikeCell, ts, pos, csTimes, ratemask, varargin)
% Per-cell paired test (TRIAL vs speed-matched NON-TRIAL) with FDR across cells.
% CompareMode:
%   'seconds' — contiguous 2 s window outside trials, matched by 2 s mean speed
%   'bins'    — select K non-trial bins (binSize), **greedy** to make the mean of
%               selected bins' mean speeds closest to the trial's mean speed.

ip = inputParser();
ip.addParameter('win',              [0 2]);
ip.addParameter('binSize',           1/7.5);
ip.addParameter('minPairs',          []);
ip.addParameter('alpha',             0.05);
ip.addParameter('test',              'ttest');
ip.addParameter('CompareMode',       'seconds');
ip.addParameter('MinPairsFrac',      0.2);
ip.addParameter('AllowReuseControl', true, @islogical);
ip.parse(varargin{:});

win         = ip.Results.win;
binSize     = ip.Results.binSize;
alpha       = ip.Results.alpha;
testType    = lower(ip.Results.test);
compareMode = lower(ip.Results.CompareMode);
allowReuse  = ip.Results.AllowReuseControl;

% normalize spikes
if ~iscell(spikeCell) && isnumeric(spikeCell)
  spikeCell = mat2cell(spikeCell, ones(size(spikeCell,1),1), size(spikeCell,2));
end

nTrials  = numel(csTimes);
if isempty(ip.Results.minPairs)
    minPairs = max(1, ceil(ip.Results.MinPairsFrac * nTrials));
else
    minPairs = ip.Results.minPairs;
end

% speed & dt
ts   = ts(:);
dts  = diff(ts); dts = dts(isfinite(dts) & dts>0);
dtMed = isempty(dts)*binSize + ~isempty(dts)*median(dts);
dXY  = [diff(pos); [0 0]];
spd  = hypot(dXY(:,1),dXY(:,2)) ./ [diff(ts); dtMed];
spd(~isfinite(spd)) = 0;

% trial vs non-trial masks
inTrial = false(size(ts));
for t=1:nTrials
  inTrial = inTrial | (ts>=csTimes(t)+win(1) & ts<csTimes(t)+win(2));
end
outMask  = ~inTrial;

% per-trial mean speed
edgeMat  = csTimes(:) + win(:).';
winLen   = diff(win);
Kbins    = max(1, round(winLen/binSize));
trialSpd = nan(nTrials,1);
for t=1:nTrials
  idxCS       = ts>=edgeMat(t,1) & ts<edgeMat(t,2);
  trialSpd(t) = mean(spd(idxCS),'omitnan');
end

% ===== prepare candidate control bins =====
L = max(1, round(binSize/dtMed));        % samples per bin
validBinStart = conv(double(outMask), ones(L,1), 'valid') == L;
startIdx_all  = find(validBinStart);
startTime_all = ts(startIdx_all); %#ok<NASGU>
spdBinMean_all = conv(spd, ones(L,1)/L, 'valid');
spdBinMean_all = spdBinMean_all(validBinStart);

% ===== state for no-reuse across whole ts grid =====
if ~allowReuse
    usedFrames = false(numel(ts),1);
end

% outputs
nCells   = numel(spikeCell);
deltaFR  = nan(nCells,1);
pVal     = nan(nCells,1);
isTested = false(nCells,1);
meanTrialRate   = nan(nCells,1);
meanMatchedRate = nan(nCells,1);
deltaRate       = nan(nCells,1);

for c=1:nCells
  if ~isempty(ratemask) && ~ratemask(c), continue; end
  st = spikeCell{c}(:);
  if isempty(st), continue; end

  FRt = nan(nTrials,1);
  FRm = nan(nTrials,1);

  switch compareMode
    case 'bins'
      for t=1:nTrials
          t0 = edgeMat(t,1); t1 = edgeMat(t,2);
          FRt(t) = sum(st>=t0 & st<t1);

          if allowReuse
              candIdx = startIdx_all; candSpd = spdBinMean_all;
          else
              keep = true(numel(startIdx_all),1);
              for ii=1:numel(startIdx_all)
                  si = startIdx_all(ii);
                  if any(usedFrames(si:si+L-1)), keep(ii) = false; end
              end
              candIdx = startIdx_all(keep);
              candSpd = spdBinMean_all(keep);
          end
          if isempty(candIdx), FRm(t) = NaN; continue; end

          % greedy: choose K bins so mean(selSpd) ≈ trialSpd(t)
          sel = []; selSpd = [];
          [~,i1] = min(abs(candSpd - trialSpd(t)));
          sel(1,1)    = candIdx(i1);
          selSpd(1,1) = candSpd(i1);
          remaining = true(numel(candIdx),1); remaining(i1)=false;
          while numel(sel) < Kbins
              bestJ = 0; bestScore = Inf;
              for jj = find(remaining).'
                  mu = mean([selSpd; candSpd(jj)], 'omitnan');
                  sc = abs(mu - trialSpd(t));
                  if sc < bestScore, bestScore = sc; bestJ = jj; end
              end
              if bestJ==0, break; end
              sel(end+1,1)    = candIdx(bestJ);   %#ok<AGROW>
              selSpd(end+1,1) = candSpd(bestJ);   %#ok<AGROW>
              remaining(bestJ)= false;
          end

          if numel(sel) < Kbins
              FRm(t) = NaN;
          else
              cnt = 0;
              for k2=1:Kbins
                  sIdx = sel(k2);
                  s0   = ts(sIdx);
                  cnt  = cnt + sum(st>=s0 & st<s0+binSize);
              end
              FRm(t) = cnt;

              if ~allowReuse
                  for k2=1:Kbins
                      sIdx = sel(k2);
                      usedFrames(sIdx:sIdx+L-1) = true;
                  end
              end
          end
      end

    case 'seconds'
      L2s = max(1, round((win(2)-win(1))/dtMed));
      validStarts = conv(double(outMask), ones(L2s,1), 'valid') == L2s;
      startIdxSec = find(validStarts);
      spdWinMean  = conv(spd, ones(L2s,1)/L2s, 'valid');

      if ~allowReuse
          usedFrames = false(numel(ts),1);
      end

      for t=1:nTrials
          t0 = edgeMat(t,1); t1 = edgeMat(t,2);
          FRt(t) = sum(st>=t0 & st<t1);

          diffs = abs(spdWinMean(validStarts) - trialSpd(t));
          [~,ordLocal] = sort(diffs,'ascend');
          chosen = false;
          for kk = 1:numel(ordLocal)
              sIdx = startIdxSec(ordLocal(kk));
              if ~allowReuse && any(usedFrames(sIdx:sIdx+L2s-1)), continue; end
              s0 = ts(sIdx);
              FRm(t) = sum(st>=s0 & st<s0+(win(2)-win(1)));
              if ~allowReuse, usedFrames(sIdx:sIdx+L2s-1) = true; end
              chosen = true; break;
          end
          if ~chosen, FRm(t) = NaN; end
      end

    otherwise
      error('Unknown CompareMode: %s', compareMode);
  end

  % keep only finite, informative pairs
  ok  = isfinite(FRt) & isfinite(FRm);
  FRt = FRt(ok); FRm = FRm(ok);
  nz  = ~(FRt==0 & FRm==0);
  FRt = FRt(nz); FRm = FRm(nz);

  if numel(FRt) < max(1,minPairs), continue; end

  % paired test on equal-duration counts
  switch testType
    case 'ttest',   [~,p] = ttest(FRt, FRm);
    otherwise,      p     = signrank(FRt, FRm);
  end

  isTested(c) = true;
  pVal(c)     = p;

  % summaries (rates)
  winLen = diff(win);
  FRt_rate          = FRt / winLen;
  FRm_rate          = FRm / winLen;
  meanTrialRate(c)   = mean(FRt_rate,'omitnan');
  meanMatchedRate(c) = mean(FRm_rate,'omitnan');
  deltaRate(c)       = mean(FRt_rate - FRm_rate,'omitnan');
  deltaFR(c)         = mean(FRt) - mean(FRm);
end

fprintf('  Processed %d cells; tested %d after mask & data checks (%s mode, AllowReuse=%d)\n', ...
        nCells, sum(isTested), compareMode, allowReuse);

% FDR across tested cells
testedIdx = find(isTested);
hFDR = false(nCells,1);
if ~isempty(testedIdx)
    p_tested   = pVal(testedIdx);
    finiteMask = isfinite(p_tested);
    h_local    = false(numel(testedIdx),1);
    if any(finiteMask)
        [~,~,adjP_sub] = fdr_bh(p_tested(finiteMask), alpha);
        h_local(finiteMask) = adjP_sub < alpha;
    end
    hFDR(testedIdx) = h_local;
end

% pack
stats.deltaFR         = deltaFR;
stats.pVal            = pVal;
stats.isVelMod        = isTested;
stats.sigDiffFDR      = hFDR;
stats.pctSigDiffFDR   = 100*mean(hFDR(isTested));
stats.meanTrialRate   = meanTrialRate;
stats.meanMatchedRate = meanMatchedRate;
stats.deltaRate       = deltaRate;
end

%% Benjamini–Hochberg (simple)
function [h, crit_p, adj_p, sorted_p] = fdr_bh(pvals, q)
  p = pvals(:); p = p(isfinite(p));
  m = numel(p);
  if m==0, h=false(0,1); crit_p=0; adj_p=[]; sorted_p=[]; return; end
  [sorted_p, ids] = sort(p);
  thresh = (1:m)'/m * q;
  rej    = find(sorted_p<=thresh,1,'last');
  if isempty(rej)
    crit_p = 0; h=false(m,1);
  else
    crit_p=sorted_p(rej);
    h     = p<=crit_p;
  end
  adj_p = nan(m,1);
  for i=m:-1:1
    if i<m, adj_p(i) = min(sorted_p(i)*m/i, adj_p(i+1));
    else,   adj_p(i) = min(sorted_p(i)*m/i, 1);
    end
  end
  tmp=adj_p; adj_p(ids)=tmp;
end
