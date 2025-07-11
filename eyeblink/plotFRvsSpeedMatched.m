function plotFRvsSpeedMatched
% Compare trial-window FR to speed-matched non-trial FR, one paired test per cell, FDR-corrected.
% No inputs; expects in your workspace:
%   ratNames, rat.Ca_peaks, rat.pos, rat.CS_times, etc.

ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
nRats    = numel(ratNames);

% parameters
win      = [0 2];         % CS window (s)
binSize  = 1/7.5;         % 133 ms bins
alpha    = 0.05;          % omnibus FDR level

% preallocate
pctPerDayFDR = nan(nRats,3);
pctPerRatFDR = nan(nRats,1);
allDeltaFR   = cell(nRats,1);

for r = 1:nRats
  fprintf('Processing %s…\n', ratNames{r});
  rat   = evalin('base', ratNames{r});
  dates = autoDateList(rat);
  idx   = find(strcmp(dates, rat.An));
  days  = dates(idx-2:idx);

  tmpDel = [];
  for d=1:3
    day = days{d};
    spk   = rat.Ca_peaks.(['CA_peaks_' day]);  % [nCells x nBins]
    pos   = rat.pos.   (['pos_'      day]);    % [nBins x 3] → [time x x y]
    ts    = pos(:,1);
    xy    = pos(:,2:3);
    cs    = rat.CS_times.(['CS_'      day]);

    S = trialVsSpeedMatched( ...
        spk, ts, xy, cs, ...
        'win',win,'binSize',binSize,'alpha',alpha,'test','signrank' );

    pctPerDayFDR(r,d)  = S.pctSigDiffFDR;
    tmpDel            = [tmpDel; S.deltaFR(S.sigDiffFDR)]; %#ok<AGROW>
  end

  pctPerRatFDR(r) = mean(pctPerDayFDR(r,:),'omitnan');
  allDeltaFR{r}   = tmpDel;
end

%— Plotting —
figure('Color','w','Position',[200 300 1400 400]);

% (1) Day-by-day: now 5 rats × 3 days
subplot(1,3,1); hold on;
bar(pctPerDayFDR,'grouped');        % <— drop the transpose here
xticks(1:nRats); xticklabels(ratNames);
ylabel('% vel-mod cells (FDR)');
legend({'day-2','day-1','day-0'},'Location','northwest');
title('Per-day % sig (FDR)');

% (2) Rat means ± std over days
subplot(1,3,2); hold on;
errs = std(pctPerDayFDR,[],2,'omitnan');    % std across the 3 days, for each rat
bar(1:nRats,pctPerRatFDR,'FaceColor',[.2 .6 .8]);
errorbar(1:nRats,pctPerRatFDR,errs,'k.','LineWidth',1.5);
xticks(1:nRats); xticklabels(ratNames);
ylabel('% vel-mod cells (FDR)');
title('Mean ± STD per Rat');

% (3) ΔFR swarmchart
subplot(1,3,3); hold on;
for r=1:nRats
    swarmchart(r*ones(size(allDeltaFR{r})), allDeltaFR{r}, ...
               5,'filled','MarkerFaceAlpha',.4);
end
plot(xlim,[0 0],'k--');
xticks(1:nRats); xticklabels(ratNames);
ylabel('\DeltaFR (trial – matched)');
title('\DeltaFR of FDR-sig cells');

%%%%%%% Console summary
fprintf('\n=== FR vs Speed Matched (FDR α=%.2f) ===\n',alpha);
for r=1:nRats
  m = pctPerRatFDR(r);
  s = std(pctPerDayFDR(r,:),[],'omitnan');
  fprintf('%s: %.1f%% ± %.1f\n', ratNames{r}, m, s);
end



fprintf('Grand mean: %.1f%% ± %.1f\n\n', ...
        mean(pctPerRatFDR,'omitnan'), std(pctPerRatFDR,'omitnan'));
end


%%===========================================================================%%
function stats = trialVsSpeedMatched(spikeCell, ts, pos, csTimes, varargin)
% Compare trial FR to speed-matched non-trial FR, one paired test per cell, FDR-corrected.
% spikeCell: either a {nCells×1} cell array of spike-time vectors OR an [nCells×nTime]
%             numeric matrix of counts-per-frame.
% ts:        [nTime×1] time stamps.
% pos:       [nTime×2] positions.
% csTimes:   [nTrials×1] CS onset times.

ip = inputParser();
ip.addParameter('win',     [0 2]);
ip.addParameter('binSize',  1/7.5);
ip.addParameter('minPairs', []);
ip.addParameter('alpha',    0.05);
ip.addParameter('test',     'ttest');
ip.parse(varargin{:});
win      = ip.Results.win;
binSize  = ip.Results.binSize;
alpha    = ip.Results.alpha;
testType = lower(ip.Results.test);

% if caller passed a numeric counts‐matrix, wrap it into a cell array of per‐cell rows
if ~iscell(spikeCell) && isnumeric(spikeCell)
  % each row is one cell’s bin‐counts
  spikeCell = mat2cell(spikeCell, ones(size(spikeCell,1),1), size(spikeCell,2));
end


nTrials = numel(csTimes);
minPairs = nTrials;

% instantaneous speed
dt  = [diff(ts); binSize];
dXY = [diff(pos); [0 0]];
spd = hypot(dXY(:,1),dXY(:,2)) ./ dt;
spd(~isfinite(spd)) = 0;

% trial mask
inTrial = false(size(ts));
for t=1:nTrials
  inTrial = inTrial | (ts>=csTimes(t)+win(1) & ts<csTimes(t)+win(2));
end
outMask  = ~inTrial;
timeOut  = ts(outMask);
speedOut = spd(outMask);

% precompute trial windows & speeds
edgeMat  = csTimes(:) + win(:).';
trialSpd = nan(nTrials,1);
for t=1:nTrials
  idxCS       = ts>=edgeMat(t,1) & ts<edgeMat(t,2);
  trialSpd(t) = mean(spd(idxCS));
end

% allocate outputs
nCells   = numel(spikeCell);
deltaFR  = nan(nCells,1);
pVal     = nan(nCells,1);
isTested = false(nCells,1);
rawP     = [];

for c=1:nCells
  % build st = list of spike times
  if iscell(spikeCell)
    st = spikeCell{c}(:);
  elseif isnumeric(spikeCell)
    counts     = spikeCell(c,:);
    frameTimes = ts;
    st = [];
    Ipos = find(counts>0);
    for k=1:numel(Ipos)
      t0   = frameTimes(Ipos(k));
      nrep = round(counts(Ipos(k)));
      st = [st; repmat(t0, nrep, 1)];
    end
  else
    error('spikeCell must be cell array or numeric matrix');
  end

  if numel(st)<minPairs
    continue;
  end
  isTested(c)=true;

  % build FRt & FRm arrays
  FRt = zeros(nTrials,1);
  FRm = zeros(nTrials,1);
  for t=1:nTrials
    % count in CS window
    edges   = edgeMat(t,:);
    FRt(t)  = sum(st>=edges(1) & st<edges(2));
    need    = FRt(t);
    if need>0
      v     = trialSpd(t);
      diffs = abs(speedOut - v);
      [~,I] = sort(diffs,'ascend');
      take  = min(need,numel(I));
      for k2=1:take
        start    = timeOut(I(k2));
        FRm(t)   = FRm(t) + sum(st>=start & st<start+binSize);
      end
    end
  end

  % drop NaNs (shouldn't occur)
  ok  = ~isnan(FRt)&~isnan(FRm);
  FRt = FRt(ok);
  FRm = FRm(ok);

  if numel(FRt)<minPairs
    continue;
  end

  % paired test
  switch testType
    case 'ttest'
      [~,p] = ttest(FRt,FRm);
    otherwise
      p     = signrank(FRt,FRm);
  end

  rawP(end+1,1)    = p;              %#ok<AGROW>
  pVal(c)          = p;
  deltaFR(c)       = mean(FRt)-mean(FRm);
end

% FDR correction
[~,~,~,adjP]       = fdr_bh(rawP, alpha);
hFDR               = adjP < alpha;
sigDiffFDR         = false(nCells,1);
sigDiffFDR(isTested)=hFDR;

stats.deltaFR       = deltaFR;
stats.pVal          = pVal;
stats.isVelMod      = isTested;
stats.sigDiffFDR    = sigDiffFDR;
stats.pctSigDiffFDR = 100*mean(hFDR);
end


%% Simple Benjamini–Hochberg
function [h, crit_p, adj_p, sorted_p] = fdr_bh(pvals, q)
  p      = pvals(:);
  m      = numel(p);
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
    if i<m
      adj_p(i) = min(sorted_p(i)*m/i, adj_p(i+1));
    else
      adj_p(i) = min(sorted_p(i)*m/i, 1);
    end
  end
  tmp=adj_p; adj_p(ids)=tmp;
end
