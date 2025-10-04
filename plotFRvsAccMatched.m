function plotFRvsAccMatched(matchBy)
  %DEP
% plotFRvsAccMatched  Compare trial-window FR to signal-matched (speed or accel) non-trial FR.
% FDR-corrected stats, with per-day and per-rat summaries, plus Δ-rate and log2FC violins.
%
% Usage:
%   plotFRvsAccMatched                 % defaults to 'accel' (SIGNED d|v|/dt)
%   plotFRvsAccMatched('speed')        % speed-matched
%   plotFRvsAccMatched('accel')        % signed acceleration (d|v|/dt)  [DEFAULT]
%   plotFRvsAccMatched('accelmag')     % |acceleration| magnitude (vector)
%   plotFRvsAccMatched('absaccel')     % alias of 'accelmag'
%
% Notes:
% - Position units + ts units define speed (posUnit/s) and accel (posUnit/s^2).
% - Matching uses the mean of the chosen signal in the [CS, CS+2] s window per trial.

if nargin<1, matchBy = 'accel'; end

ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
nRats    = numel(ratNames);

% Parameters
win       = [0 2];        % CS window (s)
binSize   = 1/7.5;        % 133 ms bins
alpha     = 0.05;         % FDR level
lfcLambda = 0.05;         % Hz pseudocount for log2FC stability

% Preallocate collectors
pctPerDayIncluded  = nan(nRats,3);
pctPerRatIncluded  = nan(nRats,1);
allDeltaFRIncl     = cell(nRats,1);
rateDelta_sig_byRat = cell(nRats,1);  % Δ rate (Hz), sig only
lfc_sig_byRat       = cell(nRats,1);  % log2 FC, sig only

for r = 1:nRats
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);   % day−2, day−1, day−0

    delIncl = [];
    rateD_sig_allDays = [];
    lfc_sig_allDays   = [];

    for d=1:3
        day      = days{d};
        spk      = rat.Ca_peaks.(['CA_peaks_' day]);   % cell array or counts matrix
        posFull  = rat.pos.     (['pos_'      day]);   % [time x y]
        posFull  = smoothpos(posFull);                 % your existing smoother
        ts       = posFull(:,1);
        xy       = posFull(:,2:3);
        cs       = rat.CS_times.(sprintf('CS_%s', day));
        ratemask = rat.ratemask.(sprintf('ratemask_%s', day));

        Sincl = trialVsSignalMatched( spk, ts, xy, cs, ratemask, ...
              'win',win,'binSize',binSize,'alpha',alpha,'test','ttest', ...
              'matchBy',matchBy);

        pctPerDayIncluded(r,d) = Sincl.pctSigDiffFDR;

        % Δ counts (sig only) for legacy panel
        delIncl = [delIncl; Sincl.deltaFR(Sincl.sigDiffFDR)]; %#ok<AGROW>

        % Δ rate (Hz) + log2 FC (sig only)
        if any(Sincl.sigDiffFDR)
            tRate = Sincl.meanTrialRate(  Sincl.sigDiffFDR);   % Hz
            mRate = Sincl.meanMatchedRate(Sincl.sigDiffFDR);   % Hz
            rateD = Sincl.deltaRate(     Sincl.sigDiffFDR);    % Hz
            lfc   = log2( (tRate + lfcLambda) ./ (mRate + lfcLambda) );

            rateD_sig_allDays = [rateD_sig_allDays; rateD]; %#ok<AGROW>
            lfc_sig_allDays   = [lfc_sig_allDays;   lfc];   %#ok<AGROW>
        end
    end

    pctPerRatIncluded(r) = mean(pctPerDayIncluded(r,:),'omitnan');
    allDeltaFRIncl{r}    = delIncl;

    rateDelta_sig_byRat{r} = rateD_sig_allDays;
    lfc_sig_byRat{r}       = lfc_sig_allDays;
end

%% — PLOTTING (per-day % sig, per-rat mean % sig, Δ counts) —
figure('Color','w','Position',[200 300 1400 800]);

% (1) INCLUDED — per-day
subplot(1,3,1); hold on;
bar(pctPerDayIncluded,'grouped');
xticks(1:nRats); xticklabels(ratNames);
ylabel('% sig (included)'); title(sprintf('INCLUDED — per day (match: %s)', matchBy));
legend({'day−2','day−1','day−0'},'Location','northwest');

% (2) INCLUDED — per-rat means
subplot(1,3,2); hold on;
bar(pctPerRatIncluded,'FaceColor',[.2 .6 .8]);
errorbar(1:nRats,pctPerRatIncluded, std(pctPerDayIncluded,[],2),'k.','LineWidth',1.5);
xticks(1:nRats); xticklabels(ratNames);
ylabel('mean % (included)'); title('INCLUDED — per rat');

% (3) ΔFR (counts) for INCLUDED — violin
subplot(1,3,3); hold on;
for r=1:nRats
    d = allDeltaFRIncl{r};
    if ~isempty(d)
        drawViolin(r, d);
        mu = mean(d,'omitnan');
        se = std(d,'omitnan')/sqrt(sum(~isnan(d)));
        errorbar(r, mu, se, 'k.','LineWidth',1.2,'CapSize',12);
    end
end
yline(0,'k--');
xticks(1:nRats); xticklabels(ratNames);
ylabel('\Delta FR (counts, included)'); title('\Delta FR (counts)');

%% — NEW FIGURE: Δ rate & log2 fold-change (sig cells) —
figure('Color','w','Position',[200 200 1400 450]);

% (A) Δ rate violin (sig cells)
subplot(1,2,1); hold on;
pool_rateD = vertcat(rateDelta_sig_byRat{:});
for r=1:nRats
    d = rateDelta_sig_byRat{r};
    if ~isempty(d)
        drawViolin(r, d);
        mu = mean(d,'omitnan');
        se = std(d,'omitnan')/sqrt(sum(~isnan(d)));
        errorbar(r, mu, se, 'k.','LineWidth',1.2,'CapSize',12);
    end
end
if ~isempty(pool_rateD)
    drawViolin(nRats+1, pool_rateD);
    muP = mean(pool_rateD,'omitnan'); seP = std(pool_rateD,'omitnan')/sqrt(sum(~isnan(pool_rateD)));
    errorbar(nRats+1, muP, seP, 'k.','LineWidth',1.2,'CapSize',12);
end
yline(0,'k--');
xticks(1:nRats+1); xticklabels([ratNames, {'All'}]);
ylabel('\Delta rate (events/s), sig cells');
title(sprintf('Sig cells — trial rate − matched rate (match: %s)', matchBy));

% (B) log2 fold-change violin (sig cells)
subplot(1,2,2); hold on;
pool_lfc = vertcat(lfc_sig_byRat{:});
for r=1:nRats
    d = lfc_sig_byRat{r};
    if ~isempty(d)
        drawViolin(r, d);
        mu = mean(d,'omitnan'); se = std(d,'omitnan')/sqrt(sum(~isnan(d)));
        errorbar(r, mu, se, 'k.','LineWidth',1.2,'CapSize',12);
    end
end
if ~isempty(pool_lfc)
    drawViolin(nRats+1, pool_lfc);
    muP = mean(pool_lfc,'omitnan'); seP = std(pool_lfc,'omitnan')/sqrt(sum(~isnan(pool_lfc)));
    errorbar(nRats+1, muP, seP, 'k.','LineWidth',1.2,'CapSize',12);
end
yline(0,'k--');             % 0 => no change (trial==matched)
yticks(-2:1:2); yticklabels({'0.25×','0.5×','1×','2×','4×'});
xticks(1:nRats+1); xticklabels([ratNames, {'All'}]);
ylabel('log_2 fold change (trial / matched)');
title('Sig cells — log_2 fold-change in rate');

end % main

%%===========================================================================%%
function stats = trialVsSignalMatched(spikeCell, ts, pos, csTimes, ratemask, varargin)
% trialVsSignalMatched  Compare trial FR to matched non-trial FR using:
%   'speed'     -> speed (|v|)
%   'accel'     -> SIGNED d|v|/dt  (derivative of speed)     [DEFAULT]
%   'accelmag'  -> |acceleration|  (vector magnitude)
%   'absaccel'  -> alias for 'accelmag'

ip = inputParser();
ip.addParameter('win',      [0 2]);
ip.addParameter('binSize',   1/7.5);
ip.addParameter('minPairs',  []);
ip.addParameter('alpha',     0.05);
ip.addParameter('test',      'ttest', @(s) any(validatestring(s,{'ttest','signrank'})));
ip.addParameter('matchBy',   'accel', @(s) any(validatestring(s,{'speed','accel','accelmag','absaccel'})));
ip.parse(varargin{:});
win      = ip.Results.win;
binSize  = ip.Results.binSize;
alpha    = ip.Results.alpha;
testType = lower(ip.Results.test);
matchBy  = lower(ip.Results.matchBy);

% normalize spikes to cell array (same behavior as your working code)
if ~iscell(spikeCell) && isnumeric(spikeCell)
  spikeCell = mat2cell(spikeCell, ones(size(spikeCell,1),1), size(spikeCell,2));
end

nTrials  = numel(csTimes);
minPairs = nTrials;
winLen   = diff(win);

% ---------- matching signals ----------
% SPEED (identical style as working code)
dt  = [diff(ts); binSize];           % s per sample (last padded with binSize)
dx  = [diff(pos(:,1)); 0];
dy  = [diff(pos(:,2)); 0];
vx  = dx ./ dt;
vy  = dy ./ dt;
speed = hypot(vx, vy);                % |v| (posUnit/s)
speed(~isfinite(speed)) = 0;

% SIGNED d|v|/dt (derivative of speed)
acc_signed = [diff(speed)./dt(1:end-1); 0];   % posUnit/s^2, can be ±
acc_signed(~isfinite(acc_signed)) = 0;

% VECTOR acceleration magnitude |a|
ax = [diff(vx)./dt(1:end-1); 0];
ay = [diff(vy)./dt(1:end-1); 0];
acc_mag = hypot(ax, ay);                  % ≥0, posUnit/s^2
acc_mag(~isfinite(acc_mag)) = 0;

% (optional smoothing)
% speed    = movmean(speed,    3);
% acc_signed = movmean(acc_signed, 3);
% acc_mag  = movmean(acc_mag,  3);

switch matchBy
  case 'speed'
    sig = speed;
  case 'accel'                 % DEFAULT: signed d|v|/dt
    sig = acc_signed;
  case {'accelmag','absaccel'} % magnitude
    sig = acc_mag;
end

% ---------- trial / non-trial masks ----------
inTrial = false(size(ts));
for t=1:nTrials
  inTrial = inTrial | (ts>=csTimes(t)+win(1) & ts<csTimes(t)+win(2));
end
outMask   = ~inTrial;
timeOut   = ts(outMask);
sigOut    = sig(outMask);

% per-trial mean of chosen signal
edgeMat   = csTimes(:) + win(:).';
trialSig  = nan(nTrials,1);
for t=1:nTrials
  idxCS        = ts>=edgeMat(t,1) & ts<edgeMat(t,2);
  trialSig(t)  = mean(sig(idxCS), 'omitnan');
end

% ---------- allocate outputs ----------
nCells   = numel(spikeCell);
deltaFR  = nan(nCells,1);
pVal     = nan(nCells,1);
isTested = false(nCells,1);
rawP     = [];

meanTrialRate   = nan(nCells,1);   % Hz
meanMatchedRate = nan(nCells,1);   % Hz
deltaRate       = nan(nCells,1);   % Hz

% ---------- per cell ----------
for c=1:nCells
  % build st = list of spike times
  if iscell(spikeCell)
    st = spikeCell{c}(:);
  else
    counts     = spikeCell(c,:);
    frameTimes = ts;
    st = [];
    Ipos = find(counts>0);
    for k=1:numel(Ipos)
      t0   = frameTimes(Ipos(k));
      nrep = round(counts(Ipos(k)));
      st = [st; repmat(t0, nrep, 1)];
    end
  end

  if sum(~isnan(st))<minPairs, continue; end
  if ratemask(c) == 0, continue; end
  isTested(c)=true;

  FRt = zeros(nTrials,1);
  FRm = zeros(nTrials,1);
  mDur= zeros(nTrials,1);

  for t=1:nTrials
    edges   = edgeMat(t,:);
    FRt(t)  = sum(st>=edges(1) & st<edges(2));
    need    = FRt(t);

    if need>0
      v       = trialSig(t);
      diffs   = abs(sigOut - v);           % 1-D nearest-neighbor on chosen signal
      [~,I]   = sort(diffs,'ascend');
      take    = min(need, numel(I));
      mDur(t) = take * binSize;
      for k2=1:take
        start  = timeOut(I(k2));
        FRm(t) = FRm(t) + sum(st>=start & st<start+binSize);
      end
    else
      mDur(t)=0;
    end
  end

  ok  = ~isnan(FRt) & ~isnan(FRm);
  FRt = FRt(ok); FRm = FRm(ok); mDur = mDur(ok);
  if numel(FRt)<minPairs, continue; end

  % counts-based summary + test
  deltaFR(c) = mean(FRt) - mean(FRm);
  switch testType
    case 'ttest',   [~,p] = ttest(FRt,FRm);
    otherwise,      p     = signrank(FRt,FRm);
  end
  rawP(end+1,1) = p;
  pVal(c)       = p;

  % rate-based summaries (Hz)
  FRt_rate = FRt / winLen;
  FRm_rate = zeros(size(FRm));
  nz = mDur>0;
  FRm_rate(nz) = FRm(nz) ./ mDur(nz);
  meanTrialRate(c)   = mean(FRt_rate,'omitnan');
  meanMatchedRate(c) = mean(FRm_rate,'omitnan');
  deltaRate(c)       = mean(FRt_rate - FRm_rate,'omitnan');
end

fprintf('[%s] Processed %d cells; tested %d\n', upper(matchBy), nCells, sum(isTested));

% ---------- FDR ----------
[~,~,adjP] = fdr_bh(rawP, alpha);
hFDR = adjP < alpha;
sigDiffFDR = false(nCells,1);
sigDiffFDR(isTested) = hFDR;

% pack outputs
stats.deltaFR         = deltaFR;
stats.pVal            = pVal;
stats.isVelMod        = isTested;
stats.sigDiffFDR      = sigDiffFDR;
stats.pctSigDiffFDR   = 100*mean(hFDR);
stats.meanTrialRate   = meanTrialRate;   % Hz
stats.meanMatchedRate = meanMatchedRate; % Hz
stats.deltaRate       = deltaRate;       % Hz
end

%% Benjamini–Hochberg
function [h, crit_p, adj_p, sorted_p] = fdr_bh(pvals, q)
p      = pvals(:);
m      = numel(p);
[sorted_p, ids] = sort(p);
thresh = (1:m)'/m * q;
rej    = find(sorted_p<=thresh,1,'last');
if isempty(rej)
  crit_p = 0; h=false(m,1);
else
  crit_p=sorted_p(rej); h = p<=crit_p;
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
