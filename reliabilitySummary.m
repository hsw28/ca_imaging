function f = reliabilitySummary(ratNames, varargin)
% reliabilitySummary  Compare spike‐count variability in place field vs trace
%
% PF side:
%   - PF masks built from running (v>=4 cm/s) & non-trial samples/spikes only
%     using velocity-timebase filtering (same logic as RateMaskVsTask_summary).
%   - PF values per visit = COUNT (default) or RATE (spikes/s).
% Trace side:
%   - values per CS window [CS, CS+2] s = COUNT (default) or RATE (spikes/s).
% Metrics:
%   - Fano = var/mean; split-half = random halves correlation.
% Plots:
%   - Per rat: Fano PF vs Trace; split-half PF vs Trace.
%   - Pooled row.
%
% Usage:
%   reliabilitySummary({'rat0222','rat0307'})                   % counts (default)
%   reliabilitySummary({'rat0222','rat0307'}, 'mode','rate')   % rate-normalized

% ---------- args ----------
p = inputParser;
addParameter(p,'mode','count',@(s) any(validatestring(s,{'count','rate'})));
addParameter(p,'minVisitDur',0.1,@(x) isnumeric(x) && isscalar(x) && x>=0); % seconds
parse(p,varargin{:});
modeChoice  = lower(p.Results.mode);
minVisitDur = p.Results.minVisitDur;

nRats = numel(ratNames);
figure('Color','w','Position',[100 100 1200 300*(nRats+1)]);

allFanoPF   = [];
allFanoTr   = [];
allSplitPF  = [];
allSplitTr  = [];

for r = 1:nRats
  rat   = evalin('base', ratNames{r});
  dates = autoDateList(rat);
  iAn   = find(strcmp(dates,rat.An),1);
  days  = dates(max(1,iAn-2):iAn);  % last 3 days up to An

  fanoPF = [];  fanoTr = [];
  spPF   = [];  spTr   = [];

  for d = 1:numel(days)
    D       = days{d};
    S       = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));     % nCells×T OR cell array
    good    = find(rat.ratemask.(sprintf('ratemask_%s',D))==1);
    pos     = rat.pos.(sprintf('pos_%s',D));               % [t x y]
    csTimes = rat.CS_times.(sprintf('CS_%s',D))(:);

    t_pos   = pos(:,1);
    xy      = pos(:,2:3);

    % ---- speed from ca_velocity (velocity timebase) ----
    vel        = ca_velocity(pos);         % [speed; time]
    vel_time   = vel(2,:)';
    vel_mag    = vel(1,:)';

    % interpolate speed to the POSITION timebase for visit segmentation
    runMaskPos = interp1(vel_time, vel_mag, t_pos, 'linear','extrap') >= 4;

    % ---- PF mask from running & non-trial ----
    [maskGrid, xEdges, yEdges, inTrialPos] = computePFmask_runNoTrial_vel( ...
        S, pos, csTimes, 2.5, 1, 1);
    if isempty(maskGrid), continue, end
    [nY,nX,~] = size(maskGrid);

    % Precompute pos-sample bin indices for visit detection
    [~,~,bx_all] = histcounts(xy(:,1), xEdges);
    [~,~,by_all] = histcounts(xy(:,2), yEdges);
    goodSamp = bx_all>=1 & bx_all<=nX & by_all>=1 & by_all<=nY;

    for c = good(:)'
      % spike times for this cell
      st = getCellSpikes(S,c); st = st(~isnan(st) & st>0);
      if numel(st) < 3, continue, end

      % ---- PF counts as per-visit spike counts (v>=4 & non-trial, inside PF mask) ----
      thisMask = maskGrid(:,:,c);
      if ~any(thisMask,'all'), continue, end

      % in-PF over position samples
      inPF_all = false(numel(t_pos),1);
      if any(goodSamp)
        linIdx = sub2ind([nY,nX], by_all(goodSamp), bx_all(goodSamp));
        inPF_all(goodSamp) = thisMask(linIdx);
      end

      pfRunNonTrial = inPF_all & runMaskPos & ~inTrialPos;
      visits = logicalRuns(pfRunNonTrial, t_pos, minVisitDur);   % visits >= minVisitDur
      if isempty(visits), continue, end

      % per-visit counts and durations
      nVisits = size(visits,1);
      vCounts = zeros(nVisits,1);
      vDur    = zeros(nVisits,1);
      for iv = 1:nVisits
        t0 = visits(iv,1); t1 = visits(iv,2);
        vCounts(iv) = sum(st >= t0 & st < t1);
        vDur(iv)    = max(eps, t1 - t0); % avoid 0
      end

      % ---- Trace counts per CS in [0,2] s ----
      if isempty(csTimes), continue, end
      countsTr = arrayfun(@(t0) sum(st>=t0 & st<t0+2), csTimes);
      trDur    = 2; % s

      % ---- values for Fano/split-half: counts vs rates ----
      switch modeChoice
        case 'count'
          valsPF = vCounts;
          valsTr = countsTr;
        case 'rate'
          valsPF = vCounts ./ vDur;     % spikes/s per visit
          valsTr = countsTr / trDur;    % spikes/s per trial
      end

      % ---- accumulate if valid ----
      if numel(valsPF)>1 && numel(valsTr)>1
        mPF = mean(valsPF); mTr = mean(valsTr);
        if mPF>0 && mTr>0
          fanoPF(end+1) = var(valsPF)/mPF;           %#ok<AGROW>
          fanoTr(end+1) = var(valsTr)/mTr;           %#ok<AGROW>
          spPF(end+1)   = splitHalfCorr(valsPF);     %#ok<AGROW>
          spTr(end+1)   = splitHalfCorr(valsTr);     %#ok<AGROW>
        else
          fanoPF(end+1) = NaN; fanoTr(end+1) = NaN;  %#ok<AGROW>
          spPF(end+1)   = NaN; spTr(end+1)   = NaN;  %#ok<AGROW>
        end
      end
    end
  end

  % aggregate
  allFanoPF  = [allFanoPF;  fanoPF(:)];
  allFanoTr  = [allFanoTr;  fanoTr(:)];
  allSplitPF = [allSplitPF; spPF(:)];
  allSplitTr = [allSplitTr; spTr(:)];

  % ---- plotting per rat ----
  ax1 = subplot(nRats+1,2, (r-1)*2+1); hold(ax1,'on'); axis(ax1,'square')
  scatter(ax1,fanoPF, fanoTr, 12, 'filled');
  lim = [0, max([fanoPF(:); fanoTr(:); 1])*1.05];
  xlim(ax1,lim); ylim(ax1,lim);
  plot(ax1, lim, lim,'k--');
  xlabel(ax1,sprintf('Fano PF (%s; run, non-trial)', modeChoice));
  ylabel(ax1,sprintf('Fano Trace (%s; 0–2 s)', modeChoice));
  title(ax1, ratNames{r});

  ax2 = subplot(nRats+1,2, r*2); hold(ax2,'on'); axis(ax2,'square')
  scatter(ax2,spPF, spTr, 12,'filled');
  lim2 = [min([spPF(:);spTr(:);-1]) max([spPF(:);spTr(:);1])];
  plot(ax2, [lim2(1) lim2(2)], [lim2(1) lim2(2)], 'k--');
  xlim(ax2,lim2); ylim(ax2,lim2);
  xlabel(ax2,sprintf('Split-half PF (%s)', modeChoice));
  ylabel(ax2,sprintf('Split-half Trace (%s)', modeChoice));
end

% ---- pooled row ----
ax1 = subplot(nRats+1,2, nRats*2+1); hold(ax1,'on'); axis(ax1,'square')
scatter(ax1,allFanoPF,allFanoTr,12,'r','filled');
lim = [0, max([allFanoPF(:); allFanoTr(:); 1])*1.05];
xlim(ax1,lim); ylim(ax1,lim);
plot(ax1, lim, lim,'k--');
xlabel(ax1,sprintf('Fano PF (%s; run, non-trial)', modeChoice));
ylabel(ax1,sprintf('Fano Trace (%s; 0–2 s)', modeChoice));
title(ax1,'All pooled');

f = [allFanoPF,allFanoTr]; % keep original return shape

ax2 = subplot(nRats+1,2, nRats*2+2); hold(ax2,'on'); axis(ax2,'square')
scatter(ax2,allSplitPF,allSplitTr,12,'r','filled');
lim2 = [min([allSplitPF(:);allSplitTr(:);-1]) max([allSplitPF(:);allSplitTr(:);1])];
plot(ax2, [lim2(1) lim2(2)], [lim2(1) lim2(2)], 'k--');
xlim(ax2,lim2); ylim(ax2,lim2);
xlabel(ax2,sprintf('Split-half PF (%s)', modeChoice));
ylabel(ax2,sprintf('Split-half Trace (%s)', modeChoice));
title(ax2,'All pooled');
end

% ================= helpers =================
function [maskGrid,xEdges,yEdges,inTrialPos] = computePFmask_runNoTrial_vel(S, pos, csTimes, binSize, smoothFlag, N)
t_pos = pos(:,1);
inTrialPos = false(size(t_pos));
for t = csTimes(:).', end
for t = csTimes(:).'
  inTrialPos = inTrialPos | (t_pos>=t & t_pos<t+2);
end
vel      = ca_velocity(pos);    % [speed; time]
vt       = vel(2,:)'; vmag = vel(1,:)';
xv = interp1(t_pos, pos(:,2), vt, 'linear', NaN);
yv = interp1(t_pos, pos(:,3), vt, 'linear', NaN);
inTrialVel = false(size(vt));
for t = csTimes(:).'
  inTrialVel = inTrialVel | (vt>=t & vt<t+2);
end
keep = (vmag>=4) & ~inTrialVel & isfinite(xv) & isfinite(yv);
posPF = [vt(keep), xv(keep), yv(keep)];
if size(posPF,1) < 5, maskGrid = []; xEdges=[]; yEdges=[]; return; end
xEdges = floor(min(posPF(:,2))):binSize:ceil(max(posPF(:,2)));
yEdges = floor(min(posPF(:,3))):binSize:ceil(max(posPF(:,3)));
ny = numel(yEdges)-1; nx = numel(xEdges)-1;
nCells = size(S,1); maskGrid = false(ny,nx,nCells);
for c = 1:nCells
  st = getCellSpikes(S,c); st = st(~isnan(st) & st>0);
  if numel(st)<3, continue, end
  spd_spk   = interp1(vt, vmag, st, 'linear', 'extrap');
  inTrialSpk= isInTrial_at_spike(st, csTimes, [0 2]);
  st_pf     = st((spd_spk>=4) & ~inTrialSpk);
  if isempty(st_pf), continue, end
  rate = CA_normalizePosData(st_pf, posPF, binSize, smoothFlag);
  if ~isequal(size(rate), [ny nx]), rate = imresize(rate, [ny nx], 'nearest'); end
  if any(isfinite(rate(:)) & rate(:)~=0)
    thr = nanmean(rate(:)) + N*nanstd(rate(:));
    maskGrid(:,:,c) = (rate >= thr) & isfinite(rate);
  end
end
end

function st = getCellSpikes(S, c)
  if iscell(S), st = S{c}(:); else, st = S(c,:).'; end
end

function in = isInTrial_at_spike(st, csTimes, win)
  st = st(:); in = false(size(st));
  for k = 1:numel(csTimes)
    in = in | (st>=csTimes(k)+win(1) & st<csTimes(k)+win(2));
  end
end

function visits = logicalRuns(mask, t, minDur)
  mask = mask(:)>0; t = t(:);
  d = diff([false; mask; false]);
  iStart = find(d==1); iEnd = find(d==-1)-1;
  if isempty(iStart), visits = zeros(0,2); return; end
  % guard last index for t(iEnd+1)
  tPad = [t; t(end)+(t(end)-t(end-1))];
  tStart = t(iStart); tEnd = tPad(iEnd+1);
  dur = tEnd - tStart;
  keep = dur >= minDur;
  visits = [tStart(keep) tEnd(keep)];
end

function r = splitHalfCorr(v)
  n = numel(v);
  if n < 2, r = NaN; return; end
  idx = randperm(n);
  h   = floor(n/2);
  idx = idx(1:2*h);
  r   = corr( v(idx(1:h)), v(idx(h+1:2*h)), 'Rows','pairwise' );
end
