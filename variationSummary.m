function f = variationSummary(ratNames, varargin)
% variationSummary  Compare variability in PF vs. trace windows
% CV = std(values) / mean(values)
%
% PF side:
%   - PF mask from NON-TRIAL running only (v >= 4 cm/s, CS/US excluded),
%     built using velocity-timebase filtering.
%   - Values per visit = spike COUNT (default) or spike RATE (spikes/s).
% Trace side:
%   - Values per trial = spike COUNT (default) or spike RATE in [CS, CS+2] s.
%
% Usage:
%   variationSummary({'rat0222','rat0307'}, 'Mode','count');  % default
%   variationSummary({'rat0222','rat0307'}, 'Mode','rate');   % rate-normalized
%   variationSummary({'rat0222'}, 'Mode','count','MinVisitDur',0.5);

% ---------------- args ----------------
p = inputParser;
addParameter(p,'Mode','count',@(s) any(validatestring(s,{'count','rate'})));
addParameter(p,'MinVisitDur',0.50,@(x) isnumeric(x) && isscalar(x) && x>=0); % seconds
parse(p,varargin{:});
cvMode      = lower(p.Results.Mode);
minVisitDur = p.Results.MinVisitDur;

nRats = numel(ratNames);

% pooled holders
allCV_PF   = []; allCV_Tr   = [];
allMean_PF = []; allMean_Tr = [];

figure('Color','w','Position',[100 100 1100 280*(nRats+1)]);

for r = 1:nRats
  rat   = evalin('base', ratNames{r});
  dates = autoDateList(rat);
  iAn   = find(strcmp(dates,rat.An),1);
  days  = dates(max(1,iAn-2):iAn);

  % per-rat holders
  cvPF = []; cvTr = [];
  meanPF = []; meanTr = [];

  for d = 1:numel(days)
    D       = days{d};
    S       = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));   % row-of-times or cell array per cell
    good    = find(rat.ratemask.(sprintf('ratemask_%s',D))==1);
    pos     = rat.pos.(sprintf('pos_%s',D));             % [t x y]
    csTimes = rat.CS_times.(sprintf('CS_%s',D))(:);

    t_pos = pos(:,1); xy = pos(:,2:3);

    % -------- speed from ca_velocity (velocity timebase) --------
    vel      = ca_velocity(pos);         % [speed; time]
    vt       = vel(2,:)';
    vmag     = vel(1,:)';

    % interpolate speed onto POSITION timebase to gate visits later
    runMaskPos = interp1(vt, vmag, t_pos, 'linear','extrap') >= 4;

    % -------- PF mask from NON-TRIAL running only (velocity timebase) --------
    [maskGrid, xEdges, yEdges, inTrialPos] = computePFmask_runNoTrial_vel( ...
        S, pos, csTimes, 2.5, 1, 1);
    if isempty(maskGrid), continue; end
    [nY,nX,~] = size(maskGrid);

    % Precompute sample-wise bin indices for position stream
    [~,~,bx_all] = histcounts(xy(:,1), xEdges);
    [~,~,by_all] = histcounts(xy(:,2), yEdges);
    goodSamp = bx_all>=1 & bx_all<=nX & by_all>=1 & by_all<=nY;

    for c = good(:)'
      st = getCellSpikes(S,c); st = st(~isnan(st) & st>0);
      if numel(st) < 3
        cvPF(end+1) = NaN;  cvTr(end+1) = NaN;
        meanPF(end+1) = NaN; meanTr(end+1) = NaN;
        continue
      end

      % ---- PF visits: in-PF & running & non-trial (position timebase) ----
      thisMask = maskGrid(:,:,c);
      if ~any(thisMask,'all')
        cvPF(end+1) = NaN;  cvTr(end+1) = NaN;
        meanPF(end+1) = NaN; meanTr(end+1) = NaN;
        continue
      end

      inPF_all = false(numel(t_pos),1);
      if any(goodSamp)
        linIdxSamp = sub2ind([nY,nX], by_all(goodSamp), bx_all(goodSamp));
        inPF_all(goodSamp) = thisMask(linIdxSamp);
      end

      pfRunNonTrial = inPF_all & runMaskPos & ~inTrialPos;

      % Segment contiguous bouts (visits)
      visits = logicalRuns(pfRunNonTrial, t_pos, minVisitDur);
      if isempty(visits)
        cvPF(end+1) = NaN;  cvTr(end+1) = NaN;
        meanPF(end+1) = NaN; meanTr(end+1) = NaN;
        continue
      end

      % Count spikes per visit and visit durations
      nVisits  = size(visits,1);
      vCounts  = zeros(nVisits,1);
      vDur     = zeros(nVisits,1);
      for iv = 1:nVisits
        t0 = visits(iv,1); t1 = visits(iv,2);
        vCounts(iv) = sum(st >= t0 & st < t1);
        vDur(iv)    = max(eps, t1 - t0);   % avoid divide-by-zero
      end

      % ---- Trace counts: spikes per CS in [0,2] s (no speed filter) ----
      if isempty(csTimes)
        cvPF(end+1) = NaN;  cvTr(end+1) = NaN;
        meanPF(end+1) = NaN; meanTr(end+1) = NaN;
        continue
      end
      countsTr = arrayfun(@(t0) sum(st>=t0 & st<t0+2), csTimes);
      trDur    = 2; % seconds

      % ---- choose values for CV: counts vs rates ----
      switch cvMode
        case 'count'
          valsPF = vCounts;
          valsTr = countsTr;
        case 'rate'
          valsPF = vCounts ./ vDur;     % spikes/s per visit
          valsTr = countsTr / trDur;    % spikes/s per trial
      end

      mPF = nanmean(valsPF); mTr = nanmean(valsTr);
      sPF = nanstd(valsPF);  sTr = nanstd(valsTr);

      if numel(valsPF)>1 && mPF>0 && numel(valsTr)>1 && mTr>0
          cvPF(end+1)   = sPF / mPF;
          cvTr(end+1)   = sTr / mTr;
          meanPF(end+1) = mPF;
          meanTr(end+1) = mTr;
      else
          cvPF(end+1)   = NaN;
          cvTr(end+1)   = NaN;
          meanPF(end+1) = NaN;
          meanTr(end+1) = NaN;
      end
    end
  end

  % aggregate + plot per rat
  allCV_PF   = [allCV_PF;   cvPF(:)];
  allCV_Tr   = [allCV_Tr;   cvTr(:)];
  allMean_PF = [allMean_PF; meanPF(:)];
  allMean_Tr = [allMean_Tr; meanTr(:)];

  ax = subplot(nRats+1,2,(r-1)*2+1); hold(ax,'on'); axis(ax,'square')
  scatter(ax,cvPF, cvTr, 12,'filled');
  lim = [0 max([cvPF(:);cvTr(:);1],[],'omitnan')*1.05];
  xlim(ax,lim); ylim(ax,lim);
  plot(ax, lim, lim,'k--');
  xlabel(ax, sprintf('CV PF (%s, run, non-trial)', cvMode));
  ylabel(ax, sprintf('CV Trace (%s, 0–2 s)', cvMode));
  title(ax, sprintf('%s', ratNames{r}));
end

% -------- pooled CV plot --------
ax = subplot(nRats+1,2,nRats*2+1); hold(ax,'on'); axis(ax,'square')
scatter(ax,allCV_PF,allCV_Tr,12,'r','filled');
lim = [0 max([allCV_PF(:);allCV_Tr(:);1],[],'omitnan')*1.05];
xlim(ax,lim); ylim(ax,lim);
plot(ax, lim, lim,'k--');
xlabel(ax, sprintf('CV PF (%s, run, non-trial)', cvMode));
ylabel(ax, sprintf('CV Trace (%s, 0–2 s)', cvMode));
title(ax,'All pooled');

% -------- CV vs CV and Fano-from-CV (paired, same cells) --------
figure;
subplot(1,2,1);
scatter(allCV_PF, allCV_Tr, 10, 'filled'); axis square
xlabel('CV PF'); ylabel('CV Trace'); hold on
lim = [0 max([allCV_PF(:);allCV_Tr(:);1],[],'omitnan')*1.05];
xlim(lim); ylim(lim); plot(lim,lim,'k--');

subplot(1,2,2);
% require validity for BOTH sides so points correspond cell-by-cell
validBoth = ~isnan(allCV_PF) & ~isnan(allCV_Tr) & ...
            ~isnan(allMean_PF) & ~isnan(allMean_Tr) & ...
            (allMean_PF>0) & (allMean_Tr>0);

fano_from_cv_PF = (allCV_PF(validBoth).^2) .* allMean_PF(validBoth);
fano_from_cv_Tr = (allCV_Tr(validBoth).^2) .* allMean_Tr(validBoth);
scatter(fano_from_cv_PF, fano_from_cv_Tr, 10, 'filled'); axis square
lim = [0 max([fano_from_cv_PF(:); fano_from_cv_Tr(:); 1],[],'omitnan')*1.05];
xlim(lim); ylim(lim); hold on; plot(lim,lim,'k--');
xlabel('Fano from CV (PF)'); ylabel('Fano from CV (Trace)');

% return pooled CVs (PF, Trace)
f = [allCV_PF, allCV_Tr];
end

% ================= helpers =================
function [maskGrid,xEdges,yEdges,inTrialPos] = computePFmask_runNoTrial_vel(S, pos, csTimes, binSize, smoothFlag, N)
% PF mask from NON-TRIAL running only (v>=4 & outside [CS,CS+2]),
% using velocity-timebase filtering.

t_pos = pos(:,1);

% trial membership at POS samples (for visit gating later)
inTrialPos = false(size(t_pos));
for t = csTimes(:).'
  inTrialPos = inTrialPos | (t_pos>=t & t_pos<t+2);
end

% Velocity stream
vel  = ca_velocity(pos);          % [speed; time]
vt   = vel(2,:)';
vmag = vel(1,:)';

% Interpolate x,y to velocity timestamps to form filtered position on vt
xv = interp1(t_pos, pos(:,2), vt, 'linear', NaN);
yv = interp1(t_pos, pos(:,3), vt, 'linear', NaN);

% trial membership at VELOCITY samples
inTrialVel = false(size(vt));
for t = csTimes(:).'
  inTrialVel = inTrialVel | (vt>=t & vt<t+2);
end

% Keep only running & non-trial velocity samples with valid position
keep = (vmag>=4) & ~inTrialVel & isfinite(xv) & isfinite(yv);
posPF = [vt(keep), xv(keep), yv(keep)];
if size(posPF,1) < 5
  maskGrid = []; xEdges=[]; yEdges=[]; return;
end

% Edges from kept pos (consistent with discretize later)
xEdges = floor(min(posPF(:,2))):binSize:ceil(max(posPF(:,2)));
yEdges = floor(min(posPF(:,3))):binSize:ceil(max(posPF(:,3)));
ny = numel(yEdges)-1; nx = numel(xEdges)-1;

nCells = size(S,1);
maskGrid = false(ny,nx,nCells);

for c = 1:nCells
  st = getCellSpikes(S,c); st = st(~isnan(st) & st>0);
  if numel(st)<3, continue, end

  % restrict spikes to running & non-trial (velocity timebase)
  spd_spk    = interp1(vt, vmag, st, 'linear', 'extrap');
  inTrialSpk = false(size(st));
  for k = 1:numel(csTimes)
    inTrialSpk = inTrialSpk | (st>=csTimes(k) & st<csTimes(k)+2);
  end
  st_pf = st((spd_spk>=4) & ~inTrialSpk);
  if isempty(st_pf), continue, end

  % rate map from running/non-trial spikes on running/non-trial pos
  rate = CA_normalizePosData(st_pf, posPF, binSize, smoothFlag);
  if ~isequal(size(rate), [ny nx])
    rate = imresize(rate, [ny nx], 'nearest');
  end

  if any(isfinite(rate(:)) & rate(:)~=0)
    thr = nanmean(rate(:)) + N*nanstd(rate(:));
    maskGrid(:,:,c) = (rate >= thr) & isfinite(rate);
  end
end
end

function visits = logicalRuns(mask, t, minDur)
% Convert a logical mask over samples to [t_start t_end] rows for each bout
% Only keep bouts with duration >= minDur (seconds).
mask = mask(:)>0; t = t(:);
d = diff([false; mask; false]);
starts = find(d==1);
ends   = find(d==-1)-1;
if isempty(starts), visits = []; return; end
tStart = t(starts);
% guard last index for t(ends+1)
tPad = [t; t(end)+(t(end)-t(end-1))];
tEnd  = tPad(ends+1);
dur   = tEnd - tStart;
keep  = dur >= minDur;
visits = [tStart(keep) tEnd(keep)];
end

function st = getCellSpikes(S, c)
% Accept row-of-times numeric or cell array per cell
  if iscell(S), st = S{c}(:);
  else,         st = S(c,:).';
  end
end
