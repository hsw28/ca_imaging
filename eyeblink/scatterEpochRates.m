function [rMat,pMat] = scatterEpochRates(ratNames, useVelFilter, velThresh)
% scatterEpochRates   Scatterplots of normalized epoch‐rate correlations
%
%   [rMat,pMat] = scatterEpochRates({'rat0222','rat0307',...})
%   [rMat,pMat] = scatterEpochRates(ratNames, true)            % speed filter ON (>=4 cm/s)
%   [rMat,pMat] = scatterEpochRates(ratNames, true, 6)         % speed filter with 6 cm/s
%
% Epochs relative to each CS:
%   Pre       = [–2  0] s
%   CS+Trace  = [ 0 0.75] s
%   US+Post   = [0.75  2] s
%   Whole     = [ 0   2] s
% Each cell’s epoch rate is normalized by its session‐wide mean rate
% (total spikes / total recording time across the 3 days).
%
% If useVelFilter==true, spikes are counted only when speed>=velThresh and
% the denominator is the *running* time inside the window.

if nargin < 2 || isempty(useVelFilter), useVelFilter = false; end
if nargin < 3 || isempty(velThresh),    velThresh    = 4;     end

% Define epochs & labels
epochWindows = [ -2.00   0.00;    ... % Pre
                  0.00   0.75;    ... % CS+Trace
                  0.75   2.00;    ... % US+Post
                  0.00   2.00 ];      % Whole
epochLabels  = {'Pre','CS+Trace','US+Post','Whole'};
comps        = [1 2; 1 3; 1 4; 2 3];  % which four pairs to plot

nRats = numel(ratNames);
% Prepare output
rMat = nan(nRats+1,4);
pMat = nan(nRats+1,4);

% Prepare pooled storage
allNormRates = cell(4,1);
for e = 1:4, allNormRates{e} = []; end

figure('Color','w','Position',[100 100 1200 300*(nRats+1)]);
panel = 0;

for r = 1:nRats
  rat   = evalin('base', ratNames{r});
  days  = autoDateList(rat);
  idx   = find(strcmp(days,rat.An),1);
  days  = days(max(1,idx-2):idx);

  % Per-day position & velocity (used if useVelFilter)
  dayInfo = struct('tpos',[],'tmin',[],'tmax',[], ...
                   'vt',[],'vmag',[], ...
                   'cs',[]);

  % Collect all CS times (keep their day index too)
  csAll = []; csDay = [];
  for d = 1:3
    posMat     = rat.pos.(sprintf('pos_%s',days{d}));
    tm         = posMat(:,1);
    dayInfo(d).tpos = tm;
    dayInfo(d).tmin = tm(1);
    dayInfo(d).tmax = tm(end);
    dayInfo(d).cs   = rat.CS_times.(sprintf('CS_%s',days{d}))(:);

    % velocity stream for this day
    if useVelFilter
      vel    = ca_velocity(posMat);          % [speed; time]
      dayInfo(d).vt   = vel(2,:)';
      dayInfo(d).vmag = vel(1,:)';
    end

    csAll  = [csAll;  dayInfo(d).cs];
    csDay  = [csDay;  d*ones(numel(dayInfo(d).cs),1)];
  end
  nTrials = numel(csAll);

  % First pass: build each cell's full spike list & mask (concatenate days)
  spikeTimesPerCell = {};
  maskCell = [];
  for d = 1:3
    spk  = rat.Ca_peaks.(sprintf('CA_peaks_%s',days{d}));
    keep = rat.ratemask.(sprintf('ratemask_%s',days{d}))==1;
    keep = keep(:);
    n_d  = numel(keep);
    if d == 1
      spikeTimesPerCell = cell(n_d,1);
      maskCell = keep;
    else
      % align sizes
      n_m = numel(maskCell);
      if n_d < n_m
        keep = [keep; false(n_m-n_d,1)];
      elseif n_d > n_m
        maskCell = [maskCell; false(n_d-n_m,1)];
        spikeTimesPerCell(n_m+1:n_d) = {[]};
      end
      maskCell = maskCell & keep;
    end
    for c = 1:numel(keep)
      if keep(c)
        st = spk(c,:); st = st(~isnan(st)&st>0);
        spikeTimesPerCell{c}(end+1:end+numel(st)) = st;
      end
    end
  end
  nCells = numel(spikeTimesPerCell);

  % Compute session duration (sum over days)
  sessionDur = 0;
  for d = 1:3
    tm = dayInfo(d).tpos;
    sessionDur = sessionDur + (tm(end)-tm(1));
  end

  % Session‐mean rate per cell
  totalSpikes = cellfun(@numel, spikeTimesPerCell);
  sessionRate = totalSpikes / sessionDur;   % nCells×1

  % Compute normalized epoch rates
  normRates = nan(nCells,4);
  for c = 1:nCells
    if ~maskCell(c), continue; end
    st_all = spikeTimesPerCell{c};
    if numel(st_all) < 1, continue; end

    % For each epoch, accumulate counts and time over trials & days
    for e = 1:4
      w0 = epochWindows(e,1);
      w1 = epochWindows(e,2);

      totalCount = 0;
      totalTime  = 0;

      % iterate trials with their day association
      for tIdx = 1:nTrials
        d  = csDay(tIdx);
        t0 = csAll(tIdx) + w0;
        t1 = csAll(tIdx) + w1;

        % Skip trials outside day bounds (shouldn't happen)
        if t1 <= dayInfo(d).tmin || t0 >= dayInfo(d).tmax, continue; end

        if ~useVelFilter
          % No speed filter: count all spikes in [t0,t1); denom is full window length
          totalCount = totalCount + sum(st_all>=t0 & st_all<t1);
          totalTime  = totalTime  + (w1 - w0);
        else
          % Speed filter: count only spikes when v>=velThresh; denom = running time
          vt   = dayInfo(d).vt;   vmag = dayInfo(d).vmag;

          % spikes in the window with speed >= velThresh
          if ~isempty(vt)
            spd_at_spk = interp1(vt, vmag, st_all, 'linear', 'extrap');
          else
            spd_at_spk = -inf(size(st_all)); % no velocity => exclude
          end
          inWin = (st_all>=t0 & st_all<t1);
          totalCount = totalCount + sum(inWin & (spd_at_spk>=velThresh));

          % running time inside the window (integrate on velocity timebase)
          if numel(vt) >= 2
            % find samples within [t0,t1)
            use = (vt>=t0 & vt<t1) & isfinite(vmag);
            if any(use)
              vt_use = vt(use);
              vm_use = vmag(use) >= velThresh;
              % duration per sample: forward diff with last dt replicated
              dt = [diff(vt_use); median(diff(vt_use),'omitnan')];
              dt(~isfinite(dt) | dt<=0) = 0;
              totalTime = totalTime + sum(dt(vm_use));
            end
          end
        end
      end

      if totalTime > 0
        rate = totalCount / totalTime;
      else
        % fall back to original denominator if no vel samples (rare)
        if ~useVelFilter
          rate = totalCount / ((w1 - w0) * nTrials);
        else
          rate = NaN; % no running time observed in these windows
        end
      end

      if sessionRate(c) > 0 && isfinite(rate)
        normRates(c,e) = rate / sessionRate(c);
      end
    end
  end

  % Append to pooled
  for e = 1:4
    allNormRates{e} = [allNormRates{e}; normRates(maskCell,e)];
  end

  % Plot each rat
  for cc = 1:4
    panel = panel + 1;
    subplot(nRats+1,4,panel); hold on

    x = normRates(maskCell, comps(cc,1));
    y = normRates(maskCell, comps(cc,2));

    scatter(x,y,20,[0 0.4470 0.7410],'filled');

    % best‐fit line
    if all(~isnan(x)) && all(~isnan(y)) && numel(x)>1
      b = polyfit(x,y,1);
      xx = linspace(0,max([x;y],[],'omitnan'),100);
      yy = polyval(b,xx);
      plot(xx,yy,'--k','LineWidth',1);
    end

    % correlation
    [R,P] = corr(x,y,'Rows','complete');
    rMat(r,cc) = R;
    pMat(r,cc) = P;

    xlabel(epochLabels{comps(cc,1)}); ylabel(epochLabels{comps(cc,2)});
    ttl = sprintf('%s (r=%.2f, p=%.3f)', ratNames{r}, R, P);
    if useVelFilter, ttl = [ttl '  | v\geq' num2str(velThresh) ' cm/s']; end
    title(ttl); axis square tight
  end
end

% Pooled "All rats" row
for cc = 1:4
  panel = panel + 1;
  subplot(nRats+1,4,panel); hold on

  x = allNormRates{comps(cc,1)};
  y = allNormRates{comps(cc,2)};

  scatter(x,y,15,[0.8500 0.3250 0.0980],'filled');

  if all(~isnan(x)) && all(~isnan(y)) && numel(x)>1
    b = polyfit(x,y,1);
    xx = linspace(0,max([x;y],[],'omitnan'),100);
    yy = polyval(b,xx);
    plot(xx,yy,'--k','LineWidth',1);
  end

  [R,P] = corr(x,y,'Rows','complete');
  rMat(nRats+1,cc) = R;
  pMat(nRats+1,cc) = P;

  xlabel(epochLabels{comps(cc,1)}); ylabel(epochLabels{comps(cc,2)});
  ttl = sprintf('All rats (r=%.2f, p=%.3f)', R, P);
  if useVelFilter, ttl = [ttl '  | v\geq' num2str(velThresh) ' cm/s']; end
  title(ttl); axis square tight
end

if useVelFilter
  sgtitle(sprintf('Normalized epoch-rate correlations (speed-filtered, v\\geq %.1f cm/s)', velThresh));
else
  sgtitle('Normalized epoch-rate correlations (no speed filter)');
end
end
