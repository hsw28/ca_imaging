function [rMat,pMat] = scatterEpochRates(ratNames)
% scatterEpochRates   Scatterplots of normalized epoch‐rate correlations
%
%   [rMat,pMat] = scatterEpochRates({'rat0222','rat0307',...})
%
% For each cell we compute its CS‐aligned mean rate in four windows:
%   Pre       = [–2 0] s before CS
%   CS+Trace  = [ 0 0.75] s after CS
%   US+Post   = [0.75 2] s after CS
%   Whole     = [ 0 2] s after CS
% Each cell’s epoch rate is normalized by its session‐wide mean rate
% (total spikes / total recording time over all 3 days).  Panels 1–4
% show the four pairwise comparisons (Pre vs CS+Trace, Pre vs US+Post,
% Pre vs Whole, CS+Trace vs US+Post) for each rat; the bottom row pools
% all rats.  Returns rMat and pMat of Pearson’s r and p‐values.

% Define epochs & labels
epochWindows = [ -2.00   0.00;    ... % Pre
                  0.00   0.75;    ... % CS+Trace
                  0.75   2.00;    ... % US+Post
                  0.00   2.00 ];     % Whole
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
  days  = days(idx-2:idx);

  % Collect all CS times
  csAll = [];
  for d = 1:3
    csAll = [csAll; rat.CS_times.(sprintf('CS_%s',days{d}))(:)];
  end
  nTrials = numel(csAll);

  % First pass: build each cell's full spike list & mask
  spikeTimesPerCell = {};
  maskCell = [];
  for d = 1:3
    spk = rat.Ca_peaks.(sprintf('CA_peaks_%s',days{d}));
    keep = rat.ratemask.(sprintf('ratemask_%s',days{d}))==1;
    keep = keep(:);
    n_d = numel(keep);
    if d == 1
      % initialize on first day
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
    for c = 1:n_d
      if keep(c)
        st = spk(c,:);
        st = st(~isnan(st)&st>0);
        spikeTimesPerCell{c}(end+1:end+numel(st)) = st;
      end
    end
  end
  nCells = numel(spikeTimesPerCell);

  % Compute session duration
  sessionDur = 0;
  for d = 1:3
    posMat = rat.pos.(sprintf('pos_%s',days{d}));
    tm     = posMat(:,1);
    sessionDur = sessionDur + (tm(end)-tm(1));
  end

  % Session‐mean rate per cell
  totalSpikes = cellfun(@numel, spikeTimesPerCell);
  sessionRate = totalSpikes / sessionDur;   % nCells×1

  % Compute normalized epoch rates
  normRates = nan(nCells,4);
  for c = 1:nCells
    if ~maskCell(c), continue; end
    st = spikeTimesPerCell{c};
    for e = 1:4
      w0 = epochWindows(e,1);
      w1 = epochWindows(e,2);
      cnt = 0;
      for t = 1:nTrials
        t0 = csAll(t) + w0;
        t1 = csAll(t) + w1;
        cnt = cnt + sum(st>=t0 & st<t1);
      end
      rate = cnt / ((w1-w0)*nTrials);
      normRates(c,e) = rate / sessionRate(c);
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
    b = polyfit(x,y,1);
    xx = linspace(0,max([x;y]),100);
    yy = polyval(b,xx);
    plot(xx,yy,'--k','LineWidth',1);

    % correlation
    [R,P] = corr(x,y,'Rows','complete');
    rMat(r,cc) = R;
    pMat(r,cc) = P;

    xlabel(epochLabels{comps(cc,1)});
    ylabel(epochLabels{comps(cc,2)});
    title(sprintf('%s (r=%.2f, p=%.3f)', ratNames{r},R,P));
    axis square tight
  end
end

% Pooled "All rats" row
for cc = 1:4
  panel = panel + 1;
  subplot(nRats+1,4,panel); hold on

  x = allNormRates{comps(cc,1)};
  y = allNormRates{comps(cc,2)};

  scatter(x,y,15,[0.8500 0.3250 0.0980],'filled');

  b = polyfit(x,y,1);
  xx = linspace(0,max([x;y]),100);
  yy = polyval(b,xx);
  plot(xx,yy,'--k','LineWidth',1);

  [R,P] = corr(x,y,'Rows','complete');
  rMat(nRats+1,cc) = R;
  pMat(nRats+1,cc) = P;

  xlabel(epochLabels{comps(cc,1)});
  ylabel(epochLabels{comps(cc,2)});
  title(sprintf('All rats (r=%.2f, p=%.3f)',R,P));
  axis square tight
end

sgtitle('Normalized epoch‐rate correlations');
end
