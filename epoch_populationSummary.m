function R = epoch_populationSummary(ratNames, varargin)

% epoch_populationSummary  Per-epoch delta histograms (no clustering) + summary
%
% Usage:
%   R = epoch_populationSummary({'rat0222','rat0307','rat0313','rat0314','rat0816'});
%
% Options:
%   'EpochEdges'  : [-15 0 0.25 0.75 0.85 2.00]   % [baseline_start baseline_end epoch2_end ...]
%   'MinTrialSpk' : 5     % min total spikes per cell in [0 2] across trials
%   'MinBaseSpk'  : 5     % min total spikes per cell in [-15 0] across trials
%   'Eps'         : 1e-3  % stabilizer in delta = (f-b)/(f+b+eps)
%   'Bins'        : -1:(1/7.5):1
%   'Plot'        : true
%   'Labels'      : []    % optional cell array of epoch labels, length must equal (#EpochEdges - 2)
%
% Notes:
%   - The code treats the first bin in EpochEdges as baseline [E(1) E(2)],
%     and all subsequent edges define the trial epochs.

  p = inputParser;
  %addParameter(p,'EpochEdges',[-2 0 0.25 0.75 0.85 2.00]);
  addParameter(p,'EpochEdges',[-20 0 0.85 2.00]);
  addParameter(p,'MinTrialSpk',0);
  addParameter(p,'MinBaseSpk',0);
  addParameter(p,'Eps',1e-5);
  addParameter(p,'Bins',-1:0.1:1);
  addParameter(p,'Plot',true);
  addParameter(p,'Labels',[],@(x) isempty(x) || iscellstr(x) || all(cellfun(@isstring,x))); % validate later against nE
  parse(p,varargin{:});

  E     = p.Results.EpochEdges;
  Ktr   = p.Results.MinTrialSpk;
  Kba   = p.Results.MinBaseSpk;
  eps0  = p.Results.Eps;
  bins  = p.Results.Bins;
  dop   = p.Results.Plot;

  % ----- Change #1: dynamic number of trial epochs -----
  % First interval is baseline [E(1), E(2)]; all others are trial epochs.
  nE = numel(E) - 2;   % number of trial epochs

  % ----- Change #2: dynamic labels with optional override -----
  labels = p.Results.Labels;
  if isempty(labels)
      if nE == 4
          labels = {'CS','Trace','US','Post'};
      else
          labels = arrayfun(@(k) sprintf('Epoch %d',k), 1:nE, 'uni', 0);
      end
  else
      if numel(labels) ~= nE
          error('Length of Labels (%d) must match number of epochs inferred from EpochEdges (%d).', numel(labels), nE);
      end
      labels = cellfun(@char, labels, 'uni', 0); % normalize possible string -> char
  end

  % ---- per-rat delta arrays (for rat-weighted summaries) ----
  Delta_by_rat = cell(numel(ratNames),1);

  for r = 1:numel(ratNames)
      rat   = evalin('base', ratNames{r});
      dates = autoDateList(rat);
      iAn   = find(strcmp(dates,rat.An),1);
      days  = dates(max(1,iAn-2):iAn);

      Delta = [];   % cells × nE
      for d = 1:numel(days)
          D  = days{d};
          S  = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));
          CS = rat.CS_times.(sprintf('CS_%s',D));
          M  = rat.ratemask.(sprintf('ratemask_%s',D))==1;

          nC = size(S,1);
          for c = 1:nC
              if ~M(c), continue, end
              st = S(c,:); st = st(~isnan(st) & st>0);
              if isempty(st), continue, end

              fr = zeros(1,numel(E)-1);
              cntBaseline = 0; cntTrial = 0;

              for e = 1:numel(E)-1
                  cnt = 0;
                  for t = 1:numel(CS)
                      t0 = CS(t)+E(e);
                      t1 = CS(t)+E(e+1);
                      cnt = cnt + sum(st>=t0 & st<t1);
                  end
                  dur   = (E(e+1)-E(e)) * numel(CS);
                  fr(e) = cnt / max(dur, eps); % use machine eps to avoid /0
                  if e==1
                      cntBaseline = cntBaseline + cnt;
                  else
                      cntTrial    = cntTrial + cnt;
                  end
              end

              % original gate logic preserved (no change #3)
              if cntTrial < Ktr && cntBaseline < Kba
                  continue
              end

              b = fr(1);
              dlt = (fr(2:end) - b) ./ (fr(2:end) + b + eps0);  % 1×nE, bounded in [-1,1]
              dlt = max(min(dlt,1),-1);                          % numerical clip
              Delta = [Delta; dlt]; %#ok<AGROW>
          end
      end
      Delta_by_rat{r} = Delta;
  end

  % ---- pooled & rat-weighted summaries ----
  % pooled cells (for a reference histogram)
  Delta_all = cat(1, Delta_by_rat{:});

  % fraction pos/neg within each epoch, rat-weighted
  fracPos = zeros(1,nE); fracNeg = zeros(1,nE);
  medDelta= zeros(1,nE);

  for e = 1:nE
      medDelta(e) = median(Delta_all(:,e),'omitnan');

      pos_r = nan(numel(ratNames),1);
      neg_r = nan(numel(ratNames),1);
      for r = 1:numel(ratNames)
          v = Delta_by_rat{r};
          if isempty(v), continue, end
          pos_r(r) = mean(v(:,e) > 0,'omitnan');
          neg_r(r) = mean(v(:,e) < 0,'omitnan');
      end
      fracPos(e) = mean(pos_r,'omitnan');
      fracNeg(e) = mean(neg_r,'omitnan');
  end

  % ---- optional plotting (nE hist panels + pos/neg summary) ----
  if dop
      figure('Color','w','Position',[100 100 1200 420]);
      tl = tiledlayout(1, nE+1, 'Padding','compact','TileSpacing','compact');

      for e = 1:nE
          nexttile; hold on
          histogram(Delta_all(:,e), bins, 'Normalization','probability');
          xline(0,'k-'); yline(0,'k-');
          xlim([bins(1) bins(end)]); ylim([0 0.30]);
          xlabel(sprintf('%s delta',labels{e})); ylabel('Fraction of cells');
          title(labels{e});
          box off
      end

      nexttile; hold on
      bh1 = bar(1:nE,  fracPos, 0.6, 'FaceColor',[0.2 0.5 0.9], 'EdgeColor','none');
      bh2 = bar(1:nE, -fracNeg, 0.6, 'FaceColor',[0.2 0.5 0.9], 'EdgeColor','none', 'FaceAlpha',0.35);
      yline(0,'k'); xlim([0.5 nE+0.5]); ylim([-1 1]);
      xticks(1:nE); xticklabels(labels);
      ylabel('Fraction (↑ pos, ↓ neg)');
      box off; legend([bh1 bh2],{'Positive','Negative'},'Location','northoutside');
      title(tl,'Per-epoch modulation (delta): distributions + summary');
  end

  % ---- return struct for caption stats ----
  R.labels      = labels;
  R.Delta       = Delta_all;        % pooled cells
  R.fracPos     = fracPos;
  R.fracNeg     = fracNeg;
  R.mediandelta = medDelta;
  R.nCells      = size(Delta_all,1);
  R.nRats       = numel(ratNames);
end
