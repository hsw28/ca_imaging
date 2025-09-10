function epoch_populationSummary(ratNames, varargin)

  % epoch_populationSummary  Per-epoch delta histograms (no clustering) + summary
  %
  % Usage:
  %   R = epochModulationHists({'rat0222','rat0307','rat0313','rat0314','rat0816'});
  %
  % Options:
  %   'EpochEdges'  : [-15 0 0.25 0.75 0.85 2.00]
  %   'MinTrialSpk' : 5     % min total spikes per cell in [0 2] across trials
  %   'MinBaseSpk'  : 5     % min total spikes per cell in [-15 0] across trials
  %   'Eps'         : 1e-3  % stabilizer in delta = (f-b)/(f+b+eps)
  %   'Bins'        : -1:(1/7.5):1
  %   'Plot'        : true

  p = inputParser;
  addParameter(p,'EpochEdges',[-2 0 0.25 0.75 0.85 2.00]);
  addParameter(p,'MinTrialSpk',0);
  addParameter(p,'MinBaseSpk',0);
  addParameter(p,'Eps',1e-5);
  addParameter(p,'Bins',-1:0.1:1);
  addParameter(p,'Plot',true);
  parse(p,varargin{:});
  E   = p.Results.EpochEdges;
  Ktr = p.Results.MinTrialSpk;
  Kba = p.Results.MinBaseSpk;
  eps0= p.Results.Eps;
  bins= p.Results.Bins;
  dop = p.Results.Plot;

  labels = {'CS','Trace','US','Post'};
  nE = 4;

  % ---- per-rat delta arrays (for rat-weighted summaries) ----
  Delta_by_rat = cell(numel(ratNames),1);

  for r = 1:numel(ratNames)
      rat   = evalin('base', ratNames{r});
      dates = autoDateList(rat);
      iAn   = find(strcmp(dates,rat.An),1);
      days  = dates(max(1,iAn-2):iAn);

      Delta = [];   % cells × 4
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
                  fr(e) = cnt / max(dur, eps);
                  if e==1, cntBaseline = cntBaseline + cnt; else, cntTrial = cntTrial + cnt; end
              end

              if cntTrial < Ktr && cntBaseline < Kba, continue, end %make this and or or

              b = fr(1);
              d = (fr(2:end) - b) ./ (fr(2:end) + b + eps0);  % 1×4, bounded in [-1,1]
              d = max(min(d,1),-1);                            % numerical clip
              Delta = [Delta; d]; %#ok<AGROW>
          end
      end
      Delta_by_rat{r} = Delta;
        sum(~isnan(Delta))
  end

  sum(~isnan(Delta))

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

  % ---- optional plotting (4 hist panels + pos/neg summary) ----
  if dop
      figure('Color','w','Position',[100 100 1200 420]);
      tl = tiledlayout(1,5,'Padding','compact','TileSpacing','compact');

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
      yline(0,'k'); xlim([0.5 4.5]); ylim([-1 1]);
      xticks(1:nE); xticklabels(labels);
      ylabel('Fraction (↑ pos, ↓ neg)'); title('Dominant-sign summary');
      box off; legend([bh1 bh2],{'Positive','Negative'},'Location','northoutside');
      title(tl,'Per-epoch modulation (delta): distributions + summary');
  end

  % ---- return struct for caption stats ----
  R.labels    = labels;
  R.Delta     = Delta_all;        % pooled cells
  R.fracPos   = fracPos;
  R.fracNeg   = fracNeg;
  R.mediandelta   = medDelta;
  R.nCells    = size(Delta_all,1);
  R.nRats     = numel(ratNames);
  end
