function plotTraceFRDist(useVelFilter, velThresh, applyTo)
% plotTraceFRDist  Compare trace‐ vs non‐trial‐FR across cells, plus corr & AUC
%
% Velocity filter (optional):
%   useVelFilter (default false) — if true, only count spikes when speed ≥ velThresh
%   velThresh    (default 4 cm/s)
%   applyTo      (default 'both') — 'trace' | 'nontrial' | 'both'
%
%  • Histograms + rank‐sum test
%  • Across‐cell Pearson & Spearman R (last dataset computed)
%  • ROC curve & AUC

if nargin < 1 || isempty(useVelFilter), useVelFilter = false; end
if nargin < 2 || isempty(velThresh),    velThresh    = 4;     end
if nargin < 3 || isempty(applyTo),      applyTo      = 'nontrial'; end
applyTo = validatestring(applyTo, {'trace','nontrial','both'});

ratNames  = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
win       = [0 2];        % CS trace window (s)
minSpikes = 0;            % min total spikes in CS windows per cell

%— aggregate per‐cell rates —%
FRt_all = [];
FRr_all = [];

for r = 1:numel(ratNames)
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);
    for d = 1:3
        spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',days{d}));
        posMat   = rat.pos.(sprintf('pos_%s',days{d}));
        csTimes  = rat.CS_times.(sprintf('CS_%s',days{d}));
        ratemask = rat.ratemask.(sprintf('ratemask_%s',days{d}));

        keep = (ratemask == 1);
        spk  = spk(keep,:);

        % flags for this run
        useTrace   = useVelFilter && (strcmp(applyTo,'trace')   || strcmp(applyTo,'both'));
        useNontrial= useVelFilter && (strcmp(applyTo,'nontrial')|| strcmp(applyTo,'both'));

        [FRt, FRr] = popVecSim(spk, posMat, csTimes, win, minSpikes, useTrace, useNontrial, velThresh);

        % collapse to **per‐cell mean** across trials (omit NaNs from vel gating)
        FRt_mean = mean(FRt, 1, 'omitnan');
        FRr_mean = mean(FRr, 1, 'omitnan');

        FRt_all = [FRt_all, FRt_mean(:)'];   %#ok<AGROW>
        FRr_all = [FRr_all, FRr_mean(:)'];   %#ok<AGROW>

        % Keep FRt/FRr from the **last dataset** for similarity distributions below
        FRt_last = FRt; %#ok<NASGU>
        FRr_last = FRr; %#ok<NASGU>
    end
end

% ==================== Trial / Random similarity (last dataset) ====================
% (kept as‐is stylistically, but robust to NaNs from vel gating)
FRt = FRt_last; FRr = FRr_last; %#ok<NASGU,ASGLU>  % for clarity with original code

Ctt = corr(FRt','Rows','pairwise');   % nTrials×nTrials
idxUT = triu(true(size(Ctt)),1);
TT_vals = Ctt(idxUT);

Crr = corr(FRr','Rows','pairwise');
RR_vals = Crr(idxUT);

Ctr = corr(FRt', FRr','Rows','pairwise');  % nTrials×nTrials
TR_vals = Ctr(:);

mean_TT = mean(TT_vals,'omitnan');
mean_TR = mean(TR_vals,'omitnan');
mean_RR = mean(RR_vals,'omitnan');
sem_TT  = std(TT_vals,'omitnan')/sqrt(nnz(~isnan(TT_vals)));
sem_TR  = std(TR_vals,'omitnan')/sqrt(nnz(~isnan(TR_vals)));
sem_RR  = std(RR_vals,'omitnan')/sqrt(nnz(~isnan(RR_vals)));
means = [mean_TT, mean_TR, mean_RR];
sems  = [sem_TT,  sem_TR,  sem_RR];
stds  = [std(TT_vals,'omitnan'), std(TR_vals,'omitnan'), std(RR_vals,'omitnan')];

fprintf('Mean TT = %.3f, TR = %.3f, RR = %.3f\n', mean_TT, mean_TR, mean_RR)

cols = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

figure('Color','w','Position',[300 300 900 400]);
subplot(1,2,1); hold on;
edges = linspace(-1,1,50);
histogram(TT_vals, edges, 'Normalization','probability','FaceColor',cols(1,:), 'FaceAlpha',.6);
histogram(TR_vals, edges, 'Normalization','probability','FaceColor',cols(2,:), 'FaceAlpha',.6);
histogram(RR_vals, edges, 'Normalization','probability','FaceColor',cols(3,:), 'FaceAlpha',.6);
legend('TT','TR','RR','Location','Best');
xlabel('Pearson similarity'); ylabel('Probability');
title('Trace vs Random pop‐vector similarities');

subplot(1,2,2); hold on;
b = bar(1:3, means, 'FaceColor','flat'); b.CData = cols;
errorbar(1:3, means, sems, 'k.', 'LineWidth',1.5);
xticks(1:3); xticklabels({'TT','TR','RR'}); ylabel('Mean similarity');
title('Mean ± SEM pop‐vector similarity'); xlim([0.5 3.5]);

% Unpaired tests
p_tt_tr = ranksum(TT_vals, TR_vals);
p_tt_rr = ranksum(TT_vals, RR_vals);
ymax = max(means + sems); dy = 0.05*ymax; y1 = ymax + dy; y2 = ymax + 2*dy;
plot([1 2],[y1 y1],'-k','LineWidth',1.5); if p_tt_tr<.05, text(1.5,y1+0.02,'*','horiz','center','FontSize',16); end
plot([1 3],[y2 y2],'-k','LineWidth',1.5); if p_tt_rr<.05, text(2.0,y2+0.02,'*','horiz','center','FontSize',16); end
fprintf('TT vs TR ranksum p = %.3g\n', p_tt_tr);
fprintf('TT vs RR ranksum p = %.3g\n', p_tt_rr);

% ==================== FR distributions & ROC ====================
figure('Color','w','Position',[300 300 900 400]);
subplot(1,2,1); hold on;
finiteAll = [FRt_all(isfinite(FRt_all)), FRr_all(isfinite(FRr_all))];
if isempty(finiteAll), finiteAll = 0; end
edges = linspace(0, max(finiteAll), 50);
histogram(FRt_all, 'Normalization','probability','FaceAlpha',.6,'BinEdges',edges);
histogram(FRr_all, 'Normalization','probability','FaceAlpha',.6,'BinEdges',edges);
set(gca,'YScale','log')
xlabel('Firing rate (Hz)'); ylabel('Probability'); legend('Trace','Non-trial','Location','Best');
title('Trace vs Non-trial FR Distributions');

nT = nnz(isfinite(FRt_all)); nR = nnz(isfinite(FRr_all));
mT = mean(FRt_all,'omitnan'); mR = mean(FRr_all,'omitnan');
semT = std(FRt_all,'omitnan')/sqrt(max(nT,1));
semR = std(FRr_all,'omitnan')/sqrt(max(nR,1));
fprintf('\nTrace FR:     mean=%.2f±%.2f Hz\n', mT, std(FRt_all,'omitnan'));
fprintf('Non-trial FR: mean=%.2f±%.2f Hz\n', mR, std(FRr_all,'omitnan'));

% two-sample (unpaired), paired (lengths may differ so show unpaired only if needed)
[~,p_unpaired] = ttest2(FRt_all(~isnan(FRt_all)), FRr_all(~isnan(FRr_all)));
fprintf('unpaired t test p = %.3g\n', p_unpaired);
try
  [~,p_paired] = ttest(FRt_all, FRr_all);
  fprintf('paired t test p = %.3g\n', p_paired);
catch
  fprintf('paired t test skipped (length mismatch)\n');
end
[~,p_ks] = kstest2(FRt_all(~isnan(FRt_all)), FRr_all(~isnan(FRr_all)));
fprintf('kstest2 p = %.3g\n', p_ks);
p_rank = ranksum(FRt_all(~isnan(FRt_all)), FRr_all(~isnan(FRr_all)));
fprintf('Wilcoxon rank-sum p = %.3g\n', p_rank);

subplot(1,2,2); hold on;
bar([1 2], [mT mR], 'FaceColor','flat'); set(gca,'XTick',[1 2],'XTickLabel',{'Trace','Non-trial'});
errorbar([1 2], [mT mR], [semT semR], 'k.', 'LineWidth',1.5);
ylabel('Mean firing rate (Hz)'); title('Mean FR ± SEM'); xlim([0.5 2.5]);
if p_unpaired < 0.05
  yl = ylim; plot([1 2],[yl(2)*0.95 yl(2)*0.95],'k-','LineWidth',1); text(1.5,yl(2)*0.97,'*','FontSize',20,'HorizontalAlignment','center');
end

% ROC & AUC (drop NaNs)
labels = [ ones(1, nnz(isfinite(FRt_all))) , zeros(1, nnz(isfinite(FRr_all))) ];
scores = [ FRt_all(isfinite(FRt_all)) ,      FRr_all(isfinite(FRr_all))      ];
[Xroc,Yroc,~,AUC] = perfcurve(labels(:), scores(:), 1);
fprintf('Trace vs non-trial FR AUC = %.3f\n', AUC);

figure('Color','w','Position',[350 350 500 400]);
plot(Xroc, Yroc, 'LineWidth',2); hold on; plot([0 1],[0 1],'k--');
xlabel('False positive rate'); ylabel('True positive rate');
title(sprintf('ROC curve (AUC = %.3f)', AUC)); axis square;

% Permutation test for AUC
nShuff = 500; permAUCs = nan(nShuff,1);
fprintf('\nRunning permutation test for AUC (n=%d)...\n', nShuff);
for s = 1:nShuff
    permLabels = labels(randperm(length(labels)));
    [~,~,~,permAUCs(s)] = perfcurve(permLabels, scores, 1);
end
figure('Color','w','Position',[400 400 500 400]);
histogram(permAUCs, 30, 'FaceColor',[.7 .7 .7], 'EdgeColor','k', 'Normalization', 'probability');
hold on; xline(AUC, 'r--', 'LineWidth', 2);
xlabel('AUC'); ylabel('Frequency');
title(sprintf('Permutation Test for AUC (Actual = %.3f)', AUC));
permP = mean(permAUCs >= AUC);
fprintf('Permutation p-value for AUC = %.4f\n', permP);
end


function [FRt,FRr] = popVecSim(spikeCell, pos, csTimes, win, minSpikes, useTraceFilter, useNontrialFilter, velThresh)
% popVecSim  Build FR matrices for CS vs. random windows
%   FRt: per-trial trace window FRs; FRr: per-trial random-window FRs.
%   If useTraceFilter or useNontrialFilter are true, only count spikes
%   when speed ≥ velThresh and divide by running time in the window.

  if ~iscell(spikeCell)
    [nCells, ~] = size(spikeCell);
    tmp = cell(nCells,1);
    for c = 1:nCells
      st = spikeCell(c,:);  st = st(~isnan(st) & st>0);
      tmp{c} = st(:);
    end
    spikeCell = tmp;
  end

  ts        = pos(:,1);
  nTrials   = numel(csTimes);
  nCells    = numel(spikeCell);
  windowLen = diff(win);
  dt        = median(diff(ts));

  % velocity (on its own timebase)
  vel  = ca_velocity(pos);    % [speed; time]
  vt   = vel(2,:)';
  vmag = vel(1,:)';

  % mask non-CS samples on ts (for choosing random windows)
  maskCS = false(size(ts));
  for t = 1:nTrials
    maskCS = maskCS | (ts >= csTimes(t)+win(1) & ts < csTimes(t)+win(2));
  end
  outs = ts(~maskCS);
  bins = max(1, round(windowLen / dt));
  maxStart = numel(outs) - bins + 1;
  if maxStart < 1
      error('Not enough non-CS timepoints to place random windows.');
  end

  FRt = nan(nTrials, nCells);
  FRr = nan(nTrials, nCells);

  for t = 1:nTrials
    t0 = csTimes(t) + win(1);
    t1 = t0 + windowLen;

    % Precompute running time in trace window (if needed)
    if useTraceFilter
      [~, runT_trace] = countSpikesAndRunTime([], vt, vmag, t0, t1, velThresh);
    end

    for c = 1:nCells
      st = spikeCell{c};

      % ---- Trace window ----
      if ~useTraceFilter
        FRt(t,c) = sum(st>=t0 & st<t1) / windowLen;
      else
        [cntT, runT] = countSpikesAndRunTime(st, vt, vmag, t0, t1, velThresh);
        FRt(t,c) = (runT>0) * (cntT / runT);  % NaN if runT==0
        if runT==0, FRt(t,c) = NaN; end
      end

      % ---- Random window ----
      sIdx = randi(maxStart);  % random start in non-CS time
      r0   = outs(sIdx);
      r1   = r0 + windowLen;

      if ~useNontrialFilter
        FRr(t,c) = sum(st>=r0 & st<r1) / windowLen;
      else
        [cntR, runR] = countSpikesAndRunTime(st, vt, vmag, r0, r1, velThresh);
        FRr(t,c) = (runR>0) * (cntR / runR);
        if runR==0, FRr(t,c) = NaN; end
      end
    end
  end

  % filter out cells with too few spikes in CS (based on counts, not rate)
  totalCounts = nansum(FRt,1) * windowLen;   % approximate counts over trials
  keepCells   = totalCounts >= minSpikes;
  FRt = FRt(:, keepCells);
  FRr = FRr(:, keepCells);
end


function [cnt, runT] = countSpikesAndRunTime(st, vt, vmag, t0, t1, velThresh)
% Count spikes with speed ≥ velThresh and compute running time within [t0,t1)

  % spikes
  if isempty(st)
    cnt = NaN;  % if asked only for runT (caller may ignore cnt)
  else
    spd_at_spk = interp1(vt, vmag, st, 'linear', 'extrap');
    inWin = (st>=t0 & st<t1) & (spd_at_spk>=velThresh);
    cnt = sum(inWin);
  end

  % running time
  use = (vt>=t0 & vt<t1) & isfinite(vmag);
  if ~any(use)
    runT = 0;
    return;
  end
  vt_use = vt(use);
  vm_use = vmag(use) >= velThresh;
  if numel(vt_use) < 2
    runT = 0;
    return;
  end
  dt = [diff(vt_use); median(diff(vt_use),'omitnan')];
  dt(~isfinite(dt) | dt<=0) = 0;
  runT = sum(dt(vm_use));
end
