function plotTraceFRDist(useVelFilter, velThresh, applyTo)
% plotTraceFRDist  Compare trace vs non-trial FR across cells + ROC/AUC
%
% Figure 2 (2 cols x 6 rows): each rat (5 rows) + pooled (last row)
%   Col 1: bar (Trace, Nontrial, Nontrial vel>=thr) with SEM + significance
%   Col 2: histogram (Trace vs ONLY one Nontrial):
%          - if useVelFilter==false: Trace vs Nontrial (no filter)
%          - if useVelFilter==true : Trace vs Nontrial (vel>=thr)
%
% Similarity/ROC sections above Figure 2 still honor useVelFilter/applyTo.

if nargin < 1 || isempty(useVelFilter), useVelFilter = false; end
if nargin < 2 || isempty(velThresh),    velThresh    = 4;     end
if nargin < 3 || isempty(applyTo),      applyTo      = 'nontrial'; end
applyTo = validatestring(applyTo, {'trace','nontrial','both'});

ratNames  = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
win       = [0 2];        % CS trace window (s)
minSpikes = 0;

FRt_all = [];
FRr_all = [];

% Per-rat storage for Figure 2
nR = numel(ratNames);
FRt_meancell  = cell(nR,1);   % Trace (no filter)
FRr0_meancell = cell(nR,1);   % Nontrial (no filter)
FRrV_meancell = cell(nR,1);   % Nontrial (vel filter)

cols = [0    0.4470 0.7410;   % Trace
        0.8500 0.3250 0.0980; % Nontrial (no filter)
        0.4660 0.6740 0.1880];% Nontrial (vel filter)

for r = 1:nR
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);

    trace_means   = [];
    nontrial0_means = [];
    nontrialV_means = [];

    for d = 1:3
        spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',days{d}));
        posMat   = rat.pos.(sprintf('pos_%s',days{d}));
        csTimes  = rat.CS_times.(sprintf('CS_%s',days{d}));
        ratemask = rat.ratemask.(sprintf('ratemask_%s',days{d}));
        keep     = (ratemask == 1);
        spk      = spk(keep,:);

        % A) Path used by similarity/ROC — unchanged logic
        useTrace    = useVelFilter && (strcmp(applyTo,'trace')   || strcmp(applyTo,'both'));
        useNontrial = useVelFilter && (strcmp(applyTo,'nontrial')|| strcmp(applyTo,'both'));
        [FRt_A, FRr_A] = popVecSim(spk, posMat, csTimes, win, minSpikes, useTrace, useNontrial, velThresh);
        FRt_all = [FRt_all, mean(FRt_A,1,'omitnan')]; %#ok<AGROW>
        FRr_all = [FRr_all, mean(FRr_A,1,'omitnan')]; %#ok<AGROW>
        FRt_last = FRt_A; %#ok<NASGU>
        FRr_last = FRr_A; %#ok<NASGU>

        % B) Figure 2 data: always compute all three series
        [FRt_traceOnly, ~] = popVecSim(spk, posMat, csTimes, win, minSpikes, false, false, velThresh);
        [~, FRr_noVF]       = popVecSim(spk, posMat, csTimes, win, minSpikes, false, false, velThresh);
        [~, FRr_velF]       = popVecSim(spk, posMat, csTimes, win, minSpikes, false, true,  velThresh);

        trace_means     = [trace_means,   mean(FRt_traceOnly,1,'omitnan')]; %#ok<AGROW>
        nontrial0_means = [nontrial0_means,mean(FRr_noVF,   1,'omitnan')]; %#ok<AGROW>
        nontrialV_means = [nontrialV_means,mean(FRr_velF,   1,'omitnan')]; %#ok<AGROW>
    end

    FRt_meancell{r}  = trace_means(:);
    FRr0_meancell{r} = nontrial0_means(:);
    FRrV_meancell{r} = nontrialV_means(:);
end

% ===== Similarity summary (unchanged) =====
FRt = FRt_last; FRr = FRr_last; %#ok<NASGU,ASGLU>
Ctt = corr(FRt','Rows','pairwise'); idxUT = triu(true(size(Ctt)),1); TT_vals = Ctt(idxUT);
Crr = corr(FRr','Rows','pairwise'); RR_vals = Crr(idxUT);
Ctr = corr(FRt', FRr','Rows','pairwise');  TR_vals = Ctr(:);

mean_TT = mean(TT_vals,'omitnan'); mean_TR = mean(TR_vals,'omitnan'); mean_RR = mean(RR_vals,'omitnan');
sem_TT  = std(TT_vals,'omitnan')/sqrt(nnz(~isnan(TT_vals)));
sem_TR  = std(TR_vals,'omitnan')/sqrt(nnz(~isnan(TR_vals)));
sem_RR  = std(RR_vals,'omitnan')/sqrt(nnz(~isnan(RR_vals)));
means = [mean_TT, mean_TR, mean_RR]; sems = [sem_TT, sem_TR, sem_RR];

fprintf('Mean TT = %.3f, TR = %.3f, RR = %.3f\n', mean_TT, mean_TR, mean_RR)

figure('Color','w','Position',[300 300 900 400]);
subplot(1,2,1); hold on;
edges = linspace(-1,1,50);
histogram(TT_vals, edges, 'Normalization','probability','FaceColor',cols(1,:), 'FaceAlpha',.6);
histogram(TR_vals, edges, 'Normalization','probability','FaceColor',cols(2,:), 'FaceAlpha',.6);
histogram(RR_vals, edges, 'Normalization','probability','FaceColor',cols(3,:), 'FaceAlpha',.6);
legend('TT','TR','RR','Location','Best');
xlabel('Pearson similarity'); ylabel('Probability'); title('Trace vs Random pop-vector similarities');

subplot(1,2,2); hold on;
b = bar(1:3, means, 'FaceColor','flat'); b.CData = cols;
errorbar(1:3, means, sems, 'k.', 'LineWidth',1.5);
xticks(1:3); xticklabels({'TT','TR','RR'}); ylabel('Mean similarity');
title('Mean ± SEM pop-vector similarity'); xlim([0.5 3.5]);

p_tt_tr = ranksum(TT_vals, TR_vals); p_tt_rr = ranksum(TT_vals, RR_vals);
ymax = max(means + sems); dy = 0.05*ymax; y1 = ymax + dy; y2 = ymax + 2*dy;
plot([1 2],[y1 y1],'-k','LineWidth',1.2); if p_tt_tr<.05, text(1.5,y1+0.02,'*','horiz','center','FontSize',14); end
plot([1 3],[y2 y2],'-k','LineWidth',1.2); if p_tt_rr<.05, text(2.0,y2+0.02,'*','horiz','center','FontSize',14); end
fprintf('TT vs TR ranksum p = %.3g\n', p_tt_tr);
fprintf('TT vs RR ranksum p = %.3g\n', p_tt_rr);

% ===== Figure 2: per-rat + pooled (2 columns × 6 rows) =====
FRt_pool  = vertcat(FRt_meancell{:});
FRr0_pool = vertcat(FRr0_meancell{:});
FRrV_pool = vertcat(FRrV_meancell{:});

figure('Color','w','Position',[200 100 1100 1600]);
tiledlayout(6,2,'Padding','compact','TileSpacing','compact');

for rr = 1:6
    if rr<=nR
        tdat  = FRt_meancell{rr};
        r0dat = FRr0_meancell{rr};
        rvdat = FRrV_meancell{rr};
        rowTitle = ratNames{rr};
    else
        tdat  = FRt_pool;
        r0dat = FRr0_pool;
        rvdat = FRrV_pool;
        rowTitle = 'All rats (pooled)';
    end

    % -------- Column 1: BAR with SEM + significance lines ----------
    nexttile; hold on;
    m  = [mean(tdat,'omitnan'), mean(r0dat,'omitnan'), mean(rvdat,'omitnan')];
    n  = [nnz(isfinite(tdat)),  nnz(isfinite(r0dat)),  nnz(isfinite(rvdat))];
    se = [std(tdat,'omitnan')/max(sqrt(n(1)),1), ...
          std(r0dat,'omitnan')/max(sqrt(n(2)),1), ...
          std(rvdat,'omitnan')/max(sqrt(n(3)),1)];
    bb = bar(1:3, m, 'FaceColor','flat'); bb.CData = cols;
    errorbar(1:3, m, se, 'k.', 'LineWidth',1.2);
    xlim([0.5 3.5]); set(gca,'XTick',1:3,'XTickLabel',{'Trace','Nontrial','Nontrial (vel≥thr)'});
    ylabel('Mean FR (Hz)'); title(sprintf('%s — mean FR ± SEM', rowTitle)); box on;

    % significance: Trace vs Nontrial (no filter) and Trace vs Nontrial (vel)
    p12 = paired_or_not(tdat, r0dat);  % Trace vs Nontrial
    p13 = paired_or_not(tdat, rvdat);  % Trace vs Nontrial (vel)
    yl = ylim; ybase = yl(2);
    yA = ybase*0.94; yB = ybase*0.98;
    draw_sigline(gca, 1, 2, yA, p12);
    draw_sigline(gca, 1, 3, yB, p13);

    % -------- Column 2: HISTOGRAM (Trace vs ONE Nontrial condition) --------
    nexttile; hold on;
    if useVelFilter
        ndat = rvdat; lab = 'Nontrial (vel≥thr)'; c2 = cols(3,:);
    else
        ndat = r0dat; lab = 'Nontrial';          c2 = cols(2,:);
    end
    finiteAll = [tdat(isfinite(tdat)); ndat(isfinite(ndat))];
    if isempty(finiteAll), finiteAll = 0; end
    edges = linspace(0, max(finiteAll), 50);
    histogram(tdat, 'BinEdges',edges,'Normalization','probability','FaceColor',cols(1,:), 'FaceAlpha',.6);
    histogram(ndat, 'BinEdges',edges,'Normalization','probability','FaceColor',c2,         'FaceAlpha',.6);
    ax = gca; ax.YAxis.Exponent = 0; ytickformat('%.3f');  % no scientific notation
    xlabel('Firing rate (Hz)'); ylabel('Probability');
    legend('Trace', lab, 'Location','best'); box on;
    title(sprintf('%s — FR distributions', rowTitle));
end

% ===== ROC & AUC (unchanged across-all) =====
labels = [ ones(1, nnz(isfinite(FRt_all))) , zeros(1, nnz(isfinite(FRr_all))) ];
scores = [ FRt_all(isfinite(FRt_all)) ,      FRr_all(isfinite(FRr_all))      ];
[Xroc,Yroc,~,AUC] = perfcurve(labels(:), scores(:), 1);
fprintf('\nTrace vs non-trial FR AUC = %.3f\n', AUC);

figure('Color','w','Position',[350 350 500 400]);
plot(Xroc, Yroc, 'LineWidth',2); hold on; plot([0 1],[0 1],'k--');
xlabel('False positive rate'); ylabel('True positive rate');
title(sprintf('ROC curve (AUC = %.3f)', AUC)); axis square;

% Permutation test for AUC
nShuff = 500; permAUCs = nan(nShuff,1);
fprintf('Running permutation test for AUC (n=%d)...\n', nShuff);
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


function p = paired_or_not(x, y)
% Try paired tests when overlapping cells exist; robust fallbacks.
mask = isfinite(x) & isfinite(y);
x = x(mask); y = y(mask);
if numel(x) >= 3 && numel(y) >= 3
    try
        [~,p] = ttest(x, y);
    catch
        try
            p = signrank(x, y);
        catch
            [~,p] = ttest2(x, y);
            if isnan(p), p = ranksum(x, y); end
        end
    end
else
    % If almost no overlap, fall back to unpaired
    x = x(isfinite(x)); y = y(isfinite(y));
    try
        [~,p] = ttest2(x, y);
    catch
        p = ranksum(x, y);
    end
end
if isnan(p), p = 1; end
end

function draw_sigline(ax, x1, x2, y, p)
% Draw horizontal line between bars x1 and x2 at height y with stars.
if ~isfinite(p), p = 1; end
stars = p_to_stars(p);
plot(ax, [x1 x2], [y y], 'k-', 'LineWidth', 1.2);
text(mean([x1 x2]), y*1.01, stars, 'HorizontalAlignment','center','FontSize',12);
end

function s = p_to_stars(p)
if p < 1e-3, s = '***';
elseif p < 1e-2, s = '**';
elseif p < 0.05, s = '*';
else, s = 'n.s.';
end
end

function [FRt,FRr] = popVecSim(spikeCell, pos, csTimes, win, minSpikes, useTraceFilter, useNontrialFilter, velThresh)
% (unchanged from your version)
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

  vel  = ca_velocity(pos);    % [speed; time]
  vt   = vel(2,:)';
  vmag = vel(1,:)';

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

    for c = 1:nCells
      st = spikeCell{c};

      % Trace window
      if ~useTraceFilter
        FRt(t,c) = sum(st>=t0 & st<t1) / windowLen;
      else
        [cntT, runT] = countSpikesAndRunTime(st, vt, vmag, t0, t1, velThresh);
        FRt(t,c) = (runT>0) * (cntT / runT);
        if runT==0, FRt(t,c) = NaN; end
      end

      % Random window
      sIdx = randi(maxStart);
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

  totalCounts = nansum(FRt,1) * windowLen;   % approximate counts over trials
  keepCells   = totalCounts >= minSpikes;
  FRt = FRt(:, keepCells);
  FRr = FRr(:, keepCells);
end


function [cnt, runT] = countSpikesAndRunTime(st, vt, vmag, t0, t1, velThresh)
% Count spikes with speed ≥ velThresh and compute running time within [t0,t1)
  if isempty(st)
    cnt = NaN;
  else
    spd_at_spk = interp1(vt, vmag, st, 'linear', 'extrap');
    inWin = (st>=t0 & st<t1) & (spd_at_spk>=velThresh);
    cnt = sum(inWin);
  end

  use = (vt>=t0 & vt<t1) & isfinite(vmag);
  if ~any(use)
    runT = 0; return;
  end
  vt_use = vt(use);
  vm_use = vmag(use) >= velThresh;
  if numel(vt_use) < 2
    runT = 0; return;
  end
  dt = [diff(vt_use); median(diff(vt_use),'omitnan')];
  dt(~isfinite(dt) | dt<=0) = 0;
  runT = sum(dt(vm_use));
end
