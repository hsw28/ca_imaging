function plotTraceFRDist()
% plotTraceFRDist  Compare trace‐ vs non‐trial‐FR across cells, plus corr & AUC
%
%  • Histograms + rank‐sum test
%  • Across‐cell Pearson & Spearman R
%  • ROC curve & AUC

% USER PARAMETERS
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
        ts       = posMat(:,1);
        csTimes  = rat.CS_times.(sprintf('CS_%s',days{d}));
        ratemask = rat.ratemask.(sprintf('ratemask_%s',days{d}));

        % --- apply ratemask: keep only cells with ratemask==1 -------------
        keep = (ratemask == 1);
        spk  = spk(keep,:);

        [FRt, FRr] = popVecSim(spk, ts, csTimes, win, minSpikes);

        % collapse FRt/FRr to **per‐cell mean** across trials
        FRt_mean = mean(FRt,1);
        FRr_null = mean(FRr,1);

        FRt_all = [FRt_all, FRt_mean(:)'];   %#ok<AGROW>
        FRr_all = [FRr_all, FRr_null(:)'];  %#ok<AGROW>
    end
end

% 1) compute trial–trial correlations
Ctt = corr(FRt');   % nTrials×nTrials
idxUT = triu(true(size(Ctt)),1);
TT_vals = Ctt(idxUT);

% 2) compute random–random correlations
Crr = corr(FRr');
RR_vals = Crr(idxUT);

% 3) compute trial–random correlations
Ctr = corr(FRt', FRr');  % nTrials×nTrials
TR_vals = Ctr(:);        % all pairwise trial vs random

% Now get means (and distributions) of each:
mean_TT = mean(TT_vals);
mean_TR = mean(TR_vals);
mean_RR = mean(RR_vals);

fprintf('Mean TT = %.3f, TR = %.3f, RR = %.3f\n', mean_TT, mean_TR, mean_RR)
sem_TT = std(TT_vals)/sqrt(numel(TT_vals));
sem_TR = std(TR_vals)/sqrt(numel(TR_vals));
sem_RR = std(RR_vals)/sqrt(numel(RR_vals));
means = [mean_TT, mean_TR, mean_RR];
sems  = [sem_TT, sem_TR, sem_RR];
stds = [std(TT_vals), std(TR_vals), std(RR_vals)];

% Define the same color‐map for TT, TR, RR
cols = [
    0     0.4470 0.7410;   % blue  for TT
    0.8500 0.3250 0.0980;  % red   for TR
    0.9290 0.6940 0.1250;  % yellow for RR
];

figure('Color','w','Position',[300 300 900 400]);

% --- Subplot 1: Histograms with explicit FaceColor ---
subplot(1,2,1); hold on;
edges = linspace(-1,1,50);
histogram(TT_vals, edges, ...
    'Normalization','probability', ...
    'FaceColor',cols(1,:), 'FaceAlpha',.6);
histogram(TR_vals, edges, ...
    'Normalization','probability', ...
    'FaceColor',cols(2,:), 'FaceAlpha',.6);
histogram(RR_vals, edges, ...
    'Normalization','probability', ...
    'FaceColor',cols(3,:), 'FaceAlpha',.6);
legend('TT','TR','RR','Location','Best');
xlabel('Cosine similarity');
ylabel('Probability');
title('Trace vs Random pop‐vector similarities');

% --- Subplot 2: Bar graph with matching colors and SEM bars ---
subplot(1,2,2); hold on;
means = [mean_TT, mean_TR, mean_RR];
sems  = [sem_TT,  sem_TR,  sem_RR];
b = bar(1:3, means, 'FaceColor','flat');
b.CData = cols;  % assign each bar its color
errorbar(1:3, means, sems, 'k.', 'LineWidth',1.5);
xticks(1:3);
xticklabels({'TT','TR','RR'});
ylabel('Mean cosine similarity');
title('Mean ± SEM pop‐vector similarity');
xlim([0.5 3.5]);

means
stds

% --- Statistical tests and stars ---
% unpaired tests (you can also use ttest2 instead of ranksum if you prefer)
[p_tt_tr, ~] = ranksum(TT_vals, TR_vals);
[p_tt_rr, ~] = ranksum(TT_vals, RR_vals);

% determine star positions
ymax = max(means + sems);
dy   = 0.05 * (ymax);      % 5% of the range
y1   = ymax + dy;
y2   = ymax + 2*dy;

% TT vs TR
plot([1 2], [y1 y1], '-k', 'LineWidth',1.5);
if p_tt_tr < .05
    text(1.5, y1 + dy*0.1, '*', 'HorizontalAlignment','center','FontSize',16);
end

% TT vs RR
plot([1 3], [y2 y2], '-k', 'LineWidth',1.5);
if p_tt_rr < .05
    text(2, y2 + dy*0.1, '*', 'HorizontalAlignment','center','FontSize',16);
end

% OPTIONAL: display p‐values in console
fprintf('TT vs TR ranksum p = %.3g\n', p_tt_tr);
fprintf('TT vs RR ranksum p = %.3g\n', p_tt_rr);



% OPTIONAL: add stars for signif TT vs TR and TT vs RR
% run your stats below and then for any p<0.05 add a star


% STATISTICAL TESTS:
% Independent two‐sample t‐tests:
[p_TT_TR,~] = ttest2(TT_vals, TR_vals);
[p_TT_RR,~] = ttest2(TT_vals, RR_vals);
[p_TR_RR,~] = ttest2(TR_vals, RR_vals);

% or nonparametric rank‐sum tests:
[p2_TT_TR,~] = ranksum(TT_vals, TR_vals);
[p2_TT_RR,~] = ranksum(TT_vals, RR_vals);
[p2_TR_RR,~] = ranksum(TR_vals, RR_vals);

[a b] = kstest2(TT_vals, TR_vals);

fprintf('\nTT vs TR: ttest2 p = %.3g, ranksum p = %.3g\n', p_TT_TR, p2_TT_TR);
fprintf('TT vs RR: ttest2 p = %.3g, ranksum p = %.3g\n', p_TT_RR, p2_TT_RR);
fprintf('TR vs RR: ttest2 p = %.3g, ranksum p = %.3g\n', p_TR_RR, p2_TR_RR);
%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––%






%— HISTOGRAM & RANK‐SUM + BAR±SEM —%
figure('Color','w','Position',[300 300 900 400]);

% 1) Histogram subplot
subplot(1,2,1); hold on;
edges = linspace(0, max([FRt_all,FRr_all]), 50);
histogram(FRt_all, 'Normalization','probability', ...
          'FaceAlpha',.6, 'BinWidth', .021);
histogram(FRr_all, 'Normalization','probability', ...
          'FaceAlpha',.6, 'BinWidth', .021);
xlabel('Firing rate (Hz)');
ylabel('Probability');
title('Trace vs Non-trial FR Distributions');
legend('Trace','Non-trial','Location','Best');

% compute means & SEMs
nT    = numel(FRt_all);
nR    = numel(FRr_all);
mT    = mean(FRt_all);
mR    = mean(FRr_all);
semT  = std(FRt_all) / sqrt(nT);
semR  = std(FRr_all) / sqrt(nR);

% print stats
fprintf('\nTrace FR:     mean=%.2f±%.2f Hz\n', mT, std(FRt_all));
fprintf('Non-trial FR: mean=%.2f±%.2f Hz\n', mR, std(FRr_all));
fprintf('unpaired t test')
[h,p,ci,stats] = ttest2(FRt_all, FRr_all)
fprintf('paired t test')
[h,p,ci,stats] = ttest(FRt_all, FRr_all)
fprintf('kstest2')
[h,p] = kstest2(FRt_all, FRr_all)
fprintf('Wilcoxon rank-sum p = %.3g\n', p);

% 2) Bar + SEM subplot
subplot(1,2,2); hold on;
barHandles = bar([1 2], [mT mR], 'FaceColor','flat');
barHandles.CData(1,:) = [0 .447 .741];  % MATLAB default blue
barHandles.CData(2,:) = [0.85 .325 .098]; % MATLAB default rust
errorbar([1 2], [mT mR], [semT semR], 'k.', 'LineWidth',1.5);
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'Trace','Non-trial'});
ylabel('Mean firing rate (Hz)');
title('Mean FR ± SEM');

% add significance star if you like
if p < 0.05
    yl = ylim;
    plot([1 2], [yl(2)*0.95 yl(2)*0.95], 'k-', 'LineWidth',1);
    text(1.5, yl(2)*0.97, '*', 'FontSize',20,'HorizontalAlignment','center');
end


%— ROC & AUC —%
labels = [ ones(size(FRt_all)), zeros(size(FRr_all)) ];
scores = [ FRt_all,            FRr_all           ];
% perfcurve requires Statistics and ML Toolbox
[Xroc,Yroc,~,AUC] = perfcurve(labels(:), scores(:), 1);
fprintf('Trace vs non-trial FR AUC = %.3f\n', AUC);

% plot ROC
figure('Color','w','Position',[350 350 500 400]);
plot(Xroc, Yroc, 'LineWidth',2);
hold on; plot([0 1],[0 1],'k--');
xlabel('False positive rate'); ylabel('True positive rate');
title(sprintf('ROC curve (AUC = %.3f)', AUC));
axis square;

% --- Permutation Test for AUC ---
nShuff = 500;
permAUCs = nan(nShuff,1);

fprintf('\nRunning permutation test for AUC (n=%d)...\n', nShuff);

for s = 1:nShuff
    permLabels = labels(randperm(length(labels)));
    [~,~,~,permAUCs(s)] = perfcurve(permLabels, scores, 1);
end

% Plot histogram of shuffled AUCs
figure('Color','w','Position',[400 400 500 400]);
histogram(permAUCs, 30, 'FaceColor',[.7 .7 .7], 'EdgeColor','k', 'Normalization', 'probability');
hold on;
xline(AUC, 'r--', 'LineWidth', 2);
xlabel('AUC'); ylabel('Frequency');
title(sprintf('Permutation Test for AUC (Actual = %.3f)', AUC));
legend('Shuffled AUCs', 'Actual AUC', 'Location','best');

mean(permAUCs)
% Compute p-value
permP = mean(permAUCs >= AUC);
fprintf('Permutation p-value for AUC = %.4f\n', permP);


end


function [FRt,FRr] = popVecSim(spikeCell, ts, csTimes, win, minSpikes)
% popVecSim  Build FR matrices for CS vs. random windows
%
%   [FRt,FRr] = popVecSim(spikeCell, ts, csTimes, win, minSpikes)
%   - spikeCell: either {nCells×1} cell array of spike-time vectors,
%                or numeric [nCells×nSpikesMax] array of spike times per row
%   - ts:        [nBins×1] timestamp vector (e.g. 7.5 Hz sampling)
%   - csTimes:   [nTrials×1] CS onset times
%   - win:       [t0 t1] CS window relative to CS onset (s)
%   - minSpikes: filter threshold on total spikes per cell in CS

  %–– 1) Convert numeric matrix to cell array of true spike times ––
  if ~iscell(spikeCell)
    [nCells, nS] = size(spikeCell);
    tmp = cell(nCells,1);
    for c = 1:nCells
      st = spikeCell(c,:);          % row of spike times (zeros/NaNs padded)
      st = st(~isnan(st) & st>0);   % drop NaN/zero padding
      tmp{c} = st(:);               % ensure column
    end
    spikeCell = tmp;
  end

  %–– 2) Set up parameters ––
  nTrials   = numel(csTimes);
  nCells    = numel(spikeCell);
  windowLen = diff(win);
  dt        = median(diff(ts));

  %–– 3) Build one mask over ts for all CS windows ––
  maskCS = false(size(ts));
  for t = 1:nTrials
    maskCS = maskCS | (ts >= csTimes(t)+win(1) & ts < csTimes(t)+win(2));
  end
  outs = ts(~maskCS);                % non‐CS timepoints

  %–– 4) FR in each CS window for each cell ––
  FRt = nan(nTrials, nCells);
  for t = 1:nTrials
    t0 = csTimes(t) + win(1);
    t1 = t0 + windowLen;
    for c = 1:nCells
      st = spikeCell{c};
      FRt(t,c) = sum(st>=t0 & st<t1) / windowLen;
    end
  end

  %–– 5) FR in randomly‐sampled windows for each cell ––
  bins    = round(windowLen / dt);
  maxStart= numel(outs) - bins + 1;
  starts  = randi(maxStart, nTrials, 1);

  FRr = nan(nTrials, nCells);
  for t = 1:nTrials
    t0 = outs(starts(t));
    t1 = t0 + windowLen;
    for c = 1:nCells
      st = spikeCell{c};
      FRr(t,c) = sum(st>=t0 & st<t1) / windowLen;
    end
  end

  %–– 6) Filter out cells with too few spikes in CS ––
  totalSpikes = sum(FRt,1) * windowLen;  % back to counts per cell
  keepCells   = totalSpikes >= minSpikes;
  FRt = FRt(:, keepCells);
  FRr = FRr(:, keepCells);
end
