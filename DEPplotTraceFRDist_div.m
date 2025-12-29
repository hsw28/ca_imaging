function DEPplotTraceFRDist_div(useVelFilter, velThresh, applyTo, groupSize, nBoot, demeanFlag)
%%DEP BC CHOOSES RANDOM TIME POINTS INSTEAD OF SPLIT HALF

% plotTraceFRDist_div  Compare trace vs non-trial FR across cells + ROC/AUC
% and add spatial (random-window) TT similarity readouts.
%
% Temporal metrics:
%   - TT: trace vs trace, disjoint groupSize-vs-groupSize grouped PVs
%   - TR: trace vs nontrial, independent groupSize-vs-groupSize grouped PVs
%   - RR: nontrial vs nontrial, disjoint groupSize-vs-groupSize grouped PVs
%
% Spatial metrics (TT-only, using random 2-s windows as pseudo-trials):
%   - TT_space_nonTask: random 2-s windows drawn from NON-TASK periods only
%   - TT_space_all    : random 2-s windows drawn from ALL frames (task+non-task)
%
% New options:
%   groupSize  : # trials/windows per grouped PV (default 25)
%   nBoot      : # bootstrap samples (default 500)
%   demeanFlag : if true, subtract per-cell mean (temporal) or per-feature
%                mean (spatial) before similarity.

if nargin < 1 || isempty(useVelFilter), useVelFilter = false; end
if nargin < 2 || isempty(velThresh),    velThresh    = 4;     end
if nargin < 3 || isempty(applyTo),      applyTo      = 'nontrial'; end
if nargin < 4 || isempty(groupSize),    groupSize    = 25;    end
if nargin < 5 || isempty(nBoot),        nBoot        = 500;   end
if nargin < 6 || isempty(demeanFlag),   demeanFlag   = false; end

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

cols = [0      0.4470 0.7410;   % Trace
        0.8500 0.3250 0.0980;   % Nontrial (no filter)
        0.4660 0.6740 0.1880];  % Nontrial (vel filter)

% keep handle to last session for spatial analysis
spk_last   = [];
pos_last   = [];
cs_last    = [];

for r = 1:nR
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);

    trace_means     = [];
    nontrial0_means = [];
    nontrialV_means = [];

    for d = 1:3
        spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',days{d}));
        posMat   = rat.pos.(sprintf('pos_%s',days{d}));
        csTimes  = rat.CS_times.(sprintf('CS_%s',days{d}));
        ratemask = rat.ratemask.(sprintf('ratemask_%s',days{d}));
        keep     = (ratemask == 1);
        spk      = spk(keep,:);

        % keep last session for spatial analysis
        spk_last = spk;
        pos_last = posMat;
        cs_last  = csTimes;

        % A) Path used by similarity/ROC (temporal FRs)
        useTrace    = useVelFilter && (strcmp(applyTo,'trace')   || strcmp(applyTo,'both'));
        useNontrial = useVelFilter && (strcmp(applyTo,'nontrial')|| strcmp(applyTo,'both'));
        [FRt_A, FRr_A] = popVecSim(spk, posMat, csTimes, win, minSpikes, useTrace, useNontrial, velThresh);
        FRt_all = [FRt_all, mean(FRt_A,1,'omitnan')]; %#ok<AGROW>
        FRr_all = [FRr_all, mean(FRr_A,1,'omitnan')]; %#ok<AGROW>
        FRt_last = FRt_A; %#ok<NASGU>
        FRr_last = FRr_A; %#ok<NASGU>

        % B) Figure 2 data: always compute all three series (per-cell FR means)
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

%% ===== Temporal similarity: grouped PVs (trace vs random) =====
FRt = FRt_last;  % [trials x cells] trace
FRr = FRr_last;  % [trials x cells] nontrial

[TT_vals, TR_vals, RR_vals] = pvSim_bootstrap(FRt, FRr, groupSize, nBoot, demeanFlag);

mean_TT = mean(TT_vals,'omitnan');
mean_TR = mean(TR_vals,'omitnan');
mean_RR = mean(RR_vals,'omitnan');
sd_TT   = std(TT_vals,'omitnan');
sd_TR   = std(TR_vals,'omitnan');
sd_RR   = std(RR_vals,'omitnan');
means_t = [mean_TT, mean_TR, mean_RR];
sds_t   = [sd_TT,   sd_TR,   sd_RR];

fprintf('TEMPORAL grouped PV similarities (groupSize=%d, nBoot=%d, demean=%d):\n', ...
        groupSize, nBoot, demeanFlag);
fprintf('  Mean TT = %.3f (SD=%.3f), TR = %.3f (SD=%.3f), RR = %.3f (SD=%.3f)\n', ...
    mean_TT, sd_TT, mean_TR, sd_TR, mean_RR, sd_RR);

%% ===== Spatial similarity: TT-only (non-task vs all-time, random 2-s windows) =====
TT_space_nonTask = [];
TT_space_all     = [];
mean_TT_s_non    = NaN; sd_TT_s_non = NaN;
mean_TT_s_all    = NaN; sd_TT_s_all = NaN;

if ~isempty(spk_last)
    % ----- spatial random-window params -----
    nWindows      = 1000;         % total random windows per condition
    winDur = 10
  %  winDur        = diff(win);    % match trace window length (e.g., 2 s)
    nXBins        = 10;
    nYBins        = 10;
    minOcc        = 0.5;          % minimum occupancy per bin (s) within a window

    % Non-task only random windows
    FR_win_non = spatialPV_randomWindows(spk_last, pos_last, cs_last, ...
                                         winDur, nWindows, velThresh, true, ...
                                         minOcc, nXBins, nYBins);

    % All frames (task + non-task) random windows
    FR_win_all = spatialPV_randomWindows(spk_last, pos_last, cs_last, ...
                                         winDur, nWindows, velThresh, false, ...
                                         minOcc, nXBins, nYBins);

    if ~isempty(FR_win_non)
        TT_space_nonTask = spatialTT_randomBootstrap(FR_win_non, groupSize, ...
                                                     nBoot, demeanFlag);
        mean_TT_s_non = mean(TT_space_nonTask,'omitnan');
        sd_TT_s_non   = std(TT_space_nonTask,'omitnan');

        fprintf('SPATIAL TT (non-task windows) grouped PV similarities (groupSize=%d, nBoot=%d, demean=%d):\n', ...
                groupSize, nBoot, demeanFlag);
        fprintf('  Mean TT_space_nonTask = %.3f (SD=%.3f)\n', mean_TT_s_non, sd_TT_s_non);
    end

    if ~isempty(FR_win_all)
        TT_space_all = spatialTT_randomBootstrap(FR_win_all, groupSize, ...
                                                 nBoot, demeanFlag);
        mean_TT_s_all = mean(TT_space_all,'omitnan');
        sd_TT_s_all   = std(TT_space_all,'omitnan');

        fprintf('SPATIAL TT (all-time windows) grouped PV similarities (groupSize=%d, nBoot=%d, demean=%d):\n', ...
                groupSize, nBoot, demeanFlag);
        fprintf('  Mean TT_space_all = %.3f (SD=%.3f)\n', mean_TT_s_all, sd_TT_s_all);
    end
end

%% ===== Combined figure: temporal + spatial TT overlay =====
figure('Color','w','Position',[200 200 900 350]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% --- Temporal hist + spatial TT overlay ---
nexttile; hold on;
edges = linspace(-1,1,50);
histogram(TT_vals, edges, 'Normalization','probability', ...
          'FaceColor',cols(1,:), 'FaceAlpha',.6);
histogram(TR_vals, edges, 'Normalization','probability', ...
          'FaceColor',cols(2,:), 'FaceAlpha',.6);
histogram(RR_vals, edges, 'Normalization','probability', ...
          'FaceColor',cols(3,:), 'FaceAlpha',.6);

legStr = {'TT','TR','RR'};

if ~isempty(TT_space_nonTask)
    histogram(TT_space_nonTask, edges, 'Normalization','probability', ...
              'FaceColor',[0.3 0.3 0.3], 'FaceAlpha',.5);
    legStr{end+1} = 'Spatial TT (non-task)';
end
if ~isempty(TT_space_all)
    histogram(TT_space_all, edges, 'Normalization','probability', ...
              'FaceColor',[0.7 0.7 0.7], 'FaceAlpha',.5);
    legStr{end+1} = 'Spatial TT (all time)';
end

legend(legStr,'Location','best');
xlabel('Pearson similarity'); ylabel('Probability');
title(sprintf('Grouped PV similarities (%d-vs-%d)', groupSize, groupSize));

% --- Temporal bar + extra spatial TT bars ---
nexttile; hold on;
means_all = means_t;
sds_all   = sds_t;
xtlbls    = {'TT','TR','RR'};

if ~isnan(mean_TT_s_non)
    means_all(end+1) = mean_TT_s_non;
    sds_all(end+1)   = sd_TT_s_non;
    xtlbls{end+1}    = 'Spatial TT (non-task)';
end
if ~isnan(mean_TT_s_all)
    means_all(end+1) = mean_TT_s_all;
    sds_all(end+1)   = sd_TT_s_all;
    xtlbls{end+1}    = 'Spatial TT (all time)';
end

nBars = numel(means_all);

b = bar(1:nBars, means_all, 'FaceColor','flat');
% colors for temporal bars
b.CData(1:min(3,nBars),:) = cols(1:min(3,nBars),:);
% spatial bars (if present)
if nBars >= 4
    b.CData(4,:) = [0.3 0.3 0.3];  % dark gray
end
if nBars >= 5
    b.CData(5,:) = [0.7 0.7 0.7];  % light gray
end

errorbar(1:nBars, means_all, sds_all, 'k.', 'LineWidth',1.5);
xticks(1:nBars); xticklabels(xtlbls);
ylabel('Mean similarity');
title('Mean ± SD grouped PV similarity'); xlim([0.5 nBars+0.5]);

% temporal significance lines only (TT vs TR, TT vs RR)
p_tt_tr = ranksum(TT_vals, TR_vals);
p_tt_rr = ranksum(TT_vals, RR_vals);
ymax = max(means_t + sds_t); dy = 0.05*ymax; y1 = ymax + dy; y2 = ymax + 2*dy;
plot([1 2],[y1 y1],'-k','LineWidth',1.2);
if p_tt_tr<.05, text(1.5,y1+0.02,'*','horiz','center','FontSize',12); end
plot([1 3],[y2 y2],'-k','LineWidth',1.2);
if p_tt_rr<.05, text(2.0,y2+0.02,'*','horiz','center','FontSize',12); end

%% ===== Figure 2: per-rat + pooled (unchanged) =====
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

    % -------- Column 1: BAR with SD + significance lines ----------
    nexttile; hold on;
    m  = [mean(tdat,'omitnan'), mean(r0dat,'omitnan'), mean(rvdat,'omitnan')];
    sd = [std(tdat,'omitnan'), std(r0dat,'omitnan'), std(rvdat,'omitnan')];

    bb = bar(1:3, m, 'FaceColor','flat'); bb.CData = cols;
    errorbar(1:3, m, sd, 'k.', 'LineWidth',1.2);
    xlim([0.5 3.5]); set(gca,'XTick',1:3,'XTickLabel',{'Trace','Nontrial','Nontrial (vel≥thr)'});
    ylabel('Mean FR (Hz)'); title(sprintf('%s — mean FR ± SD', rowTitle)); box on;

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
    edges2 = linspace(0, max(finiteAll), 50);
    histogram(tdat, 'BinEdges',edges2,'Normalization','probability', ...
              'FaceColor',cols(1,:), 'FaceAlpha',.6);
    histogram(ndat, 'BinEdges',edges2,'Normalization','probability', ...
              'FaceColor',c2,         'FaceAlpha',.6);
    ax = gca; ax.YAxis.Exponent = 0; ytickformat('%.3f');
    xlabel('Firing rate (Hz)'); ylabel('Probability');
    legend('Trace', lab, 'Location','best'); box on;
    title(sprintf('%s — FR distributions', rowTitle));
end

%% ===== ROC & AUC (unchanged across-all) =====
labels = [ ones(1, nnz(isfinite(FRt_all))) , ...
           zeros(1, nnz(isfinite(FRr_all))) ];
scores = [ FRt_all(isfinite(FRt_all)) , ...
           FRr_all(isfinite(FRr_all))      ];
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
histogram(permAUCs, 30, 'FaceColor',[.7 .7 .7], 'EdgeColor','k', ...
          'Normalization', 'probability');
hold on; xline(AUC, 'r--', 'LineWidth', 2);
xlabel('AUC'); ylabel('Frequency');
title(sprintf('Permutation Test for AUC (Actual = %.3f)', AUC));
permP = mean(permAUCs >= AUC);
fprintf('Permutation p-value for AUC = %.4f\n', permP);
end

%% ===== Helper: temporal FR computation =====

function [FRt,FRr] = popVecSim(spikeCell, pos, csTimes, win, minSpikes, useTraceFilter, useNontrialFilter, velThresh)
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

      % Random nontrial window
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

totalCounts = nansum(FRt,1) * windowLen;
keepCells   = totalCounts >= minSpikes;
FRt = FRt(:, keepCells);
FRr = FRr(:, keepCells);
end

function [cnt, runT] = countSpikesAndRunTime(st, vt, vmag, t0, t1, velThresh)
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

%% ===== Spatial helpers: random-window PVs with spatial-bin overlap =====

function FR_windows = spatialPV_randomWindows(spikeCell,pos,csTimes,winDur, ...
                                              nWindows,velThresh,excludeTask, ...
                                              minOcc,nXBins,nYBins)
% spatialPV_randomWindows  Build random-window spatial population vectors.
%
% Each window:
%   - picks a random contiguous block of frames (~winDur)
%   - accumulates occupancy in 2D spatial bins
%   - builds per-neuron rate maps (Hz) per bin
%
% Output FR_windows: [nWindows x nBins x nCells]

if nargin < 5 || isempty(nWindows),    nWindows    = 50; end
if nargin < 6,                         velThresh   = []; end
if nargin < 7 || isempty(excludeTask), excludeTask = true; end
if nargin < 8 || isempty(minOcc),      minOcc      = 1; end   % seconds
if nargin < 9 || isempty(nXBins),      nXBins      = 7;   end
if nargin < 10|| isempty(nYBins),      nYBins      = 7;   end

% make sure spikes are in cell array
if ~iscell(spikeCell)
    [nCells, ~] = size(spikeCell);
    tmp = cell(nCells,1);
    for c = 1:nCells
        st = spikeCell(c,:); st = st(~isnan(st) & st>0);
        tmp{c} = st(:);
    end
    spikeCell = tmp;
end
nCells = numel(spikeCell);

t  = pos(:,1);
x  = pos(:,2);
y  = pos(:,3);
dt = median(diff(t),'omitnan');

includeMask = true(size(t));

% optionally exclude trace windows around CS
if excludeTask && ~isempty(csTimes)
    maskCS = false(size(t));
    for k = 1:numel(csTimes)
        maskCS = maskCS | (t >= csTimes(k) & t < csTimes(k) + winDur);
    end
    includeMask = includeMask & ~maskCS;
end

% optional running-only mask
if ~isempty(velThresh)
    vel  = ca_velocity(pos);          % [speed; time]
    vt   = vel(2,:)';
    vmag = vel(1,:)';
    spd  = interp1(vt, vmag, t, 'linear','extrap');
    includeMask = includeMask & (spd >= velThresh);
end

t_inc = t(includeMask);
x_inc = x(includeMask);
y_inc = y(includeMask);

if numel(t_inc) < 2
    FR_windows = [];
    warning('Not enough included samples for spatial random windows.');
    return;
end

% spatial bins from included positions
xEdges = linspace(min(x_inc), max(x_inc), nXBins+1);
yEdges = linspace(min(y_inc), max(y_inc), nYBins+1);
nBins  = nXBins * nYBins;

% For all INCLUDED frames, precompute bin index
[binX,~] = discretize(x_inc,xEdges);
[binY,~] = discretize(y_inc,yEdges);
valid   = binX>=1 & binX<=nXBins & binY>=1 & binY<=nYBins;
binX    = binX(valid);
binY    = binY(valid);
t_inc   = t_inc(valid);

binIdx_all = sub2ind([nXBins nYBins], binX, binY);  % indices 1..nBins
nFramesInc = numel(t_inc);

W = max(1, round(winDur / dt));     % # frames ~ window duration
maxStart = nFramesInc - W + 1;
if maxStart < 1
    FR_windows = [];
    warning('Not enough continuous included frames for windowDur.');
    return;
end

FR_windows = nan(nWindows, nBins, nCells);

for s = 1:nWindows
    sIdx = randi(maxStart);
    idxWin = sIdx:(sIdx+W-1);
    binIdx_win = binIdx_all(idxWin);
    t0   = t_inc(idxWin(1));
    t1   = t0 + winDur;  % used only for spike windowing

    % occupancy per bin (seconds) within this window
    occCounts = accumarray(binIdx_win, dt, [nBins 1], @sum, 0);
    occCounts(occCounts < minOcc) = NaN;   % mark low-occupancy bins as invalid

    % per-neuron rate maps
    for c = 1:nCells
        st = spikeCell{c};
        if isempty(st)
            FR_windows(s,:,c) = NaN;
            continue;
        end

        stWin = st(st>=t0 & st<t1);
        if isempty(stWin)
            FR_windows(s,:,c) = 0;  % no spikes but occupancy is info
            continue;
        end

        % position at spike times via interpolation
        x_spk = interp1(t, x, stWin, 'linear','extrap');
        y_spk = interp1(t, y, stWin, 'linear','extrap');

        [bx_spk,~] = discretize(x_spk, xEdges);
        [by_spk,~] = discretize(y_spk, yEdges);

        validSpk = bx_spk>=1 & bx_spk<=nXBins & by_spk>=1 & by_spk<=nYBins;
        bx_spk   = bx_spk(validSpk);
        by_spk   = by_spk(validSpk);

        if isempty(bx_spk)
            FR_windows(s,:,c) = 0;
            continue;
        end

        bIdx_spk = sub2ind([nXBins nYBins], bx_spk, by_spk);
        spkCounts = accumarray(bIdx_spk, 1, [nBins 1], @sum, 0);

        rate_c = spkCounts ./ occCounts;      % Hz; NaN where occCounts is NaN
        FR_windows(s,:,c) = rate_c;
    end
end

% drop cells that are NaN everywhere across all windows/bins
tmp = reshape(FR_windows, [], nCells);  % [nWindows*nBins x nCells]
keepCells = any(isfinite(tmp),1);
FR_windows = FR_windows(:,:,keepCells);
end

function TT_vals = spatialTT_randomBootstrap(FR_windows, groupSize, nBoot, demeanFlag)
% spatialTT_randomBootstrap
% FR_windows: [nWindows x nBins x nCells]
% groupSize : windows per group
% demeanFlag: if true, demean flattened vectors before corr
%
% Overlap rule:
%   Let visited bins = any finite value across cells.
%   Require:
%     nShared >= max(minBinsAbs, fracSmall * min(nVisited1, nVisited2))

FR_windows = double(FR_windows);
[nWindows, nBins, nCells] = size(FR_windows);

if 2*groupSize > nWindows
    error('groupSize=%d too large for nWindows=%d (need 2*groupSize <= nWindows).', ...
          groupSize, nWindows);
end

% parameters for spatial overlap
fracSmall  = 0.75;  % require at least 50% overlap relative to smaller map
minBinsAbs = 0;    % and at least 5 bins shared in absolute terms

TT_vals = nan(nBoot,1);

for b = 1:nBoot
    idx = randperm(nWindows);
    g1 = idx(1:groupSize);
    g2 = idx(groupSize+1:2*groupSize);

    % mean rate maps per group: [nBins x nCells]
    pv1 = squeeze(mean(FR_windows(g1,:,:), 1, 'omitnan'));  % [nBins x nCells]
    pv2 = squeeze(mean(FR_windows(g2,:,:), 1, 'omitnan'));  % [nBins x nCells]

    if isempty(pv1) || isempty(pv2)
        TT_vals(b) = NaN;
        continue;
    end

    % visited bins = any finite cell
    visited1 = any(isfinite(pv1), 2);   % [nBins x 1]
    visited2 = any(isfinite(pv2), 2);

    nVisited1 = nnz(visited1);
    nVisited2 = nnz(visited2);

    if nVisited1 == 0 || nVisited2 == 0
        TT_vals(b) = NaN;
        continue;
    end

    sharedBins = visited1 & visited2;
    nShared    = nnz(sharedBins);

    minSharedThisPair = max(minBinsAbs, round(fracSmall * min(nVisited1, nVisited2)));

    if nShared < minSharedThisPair
        TT_vals(b) = NaN;
        continue;
    end

    pv1_use = pv1(sharedBins,:);  % [nShared x nCells]
    pv2_use = pv2(sharedBins,:);

    v1 = pv1_use(:);
    v2 = pv2_use(:);

    if demeanFlag
        X  = [v1; v2];
        mu = mean(X,'omitnan');
        v1 = v1 - mu;
        v2 = v2 - mu;
    end

    TT_vals(b) = corr(v1, v2, 'Rows','pairwise');
end
end

%% ===== Generic grouped-similarity helper (temporal) =====

function [TT_vals, TR_vals, RR_vals] = pvSim_bootstrap(FRt, FRr, groupSize, nBoot, demeanFlag)
FRt = double(FRt);
FRr = double(FRr);

keepCells = any(isfinite(FRt) | isfinite(FRr),1);
FRt = FRt(:,keepCells);
FRr = FRr(:,keepCells);

[nTrials_t, nCells]   = size(FRt);
[nTrials_r, nCells_r] = size(FRr);
if nCells_r ~= nCells
    error('FRt and FRr must have the same number of cells.');
end
if nTrials_t ~= nTrials_r
    warning('FRt and FRr have different # of trials; using min for grouping.');
end
nTrials = min(nTrials_t, nTrials_r);

FRt = FRt(1:nTrials,:);
FRr = FRr(1:nTrials,:);

if 2*groupSize > nTrials
    error('groupSize=%d too large for nTrials=%d (need 2*groupSize <= nTrials).', ...
          groupSize, nTrials);
end

if demeanFlag
    X = [FRt; FRr];
    mu = mean(X,1,'omitnan');
    FRt = FRt - mu;
    FRr = FRr - mu;
end

TT_vals = nan(nBoot,1);
TR_vals = nan(nBoot,1);
RR_vals = nan(nBoot,1);

for b = 1:nBoot
    % TT
    idx = randperm(nTrials);
    g1 = idx(1:groupSize);
    g2 = idx(groupSize+1:2*groupSize);
    pv_t1 = mean(FRt(g1,:),1,'omitnan')';
    pv_t2 = mean(FRt(g2,:),1,'omitnan')';
    TT_vals(b) = corr(pv_t1, pv_t2, 'Rows','pairwise');

    % RR
    idx = randperm(nTrials);
    g1 = idx(1:groupSize);
    g2 = idx(groupSize+1:2*groupSize);
    pv_r1 = mean(FRr(g1,:),1,'omitnan')';
    pv_r2 = mean(FRr(g2,:),1,'omitnan')';
    RR_vals(b) = corr(pv_r1, pv_r2, 'Rows','pairwise');

    % TR
    gT = randperm(nTrials, groupSize);
    gR = randperm(nTrials, groupSize);
    pv_t = mean(FRt(gT,:),1,'omitnan')';
    pv_r = mean(FRr(gR,:),1,'omitnan')';
    TR_vals(b) = corr(pv_t, pv_r, 'Rows','pairwise');
end
end

%% ===== Stats/plot helpers =====

function p = paired_or_not(x, y)
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
