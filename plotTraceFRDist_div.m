function plotTraceFRDist_div(useVelFilter, velThresh, applyTo, groupSize, nBoot, demeanFlag)
% plotTraceFRDist_div  Compare trace vs non-trial FR across cells + ROC/AUC
% + within-session temporal TT/TR/RR pooled across 15 sessions (no spatial metrics).
%
% Temporal metrics (within-session):
%   - TT: trace vs trace, disjoint groupSize-vs-groupSize grouped PVs
%   - TR: trace vs nontrial, grouped PVs from each (pairing randomized)
%   - RR: nontrial vs nontrial, disjoint groupSize-vs-groupSize grouped PVs
%
% Pooling:
%   - compute TT/TR/RR per session (within-session only)
%   - bar plot shows mean ± SD across the 15 session means (equal weight per session)
%
% Diagnostics:
%   - per-session labeled "Non-task running diagnostics" computed from the exact
%     random non-task windows used for FRr (printed once per session)

if nargin < 1 || isempty(useVelFilter), useVelFilter = true; end
if nargin < 2 || isempty(velThresh),    velThresh    = 4;     end
if nargin < 3 || isempty(applyTo),      applyTo      = 'nontrial'; end
if nargin < 4 || isempty(groupSize),    groupSize    = 25;    end
if nargin < 5 || isempty(nBoot),        nBoot        = 500;   end
if nargin < 6 || isempty(demeanFlag),   demeanFlag   = false; end

applyTo = validatestring(applyTo, {'trace','nontrial','both'});

% ---- sanitize numeric inputs ----
if ischar(nBoot) || isstring(nBoot), nBoot = str2double(nBoot); end
if isempty(nBoot) || ~isfinite(nBoot) || ~isscalar(nBoot), nBoot = 500; end
nBoot = max(1, round(nBoot));

if ischar(groupSize) || isstring(groupSize), groupSize = str2double(groupSize); end
if isempty(groupSize) || ~isfinite(groupSize) || ~isscalar(groupSize), groupSize = 25; end
groupSize = max(1, round(groupSize));

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

% ---- SESSION-LEVEL pooling (equal weight per session) ----
nSessions = nR * 3;
sess_TT_mean = nan(nSessions,1);
sess_TR_mean = nan(nSessions,1);
sess_RR_mean = nan(nSessions,1);

showHist = true;   % set false if you only want the summary bars
TT_draws = [];
TR_draws = [];
RR_draws = [];

sessIdx = 0;

for r = 1:nR
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);

    trace_means     = [];
    nontrial0_means = [];
    nontrialV_means = [];

    for d = 1:3
        sessIdx = sessIdx + 1;

        spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',days{d}));
        posMat   = rat.pos.(sprintf('pos_%s',days{d}));
        csTimes  = rat.CS_times.(sprintf('CS_%s',days{d}));
        ratemask = rat.ratemask.(sprintf('ratemask_%s',days{d}));
        keep     = (ratemask == 1);
        spk      = spk(keep,:);

        sessionLabel = sprintf('%s %s', ratNames{r}, days{d});

        % A) Path used by similarity/ROC (temporal FRs)
        useTrace    = useVelFilter && (strcmp(applyTo,'trace')   || strcmp(applyTo,'both'));
        useNontrial = useVelFilter && (strcmp(applyTo,'nontrial')|| strcmp(applyTo,'both'));

        [FRt_A, FRr_A] = popVecSim(spk, posMat, csTimes, win, minSpikes, ...
                                   useTrace, useNontrial, velThresh, sessionLabel);

        FRt_all = [FRt_all, mean(FRt_A,1,'omitnan')]; %#ok<AGROW>
        FRr_all = [FRr_all, mean(FRr_A,1,'omitnan')]; %#ok<AGROW>

        % ---- temporal TT/TR/RR: WITHIN THIS SESSION ONLY ----
        [TT_s, TR_s, RR_s] = pvSim_bootstrap(FRt_A, FRr_A, groupSize, nBoot, demeanFlag);
        sess_TT_mean(sessIdx) = mean(TT_s,'omitnan');
        sess_TR_mean(sessIdx) = mean(TR_s,'omitnan');
        sess_RR_mean(sessIdx) = mean(RR_s,'omitnan');

        if showHist
            TT_draws = [TT_draws; TT_s(:)]; %#ok<AGROW>
            TR_draws = [TR_draws; TR_s(:)]; %#ok<AGROW>
            RR_draws = [RR_draws; RR_s(:)]; %#ok<AGROW>
        end

        % B) Figure 2 data: always compute all three series (per-cell FR means)
        [FRt_traceOnly, ~] = popVecSim(spk, posMat, csTimes, win, minSpikes, false, false, velThresh, '');
        [~, FRr_noVF]       = popVecSim(spk, posMat, csTimes, win, minSpikes, false, false, velThresh, '');
        [~, FRr_velF]       = popVecSim(spk, posMat, csTimes, win, minSpikes, false, true,  velThresh, '');

        trace_means     = [trace_means,     mean(FRt_traceOnly,1,'omitnan')]; %#ok<AGROW>
        nontrial0_means = [nontrial0_means, mean(FRr_noVF,     1,'omitnan')]; %#ok<AGROW>
        nontrialV_means = [nontrialV_means, mean(FRr_velF,     1,'omitnan')]; %#ok<AGROW>
    end

    FRt_meancell{r}  = trace_means(:);
    FRr0_meancell{r} = nontrial0_means(:);
    FRrV_meancell{r} = nontrialV_means(:);
end

%% ===== POOLED (session-averaged) summaries for plotting =====
mean_TT = mean(sess_TT_mean,'omitnan');
mean_TR = mean(sess_TR_mean,'omitnan');
mean_RR = mean(sess_RR_mean,'omitnan');

sd_TT = std(sess_TT_mean,'omitnan');
sd_TR = std(sess_TR_mean,'omitnan');
sd_RR = std(sess_RR_mean,'omitnan');

fprintf('ttest for TT vs TR')
[a b c d] = ttest(sess_TT_mean,sess_TR_mean)
fprintf('ttest for TT vs RR')
[a b c d] = ttest(sess_TT_mean,sess_RR_mean)

means_t = [mean_TT, mean_TR, mean_RR];
sds_t   = [sd_TT,   sd_TR,   sd_RR];

fprintf('TEMPORAL grouped PV similarities (SESSION-AVERAGED; nSess=%d; groupSize=%d; nBoot/session=%d; demean=%d):\n', ...
        nSessions, groupSize, nBoot, demeanFlag);
fprintf('  TT: mean=%.3f SD_acrossSess=%.3f | TR: mean=%.3f SD=%.3f | RR: mean=%.3f SD=%.3f\n', ...
        mean_TT, sd_TT, mean_TR, sd_TR, mean_RR, sd_RR);

%% ===== Figure: temporal only (hist + session-mean bars) =====
figure('Color','w','Position',[200 200 900 350]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Histogram of pooled draws (visualization only)
nexttile; hold on;
edges = linspace(-1,1,50);
if showHist
    histogram(TT_draws, edges, 'Normalization','probability', 'FaceColor',cols(1,:), 'FaceAlpha',.6);
    histogram(TR_draws, edges, 'Normalization','probability', 'FaceColor',cols(2,:), 'FaceAlpha',.6);
    histogram(RR_draws, edges, 'Normalization','probability', 'FaceColor',cols(3,:), 'FaceAlpha',.6);
    legend({'TT draws','TR draws','RR draws'},'Location','best');
else
    text(0.1,0.5,'showHist=false','Units','normalized');
end
xlabel('Pearson similarity'); ylabel('Probability');
title(sprintf('Within-session bootstrap draws (nSess=%d)', nSessions));

% Bar plot: session-mean TT/TR/RR mean ± SD(across sessions)
nexttile; hold on;
b = bar(1:3, means_t, 'FaceColor','flat');
b.CData = cols(1:3,:);
errorbar(1:3, means_t, sds_t, 'k.', 'LineWidth',1.5);
xticks(1:3); xticklabels({'TT','TR','RR'});
ylabel('Mean similarity');
title(sprintf('Session-averaged mean ± SD (nSess=%d)', nSessions));
xlim([0.5 3.5]);

% significance on SESSION MEANS (paired across sessions)
p_tt_tr = signrank(sess_TT_mean, sess_TR_mean);
p_tt_rr = signrank(sess_TT_mean, sess_RR_mean);
ymax = max(means_t + sds_t); if ~isfinite(ymax) || ymax==0, ymax = 1; end
dy = 0.05*ymax; y1 = ymax + dy; y2 = ymax + 2*dy;
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

    nexttile; hold on;
    m  = [mean(tdat,'omitnan'), mean(r0dat,'omitnan'), mean(rvdat,'omitnan')];
    sd = [std(tdat,'omitnan'),  std(r0dat,'omitnan'),  std(rvdat,'omitnan')];

    bb = bar(1:3, m, 'FaceColor','flat'); bb.CData = cols;
    errorbar(1:3, m, sd, 'k.', 'LineWidth',1.2);
    xlim([0.5 3.5]);
    set(gca,'XTick',1:3,'XTickLabel',{'Trace-task','Non-trace-task','Non-trace-task (vel≥thr)'});
    ylabel('Mean FR (Hz)'); title(sprintf('%s — mean FR ± SD', rowTitle)); box on;

    p12 = paired_or_not(tdat, r0dat);
    p13 = paired_or_not(tdat, rvdat);
    yl = ylim; ybase = yl(2);
    draw_sigline(gca, 1, 2, ybase*0.94, p12);
    draw_sigline(gca, 1, 3, ybase*0.98, p13);

    nexttile; hold on;
    if useVelFilter
        ndat = rvdat; lab = 'Non-trace-task (vel≥thr)'; c2 = cols(3,:);
    else
        ndat = r0dat; lab = 'Non-trace-task';          c2 = cols(2,:);
    end
    finiteAll = [tdat(isfinite(tdat)); ndat(isfinite(ndat))];
    if isempty(finiteAll), finiteAll = 0; end
    edges2 = linspace(0, max(finiteAll), 50);
    histogram(tdat, 'BinEdges',edges2,'Normalization','probability', 'FaceColor',cols(1,:), 'FaceAlpha',.6);
    histogram(ndat, 'BinEdges',edges2,'Normalization','probability', 'FaceColor',c2,         'FaceAlpha',.6);
    ax = gca; ax.YAxis.Exponent = 0; ytickformat('%.3f');
    xlabel('Firing rate (Hz)'); ylabel('Probability');
    legend('Trace', lab, 'Location','best'); box on;
    title(sprintf('%s — FR distributions', rowTitle));
end

%% ===== ROC & AUC (pooled across cells) =====
labels = [ ones(1, nnz(isfinite(FRt_all))) , zeros(1, nnz(isfinite(FRr_all))) ];
scores = [ FRt_all(isfinite(FRt_all)) , FRr_all(isfinite(FRr_all)) ];
[Xroc,Yroc,~,AUC] = perfcurve(labels(:), scores(:), 1);
fprintf('\nTrace vs non-trial FR AUC = %.3f\n', AUC);

figure('Color','w','Position',[350 350 500 400]);
plot(Xroc, Yroc, 'LineWidth',2); hold on; plot([0 1],[0 1],'k--');
xlabel('False positive rate'); ylabel('True positive rate');
title(sprintf('ROC curve (AUC = %.3f)', AUC)); axis square;

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

%% ===== Helper: temporal FR computation (with labeled diagnostics) =====
function [FRt,FRr] = popVecSim(spikeCell, pos, csTimes, win, minSpikes, useTraceFilter, useNontrialFilter, velThresh, sessionLabel)

if nargin < 9 || isempty(sessionLabel), sessionLabel = ''; end

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

vel  = ca_velocity(pos);    % [speed; time]
vt   = vel(2,:)';
vmag = vel(1,:)';

maskCS = false(size(ts));
for t = 1:nTrials
    maskCS = maskCS | (ts >= csTimes(t)+win(1) & ts < csTimes(t)+win(2));
end
outs = ts(~maskCS);

dt   = median(diff(ts),'omitnan');
bins = max(1, round(windowLen / dt));
maxStart = numel(outs) - bins + 1;
if maxStart < 1
    error('Not enough non-CS timepoints to place random windows.');
end

FRt = nan(nTrials, nCells);
FRr = nan(nTrials, nCells);

r0_used = nan(nTrials,1);
r1_used = nan(nTrials,1);

for t = 1:nTrials
    t0 = csTimes(t) + win(1);

    % choose ONE non-task window per trial, shared across cells
    sIdx = randi(maxStart);
    r0   = outs(sIdx);
    r1   = r0 + windowLen;
    r0_used(t) = r0;
    r1_used(t) = r1;

    for c = 1:nCells
        st = spikeCell{c};

        if ~useTraceFilter
            FRt(t,c) = sum(st>=t0 & st<(t0+windowLen)) / windowLen;
        else
            [cntT, runT] = countSpikesAndRunTime(st, vt, vmag, t0, t0+windowLen, velThresh);
            if runT>0, FRt(t,c) = cntT / runT; else, FRt(t,c) = NaN; end
        end

        if ~useNontrialFilter
            FRr(t,c) = sum(st>=r0 & st<r1) / windowLen;
        else
            [cntR, runR] = countSpikesAndRunTime(st, vt, vmag, r0, r1, velThresh);
            if runR>0, FRr(t,c) = cntR / runR; else, FRr(t,c) = NaN; end
        end
    end
end

if useNontrialFilter && ~isempty(sessionLabel)
    runDur = nan(nTrials,1);
    for tt = 1:nTrials
        r0_local = r0_used(tt);
        r1_local = r1_used(tt);
        if ~isfinite(r0_local) || ~isfinite(r1_local)
            runDur(tt) = NaN; continue;
        end

        use_local = (vt >= r0_local & vt < r1_local);
        if ~any(use_local)
            runDur(tt) = 0;
        else
            vt_use_local = vt(use_local);
            vm_use_local = (vmag(use_local) >= velThresh);
            if numel(vt_use_local) < 2
                runDur(tt) = 0;
            else
                dt_local = [diff(vt_use_local); median(diff(vt_use_local),'omitnan')];
                dt_local(~isfinite(dt_local) | dt_local<=0) = 0;
                runDur(tt) = sum(dt_local(vm_use_local));
            end
        end
    end

    fprintf('\n=== %s ===\n', sessionLabel);
    fprintf('=== Non-task running diagnostics ===\n');
    fprintf('Average usable run time per window: %.3f s\n', mean(runDur,'omitnan'));
    fprintf('Median usable run time per window:  %.3f s\n', median(runDur,'omitnan'));
    fprintf('Min / max usable run time:          %.3f / %.3f s\n', ...
            min(runDur,[],'omitnan'), max(runDur,[],'omitnan'));
    fprintf('%% windows with >0 run time:          %.1f%%\n', ...
            100 * nnz(runDur > 0) / nnz(isfinite(runDur)));
    fprintf('%% windows with >=0.5 s run time:     %.1f%%\n\n', ...
            100 * nnz(runDur >= 0.5) / nnz(isfinite(runDur)));
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
if ~any(use), runT = 0; return; end
vt_use = vt(use);
vm_use = vmag(use) >= velThresh;
if numel(vt_use) < 2, runT = 0; return; end
dt = [diff(vt_use); median(diff(vt_use),'omitnan')];
dt(~isfinite(dt) | dt<=0) = 0;
runT = sum(dt(vm_use));
end

%% ===== temporal grouped-similarity helper =====
function [TT_vals, TR_vals, RR_vals] = pvSim_bootstrap(FRt, FRr, groupSize, nBoot, demeanFlag)
FRt = double(FRt);
FRr = double(FRr);

if ischar(nBoot) || isstring(nBoot), nBoot = str2double(nBoot); end
if isempty(nBoot) || ~isfinite(nBoot) || ~isscalar(nBoot), nBoot = 500; end
nBoot = max(1, round(nBoot));

if ischar(groupSize) || isstring(groupSize), groupSize = str2double(groupSize); end
if isempty(groupSize) || ~isfinite(groupSize) || ~isscalar(groupSize), groupSize = 25; end
groupSize = max(1, round(groupSize));

keepCells = any(isfinite(FRt) | isfinite(FRr),1);
FRt = FRt(:,keepCells);
FRr = FRr(:,keepCells);

nTrials = min(size(FRt,1), size(FRr,1));
FRt = FRt(1:nTrials,:);
FRr = FRr(1:nTrials,:);

maxGroup = floor(nTrials/2);
if groupSize > maxGroup
    oldGS = groupSize;
    groupSize = maxGroup;
    fprintf('NOTE: groupSize clamped from %d -> %d (nTrials=%d)\n', oldGS, groupSize, nTrials);
end
if groupSize < 1
    TT_vals = []; TR_vals = []; RR_vals = [];
    return;
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
    idx = randperm(nTrials);
    g1 = idx(1:groupSize);
    g2 = idx(groupSize+1:2*groupSize);
    pv_t1 = mean(FRt(g1,:),1,'omitnan')';
    pv_t2 = mean(FRt(g2,:),1,'omitnan')';
    TT_vals(b) = corr(pv_t1, pv_t2, 'Rows','pairwise');

    idx = randperm(nTrials);
    g1 = idx(1:groupSize);
    g2 = idx(groupSize+1:2*groupSize);
    pv_r1 = mean(FRr(g1,:),1,'omitnan')';
    pv_r2 = mean(FRr(g2,:),1,'omitnan')';
    RR_vals(b) = corr(pv_r1, pv_r2, 'Rows','pairwise');

    xsel = randi(4);
    if xsel == 1, pv_t = pv_t1; pv_r = pv_r1;
    elseif xsel == 2, pv_t = pv_t2; pv_r = pv_r2;
    elseif xsel == 3, pv_t = pv_t2; pv_r = pv_r1;
    else,             pv_t = pv_t1; pv_r = pv_r2;
    end
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
