function plotTraceFRDist_div(useVelFilter, velThresh, applyTo, groupSize, nBoot, normMode, doSpeedMatch, nonTaskBufferSec, speedMatchMode, speedProfileNBins)
% plotTraceFRDist_div  Compare trace vs non-trial FR across cells + ROC/AUC
% + within-session temporal TT/TR/RR/TS/SS pooled across 15 sessions (no spatial metrics).
%
% Temporal metrics (within-session):
%   - TT: task vs task (trace-window FR PVs, disjoint groupSize-vs-groupSize)
%   - RR: non-task vs non-task (random non-task windows, disjoint)
%   - TR: task vs non-task (random pairing of grouped PVs)
%   - TS: task vs speed-matched non-task (random pairing)
%   - SS: speed-matched vs speed-matched (disjoint)
%
% Normalization (normMode):
%   - 'none'     : no normalization
%   - 'demean'   : subtract per-cell mean across [FRt; FRr; FRsm] within session
%   - 'meanrate' : divide each cell by its mean rate over the entire day (session)
%
% Speed matching:
%   - doSpeedMatch=true constructs a non-task window per trial whose speed matches the task window
%   - speedMatchMode:
%       'mean'    : match mean speed (original)
%       'meanstd' : (Option A) match [mean speed, std speed]
%       'profile' : (Option B) match binned speed profile across window (speedProfileNBins bins)
%
% Non-task buffer:
%   - nonTaskBufferSec excludes an additional buffer after each task window
%     when sampling non-task windows (applies to both random and speed-matched)
%
% Strict non-task windows:
%   - candidates are STRICT: the entire 2s window must lie in ~maskCS (task+buffer excluded)

if nargin < 1 || isempty(useVelFilter), useVelFilter = false; end
if nargin < 2 || isempty(velThresh),    velThresh    = 4;     end
if nargin < 3 || isempty(applyTo),      applyTo      = 'both'; end
if nargin < 4 || isempty(groupSize),    groupSize    = 25;    end
if nargin < 5 || isempty(nBoot),        nBoot        = 500;   end
if nargin < 6 || isempty(normMode),     normMode     = 'none'; end
if nargin < 7 || isempty(doSpeedMatch), doSpeedMatch = true; end
if nargin < 8 || isempty(nonTaskBufferSec), nonTaskBufferSec = 0; end
if nargin < 9 || isempty(speedMatchMode), speedMatchMode = 'profile'; end
    % meanstd: mean + std matching (default)
    % profile: binned speed-profile matching (e.g., 4 bins of 0.5s each)
    % mean: mean speed only
if nargin < 10 || isempty(speedProfileNBins), speedProfileNBins = 15; end %only if matching is profile


applyTo  = validatestring(applyTo,  {'trace','nontrial','both'});
normMode = validatestring(normMode, {'none','demean','meanrate'});
speedMatchMode = validatestring(speedMatchMode, {'mean','meanstd','profile'});

% ---- sanitize numeric inputs ----
if ischar(nBoot) || isstring(nBoot), nBoot = str2double(nBoot); end
if isempty(nBoot) || ~isfinite(nBoot) || ~isscalar(nBoot), nBoot = 500; end
nBoot = max(1, round(nBoot));

if ischar(groupSize) || isstring(groupSize), groupSize = str2double(groupSize); end
if isempty(groupSize) || ~isfinite(groupSize) || ~isscalar(groupSize), groupSize = 25; end
groupSize = max(1, round(groupSize));

if ~isscalar(nonTaskBufferSec) || ~isfinite(nonTaskBufferSec) || nonTaskBufferSec < 0
    nonTaskBufferSec = 2;
end

if ~isscalar(speedProfileNBins) || ~isfinite(speedProfileNBins) || speedProfileNBins < 2
    speedProfileNBins = 4;
end
speedProfileNBins = round(speedProfileNBins);

ratNames  = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
win       = [0 2];        % CS trace window (s)
minSpikes = 0;

FRt_all = [];
FRr_all = [];

% Per-rat storage for Figure 2 (FR distributions)
nR = numel(ratNames);
FRt_meancell  = cell(nR,1);   % Trace (no filter)
FRr0_meancell = cell(nR,1);   % Nontrial (no filter)
FRrV_meancell = cell(nR,1);   % Nontrial (vel filter)

cols = [0      0.4470 0.7410;   % Trace
        0.8500 0.3250 0.0980;   % Nontrial (no filter)
        0.4660 0.6740 0.1880;   % Nontrial (vel filter)
        0.4940 0.1840 0.5560;   % Speed-matched
        0.3010 0.7450 0.9330];  % SS

% ---- SESSION-LEVEL pooling (equal weight per session) ----
nSessions = nR * 3;
sess_TT_mean = nan(nSessions,1);
sess_TR_mean = nan(nSessions,1);
sess_RR_mean = nan(nSessions,1);
sess_TS_mean = nan(nSessions,1);
sess_SS_mean = nan(nSessions,1);

% Per-rat session means (3 sessions per rat)
rat_TT = nan(nR,3); rat_TR = nan(nR,3); rat_RR = nan(nR,3); rat_TS = nan(nR,3); rat_SS = nan(nR,3);

showHist = true;   % set false if you only want the summary bars
TT_draws = [];
TR_draws = [];
RR_draws = [];
TS_draws = [];
SS_draws = [];

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

        [FRt_A, FRr_A, FRsm_A, meanRateDay] = popVecSim( ...
            spk, posMat, csTimes, win, minSpikes, ...
            useTrace, useNontrial, velThresh, sessionLabel, ...
            doSpeedMatch, nonTaskBufferSec, speedMatchMode, speedProfileNBins);

        FRt_all = [FRt_all, mean(FRt_A,1,'omitnan')]; %#ok<AGROW>
        FRr_all = [FRr_all, mean(FRr_A,1,'omitnan')]; %#ok<AGROW>

        % ---- temporal TT/TR/RR/TS/SS: WITHIN THIS SESSION ONLY ----
        [TT_s, TR_s, RR_s, TS_s, SS_s] = pvSim_bootstrap_5( ...
            FRt_A, FRr_A, FRsm_A, meanRateDay, groupSize, nBoot, normMode);

        sess_TT_mean(sessIdx) = mean(TT_s,'omitnan');
        sess_TR_mean(sessIdx) = mean(TR_s,'omitnan');
        sess_RR_mean(sessIdx) = mean(RR_s,'omitnan');
        sess_TS_mean(sessIdx) = mean(TS_s,'omitnan');
        sess_SS_mean(sessIdx) = mean(SS_s,'omitnan');

        rat_TT(r,d) = sess_TT_mean(sessIdx);
        rat_TR(r,d) = sess_TR_mean(sessIdx);
        rat_RR(r,d) = sess_RR_mean(sessIdx);
        rat_TS(r,d) = sess_TS_mean(sessIdx);
        rat_SS(r,d) = sess_SS_mean(sessIdx);

        if showHist
            TT_draws = [TT_draws; TT_s(:)]; %#ok<AGROW>
            TR_draws = [TR_draws; TR_s(:)]; %#ok<AGROW>
            RR_draws = [RR_draws; RR_s(:)]; %#ok<AGROW>
            TS_draws = [TS_draws; TS_s(:)]; %#ok<AGROW>
            SS_draws = [SS_draws; SS_s(:)]; %#ok<AGROW>
        end

        % B) Figure 2 data: per-cell FR means (Trace; Non-task; Non-task VF)
        [FRt_traceOnly, FRr_noVF, ~, ~] = popVecSim( ...
            spk, posMat, csTimes, win, minSpikes, false, false, velThresh, '', ...
            false, nonTaskBufferSec, speedMatchMode, speedProfileNBins);
        [~,            FRr_velF, ~, ~] = popVecSim( ...
            spk, posMat, csTimes, win, minSpikes, false, true,  velThresh, '', ...
            false, nonTaskBufferSec, speedMatchMode, speedProfileNBins);

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
mean_TS = mean(sess_TS_mean,'omitnan');
mean_SS = mean(sess_SS_mean,'omitnan');

sd_TT = std(sess_TT_mean,'omitnan');
sd_TR = std(sess_TR_mean,'omitnan');
sd_RR = std(sess_RR_mean,'omitnan');
sd_TS = std(sess_TS_mean,'omitnan');
sd_SS = std(sess_SS_mean,'omitnan');

means_t = [mean_TT, mean_TR, mean_RR, mean_TS, mean_SS];
sds_t   = [sd_TT,   sd_TR,   sd_RR,   sd_TS,   sd_SS];

fprintf('TEMPORAL grouped PV similarities (SESSION-AVERAGED; nSess=%d; groupSize=%d; nBoot/session=%d; normMode=%s; speedmatch=%d; buffer=%.2fs; matchMode=%s):\n', ...
        nSessions, groupSize, nBoot, normMode, doSpeedMatch, nonTaskBufferSec, speedMatchMode);
fprintf('  TT %.3f±%.3f | TR %.3f±%.3f | RR %.3f±%.3f | TS %.3f±%.3f | SS %.3f±%.3f (mean±SD across sessions)\n', ...
        mean_TT, sd_TT, mean_TR, sd_TR, mean_RR, sd_RR, mean_TS, sd_TS, mean_SS, sd_SS);

[a p_tt_tr c d] = ttest(sess_TT_mean, sess_TR_mean)
[a p_tt_rr c d]= ttest(sess_TT_mean, sess_RR_mean)
[a p_tt_ts c d]= ttest(sess_TT_mean, sess_TS_mean)
[a p_tt_ss c d]= ttest(sess_TT_mean, sess_SS_mean)

fprintf('ttest (paired across sessions): TT vs TR p=%.4g | TT vs RR p=%.4g | TT vs TS p=%.4g | TT vs SS p=%.4g\n', ...
    p_tt_tr, p_tt_rr, p_tt_ts, p_tt_ss);

%% ===== Figure: temporal only (hist + session-mean bars) =====
figure('Color','w','Position',[200 200 1100 350]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
edges = linspace(-1,1,50);
if showHist
    histogram(TT_draws, edges, 'Normalization','probability', 'FaceAlpha',.35);
    histogram(TR_draws, edges, 'Normalization','probability', 'FaceAlpha',.35);
    histogram(RR_draws, edges, 'Normalization','probability', 'FaceAlpha',.35);
    histogram(TS_draws, edges, 'Normalization','probability', 'FaceAlpha',.35);
    histogram(SS_draws, edges, 'Normalization','probability', 'FaceAlpha',.35);
    legend({'TT','TR','RR','TS','SS'},'Location','best');
else
    text(0.1,0.5,'showHist=false','Units','normalized');
end
xlabel('Pearson similarity'); ylabel('Probability');
title(sprintf('Within-session bootstrap draws (nSess=%d)', nSessions));

nexttile; hold on;
bar(1:5, means_t);
errorbar(1:5, means_t, sds_t, 'k.', 'LineWidth',1.5);
xticks(1:5); xticklabels({'TT','TR','RR','TS','SS'});
ylabel('Mean similarity');
title(sprintf('Session-averaged mean ± SD (nSess=%d)', nSessions));
xlim([0.5 5.5]);

ymax = max(means_t + sds_t); if ~isfinite(ymax) || ymax==0, ymax = 1; end
dy = 0.05*ymax;
pairs = [1 2; 1 3; 1 4; 1 5];
ps    = [p_tt_tr; p_tt_rr; p_tt_ts; p_tt_ss];
for k = 1:size(pairs,1)
    x1 = pairs(k,1); x2 = pairs(k,2);
    y  = ymax + k*dy;
    plot([x1 x2],[y y],'-k','LineWidth',1.2);
    if ps(k) < .05
        text(mean([x1 x2]), y + 0.02*dy, '*', 'horiz','center','FontSize',12);
    end
end

%% ===== Figure 2: per-rat + pooled FR distributions (unchanged) =====
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

    bb = bar(1:3, m, 'FaceColor','flat'); bb.CData = cols(1:3,:);
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

%% ===== per-rat page with 5 bars =====
figure('Color','w','Position',[200 100 1100 1600]);
tiledlayout(6,1,'Padding','compact','TileSpacing','compact');

condLabs = {'TT','TR','RR','TS','SS'};

for rr = 1:6
    if rr <= nR
        vTT = rat_TT(rr,:); vTR = rat_TR(rr,:); vRR = rat_RR(rr,:); vTS = rat_TS(rr,:); vSS = rat_SS(rr,:);
        rowTitle = ratNames{rr};
    else
        vTT = sess_TT_mean(:)'; vTR = sess_TR_mean(:)'; vRR = sess_RR_mean(:)'; vTS = sess_TS_mean(:)'; vSS = sess_SS_mean(:)';
        rowTitle = 'All rats (pooled sessions)';
    end

    vals = {vTT, vTR, vRR, vTS, vSS};
    m = nan(1,5); sd = nan(1,5);
    for k = 1:5
        m(k)  = mean(vals{k},'omitnan');
        sd(k) = std(vals{k}, 'omitnan');
    end

    nexttile; hold on;
    bb = bar(1:5, m, 'FaceColor','flat');
    bb.CData = [cols(1,:); cols(2,:); cols(3,:); cols(4,:); cols(5,:)];
    errorbar(1:5, m, sd, 'k.', 'LineWidth',1.2);

    nPts = numel(vTT);
    if rr <= nR
        xJ = (-0.10:0.10:0.10);
        for j = 1:nPts
            yj = [vTT(j) vTR(j) vRR(j) vTS(j) vSS(j)];
            plot((1:5) + xJ(j), yj, 'k.', 'MarkerSize',10);
        end
    else
        for j = 1:numel(sess_TT_mean)
            yj = [sess_TT_mean(j) sess_TR_mean(j) sess_RR_mean(j) sess_TS_mean(j) sess_SS_mean(j)];
            plot((1:5) + (rand(1,5)-0.5)*0.06, yj, 'k.', 'MarkerSize',7);
        end
    end

    p12 = signrank(vTT, vTR);
    p13 = signrank(vTT, vRR);
    p14 = signrank(vTT, vTS);
    p15 = signrank(vTT, vSS);

    yl = ylim; ymax = yl(2);
    dy = 0.06*(yl(2)-yl(1));
    y1 = ymax + 1*dy; y2 = ymax + 2*dy; y3 = ymax + 3*dy; y4 = ymax + 4*dy;
    plot([1 2],[y1 y1],'-k','LineWidth',1.0); if p12<.05, text(1.5,y1,'*','horiz','center'); end
    plot([1 3],[y2 y2],'-k','LineWidth',1.0); if p13<.05, text(2.0,y2,'*','horiz','center'); end
    plot([1 4],[y3 y3],'-k','LineWidth',1.0); if p14<.05, text(2.5,y3,'*','horiz','center'); end
    plot([1 5],[y4 y4],'-k','LineWidth',1.0); if p15<.05, text(3.0,y4,'*','horiz','center'); end
    ylim([yl(1), y4 + dy]);

    xticks(1:5); xticklabels(condLabs);
    ylabel('Mean similarity');
    title(sprintf('%s — session means ± SD', rowTitle));
    xlim([0.5 5.5]); box on;
end

end % end main function


%% ===== Helper: temporal FR computation (STRICT non-task candidates + speed match options) =====
function [FRt,FRr,FRsm,meanRateDay] = popVecSim( ...
    spikeCell, pos, csTimes, win, minSpikes, useTraceFilter, useNontrialFilter, ...
    velThresh, sessionLabel, doSpeedMatch, nonTaskBufferSec, speedMatchMode, speedProfileNBins)

if nargin < 9  || isempty(sessionLabel),    sessionLabel = ''; end
if nargin < 10 || isempty(doSpeedMatch),    doSpeedMatch = true; end
if nargin < 11 || isempty(nonTaskBufferSec), nonTaskBufferSec = 0; end
if nargin < 12 || isempty(speedMatchMode), speedMatchMode = 'meanstd'; end
if nargin < 13 || isempty(speedProfileNBins), speedProfileNBins = 4; end

speedMatchMode = validatestring(speedMatchMode, {'mean','meanstd','profile'});
speedProfileNBins = max(2, round(speedProfileNBins));

if ~isscalar(nonTaskBufferSec) || ~isfinite(nonTaskBufferSec) || nonTaskBufferSec < 0
    nonTaskBufferSec = 0;
end

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

% Mean rate per cell over the whole session (Hz)
sessDur = ts(end) - ts(1);
meanRateDay = nan(1, nCells);
for c = 1:nCells
    st = spikeCell{c};
    st = st(isfinite(st) & st>0);
    meanRateDay(c) = numel(st) / max(eps, sessDur);
end

% ---- Build STRICT non-task candidate windows: entire window must be outside task+buffer ----
maskCS = false(size(ts));
for t = 1:nTrials
    maskCS = maskCS | (ts >= csTimes(t)+win(1) & ts < csTimes(t)+win(2) + nonTaskBufferSec);
end

dt   = median(diff(ts),'omitnan');
bins = max(1, round(windowLen / dt));

isNonTask = ~maskCS(:);
validStart = false(numel(ts),1);
if numel(ts) >= bins
    nonTaskCount = conv(double(isNonTask), ones(bins,1), 'valid');  % length = numel(ts)-bins+1
    validStart(1:numel(nonTaskCount)) = (nonTaskCount == bins);
end

candStart = ts(validStart);
maxStart  = numel(candStart);
if maxStart < 1
    error('Not enough strictly non-task timepoints to place windows (after buffer).');
end

% ---- Precompute candidate speed features (once per session) ----
candMean = nan(maxStart,1);
candStd  = nan(maxStart,1);
candProf = nan(maxStart, speedProfileNBins);

for i = 1:maxStart
    r0i = candStart(i);
    r1i = r0i + windowLen;
    use_i = (vt >= r0i & vt < r1i);
    if any(use_i)
        vseg = vmag(use_i);
        candMean(i) = mean(vseg,'omitnan');
        candStd(i)  = std(vseg,'omitnan');

        % profile: mean speed in equal-time bins
        if strcmp(speedMatchMode,'profile')
            edges = linspace(r0i, r1i, speedProfileNBins+1);
            for b = 1:speedProfileNBins
                ub = (vt >= edges(b) & vt < edges(b+1));
                candProf(i,b) = mean(vmag(ub), 'omitnan');
            end
        end
    end
end

% For meanstd matching: z-score candidate features so mean and std are comparable
if strcmp(speedMatchMode,'meanstd')
    mu1 = mean(candMean,'omitnan'); sd1 = std(candMean,'omitnan');
    mu2 = mean(candStd ,'omitnan'); sd2 = std(candStd ,'omitnan');
    if ~isfinite(sd1) || sd1==0, sd1 = 1; end
    if ~isfinite(sd2) || sd2==0, sd2 = 1; end
    candMean_z = (candMean - mu1) ./ sd1;
    candStd_z  = (candStd  - mu2) ./ sd2;
end

FRt  = nan(nTrials, nCells);
FRr  = nan(nTrials, nCells);
FRsm = nan(nTrials, nCells);

r0_used = nan(nTrials,1);
r1_used = nan(nTrials,1);

for t = 1:nTrials
    t0 = csTimes(t) + win(1);
    t1 = t0 + windowLen;

    % --- Random non-task window (FRr) ---
    sIdx = randi(maxStart);
    r0   = candStart(sIdx);
    r1   = r0 + windowLen;
    r0_used(t) = r0;
    r1_used(t) = r1;

    % --- Speed-matched non-task window (FRsm) ---
    sm0 = NaN; sm1 = NaN;

    if doSpeedMatch
        use_t = (vt >= t0 & vt < t1);
        if any(use_t)
            vtask = vmag(use_t);
        else
            vtask = [];
        end

        if ~isempty(vtask) && any(isfinite(candMean))
            switch speedMatchMode
                case 'mean'
                    trialMean = mean(vtask,'omitnan');
                    diffs = abs(candMean - trialMean);

                case 'meanstd'
                    trialMean = mean(vtask,'omitnan');
                    trialStd  = std(vtask,'omitnan');
                    % z-score using candidate stats
                    trialMean_z = (trialMean - mu1) ./ sd1;
                    trialStd_z  = (trialStd  - mu2) ./ sd2;
                    diffs = hypot(candMean_z - trialMean_z, candStd_z - trialStd_z);

                case 'profile'
                    edgesT = linspace(t0, t1, speedProfileNBins+1);
                    profT = nan(1, speedProfileNBins);
                    for b = 1:speedProfileNBins
                        ub = (vt >= edgesT(b) & vt < edgesT(b+1));
                        profT(b) = mean(vmag(ub), 'omitnan');
                    end
                    diffs = sqrt(sum((candProf - profT).^2, 2, 'omitnan'));
            end

            diffs(~isfinite(diffs)) = inf;
            [~, ord] = sort(diffs, 'ascend');

            % random pick among top K closest
            K = min(200, numel(ord));
            pick = ord(randi(K));
            sm0 = candStart(pick);
            sm1 = sm0 + windowLen;
        end
    end

    for c = 1:nCells
        st = spikeCell{c};

        % Task window FR
        if ~useTraceFilter
            FRt(t,c) = sum(st>=t0 & st<t1) / windowLen;
        else
            [cntT, runT] = countSpikesAndRunTime(st, vt, vmag, t0, t1, velThresh);
            if runT>0, FRt(t,c) = cntT / runT; else, FRt(t,c) = NaN; end
        end

        % Random non-task FR
        if ~useNontrialFilter
            FRr(t,c) = sum(st>=r0 & st<r1) / windowLen;
        else
            [cntR, runR] = countSpikesAndRunTime(st, vt, vmag, r0, r1, velThresh);
            if runR>0, FRr(t,c) = cntR / runR; else, FRr(t,c) = NaN; end
        end

        % Speed-matched non-task FR
        if ~isfinite(sm0) || ~isfinite(sm1)
            FRsm(t,c) = NaN;
        else
            if ~useNontrialFilter
                FRsm(t,c) = sum(st>=sm0 & st<sm1) / windowLen;
            else
                [cntS, runS] = countSpikesAndRunTime(st, vt, vmag, sm0, sm1, velThresh);
                if runS>0, FRsm(t,c) = cntS / runS; else, FRsm(t,c) = NaN; end
            end
        end
    end
end

% Diagnostics on the random non-task windows used for FRr
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

% minSpikes filter based on task window counts
totalCounts = nansum(FRt,1) * windowLen;
keepCells   = totalCounts >= minSpikes;

FRt        = FRt(:,  keepCells);
FRr        = FRr(:,  keepCells);
FRsm       = FRsm(:, keepCells);
meanRateDay = meanRateDay(keepCells);
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


%% ===== temporal grouped-similarity helper (5 comparisons) =====
function [TT_vals, TR_vals, RR_vals, TS_vals, SS_vals] = pvSim_bootstrap_5( ...
    FRt, FRr, FRsm, meanRateDay, groupSize, nBoot, normMode)

FRt  = double(FRt);
FRr  = double(FRr);

if isempty(FRsm)
    FRsm = nan(size(FRt));
else
    FRsm = double(FRsm);
end

keepCells = any(isfinite(FRt) | isfinite(FRr) | isfinite(FRsm), 1);
FRt  = FRt(:, keepCells);
FRr  = FRr(:, keepCells);
FRsm = FRsm(:,keepCells);

if nargin < 4 || isempty(meanRateDay)
    meanRateDay = nan(1, size(FRt,2));
else
    meanRateDay = meanRateDay(keepCells);
end

nTrials = min([size(FRt,1), size(FRr,1), size(FRsm,1)]);
FRt  = FRt(1:nTrials,:);
FRr  = FRr(1:nTrials,:);
FRsm = FRsm(1:nTrials,:);

maxGroup = floor(nTrials/2);
if groupSize > maxGroup
    oldGS = groupSize;
    groupSize = maxGroup;
    fprintf('NOTE: groupSize clamped from %d -> %d (nTrials=%d)\n', oldGS, groupSize, nTrials);
end
if groupSize < 1
    TT_vals = []; TR_vals = []; RR_vals = []; TS_vals = []; SS_vals = [];
    return;
end

switch normMode
    case 'demean'
        X  = [FRt; FRr; FRsm];
        mu = mean(X,1,'omitnan');
        FRt  = FRt  - mu;
        FRr  = FRr  - mu;
        FRsm = FRsm - mu;

    case 'meanrate'
        mr = meanRateDay(:)'; % 1 x nCells
        mr(~isfinite(mr) | mr<=0) = NaN;
        FRt  = FRt  ./ mr;
        FRr  = FRr  ./ mr;
        FRsm = FRsm ./ mr;

    case 'none'
end

TT_vals = nan(nBoot,1);
TR_vals = nan(nBoot,1);
RR_vals = nan(nBoot,1);
TS_vals = nan(nBoot,1);
SS_vals = nan(nBoot,1);

parfor b = 1:nBoot
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

    idx = randperm(nTrials);
    g1 = idx(1:groupSize);
    g2 = idx(groupSize+1:2*groupSize);
    pv_s1 = mean(FRsm(g1,:),1,'omitnan')';
    pv_s2 = mean(FRsm(g2,:),1,'omitnan')';
    SS_vals(b) = corr(pv_s1, pv_s2, 'Rows','pairwise');

    xsel = randi(4);
    if xsel == 1, pv_t = pv_t1; pv_r = pv_r1;
    elseif xsel == 2, pv_t = pv_t2; pv_r = pv_r2;
    elseif xsel == 3, pv_t = pv_t2; pv_r = pv_r1;
    else,             pv_t = pv_t1; pv_r = pv_r2;
    end
    TR_vals(b) = corr(pv_t, pv_r, 'Rows','pairwise');

    xsel = randi(4);
    if xsel == 1, pv_t = pv_t1; pv_s = pv_s1;
    elseif xsel == 2, pv_t = pv_t2; pv_s = pv_s2;
    elseif xsel == 3, pv_t = pv_t2; pv_s = pv_s1;
    else,             pv_t = pv_t1; pv_s = pv_s2;
    end
    TS_vals(b) = corr(pv_t, pv_s, 'Rows','pairwise');
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
