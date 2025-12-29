function OUT = stateShift_postTask(ratNames, varargin)
% stateShift_postTask
% Quantify post-trial state shifts vs generic non-task running.
%
% For each rat (last pre-An day):
%   1) Build a running-only non-task template PV.
%   2) Compute Pre/Post state similarity to this template:
%        simPre_real, simPost_real, dSim_real = simPost - simPre.
%   3) Optional: do the same for pseudo-CS times (non-task control).
%   4) Compute sliding state index S(t) = cos(PV(t), template) around CS.
%   5) Plot across-rat summary for Pre vs Post and dSim.
%
% Example:
%   OUT = stateShift_postTask({'rat0222','rat0307',...}, ...
%                             'VelMatch',true,'PlotPerRat',false);

if nargin < 1 || isempty(ratNames)
    ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end

% ---------------- options ----------------
p = inputParser;
p.addParameter('preWin',[0 .85]);           % relative to CS
p.addParameter('postWin',[0.85 2.0]);      % relative to CS
p.addParameter('trialWin',[0 4]);          % for excluding from pool
p.addParameter('VelMatch',true,@islogical);
p.addParameter('velThresh',4);             % cm/s threshold for "running"
p.addParameter('minRunPre',.85/4);         % minimum run time in Pre (s)
p.addParameter('minRunPost',1.15/4);       % minimum run time in Post (s)
p.addParameter('nPseudo', 200);            % # pseudo CS times for control
p.addParameter('UsePseudo',false,@islogical);   % NEW
p.addParameter('PlotPerRat',false,@islogical);
p.addParameter('PlotSummary',true,@islogical);
p.addParameter('SimMetric','cosine', @(s) any(strcmpi(s,{'cosine','pearson'})));

% NEW: PV normalization
p.addParameter('PVNorm','demean', @(s) any(strcmpi(s,{'none','demean'})));

% sliding state index options
p.addParameter('stateWin',[-2 4]);         % sliding window relative to CS
p.addParameter('nStateBins',25);           % # bins in sliding index
p.addParameter('minRunBin',0.20);          % min run per sliding bin (s)

p.parse(varargin{:});
opt = p.Results;

nR = numel(ratNames);
OUT = struct([]);

% For across-rat summary
all_simPre_real  = nan(nR,1);
all_simPost_real = nan(nR,1);
all_dSim_real    = nan(nR,1);

all_dSim_pseudo  = nan(nR,1);  % optional
all_ratLabels    = cell(nR,1);

all_slide_real   = nan(nR, opt.nStateBins);
all_slide_pseudo = nan(nR, opt.nStateBins);  % optional

for r = 1:nR
    ratName = ratNames{r};
    rat     = evalin('base', ratName);

    % ---- choose days: last up-to-3 pre-An days ----
    dates = autoDateList(rat);
    idxAn = find(strcmp(dates, rat.An), 1);
    if isempty(idxAn)
        warning('Rat %s has no An date in autoDateList; skipping.', ratName);
        continue;
    end
    dayIdx = max(1, idxAn-2):idxAn;   % up to last 3 days

    % accumulators across days for this rat
    simPre_days    = [];
    simPost_days   = [];
    dSim_days      = [];
    nKeep_days     = [];
    slide_stack    = [];      % [nDays x nStateBins]
    postDrift_days = struct('EarlyMid',{},'EarlyLate',{},'MidLate',{});

    for di = 1:numel(dayIdx)
        d = dayIdx(di);
        dayStr = dates{d};

        spk   = rat.Ca_peaks.(sprintf('CA_peaks_%s', dayStr));
        pos   = rat.pos.(sprintf('pos_%s', dayStr));
        csTimes = rat.CS_times.(sprintf('CS_%s', dayStr));
        ratemask = rat.ratemask.(sprintf('ratemask_%s', dayStr))==1;

        if isempty(csTimes)
            warning('No CS times for %s %s; skipping this day.', ratName, dayStr);
            continue;
        end

        % restrict to kept cells
        spk = spk(ratemask,:);

        % spike → time × cell FR matrix
        t = pos(:,1);
        FR = spikesToRate(spk, t);

        % velocity and running mask
        vel  = ca_velocity(pos);
        vt   = vel(2,:)';
        vmag = vel(1,:)';
        spd  = interp1(vt, vmag, t, 'linear','extrap');
        runMask = spd >= opt.velThresh;

        % ---- build running-only non-task pool and template PV ----
        excl = [csTimes(:)+opt.trialWin(1), csTimes(:)+opt.trialWin(2)];
        sessMin = min(csTimes) - 60;
        sessMax = max(csTimes) + 60;
        pool = complementIntervals([sessMin sessMax], excl);

        poolMask = false(size(t));
        for i = 1:size(pool,1)
            poolMask = poolMask | (t>=pool(i,1) & t<pool(i,2));
        end
        runPoolMask = runMask & poolMask;

        if ~any(runPoolMask)
            warning('No running non-task pool for %s %s; skipping this day.', ratName, dayStr);
            continue;
        end

        % ---- optional per-cell demeaning using running non-task pool baseline ----
        FRuse = FR;
        if strcmpi(opt.PVNorm,'demean')
            base = mean(FR(runPoolMask,:), 1, 'omitnan');   % 1 x nCells
            FRuse = FR - base;                              % time x cell
        end

        templatePV = mean(FRuse(runPoolMask,:),1,'omitnan');

        % 1) Pre/Post similarity for this day
        [simPre_d, simPost_d, dSim_d, nKeep_d] = ...
            prePostSimilarity(FRuse, t, runMask, runPoolMask, csTimes, templatePV, opt);

        if nKeep_d == 0
            continue;
        end

        simPre_days(end+1,1)  = simPre_d;
        simPost_days(end+1,1) = simPost_d;
        dSim_days(end+1,1)    = dSim_d;
        nKeep_days(end+1,1)   = nKeep_d;

        % 2) Within-Post drift for this day (struct)
        postDrift_days(end+1) = postDriftAnalysis(FRuse, t, runMask, csTimes, opt.postWin, opt);

        % 3) Sliding state index for this day
        slide_d = slidingStateIndex_internal(FRuse, t, runMask, runPoolMask, csTimes, templatePV, opt);
        slide_stack = [slide_stack; slide_d.simReal];   %#ok<AGROW>
        last_timeBins = slide_d.timeBins;
    end

    % if no valid days, skip rat
    if isempty(simPre_days)
        warning('No valid days for rat %s; skipping.', ratName);
        continue;
    end

    % --- combine days for this rat (trial-count-weighted average) ---
    totalTrials = sum(nKeep_days);
    w = nKeep_days / totalTrials;

    simPre_real  = sum(simPre_days  .* w);
    simPost_real = sum(simPost_days .* w);
    dSim_real    = simPost_real - simPre_real;

    % average sliding index across days
    slide.simReal   = mean(slide_stack, 1, 'omitnan');
    slide.simPseudo = nan(1, opt.nStateBins);   % unused if UsePseudo=false
    slide.timeBins  = last_timeBins;
    slide.opt       = opt;

    % average postDrift fields across days
    EM  = [postDrift_days.EarlyMid];
    EL  = [postDrift_days.EarlyLate];
    ML  = [postDrift_days.MidLate];
    postDrift.EarlyMid  = mean(EM,'omitnan');
    postDrift.EarlyLate = mean(EL,'omitnan');
    postDrift.MidLate   = mean(ML,'omitnan');

    % ---- store in OUT ----
    OUT(r).ratName        = ratName;
    OUT(r).days           = dates(dayIdx);      % cell array of 2–3 days
    OUT(r).simPre_real    = simPre_real;
    OUT(r).simPost_real   = simPost_real;
    OUT(r).dSim_real      = dSim_real;
    OUT(r).nReal_keep     = totalTrials;

    OUT(r).simPre_pseudo  = NaN;
    OUT(r).simPost_pseudo = NaN;
    OUT(r).dSim_pseudo    = NaN;
    OUT(r).nPseudo_keep   = NaN;

    OUT(r).postDrift      = postDrift;
    OUT(r).slide          = slide;

    all_simPre_real(r)  = simPre_real;
    all_simPost_real(r) = simPost_real;
    all_dSim_real(r)    = dSim_real;
    all_ratLabels{r}    = ratName;

    fprintf('\n=== %s (last %d pre-An days): stateShift_postTask (VelMatch=%d, PVNorm=%s) ===\n', ...
        ratName, numel(EM), opt.VelMatch, opt.PVNorm);
    fprintf('Real trials (combined days): simPre=%.3f, simPost=%.3f, dSim=%.3f (n=%d kept)\n', ...
        simPre_real, simPost_real, dSim_real, totalTrials);

    if opt.PlotPerRat
        plotPerRatSummary(OUT(r));
    end
end

% ============================================
% Across-rat summary
% ============================================
if opt.PlotSummary
    valid = ~isnan(all_simPre_real) & ~isnan(all_simPost_real);
    simPre = all_simPre_real(valid);
    simPost= all_simPost_real(valid);
    dSimR  = all_dSim_real(valid);
    ratLabels = all_ratLabels(valid); %#ok<NASGU>
    nGood  = numel(simPre);

    if nGood > 0
        % ----- stats -----
        [~,p_pair,~,stats_pair] = ttest(simPre, simPost);
        diffVals = simPost - simPre;
        mPre  = mean(simPre,'omitnan');
        mPost = mean(simPost,'omitnan');
        mDiff = mean(diffVals,'omitnan');
        sDiff = std(diffVals,'omitnan');
        d_pair = mDiff / sDiff;

        [~,p_diff,~,stats_diff] = ttest(dSimR,0);

        fprintf('\n=== stateShift_postTask across rats ===\n');
        fprintf('Pre vs Post similarity (paired t-test): t(%d)=%.3f, p=%.4g, meanPre=%.3f, meanPost=%.3f, d=%.3f\n', ...
            nGood-1, stats_pair.tstat, p_pair, mPre, mPost, d_pair);
        fprintf('Post-Pre difference vs 0 (one-sample t-test): t(%d)=%.3f, p=%.4g, meanΔ=%.3f\n', ...
            nGood-1, stats_diff.tstat, p_diff, mDiff);

        figure('Color','w','Position',[200 200 900 350]);
        tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

        % -------- Left: paired Pre vs Post --------
        nexttile; hold on;

        xPre  = ones(nGood,1);
        xPost = 2*ones(nGood,1);

        for i = 1:nGood
            plot([xPre(i) xPost(i)], [simPre(i) simPost(i)], ...
                '-','Color',[0.7 0.7 0.7]);
        end

        plot(xPre,  simPre,  'o', ...
            'MarkerFaceColor',[0.1 0.1 0.1], 'MarkerEdgeColor','none');
        plot(xPost, simPost, 'o', ...
            'MarkerFaceColor',[0.0 0.45 0.9], 'MarkerEdgeColor','none');

        mPre  = mean(simPre,'omitnan');
        mPost = mean(simPost,'omitnan');
        sPre  = std(simPre,'omitnan')/sqrt(nGood);
        sPost = std(simPost,'omitnan')/sqrt(nGood);

        errorbar(1, mPre,  sPre,  'k','LineWidth',1.5,'CapSize',0);
        errorbar(2, mPost, sPost, 'k','LineWidth',1.5,'CapSize',0);

        plot([1 2],[mPre mPost],'k-','LineWidth',2);

        xlim([0.5 2.5]);
        set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post'});
        ylabel('Similarity to running template');
        title('State similarity before vs after trial');

        % -------- Right: dSim distribution (Post–Pre) --------
        nexttile; hold on;

        x0 = 1;
        xj = x0 + 0.05*randn(size(dSimR));
        plot(xj, dSimR, 'o', ...
            'MarkerFaceColor',[0 0.2 0.6], 'MarkerEdgeColor','none');

        mD = mean(dSimR,'omitnan');
        sD = std(dSimR,'omitnan')/sqrt(nGood);
        plot([0.8 1.2],[mD mD],'k-','LineWidth',2);
        plot([1 1],[mD-sD mD+sD],'k-','LineWidth',1.5);

        yline(0,'k:');
        xlim([0.5 1.5]);
        set(gca,'XTick',1,'XTickLabel',{'Post - Pre'});
        ylabel('dSim (Post - Pre)');
        title(sprintf('Post state closer to running (p = %.3f)', p_diff));
    end
end

% --- sliding-state index summary (Real only by default) ---
if opt.PlotSummary
    valid = all(~isnan(all_slide_real),2);
    slideR = all_slide_real(valid,:);
    if ~isempty(slideR)
        tBins = OUT(find(valid,1,'first')).slide.timeBins;
        mR = mean(slideR,1,'omitnan');
        sR = std(slideR,[],1,'omitnan')/sqrt(size(slideR,1));

        figure('Color','w','Position',[200 600 900 300]); hold on;
        shadedErrorBar_local(tBins, mR, sR, 'b');
        xlabel('Time from CS (s)');
        ylabel('Similarity to running template');
        yline(0,'k:');
        title('Sliding state index (real CS, across rats)');
    end
end

end % main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------- Helper functions -------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FR = spikesToRate(spk, t)
dt = median(diff(t),'omitnan');
nT = numel(t);
nC = size(spk,1);
FR = zeros(nT,nC);
for c = 1:nC
    st = spk(c,:);
    st = st(~isnan(st) & st>0);
    if isempty(st), continue; end
    idx = discretize(st, t);
    idx = idx(isfinite(idx) & idx>=1 & idx<=nT);
    FR(idx,c) = FR(idx,c) + 1;
end
FR = FR ./ dt;
end

function s = pvSimilarity(a,b,metric)
% pvSimilarity  Population-vector similarity: cosine or Pearson r
a = a(:); b = b(:);

switch lower(metric)
    case 'cosine'
        s = dot(a,b) / (norm(a)*norm(b));
        if ~isfinite(s), s = NaN; end
    case 'pearson'
        if all(isnan(a)) || all(isnan(b))
            s = NaN;
        else
            s = corr(a, b, 'rows','pairwise');
        end
    otherwise
        error('Unknown similarity metric: %s', metric);
end
end

function [simPre, simPost, dSim, nKeep] = prePostSimilarity(FR, t, runMask, runPoolMask, csTimes, templatePV, opt) %#ok<INUSD>
nTrials = numel(csTimes);

pre_pv  = [];
post_pv = [];
keepCount = 0;

for k = 1:nTrials
    % Pre/Post windows
    t0_pre  = csTimes(k) + opt.preWin(1);
    t1_pre  = csTimes(k) + opt.preWin(2);
    t0_post = csTimes(k) + opt.postWin(1);
    t1_post = csTimes(k) + opt.postWin(2);

    mask_pre  = t>=t0_pre  & t<t1_pre;
    mask_post = t>=t0_post & t<t1_post;

    if opt.VelMatch
        mask_pre  = mask_pre  & runMask;
        mask_post = mask_post & runMask;
    end

    runPre  = sum(diff(t(mask_pre)));
    runPost = sum(diff(t(mask_post)));

    if runPre < opt.minRunPre || runPost < opt.minRunPost
        continue;
    end

    keepCount = keepCount + 1;
    pre_pv(keepCount,:)  = mean(FR(mask_pre,:),1,'omitnan');  %#ok<AGROW>
    post_pv(keepCount,:) = mean(FR(mask_post,:),1,'omitnan'); %#ok<AGROW>
end

if keepCount == 0
    simPre = NaN; simPost = NaN; dSim = NaN; nKeep = 0; return;
end

simPre  = mean(arrayfun(@(i) pvSimilarity(pre_pv(i,:),  templatePV, opt.SimMetric), 1:keepCount),'omitnan');
simPost = mean(arrayfun(@(i) pvSimilarity(post_pv(i,:), templatePV, opt.SimMetric), 1:keepCount),'omitnan');

dSim  = simPost - simPre;
nKeep = keepCount;
end

function [simPre, simPost, dSim, nKeep] = prePostSimilarity_pseudo(FR, t, runMask, runPoolMask, csTimes, templatePV, opt) %#ok<INUSD>
validTimes = t(runPoolMask);
if numel(validTimes) < opt.nPseudo
    pseudoCS = validTimes;
else
    pseudoCS = datasample(validTimes, opt.nPseudo, 'Replace', false);
end
nPseudo = numel(pseudoCS);

pre_pv  = [];
post_pv = [];
keepCount = 0;

for k = 1:nPseudo
    cs = pseudoCS(k);

    t0_pre  = cs + opt.preWin(1);
    t1_pre  = cs + opt.preWin(2);
    t0_post = cs + opt.postWin(1);
    t1_post = cs + opt.postWin(2);

    mask_pre  = t>=t0_pre  & t<t1_pre;
    mask_post = t>=t0_post & t<t1_post;

    if opt.VelMatch
        mask_pre  = mask_pre  & runMask;
        mask_post = mask_post & runMask;
    end

    runPre  = sum(diff(t(mask_pre)));
    runPost = sum(diff(t(mask_post)));

    if runPre < opt.minRunPre || runPost < opt.minRunPost
        continue;
    end

    keepCount = keepCount + 1;
    pre_pv(keepCount,:)  = mean(FR(mask_pre,:),1,'omitnan');  %#ok<AGROW>
    post_pv(keepCount,:) = mean(FR(mask_post,:),1,'omitnan'); %#ok<AGROW>
end

if keepCount == 0
    simPre = NaN; simPost = NaN; dSim = NaN; nKeep = 0; return;
end

simPre  = mean(arrayfun(@(i) pvSimilarity(pre_pv(i,:),  templatePV, opt.SimMetric), 1:keepCount),'omitnan');
simPost = mean(arrayfun(@(i) pvSimilarity(post_pv(i,:), templatePV, opt.SimMetric), 1:keepCount),'omitnan');
dSim    = simPost - simPre;
nKeep   = keepCount;
end

function postDrift = postDriftAnalysis(FR, t, runMask, csTimes, postWin, opt)
% divide postWin into Early/Mid/Late thirds and compare PVs
t0 = postWin(1); t1 = postWin(2);
edges = linspace(t0, t1, 4); % 3 segments
segPV = nan(3, size(FR,2));

for s = 1:3
    mask = false(size(t));
    for k = 1:numel(csTimes)
        a = csTimes(k) + edges(s);
        b = csTimes(k) + edges(s+1);
        mask = mask | (t>=a & t<b);
    end
    mask = mask & runMask;
    segPV(s,:) = mean(FR(mask,:),1,'omitnan');
end

postDrift.EarlyMid  = pvSimilarity(segPV(1,:), segPV(2,:), opt.SimMetric);
postDrift.EarlyLate = pvSimilarity(segPV(1,:), segPV(3,:), opt.SimMetric);
postDrift.MidLate   = pvSimilarity(segPV(2,:), segPV(3,:), opt.SimMetric);
end

function slide = slidingStateIndex_internal(FR, t, runMask, runPoolMask, csTimes, templatePV, opt)
% Sliding state index around real CS vs pseudo CS

% --- time bins ---
edges = linspace(opt.stateWin(1), opt.stateWin(2), opt.nStateBins+1);
bStart = edges(1:end-1);
bEnd   = edges(2:end);
timeBins = (bStart + bEnd)/2;

% pseudo CS from running pool
validTimes = t(runPoolMask);
if numel(validTimes) < opt.nPseudo
    pseudoCS = validTimes;
else
    pseudoCS = datasample(validTimes, opt.nPseudo, 'Replace', false);
end

simReal   = nan(1,opt.nStateBins);
simPseudo = nan(1,opt.nStateBins);

for bi = 1:opt.nStateBins
    % REAL
    maskR = false(size(t));
    for k = 1:numel(csTimes)
        a = csTimes(k) + bStart(bi);
        bnd = csTimes(k) + bEnd(bi);
        maskR = maskR | (t>=a & t<bnd);
    end
    maskR = maskR & runMask;
    runDurR = sum(diff(t(maskR)));
    if runDurR >= opt.minRunBin
        pvR = mean(FR(maskR,:),1,'omitnan');
        simReal(bi) = pvSimilarity(pvR, templatePV, opt.SimMetric);
    end

    % PSEUDO
    maskP = false(size(t));
    for k = 1:numel(pseudoCS)
        a = pseudoCS(k) + bStart(bi);
        bnd = pseudoCS(k) + bEnd(bi);
        maskP = maskP | (t>=a & t<bnd);
    end
    maskP = maskP & runMask;
    runDurP = sum(diff(t(maskP)));
    if runDurP >= opt.minRunBin
        pvP = mean(FR(maskP,:),1,'omitnan');
        simPseudo(bi) = pvSimilarity(pvP, templatePV, opt.SimMetric);
    end
end

slide.simReal   = simReal;
slide.simPseudo = simPseudo;
slide.timeBins  = timeBins;
slide.opt       = opt;
end

function plotPerRatSummary(R)
% Quick per-rat diagnostic figure (optional)
s  = R.slide;
tb = s.timeBins;

figure('Color','w','Position',[200 200 800 350]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
plot(tb, s.simReal,'-o','LineWidth',1.5);
plot(tb, s.simPseudo,'--s','LineWidth',1.5);
yline(0,'k:');
xlabel('Time from CS (s)');
ylabel('Similarity to template');
legend({'Real','Pseudo'},'Location','best');
title(sprintf('%s sliding state index', R.ratName));

nexttile; hold on;
bar(1,[R.dSim_real],'FaceColor',[0.2 0.6 0.9]);
bar(2,[R.dSim_pseudo],'FaceColor',[0.9 0.4 0.2]);
xlim([0.5 2.5]); set(gca,'XTick',[1 2],'XTickLabel',{'Real','Pseudo'});
yline(0,'k:');
ylabel('dSim (Post - Pre)');
title('dSim');
end

function shadedErrorBar_local(x, m, s, colorChar)
% simple shaded error helper (no toolbox)
if ischar(colorChar)
    switch colorChar
        case 'b', col = [0.2 0.4 0.9];
        case 'r', col = [0.9 0.3 0.2];
        otherwise, col = [0.3 0.3 0.3];
    end
else
    col = colorChar;
end
fill([x fliplr(x)], [m-s fliplr(m+s)], col, ...
    'FaceAlpha',0.2,'EdgeColor','none');
plot(x, m, 'Color', col, 'LineWidth',2);
end

function C = complementIntervals(full, excl)
a = full(1); b = full(2);
E = excl;
E = E(~any(isnan(E),2),:);
E = sortrows(E,1);
E(:,1) = max(E(:,1),a);
E(:,2) = min(E(:,2),b);
E = E(E(:,1)<E(:,2),:);
C = [];
t = a;
for i = 1:size(E,1)
    if E(i,1) > t
        C = [C; t E(i,1)]; %#ok<AGROW>
    end
    t = max(t, E(i,2));
end
if t < b
    C = [C; t b]; %#ok<AGROW>
end
end
