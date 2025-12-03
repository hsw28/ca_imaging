function R = run_accelBinMatched(ratNames, varargin)
% run_accelBinMatched
% R is a struct array per rat with per-day stats and pooled per-rat summary (ACCEL bins).
% Usage:
%   R = run_accelBinMatched([], 'AccelMode','signed', 'AccelBinWidth', 0.2, ...)

if nargin<1 || isempty(ratNames)
    ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end

nRats = numel(ratNames);

% defaults, override via varargin
p = inputParser;
p.addParameter('AccelMode','signed', @(s) any(validatestring(lower(s),{'signed','magnitude'}))); % 'signed' = d|v|/dt; 'magnitude' = ||a||
p.addParameter('AccelBinWidth',3);      % acceleration units (e.g., cm/s^2)
p.addParameter('AccelEdges',[]);
p.addParameter('MinDurPerBin',(1/7.5));       % seconds required per side in each accel bin
p.addParameter('MinBins',1);
p.addParameter('alpha',0.05);
p.addParameter('test','ttest');         % 'ttest' | 'signrank'
p.addParameter('win',[0 2]);            % trial window
p.addParameter('binSize',1);            % frame/bin width (s)
p.parse(varargin{:});
P = p.Results;

abort_loops = false; % Initialize a flag variable

R = struct('rat',[], 'perDay',[], 'pctSig',[], 'deltaRate_sig',[]);
for r=1:nRats
    ratVar = ratNames{r};
    rat = evalin('base', ratVar);
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    if isempty(idx) || idx<3
        warning('%s does not have >=3 days ending at An. Skipping...', ratVar);
        continue;
    end
    days  = dates(idx-2:idx);

    perDay = cell(1,3);
    pct    = nan(1,3);
    delSig = [];
    medp   = [];

    for d=1:3
        day      = days{d};
        spk      = rat.Ca_peaks.(['CA_peaks_' day]);
        posFull  = smoothpos(rat.pos.(['pos_' day]));
        ts       = posFull(:,1);
        xy       = posFull(:,2:3);
        cs       = rat.CS_times.(['CS_' day]);
        mask     = rat.ratemask.(['ratemask_' day]);

        S = accelBinMatched(spk, ts, xy, cs, mask, ...
            'win',P.win, 'binSize',P.binSize, ...
            'AccelEdges',P.AccelEdges, 'AccelBinWidth',P.AccelBinWidth, ...
            'AccelMode',P.AccelMode, ...
            'MinDurPerBin',P.MinDurPerBin, 'MinBins',P.MinBins, ...
            'alpha',P.alpha, 'test',P.test);

        perDay{d} = S;
        if any(isfinite(S.pVal))
            pct(d) = S.pctSigDiffFDR;
            delSig = [delSig; S.deltaRate(S.sigDiffFDR)]; %#ok<AGROW>
        end

        fprintf('[%s %s] tested=%d, sig=%d (%.1f%%), median p=%.3g\n', ...
            ratVar, day, sum(isfinite(S.pVal)), sum(S.sigDiffFDR), pct(d), ...
            median(S.pVal(isfinite(S.pVal)),'omitnan'));

        medp(end+1) = median(S.pVal(isfinite(S.pVal)),'omitnan');

        % Early exit criteria (commented out to match current speed version behavior)
        % if isnan(pct(d)) || pct(d)==0 || medp(end)>0.15
        %     abort_loops = true;
        %     fprintf('breaking loop bc not enough points or med p too high\n');
        %     break;
        % end
    end

    if abort_loops
        break; % Exit the outer loop
    end

    R(r).rat           = ratVar;
    R(r).perDay        = perDay;
    R(r).pctSig        = pct;
    R(r).deltaRate_sig = delSig;
end














% === Population summaries + PLOTS (ACCEL version) ===
fprintf('Running ACCEL population summaries on %d rats...\n', numel(R));

Pop      = accel_population_summary(       R, 'nPerm', 500, 'useSignrank', false, ...
                                           'titleTag', 'Accel-bin analysis'); %#ok<NASGU>
PopSplit_POP = accel_population_summary_split_POP(R, ...
    'nPerm', 500, 'titleTag', 'Accel-bin analysis'); %#ok<NASGU>

    % Early-exit check (same pattern as before; based on population summary)
    if isstruct(PopSplit_POP) && isfield(PopSplit_POP,'early_exit') && PopSplit_POP.early_exit
        fprintf('%s\n', PopSplit_POP.reason);
        R = NaN;
        return;
    end

% Per-cell sig panels
speedBin_cellwise_sigpanels(R, 'titleTag', 'Accel-bin analysis');
speedBin_population_taskNon_violin(R, 'titleTag','Accel-bin analysis', ...
    'ShowData', true, 'DataAlpha', 0.15, 'Bandwidth', 0.02);

plot_per_rat_means(R, 'Accel-bin analysis');
plot_delta_only_fromR(R, 'Accel-bin analysis');

plot_paired_means_slopegraph(R, 'Accel-bin analysis');
plot_paired_means_slopegraph_POP(R,  'Accel-bin analysis (population POP)'); % NEW

% Change bars (percent + fold), with SEM like in speed version
plot_rate_change_bars(R, 'Accel-bin analysis', 'metric','percent', 'errType','sem');
plot_rate_change_bars(R, 'Accel-bin analysis', 'metric','fold',    'errType','sem');

% POPULATION (speed-matched) bars + text
plot_rate_change_bars_POP(R, 'Accel-bin analysis', 'metric','percent', 'errType','sem');
plot_rate_change_bars_POP(R, 'Accel-bin analysis', 'metric','fold',   'errType','sem');

speedBin_dayLevel_heatmap(R, 'titleTag','Speed-bin analysis');
speedBin_dayLevel_bars(R,    'titleTag','Speed-bin analysis');

% Combo panel: Δ violin, per-rat paired violins, permutation null
speedBin_population_combo_plots(R, ...
    'nPerm', 500, ...
    'useSignrank', false, ...
    'titleTag', 'Accel-bin analysis');

% Per-animal mean %sig + across bar
speedBin_sigbar_acrossAnimals(R, 'titleTag','Accel-bin analysis');

end   % <-- end of run_accelBinMatched


% ======================================================================
function stats = accelBinMatched(spikeCell, ts, pos, csTimes, ratemask, varargin)
% accelBinMatched
% Compare TRIAL vs NON-TRIAL firing **rates** matched within ACCELERATION BINS.
% One paired test per cell using paired accel bins (rate_trial_bin vs rate_nontrial_bin).
%
% AccelMode:
%   'signed'    : d|v|/dt (change in speed; can be negative/positive)
%   'magnitude' : ||a|| = sqrt((dvx/dt)^2 + (dvy/dt)^2), always >=0

ip = inputParser();
ip.addParameter('win',            [0 2]);
ip.addParameter('binSize',         1/7.5);
ip.addParameter('AccelEdges',      []);
ip.addParameter('AccelBinWidth',   0.1);
ip.addParameter('AccelMode',       'signed', @(s) any(validatestring(lower(s),{'signed','magnitude'})));
ip.addParameter('MinDurPerBin',    0.5);
ip.addParameter('MinBins',         3);
ip.addParameter('alpha',           0.05);
ip.addParameter('test',            'ttest', @(s) any(validatestring(s,{'ttest','signrank'})));
ip.parse(varargin{:});
win           = ip.Results.win;
binSize       = ip.Results.binSize;
AccelEdges    = ip.Results.AccelEdges;
AccelBinWidth = ip.Results.AccelBinWidth;
AccelMode     = lower(ip.Results.AccelMode);
MinDurPerBin  = ip.Results.MinDurPerBin;
MinBins       = ip.Results.MinBins;
alpha         = ip.Results.alpha;
testType      = lower(ip.Results.test);

if ~iscell(spikeCell) && isnumeric(spikeCell)
    spikeCell = mat2cell(spikeCell, ones(size(spikeCell,1),1), size(spikeCell,2));
end
nCells = numel(spikeCell);
if nargin < 5 || isempty(ratemask), ratemask = true(nCells,1); end
ratemask = logical(ratemask(:));
if numel(ratemask) < nCells, ratemask(end+1:nCells) = true; end

% ---------- drop counters (for debugging/QA) ----------
drop = struct();
drop.nCells          = numel(spikeCell);
drop.nMasked         = sum(~logical(ratemask(:)));
drop.nEmpty          = 0;
drop.nBelowMinBins   = 0;
drop.nTested         = 0;
drop.pairedBinsAvail = 0;
drop.MinBins         = MinBins;
drop.MinDurPerBin    = MinDurPerBin;

% velocities from pos/ts
dt  = [diff(ts); binSize];
dx  = [diff(pos(:,1)); 0];
dy  = [diff(pos(:,2)); 0];
vx  = dx ./ dt;
vy  = dy ./ dt;

% speeds
spd = hypot(vx, vy);
spd(~isfinite(spd)) = 0;

% accelerations
switch AccelMode
    case 'signed'      % d|v|/dt (change in speed)
        acc = [diff(spd)./dt(1:end-1); 0];
    case 'magnitude'   % true vector acceleration magnitude
        ax  = [diff(vx)./dt(1:end-1); 0];
        ay  = [diff(vy)./dt(1:end-1); 0];
        acc = hypot(ax, ay);
end
acc(~isfinite(acc)) = 0;

% trial mask
inTrial = false(size(ts));
for t = 1:numel(csTimes)
    inTrial = inTrial | (ts >= csTimes(t)+win(1) & ts < csTimes(t)+win(2));
end
outMask = ~inTrial;

% robust per-sample dt (fallback to binSize on the last sample)
dts   = diff(ts);
dts   = dts(isfinite(dts) & dts > 0);
dtEff = isempty(dts) * binSize + ~isempty(dts) * median(dts);
dtVec = [diff(ts); dtEff];

dx = [diff(pos(:,1)); 0];
dy = [diff(pos(:,2)); 0];
vx = dx ./ max(dtVec, eps);
vy = dy ./ max(dtVec, eps);

spd = hypot(vx, vy); spd(~isfinite(spd)) = 0;

switch AccelMode
  case 'signed'      % d|v|/dt (change in speed)
    acc = [diff(spd) ./ max(dtVec(1:end-1), eps); 0];
  case 'magnitude'   % ||a|| from vx,vy
    ax  = [diff(vx) ./ max(dtVec(1:end-1), eps); 0];
    ay  = [diff(vy) ./ max(dtVec(1:end-1), eps); 0];
    acc = hypot(ax, ay);
end
acc(~isfinite(acc)) = 0;

% acceleration binning
if isempty(AccelEdges)
    aMin = floor(min(acc));
    aMax = ceil(max(acc));
    if aMin==aMax
        aMax = aMin + AccelBinWidth;
    end
    AccelEdges = aMin:AccelBinWidth:aMax;
    if numel(AccelEdges)<2, AccelEdges = [aMin aMin+AccelBinWidth]; end
end
[~,~,binIdx] = histcounts(acc, AccelEdges);
binCenters   = (AccelEdges(1:end-1)+AccelEdges(2:end))/2;
nBins        = numel(binCenters);

% ---------- durations per bin (real seconds; no overlap) ----------
okT = inTrial  & (binIdx > 0);
okN = outMask  & (binIdx > 0);
durTrial = accumarray(binIdx(okT), dtVec(okT), [nBins 1], @sum, 0);   % s in each speed bin during TRIAL
durNon   = accumarray(binIdx(okN), dtVec(okN), [nBins 1], @sum, 0);   % s in each speed bin during NON-TRIAL

% ---------- union-of-time intervals per bin (for spike counting) ----------
trialIntervals = cell(nBins,1);
nonIntervals   = cell(nBins,1);
dtMedLocal = isempty(dts) * binSize + ~isempty(dts) * median(dts);
for b = 1:nBins
    trialIntervals{b} = maskToIntervals(ts, okT & (binIdx==b), dtMedLocal);
    nonIntervals{b}   = maskToIntervals(ts, okN & (binIdx==b), dtMedLocal);
end

% ---------- bins that qualify on BOTH sides ----------
keepTemplate         = (durTrial >= MinDurPerBin) & (durNon >= MinDurPerBin);
drop.pairedBinsAvail = sum(keepTemplate);

% POPULATION accumulators (across cells, per accel bin)
popTrialSpikes = zeros(nBins,1);   % total spikes in trial intervals per bin
popNonSpikes   = zeros(nBins,1);   % total spikes in non-trial intervals per bin

% outputs (per-cell stats)
pVal             = nan(nCells,1);
sigDiffFDR       = false(nCells,1);
nBinsUsed        = zeros(nCells,1);
meanTrialRate    = nan(nCells,1);
meanNonTrialRate = nan(nCells,1);
deltaRate        = nan(nCells,1);
binRatesTrial    = cell(nCells,1);
binRatesNon      = cell(nCells,1);
accBinsUsed      = cell(nCells,1);

% per-cell counting
for c = 1:nCells
    if ~ratemask(c), continue; end
    st = spikeCell{c};
    if isempty(st)
        drop.nEmpty = drop.nEmpty + 1;
        continue;
    end
    st = st(:);

    % spikes in each bin, using interval counting
    spkTrial = zeros(nBins,1);
    spkNon   = zeros(nBins,1);
    for b = 1:nBins
        if durTrial(b) > 0
            spkTrial(b) = countInIntervals(st, trialIntervals{b});
        end
        if durNon(b) > 0
            spkNon(b) = countInIntervals(st, nonIntervals{b});
        end
    end

    % accumulate for POPULATION stats (only included cells)
    popTrialSpikes = popTrialSpikes + spkTrial;
    popNonSpikes   = popNonSpikes   + spkNon;

    % rates per bin (Hz) for eligible paired bins
    keep  = keepTemplate;                 % same MinDurPerBin filter for everyone
    rateT = zeros(nBins,1);
    rateM = zeros(nBins,1);
    rateT(keep) = spkTrial(keep) ./ max(durTrial(keep), eps);
    rateM(keep) = spkNon(keep)   ./ max(durNon(keep),   eps);

    x = rateT(keep);
    y = rateM(keep);
    informative = ~(x==0 & y==0);         % drop bins where both sides are 0
    x = x(informative);  y = y(informative);
    binsUsed = find(keep);
    binsUsed = binsUsed(informative);

    if numel(x) < MinBins
        drop.nBelowMinBins = drop.nBelowMinBins + 1;
        continue;                          % not enough paired bins to test this cell
    end

    % paired test across bins
    switch testType
        case 'ttest',   [~, p] = ttest(x, y);
        otherwise,      p = signrank(x, y);
    end

    pVal(c) = p;
    nBinsUsed(c) = numel(x);
    drop.nTested = drop.nTested + 1;

    % duration-unweighted means across used bins
    meanTrialRate(c)    = mean(x,'omitnan');
    meanNonTrialRate(c) = mean(y,'omitnan');
    deltaRate(c)        = mean(x - y,'omitnan');

    binRatesTrial{c}    = x;
    binRatesNon{c}      = y;
    accBinsUsed{c}      = binCenters(binsUsed);
end

% ---------- POPULATION rates per accel bin (matched) ----------
pop_keep = keepTemplate & (durTrial > 0) & (durNon > 0);
if any(pop_keep)
    % current POP rates are *sums* across cells (Hz); we'll convert to per-cell
    trialRatePerBin_sum = popTrialSpikes(pop_keep) ./ max(durTrial(pop_keep), eps);
    nonRatePerBin_sum   = popNonSpikes(pop_keep)   ./ max(durNon(pop_keep),   eps);

    nCells_pop = sum(ratemask);               % # cells included in POP
    if nCells_pop < 1
        nCells_pop = 1;
    end

    % store both, but downstream code will use the per-cell versions
    trialRatePerBin = trialRatePerBin_sum(:) / nCells_pop;   % Hz per cell
    nonRatePerBin   = nonRatePerBin_sum(:)   / nCells_pop;   % Hz per cell
    deltaPerBin     = trialRatePerBin - nonRatePerBin;

    stats.pop = struct( ...
        'trialRatePerBin', trialRatePerBin(:), ...   % Hz per cell
        'nonRatePerBin',   nonRatePerBin(:),   ...   % Hz per cell
        'deltaPerBin',     deltaPerBin(:),     ...
        'accelBinCenters', binCenters(pop_keep), ...
        'nCells',          nCells_pop,        ...   % keep this for reference
        'trialRatePerBin_sum', trialRatePerBin_sum(:), ... % optional
        'nonRatePerBin_sum',   nonRatePerBin_sum(:));      % optional
else
    stats.pop = [];
end

% FDR across tested cells
tested = isfinite(pVal);
sig = false(sum(tested),1);
if any(tested)
    [~,~,adjP] = fdr_bh_local(pVal(tested), alpha);
    sig = adjP < alpha;
end
sigDiffFDR(tested) = sig;

stats.pVal               = pVal;
stats.sigDiffFDR         = sigDiffFDR;
stats.pctSigDiffFDR      = 100*mean(sigDiffFDR(tested));
stats.nBinsUsed          = nBinsUsed;
stats.meanTrialRate      = meanTrialRate;
stats.meanNonTrialRate   = meanNonTrialRate;
stats.deltaRate          = deltaRate;
stats.binRatesTrial      = binRatesTrial;
stats.binRatesNon        = binRatesNon;
stats.accelBinCenters    = accBinsUsed;
stats.AccelEdges         = AccelEdges;
stats.params             = struct('win',win,'binSize',binSize,'AccelEdges',AccelEdges, ...
                                  'AccelBinWidth',AccelBinWidth,'AccelMode',AccelMode, ...
                                  'MinDurPerBin',MinDurPerBin,'MinBins',MinBins, ...
                                  'alpha',alpha,'test',testType);
end

% ======================================================================
function [h, crit_p, adj_p] = fdr_bh_local(pvals, q)
p = pvals(:); m = numel(p);
[sp, idx] = sort(p);
thresh = (1:m)'/m * q;
rej = find(sp<=thresh,1,'last');
if isempty(rej), crit_p=0; h=false(m,1);
else, crit_p=sp(rej); h = p<=crit_p; end
adj_p = nan(m,1);
for i=m:-1:1
    if i<m, adj_p(i) = min(sp(i)*m/i, adj_p(i+1));
    else,   adj_p(i) = min(sp(i)*m/i, 1);
    end
end
tmp = adj_p; adj_p(idx) = tmp;
end


% ======================================================================

function Pop = accel_population_summary(R, varargin)
% Population-level inference with ONE p (no FDR):
% (A) one-sample t/signrank on per-cell deltaRate
% (B) permutation (sign-flip) test on the mean deltaRate
%
% Also draws:
%   Fig A: violin+swarm of all per-cell deltaRate; per-rat means w/ CI
%   Fig B: permutation null of mean deltaRate with observed line

p = inputParser;
p.addParameter('nPerm', 500);
p.addParameter('useSignrank', false);  % if true, report signrank instead of t-test
p.addParameter('titleTag', '');
p.parse(varargin{:});
nPerm       = p.Results.nPerm;
useSignrank = p.Results.useSignrank;
titleTag    = p.Results.titleTag;

% ---- collect per-cell deltas (tested cells only) & per-rat summaries ----
delta_all    = [];
ratNames     = {R.rat};
perRat_means = nan(numel(R),1);
perRat_ns    = zeros(numel(R),1);

for rr = 1:numel(R)
    perDay      = R(rr).perDay;
    del_thisRat = [];
    for d = 1:numel(perDay)
        S = perDay{d};
        if isempty(S), continue; end
        tested = isfinite(S.pVal);
        del    = S.deltaRate(tested);
        del_thisRat = [del_thisRat; del]; %#ok<AGROW>
    end
    perRat_ns(rr)    = numel(del_thisRat);
    perRat_means(rr) = mean(del_thisRat,'omitnan');
    delta_all        = [delta_all; del_thisRat]; %#ok<AGROW>
end

% --- NEW: handle no usable deltas gracefully ---
if isempty(delta_all)
    fprintf('[accel_population_summary] No cells with finite deltaRate for %s — skipping population t-test / perms.\n', titleTag);
    Pop = struct();
    Pop.early_exit   = true;
    Pop.reason       = 'no_delta_data';
    Pop.delta_all    = [];
    Pop.perRat_means = perRat_means;
    Pop.perRat_ns    = perRat_ns;
    return;
end

% ---- A) One-sample test on per-cell deltas ----
if useSignrank
    [pA,~,statsA] = signrank(delta_all, 0);
    testA_name = 'signrank';
    Tobs_text  = sprintf('median=%.3f Hz', median(delta_all,'omitnan'));
    df_text    = 'N/A';
    t_text     = 'N/A';
else
    [~,pA,~,statsA] = ttest(delta_all, 0);
    testA_name = 'one-sample t';
    Tobs_text  = sprintf('mean=%.3f Hz, t(%d)=%.2f', ...
        mean(delta_all,'omitnan'), statsA.df, statsA.tstat);
    df_text    = sprintf('%d', statsA.df);
    t_text     = sprintf('%.3f', statsA.tstat);
end

% ---- B) Permutation (sign-flip) test on mean ----
T_obs  = mean(delta_all,'omitnan');
T_perm = zeros(nPerm,1);
for b=1:nPerm
    flips     = (rand(size(delta_all))<0.5)*2 - 1;   % ±1
    T_perm(b) = mean(delta_all .* flips,'omitnan');
end
p_perm = mean(abs(T_perm) >= abs(T_obs));

% ---- POPULATION PRINTOUT ----
fprintf('\n=== Population \\Delta-rate summary (%s) ===\n', titleTag);
fprintf('nCells_total = %d\n', numel(delta_all));
fprintf('One-sample test on per-cell \\Delta (trial−nontrial): %s\n', testA_name);
fprintf('  mean(Delta) = %.4f Hz, median = %.4f Hz\n', ...
    mean(delta_all,'omitnan'), median(delta_all,'omitnan'));
fprintf('  t/stat = %s, df = %s, p = %.4g\n', t_text, df_text, pA);
fprintf('Permutation test on mean(Delta):\n');
fprintf('  T_obs = %.4f Hz, p_perm = %.4g   (nPerm = %d)\n\n', T_obs, p_perm, nPerm);

% ---- PLOTS ----
figure('Color','w','Position',[180 220 1200 480]);

% (Left) All cells violin + swarm; per-rat means with CI on the side
subplot(1,2,1); hold on;
drawViolin(1, delta_all);
swarmchart(1*ones(size(delta_all)), delta_all, 6, 'filled', ...
           'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.15);
% mean ± SE for all cells
muAll = mean(delta_all,'omitnan');
seAll = std(delta_all,'omitnan')/sqrt(numel(delta_all));
errorbar(1.15, muAll, seAll, 'k.','LineWidth',1.2,'CapSize',12);
% per-rat means
for rr=1:numel(R)
    plot(1.35 + 0.05*randn, perRat_means(rr), 'd', 'MarkerSize',6, ...
        'MarkerFaceColor',[.2 .6 .8], 'MarkerEdgeColor','k');
end
yline(0,'k--');
xlim([0.7 1.6]); xticks([1]); xticklabels({'All cells'});
ylabel('\Delta rate (trial − nontrial), Hz');
title(sprintf('Population effect %s\n%s, p=%.3g', titleTag, Tobs_text, pA));

% (Right) Permutation null of mean deltaRate
subplot(1,2,2); hold on;
histogram(T_perm, 'Normalization','percentage','EdgeColor','none');
xline(T_obs, 'r-', 'LineWidth',2);
yL = ylim;
text(T_obs, 0.9*yL(2), sprintf('obs=%.3f Hz\np_{perm}=%.4g', T_obs, p_perm), ...
    'Color','r','FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
xlabel('Mean \Delta rate under sign-flip null'); ylabel('Density');
title('Permutation test (sign flip across cells)');
box on;

% ---- return struct ----
Pop = struct();
Pop.early_exit     = false;
Pop.delta_all      = delta_all;
Pop.perRat_means   = perRat_means;
Pop.perRat_ns      = perRat_ns;
Pop.oneSample      = struct('test',testA_name,'p',pA,'stats',statsA,'text',Tobs_text);
Pop.permutation    = struct('nPerm',nPerm,'T_obs',T_obs,'T_perm',T_perm,'p',p_perm);
end

function Pop = accel_population_summary_split(R, varargin)
% Left: 6 violins (5 rats + pooled "All") of per-cell Δrate
% Right: 3x2 permutation tests (per rat + pooled)
% EARLY EXIT: if any rat’s permutation p exceeds alphaStop, return immediately.

p = inputParser;
p.addParameter('nPerm', 500);
p.addParameter('useSignrank', false);
p.addParameter('titleTag', '');
p.addParameter('alphaStop', 0.05);   % early-exit threshold
p.parse(varargin{:});
nPerm       = p.Results.nPerm;
useSignrank = p.Results.useSignrank; %#ok<NASGU>
titleTag    = p.Results.titleTag;
alphaStop   = p.Results.alphaStop;

nR       = numel(R);
ratNames = {R.rat};

% ---- gather per-rat vectors (deltas + rates) ----
perRat_cellDeltas  = cell(nR,1);
perRat_cellTrialHz = cell(nR,1);
perRat_cellNonHz   = cell(nR,1);
perRat_meanDelta   = nan(nR,1);
perRat_meanTrialHz = nan(nR,1);
perRat_meanNonHz   = nan(nR,1);
perRat_nCells      = zeros(nR,1);

for rr=1:nR
    dDel = []; dTri = []; dNon = [];
    for d=1:numel(R(rr).perDay)
        S = R(rr).perDay{d};
        if isempty(S) || ~isfield(S,'pop') || isempty(S.pop)
            continue;
        end

        % convert population-summed rates to per-cell rates
        nC = 1;
        if isfield(S.pop,'nCells') && ~isempty(S.pop.nCells) && S.pop.nCells > 0
            nC = S.pop.nCells;
        end

        dDel = [dDel; S.pop.deltaPerBin(:)      / nC]; %#ok<AGROW>
        dTri = [dTri; S.pop.trialRatePerBin(:)  / nC]; %#ok<AGROW>
        dNon = [dNon; S.pop.nonRatePerBin(:)    / nC]; %#ok<AGROW>
    end
    dDel = dDel(isfinite(dDel));
    dTri = dTri(isfinite(dTri));
    dNon = dNon(isfinite(dNon));

    perRat_binDeltas{rr}   = dDel;
    perRat_binTrialHz{rr}  = dTri;
    perRat_binNonHz{rr}    = dNon;

    perRat_nBins(rr)       = numel(dDel);
    perRat_meanDelta(rr)   = mean(dDel,'omitnan');   % Hz per cell
    perRat_meanTrialHz(rr) = mean(dTri,'omitnan');   % Hz per cell
    perRat_meanNonHz(rr)   = mean(dNon,'omitnan');   % Hz per cell
end

delta_all = vertcat(perRat_cellDeltas{:});
trial_all = vertcat(perRat_cellTrialHz{:});
non_all   = vertcat(perRat_cellNonHz{:});

% ---- EARLY-EXIT permutation loop (per rat) ----
T_obs_rats  = nan(nR,1);
p_perm_rats = nan(nR,1);
T_perm_rats = cell(nR,1);

for rr=1:nR
    x = perRat_cellDeltas{rr};
    x = x(isfinite(x));
    if isempty(x)
        continue;
    end
    T_obs_rats(rr) = mean(x,'omitnan');

    % permutation null for THIS rat
    Tperm = zeros(nPerm,1);
    for b=1:nPerm
        flips = (rand(size(x))<0.5)*2 - 1;
        Tperm(b) = mean(x .* flips,'omitnan');
    end
    p_perm_rats(rr) = mean(abs(Tperm) >= abs(T_obs_rats(rr)));
    T_perm_rats{rr} = Tperm;

    % short-circuit
    if p_perm_rats(rr) > alphaStop
        Pop = struct();
        Pop.early_exit   = true;
        Pop.reason       = sprintf('Rat %s exceeded alphaStop (p_{perm}=%.3g > %.3g).', ...
                                   ratNames{rr}, p_perm_rats(rr), alphaStop);
        Pop.offending_rat= ratNames{rr};
        Pop.p_perm_rat   = p_perm_rats(rr);
        Pop.nPerm        = nPerm;
        Pop.titleTag     = titleTag;
        Pop.partial      = struct('T_obs_rats', T_obs_rats, ...
                                  'p_perm_rats', p_perm_rats, ...
                                  'T_perm_rats', {T_perm_rats});
        return;
    end
end

% ---- pooled permutation
T_obs_all  = mean(delta_all,'omitnan');
T_perm_all = zeros(nPerm,1);
for b=1:nPerm
    flips = (rand(size(delta_all))<0.5)*2 - 1;
    T_perm_all(b) = mean(delta_all .* flips,'omitnan');
end
p_perm_all = mean(abs(T_perm_all) >= abs(T_obs_all));

% ---- pooled one-sample test on delta_all (defines p_all/st_all/oneText) ----
if useSignrank
    [p_all,~,st_all] = signrank(delta_all, 0);
    oneText = sprintf('median=%.3f Hz', median(delta_all,'omitnan'));
    df_text = 'N/A';
    t_text  = 'N/A';
else
    [~,p_all,~,st_all] = ttest(delta_all, 0);
    oneText = sprintf('mean=%.3f Hz, t(%d)=%.2f', ...
        mean(delta_all,'omitnan'), st_all.df, st_all.tstat);
    df_text = sprintf('%d', st_all.df);
    t_text  = sprintf('%.3f', st_all.tstat);
end

% ---- POPULATION-LEVEL PRINT TABLE (per rat + pooled) ----
fprintf('\n=== Mean firing rates (per rat; tested cells) — %s ===\n', titleTag);
fprintf('%-8s  nCells  Trial_Hz  Nontrial_Hz   Delta_Hz   p_perm\n', 'Rat');
for rr=1:nR
    fprintf('%-8s  %6d  %8.4f  %10.4f  %8.4f   %8.4g\n', ...
        ratNames{rr}, perRat_nCells(rr), ...
        perRat_meanTrialHz(rr), perRat_meanNonHz(rr), ...
        perRat_meanDelta(rr), p_perm_rats(rr));
end
fprintf('Pooled:  nCells=%d  Trial=%.4f Hz  Nontrial=%.4f Hz  Δ=%.4f Hz  p_perm=%.4g  p_one=%.4g  (stat=%s, df=%s)\n\n', ...
  numel(delta_all), ...
  mean(trial_all,'omitnan'), mean(non_all,'omitnan'), mean(delta_all,'omitnan'), ...
  p_perm_all, p_all, t_text, df_text);

% ---- FIGURE: top-level tiling ----
fig = figure('Color','w','Position',[160 200 1450 560]);
TL  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

% LEFT tile: violins (per rat + pooled)
axL = nexttile(TL,1); hold(axL,'on');
for rr=1:nR
    d = perRat_cellDeltas{rr};
    if isempty(d), continue; end
    drawViolin(rr, d);
    swarmchart(rr*ones(size(d)), d, 6, 'filled', ...
        'MarkerFaceAlpha',0.20, 'MarkerEdgeAlpha',0.08);
    mu = mean(d,'omitnan'); se = std(d,'omitnan')/sqrt(numel(d));
    errorbar(rr+0.12, mu, se, 'k.','LineWidth',1.2,'CapSize',12);
end
% pooled at nR+1
drawViolin(nR+1, delta_all);
swarmchart((nR+1)*ones(size(delta_all)), delta_all, 6, 'filled', ...
    'MarkerFaceAlpha',0.20, 'MarkerEdgeAlpha',0.08);
muAll = mean(delta_all,'omitnan'); seAll = std(delta_all,'omitnan')/sqrt(numel(delta_all));
errorbar(nR+1+0.12, muAll, seAll, 'k.','LineWidth',1.2,'CapSize',12);

yline(0,'k--');
xticks(1:nR+1); xticklabels([ratNames, {'All'}]); xtickangle(15);
ylabel('\Delta rate (trial − nontrial), Hz');

% RIGHT tile: nested 3x2 permutation grid
subTL = tiledlayout(TL,3,2,'TileSpacing','compact','Padding','compact');
subTL.Layout.Tile = 2;  % put nested layout in right tile

for rr=1:nR
    ax = nexttile(subTL); hold(ax,'on');
    Tperm = T_perm_rats{rr};
    if isempty(Tperm)
        title(ax, sprintf('%s (no cells)', ratNames{rr}));
        axis(ax,'off');
        continue;
    end
    histogram(ax, Tperm, 'Normalization','percentage','EdgeColor','none');
    xline(ax, T_obs_rats(rr), 'r-', 'LineWidth',2);
    title(ax, sprintf('%s  p_{perm}=%.3g', ratNames{rr}, p_perm_rats(rr)));
    xlabel(ax,'mean \Delta rate'); ylabel(ax,'pdf');
    box(ax,'on');
end
ax = nexttile(subTL); hold(ax,'on');
histogram(ax, T_perm_all, 'Normalization','percentage','EdgeColor','none');
xline(ax, T_obs_all, 'r-', 'LineWidth',2);
title(ax, sprintf('All  p_{perm}=%.3g', p_perm_all));
xlabel(ax,'mean \Delta rate'); ylabel(ax,'pdf');
box(ax,'on');
title(subTL, 'Permutation tests (sign flip on cell-level effects)');

% ---- pack outputs ----
Pop = struct();
Pop.early_exit = false;

Pop.perRat = struct( ...
    'name',           ratNames(:), ...
    'nCells',         num2cell(perRat_nCells(:)), ...
    'cellDelta',      perRat_cellDeltas(:), ...
    'cellTrialHz',    perRat_cellTrialHz(:), ...
    'cellNontrialHz', perRat_cellNonHz(:), ...
    'meanDelta',      num2cell(perRat_meanDelta(:)), ...
    'meanTrialHz',    num2cell(perRat_meanTrialHz(:)), ...
    'meanNontrialHz', num2cell(perRat_meanNonHz(:)), ...
    'perm_Tobs',      num2cell(T_obs_rats(:)), ...
    'perm_p',         num2cell(p_perm_rats(:)) );

Pop.pooled = struct( ...
    'cellDelta',      delta_all, ...
    'cellTrialHz',    trial_all, ...
    'cellNontrialHz', non_all, ...
    'meanDelta',      mean(delta_all,'omitnan'), ...
    'meanTrialHz',    mean(trial_all,'omitnan'), ...
    'meanNontrialHz', mean(non_all,'omitnan'), ...
    'perm_Tobs',      T_obs_all, ...
    'perm_p',         p_perm_all );

if useSignrank
    Pop.oneSample = struct('test','signrank','p',p_all,'stats',st_all,'text',oneText);
else
    Pop.oneSample = struct('test','one-sample t','p',p_all,'stats',st_all,'text',oneText);
end
Pop.figure = fig;
end

function speedBin_cellwise_sigpanels(R, varargin)
% Draws:
%  • Figure 1: 3x2 grid — one subplot per animal (3 bars = day−2/−1/0),
%    bars show % significant cells; labels show "sig/tested".
%  • Figure 2: summary — per-rat mean % sig (± SD across days) + pooled line.

p = inputParser;
p.addParameter('titleTag','');
p.parse(varargin{:});
titleTag = p.Results.titleTag;

nR       = numel(R);
ratNames = {R.rat};

% ------- collect per-day tested/sig for each rat -------
tested_counts = nan(nR,3);
sig_counts    = nan(nR,3);
pct_sig       = nan(nR,3);

for rr = 1:nR
    for d = 1:3
        S = R(rr).perDay{d};
        if isempty(S), continue; end
        tested = isfinite(S.pVal);
        nTest  = sum(tested);
        nSig   = sum(S.sigDiffFDR(tested));
        tested_counts(rr,d) = nTest;
        sig_counts(rr,d)    = nSig;
        if nTest>0
            pct_sig(rr,d) = 100 * nSig / nTest;
        end
    end
end

% ------- Figure 1: per-animal 3-bar panels -------
figure('Color','w','Position',[160 160 1200 720]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
dayLabels = {'day−2','day−1','day 0'};

for rr = 1:nR
    nexttile; hold on;
    bar(pct_sig(rr,:),'FaceColor',[0.70 0.80 0.95]);
    xticks(1:3); xticklabels(dayLabels); ylim([0 100]);
    ylabel('% sig cells'); title(ratNames{rr});
    % annotate counts "sig/tested" above bars
    for d=1:3
        if ~isnan(pct_sig(rr,d))
            txt = sprintf('%d/%d', sig_counts(rr,d), tested_counts(rr,d));
            y   = pct_sig(rr,d);
            text(d, y + 2, txt, 'HorizontalAlignment','center','FontSize',9);
        else
            text(d, 5, 'n/a', 'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);
        end
    end
    box on; yline(0,'k-');
end
% optional empty tiles
if nR < 6
    for k = (nR+1):6, nexttile; axis off; end
end
sgtitle(sprintf('Per-cell significance per day — %s', titleTag));

% ------- Figure 2: summary per rat + pooled -------
mean_perRat = mean(pct_sig,2,'omitnan');
sd_perRat   = std(pct_sig,0,2,'omitnan');

% pooled across all rats/days
pooled_num = sum(sig_counts(:), 'omitnan');
pooled_den = sum(tested_counts(:), 'omitnan');
pooled_pct = 100 * pooled_num / max(pooled_den,1);

figure('Color','w','Position',[220 220 960 420]); hold on;
bar(mean_perRat,'FaceColor',[0.35 0.65 0.85]);
errorbar(1:nR, mean_perRat, sd_perRat, 'k.','LineWidth',1.3,'CapSize',12);
xline(nR+0.6,'k:');
plot([nR+1 nR+1], [pooled_pct pooled_pct], 'r.', 'MarkerSize',18);
yline(pooled_pct,'r--','LineWidth',1.5);
legend({'Per-rat mean % sig','SD across days','Pooled % sig'},'Location','best');
xticks(1:nR); xticklabels(ratNames); xtickangle(15);
ylabel('Mean % sig (across days)'); ylim([0 100]);
title(sprintf('Per-cell significance (%% sig) — pooled=%.1f%%', pooled_pct));
box on;
end

function PopTN = speedBin_population_taskNon_violin(R, varargin)
% Split violins using MATLAB's built-in violinplot (R2024b).
% Left = Non-trial ('negative'), Right = Trial ('positive').

p = inputParser;
p.addParameter('titleTag','',@(s)ischar(s)||isstring(s));
p.addParameter('ShowData', true, @islogical);           % our own option (uses scatter)
p.addParameter('DataAlpha', 0.18, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('Bandwidth', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0)); % accepted but IGNORED
p.parse(varargin{:});
titleTag  = string(p.Results.titleTag);
showData  = p.Results.ShowData;
dataAlpha = p.Results.DataAlpha;

ratNames = {R.rat};
nRats    = numel(R);

% ----- collect per-rat rates (tested cells only) -----
allNon = cell(nRats,1);
allTri = cell(nRats,1);
for rr = 1:nRats
    n = []; t = [];
    for d = 1:numel(R(rr).perDay)
        S = R(rr).perDay{d};
        if isempty(S), continue; end
        tested = isfinite(S.pVal);
        n = [n; double(S.meanNonTrialRate(tested))];
        t = [t; double(S.meanTrialRate(tested))];
    end
    allNon{rr} = n(isfinite(n));
    allTri{rr} = t(isfinite(t));
end
non_all = vertcat(allNon{:});
tri_all = vertcat(allTri{:});

% ----- long-form inputs (groups first, numeric second) -----
catOrder = [ratNames, {'All'}];

% Non-trial (left halves -> 'negative')
Y_non = zeros(0,1,'double'); G_non = {};
for rr = 1:nRats
    y = allNon{rr}; if isempty(y), continue; end
    Y_non = [Y_non; max(y(:),0)];                          % clamp tiny negatives
    G_non = [G_non; repmat(ratNames(rr), numel(y), 1)];
end
if ~isempty(non_all)
    Y_non = [Y_non; max(non_all(:),0)];
    G_non = [G_non; repmat({'All'}, numel(non_all), 1)];
end
G_non = categorical(G_non, catOrder);
maskN = isfinite(Y_non) & ~isundefined(G_non);
Y_non = Y_non(maskN);  G_non = G_non(maskN);

% Trial (right halves -> 'positive')
Y_tri = zeros(0,1,'double'); G_tri = {};
for rr = 1:nRats
    y = allTri{rr}; if isempty(y), continue; end
    Y_tri = [Y_tri; max(y(:),0)];
    G_tri = [G_tri; repmat(ratNames(rr), numel(y), 1)];
end
if ~isempty(tri_all)
    Y_tri = [Y_tri; max(tri_all(:),0)];
    G_tri = [G_tri; repmat({'All'}, numel(tri_all), 1)];
end
G_tri = categorical(G_tri, catOrder);
maskT = isfinite(Y_tri) & ~isundefined(G_tri);
Y_tri = Y_tri(maskT);  G_tri = G_tri(maskT);

% ----- plot (STRICT doc-compliant) -----
figure('Color','w','Position',[180 260 1200 480]); hold on;

vL = violinplot(G_non, Y_non, 'DensityDirection','negative');  % left half (Non)
vR = violinplot(G_tri, Y_tri, 'DensityDirection','positive');  % right half (Trial)

% style
for k = 1:numel(vL)
    vL(k).FaceColor = [0.50 0.70 0.95];
    vL(k).FaceAlpha = 0.65;
    vL(k).EdgeColor = [0.25 0.35 0.55];
end
for k = 1:numel(vR)
    vR(k).FaceColor = [0.95 0.60 0.50];
    vR(k).FaceAlpha = 0.65;
    vR(k).EdgeColor = [0.55 0.35 0.25];
end

% optional rugs (verify density vs raw points)
if showData
    xN = grp2idx(G_non);
    xT = grp2idx(G_tri);
    scatter(xN-0.02, Y_non, 6, [0 0 1], 'filled', ...
        'MarkerFaceAlpha', dataAlpha, 'MarkerEdgeAlpha', dataAlpha, 'HandleVisibility','off');
    scatter(xT+0.02, Y_tri, 6, [1 0 0], 'filled', ...
        'MarkerFaceAlpha', dataAlpha, 'MarkerEdgeAlpha', dataAlpha, 'HandleVisibility','off');
end

ylabel('Rate (Hz)');
title(sprintf('Trial vs Non-trial rates %s', titleTag));
yline(0,'k--'); box on;

% simple legend
hN = plot(nan,nan,'-','Color',[0.50 0.70 0.95],'LineWidth',8);
hT = plot(nan,nan,'-','Color',[0.95 0.60 0.50],'LineWidth',8);
legend([hN hT], {'Non-trial','Trial'}, 'Location','best');

% ---- return ----
PopTN = struct();
PopTN.non_all   = non_all;
PopTN.trial_all = tri_all;
PopTN.perRat    = struct('non', {allNon}, 'trial', {allTri});
end









function speedBin_sigbar_acrossAnimals(R, varargin)
% One figure: bars = per-animal mean % significant across days,
% plus one final bar = mean of those per-animal means ("Across animals").

p = inputParser;
p.addParameter('titleTag','');
p.parse(varargin{:});
titleTag = p.Results.titleTag;

nR = numel(R);
ratNames = {R.rat};

pct_sig = nan(nR,3);
tested_counts = nan(nR,3);
sig_counts    = nan(nR,3);

% collect per-day %sig
for rr=1:nR
    for d=1:3
        S = R(rr).perDay{d};
        if isempty(S), continue; end
        tested = isfinite(S.pVal);
        nTest  = sum(tested);
        nSig   = sum(S.sigDiffFDR(tested));
        tested_counts(rr,d) = nTest;
        sig_counts(rr,d)    = nSig;
        if nTest>0
            pct_sig(rr,d) = 100*nSig/nTest;
        end
    end
end

% per-animal mean across days
meanPct_perRat = mean(pct_sig,2,'omitnan');

% across-animals bar = mean of the per-animal means
acrossAnimals  = mean(meanPct_perRat,'omitnan');

% plot
figure('Color','w','Position',[220 220 900 420]); hold on;
x = 1:nR;
bar(x, meanPct_perRat, 'FaceColor',[0.35 0.65 0.85]);
% add last bar "Across"
bar(nR+1, acrossAnimals, 'FaceColor',[0.85 0.45 0.45]);
yline(0,'k-'); box on;
xticks([x nR+1]); xticklabels([ratNames, {'Across'}]); xtickangle(15);
ylabel('Mean % sig (across days)');
title(sprintf('Per-cell significance (%% sig) — %s', titleTag));
ylim([0 max(100, max([meanPct_perRat; acrossAnimals])+10)]);

% annotate each bar with "Σsig/Σtested" across days (so you see raw counts)
for rr=1:nR
    nSig_sum = nansum(sig_counts(rr,:));
    nTest_sum= nansum(tested_counts(rr,:));
    if nTest_sum>0
        text(rr, meanPct_perRat(rr)+3, sprintf('%d/%d', nSig_sum, nTest_sum), ...
            'HorizontalAlignment','center','FontSize',9);
    end
end
% across-animals counts pooled across all rats/days
nSig_all  = nansum(sig_counts(:));
nTest_all = nansum(tested_counts(:));
text(nR+1, acrossAnimals+3, sprintf('%d/%d', nSig_all, nTest_all), ...
    'HorizontalAlignment','center','FontSize',9);

% console recap

%fprintf('\n=== Per-animal mean %% sig (across days) ===\n');
%for rr=1:nR
%    fprintf('%s: %.1f%%   (Σ %d/%d)\n', ratNames{rr}, ...
%        meanPct_perRat(rr), nansum(sig_counts(rr,:)), nansum(tested_counts(rr,:)));
%end

%fprintf('Across animals (mean of per-animal means): %.1f%%   (pooled Σ %d/%d)\n\n', ... %%%
%    acrossAnimals, nSig_all, nTest_all);
%
end


function speedBin_population_combo_plots(R, varargin)
% Three subplots:
% (1) Δ-rate violin across all cells (your current panel A)
% (2) NEW: per-rat paired violins (Non vs Trial) + an "All" pair; stars per rat
% (3) Permutation histogram of mean Δ-rate (your current panel C)

p = inputParser;
p.addParameter('nPerm', 5000);
p.addParameter('useSignrank', false);   % paired t-test vs signrank for the violins
p.addParameter('titleTag', '');
p.parse(varargin{:});
nPerm       = p.Results.nPerm;
useSignrank = p.Results.useSignrank;
titleTag    = p.Results.titleTag;

nR = numel(R);
ratNames = {R.rat};

% ---------- gather cellwise deltas and rates ----------
delta_all = [];
perRat_non = cell(nR,1);
perRat_tri = cell(nR,1);
for rr=1:nR
    non_r = []; tri_r = [];
    for d=1:3
        S = R(rr).perDay{d}; if isempty(S), continue; end
        tested = isfinite(S.pVal);
        delta_all = [delta_all; S.deltaRate(tested)]; %#ok<AGROW>
        non_r     = [non_r;     S.meanNonTrialRate(tested)]; %#ok<AGROW>
        tri_r     = [tri_r;     S.meanTrialRate(tested)];    %#ok<AGROW>
    end
    perRat_non{rr} = non_r;
    perRat_tri{rr} = tri_r;
end
non_all = vertcat(perRat_non{:});
tri_all = vertcat(perRat_tri{:});

% ---------- basic tests for panel (1) and (3) ----------
% one-sample t on delta_all and permutation test (same as before)
[~,pA,~,stA] = ttest(delta_all,0);
T_obs  = mean(delta_all,'omitnan');
T_perm = zeros(nPerm,1);
for b=1:nPerm
    flips     = (rand(size(delta_all))<0.5)*2 - 1;
    T_perm(b) = mean(delta_all .* flips,'omitnan');
end
p_perm = mean(abs(T_perm) >= abs(T_obs));

% ---------- figure with 3 subplots ----------
figure('Color','w','Position',[160 160 1500 500]);

% (1) Δ-rate violin across all cells
subplot(1,3,1); hold on;
drawViolin(1, delta_all);
swarmchart(1*ones(size(delta_all)), delta_all, 6, 'filled', ...
           'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.15);
muAll = mean(delta_all,'omitnan'); seAll = std(delta_all,'omitnan')/sqrt(numel(delta_all));
errorbar(1.15, muAll, seAll, 'k.','LineWidth',1.2,'CapSize',12);
yline(0,'k--'); xlim([0.7 1.6]); xticks(1); xticklabels({'All cells'});
ylabel('\Delta rate (trial − nontrial), Hz');
title(sprintf('\\Delta-rate summary %s\nmean=%.3f Hz, t(%d)=%.2f, p=%.3g', ...
      titleTag, muAll, stA.df, stA.tstat, pA));

% (2) Paired violins per rat (+ All)
subplot(1,3,2); hold on;

% layout: for each rat use two x's [g, g+0.35], then gap 0.3
xPairs = [];
labels = {};
x0 = 1;
for rr=1:nR
    non = perRat_non{rr}; tri = perRat_tri{rr};
    if isempty(non) || isempty(tri), continue; end
    % violins
    drawViolin(x0,     non);
    drawViolin(x0+0.35, tri);
    % swarms
    swarmchart(x0*ones(size(non)), non, 5, 'filled', 'MarkerFaceAlpha',0.18,'MarkerEdgeAlpha',0.05);
    swarmchart((x0+0.35)*ones(size(tri)), tri, 5, 'filled', 'MarkerFaceAlpha',0.18,'MarkerEdgeAlpha',0.05);

    % paired test
    minN = min(numel(non), numel(tri));
    if useSignrank
        [pval,~,~] = signrank(tri(1:minN), non(1:minN));
    else
        [~,pval] = ttest(tri(1:minN), non(1:minN));
    end
    stars = p2stars(pval);
    % star bar
    yMax = max([non; tri],[],'omitnan');
    yMin = min([non; tri],[],'omitnan');
    yPad = 0.06 * max(1e-6, yMax - yMin);
    yBar = yMax + 2*yPad;
    plot([x0, x0+0.35], [yBar yBar], 'k-', 'LineWidth',1);
    text(x0+0.175, yBar + yPad, stars, 'HorizontalAlignment','center','FontWeight','bold');

    xPairs = [xPairs, x0, x0+0.35]; %#ok<AGROW>
    labels{end+1} = sprintf('%s\nNon', ratNames{rr}); %#ok<AGROW>
    labels{end+1} = sprintf('%s\nTrial', ratNames{rr}); %#ok<AGROW>
    x0 = x0 + 1.1; % advance group
end

% add "All" pair at the end
drawViolin(x0,     non_all);
drawViolin(x0+0.35, tri_all);
swarmchart(x0*ones(size(non_all)), non_all, 5, 'filled', 'MarkerFaceAlpha',0.18,'MarkerEdgeAlpha',0.05);
swarmchart((x0+0.35)*ones(size(tri_all)), tri_all, 5, 'filled', 'MarkerFaceAlpha',0.18,'MarkerEdgeAlpha',0.05);
% paired test for All
minN = min(numel(non_all), numel(tri_all));
if useSignrank
    [p_all,~,~] = signrank(tri_all(1:minN), non_all(1:minN));
else
    [~,p_all] = ttest(tri_all(1:minN), non_all(1:minN));
end
stars_all = p2stars(p_all);
yMax = max([non_all; tri_all],[],'omitnan');
yMin = min([non_all; tri_all],[],'omitnan');
yPad = 0.06 * max(1e-6, yMax - yMin);
yBar = yMax + 2*yPad;
plot([x0, x0+0.35], [yBar yBar], 'k-', 'LineWidth',1);
text(x0+0.175, yBar + yPad, stars_all, 'HorizontalAlignment','center','FontWeight','bold');

% axis cosmetics
xticks([1:1.1:(x0-1.1), x0, x0+0.35]);
xticklabels([labels, {'All','Trial'}]);  % the last two labels will look stacked; OK
xtickangle(0);
ylabel('Rate (Hz)');
title('Per-rat paired violins (Non vs Trial) with significance stars');
box on;

% (3) Permutation histogram of mean Δ-rate
subplot(1,3,3); hold on;
histogram(T_perm, 'Normalization','percentage','EdgeColor','none');
xline(T_obs, 'r-', 'LineWidth',2);
yL = ylim;
text(T_obs, 0.9*yL(2), sprintf('obs=%.3f Hz\np_{perm}=%.4g', T_obs, p_perm), ...
    'Color','r','FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top');
xlabel('Mean \Delta rate under sign-flip null'); ylabel('Density');
title('Permutation test (sign flip across cells)');
box on;

% console sanity check of raw means
%%%
%{
fprintf('\n=== Raw rate means (tested cells) ===\n');
for rr=1:nR
    non = perRat_non{rr}; tri = perRat_tri{rr};
    if isempty(non) || isempty(tri), continue; end
    fprintf('%s: non=%.6f Hz (n=%d), trial=%.6f Hz (n=%d), Δ=%.6f\n', ...
        ratNames{rr}, mean(non,'omitnan'), numel(non), mean(tri,'omitnan'), numel(tri), ...
        mean(tri-non,'omitnan'));
end
fprintf('All : non=%.6f Hz (n=%d), trial=%.6f Hz (n=%d), Δ=%.6f\n\n', ...
    mean(non_all,'omitnan'), numel(non_all), mean(tri_all,'omitnan'), numel(tri_all), ...
    mean(tri_all - non_all,'omitnan'));
%}
end


function s = p2stars(p)
if     p < 1e-3, s='***';
elseif p < 1e-2, s='**';
elseif p < 5e-2, s='*';
else,            s='n.s.';
end
end

function drawSplitViolin_safe(x0, leftData, rightData, halfWidth, varargin)
% Split violin (left=Non-trial, right=Trial) with robust sanity checks.
% Always plots y = rate (Hz), width ∝ KDE density.
if nargin<4 || isempty(halfWidth), halfWidth = 0.22; end
p = inputParser;
p.addParameter('Bandwidth',[]);        % empty -> Scott's rule
p.addParameter('GridN',512);
p.addParameter('CapQuantile',99);    % y-limit based on pooled percentile
p.addParameter('ShowTicks',true,@islogical);
p.parse(varargin{:});
BW   = p.Results.Bandwidth;
N    = p.Results.GridN;
capQ = p.Results.CapQuantile;
show = p.Results.ShowTicks;

L = leftData(:);  L = L(isfinite(L));
R = rightData(:); R = R(isfinite(R));
if isempty(L) && isempty(R), return; end
Xall = [L; R];

% ---- sanity prints ----
local_print('NON',L);
local_print('TRI',R);

% ---- KDE grid (reflect at 0) ----
yTop = max(prctile(Xall, capQ));
yTop = max(yTop, max(Xall));                   % never below max
yTop = yTop * 1.02;                            % tiny headroom
y    = linspace(0, yTop, N);

[fL, y] = local_kde_reflect(L, y, BW);         % ALWAYS returns (f,y)
[fR, ~] = local_kde_reflect(R, y, BW);

% add zero-width endpoints for smooth taper
y  = y(:);
fL = [0; fL(:); 0];
fR = [0; fR(:); 0];
y  = [0;  y;    yTop];

% normalize widths with a shared max so halves are comparable
mx = max([fL; fR; eps]);
WL = (fL / mx) * halfWidth;
WR = (fR / mx) * halfWidth;

% left half
patch([x0 - WL; repmat(x0,size(y))], [y; flipud(y)], [0.50 0.70 0.95], ...
      'EdgeColor','none', 'FaceAlpha',0.65);

% right half
patch([repmat(x0,size(y)); x0 + WR], [y; flipud(y)], [0.95 0.60 0.50], ...
      'EdgeColor','none', 'FaceAlpha',0.65);


% optional: mean/median ticks (so the picture matches your numbers)
if show
    muL = mean(L,'omitnan'); mdL = median(L,'omitnan');
    muR = mean(R,'omitnan'); mdR = median(R,'omitnan');
    plot([x0-0.13 x0-0.02],[muL muL],'b-','LineWidth',2);
    plot([x0-0.13 x0-0.02],[mdL mdL],'b:','LineWidth',1);
    plot([x0+0.02 x0+0.13],[muR muR],'r-','LineWidth',2);
    plot([x0+0.02 x0+0.13],[mdR mdR],'r:','LineWidth',1);
end

% keep the axes honest
ylim([0, yTop]);  % rate (Hz) on vertical axis
end

% ---- KDE that ALWAYS returns (f,y), reflected at 0 ----
function [f,y] = local_kde_reflect(x, ygrid, BW)
x = x(:); x = x(isfinite(x));
x(x<0) = 0;                      % <- clamp tiny negatives
y = ygrid(:)';

if isempty(x), f = zeros(size(y)); return; end
if isempty(BW) || ~isfinite(BW) || BW<=0
    s = std(x); if ~isfinite(s) || s==0, s = max(x)/10; if s==0, s=1; end, end
    BW = 1.06*s*max(1,numel(x))^(-1/5);
end

try
    [f,y] = ksdensity(x, y, 'Bandwidth',BW, ...
        'Support',[0 Inf], 'BoundaryCorrection','reflection');
    return;
catch
    phi = @(z) exp(-0.5*z.^2)/sqrt(2*pi);
    n   = numel(x);
    f   = zeros(size(y));
    for i=1:n
        xi = x(i);
        f = f + phi((y - xi)/BW) + phi((y + xi)/BW);
    end
    f = f / (n*BW);
end
end


function local_print(tag, x)
x = x(:); x = x(isfinite(x));
if isempty(x), fprintf('%s: EMPTY\n', tag); return; end
q = prctile(x,[0 25 50 75 90 95 99 99.9]);
fprintf(['%s  n=%d  mean=%.4f  med=%.4f  min=%.4f  p90=%.4f  p95=%.4f  p99=%.4f  p99.9=%.4f  max=%.4f\n'], ...
        tag, numel(x), mean(x), median(x), q(1), q(5), q(6), q(7), q(8), max(x));
end


function print_violin_sanity(tag, x)
x = x(:); x = x(isfinite(x));
if isempty(x), fprintf('%s: EMPTY\n', tag); return; end
x(x<0) = 0;                            % <- clamp tiny negatives

gridMax = max(x);
gridMax = max(gridMax, eps);
grid    = linspace(0, gridMax*1.10, 512);

% KDE for the mode (robust to tiny negatives now)
try
    [f,y] = ksdensity(x, grid, 'Support',[0 Inf], 'BoundaryCorrection','reflection');
catch
    [f,y] = ksdensity(x, grid);        % fallback if Statistics TB options differ
end
[~,imax] = max(f);

qs = prctile(x,[50 90 95 99 99.9]);
fprintf(['%s  n=%d  mean=%.4f  median=%.4f  mode≈%.4f  ' ...
         'p90=%.4f  p95=%.4f  p99=%.4f  p99.9=%.4f  max=%.4f\n'], ...
        tag, numel(x), mean(x), median(x), y(imax), qs(2), qs(3), qs(4), qs(5), max(x));
end

function plot_per_rat_means(R, titleTag)
names = {R.rat};
mNon = nan(numel(R),1); mTri = nan(numel(R),1);

for rr=1:numel(R)
    non=[]; tri=[];
    for d=1:numel(R(rr).perDay)
        S = R(rr).perDay{d}; if isempty(S), continue; end
        tested = isfinite(S.pVal);
        non = [non; double(S.meanNonTrialRate(tested))];
        tri = [tri; double(S.meanTrialRate(tested))];
    end
    mNon(rr) = mean(non,'omitnan');
    mTri(rr) = mean(tri,'omitnan');
end

figure('Color','w','Position',[300 300 650 420]); hold on
x = (1:numel(R))';
plot([x x]', [mNon mTri]','-','Color',[.7 .7 .7]);                 % connectors
scatter(x-0.05,mNon,60,[0.50 0.70 0.95],'filled','MarkerEdgeColor','k');
scatter(x+0.05,mTri,60,[0.95 0.60 0.50],'filled','MarkerEdgeColor','k');
xlim([0.5 numel(R)+0.5]); xticks(x); xticklabels(names); xtickangle(15)
ylabel('Mean rate per tested cell (Hz)'); yline(0,'k--'); box on
legend({'paired mean per rat','Non-trial','Trial'},'Location','best')
title(sprintf('Non vs Trial (per-rat means) %s',titleTag))
end

function plot_delta_only_fromR(R, titleTag)
% Build pooled, within-cell deltas
delta_all = [];
for rr = 1:numel(R)
    for d = 1:numel(R(rr).perDay)
        S = R(rr).perDay{d}; if isempty(S), continue; end
        tested = isfinite(S.pVal);
        tri = double(S.meanTrialRate(tested));
        non = double(S.meanNonTrialRate(tested));
        delta_all = [delta_all; tri - non]; %#ok<AGROW>
    end
end
delta_all = delta_all(isfinite(delta_all));

% Plot: histogram + KDE + mean/CI
figure('Color','w','Position',[300 300 560 420]); hold on
histogram(delta_all,'Normalization','percentage','FaceColor',[0.75 0.85 1],'EdgeColor','none');
[f,x] = ksdensity(delta_all); plot(x,f,'k-','LineWidth',1.8);
xline(0,'k--');

mu = mean(delta_all,'omitnan');
% delta_all is your vector (can include NaNs)
statfun = @(x) mean(x,'omitnan');              % scalar stat
ci = bootci(500, {statfun, delta_all}, 'type','percentile');   % 95% CI

xline(mu,'r-','LineWidth',2);
plot(ci,[0 0],'r|','MarkerSize',16,'LineWidth',2);

title(sprintf('\\Delta = Trial − Non-trial  (mean=%.3f Hz, 95%% CI [%.3f %.3f])',mu,ci(1),ci(2)));
xlabel('\Delta rate (Hz)'); ylabel('Density'); box on
end

function plot_paired_means_slopegraph(R, titleTag)
if nargin<2, titleTag = ''; end
ratNames = {R.rat};
nR = numel(R);
meanNon = nan(nR,1);
meanTri = nan(nR,1);

for rr = 1:nR
    tri = []; non = [];
    for d = 1:numel(R(rr).perDay)
        S = R(rr).perDay{d}; if isempty(S), continue; end
        tested = isfinite(S.pVal);
        tri = [tri; S.meanTrialRate(tested)];
        non = [non; S.meanNonTrialRate(tested)];
    end
    meanTri(rr) = mean(tri,'omitnan');
    meanNon(rr) = mean(non,'omitnan');
end

figure('Color','w','Position',[300 260 520 420]); hold on
for rr = 1:nR
    plot([1 2], [meanNon(rr) meanTri(rr)], '-', 'Color',[.6 .6 .6], 'LineWidth',1.6);
end
scatter(ones(nR,1), meanNon, 40, [0.35 0.65 0.95], 'filled');
scatter(2*ones(nR,1), meanTri, 40, [0.95 0.45 0.45], 'filled');
xlim([0.8 2.2]); xticks([1 2]); xticklabels({'Non-trial','Trial'});
ymin = min([meanNon; meanTri]); ymax = max([meanNon; meanTri]); pad = 0.1*max(ymax-ymin,eps);
ylim([ymin - pad, ymax + pad]); yline(0,'k--'); box on
ylabel('Mean rate per tested cell (Hz)');
title(sprintf('Non vs Trial (per-rat means) %s', titleTag));
for rr = 1:nR
    text(0.98, meanNon(rr), ratNames{rr}, 'HorizontalAlignment','right','Color',[.4 .4 .4]);
end
end


function plot_rate_change_bars(R, titleTag, varargin)
% Bar chart of change from Non-trial -> Trial per rat (+ pooled "All").
% Each bar = mean of per-cell paired change; error bars = SD (default) or SEM across cells.
%
% metric:
%   'percent' (default): 100*(trial - non)/non  (per cell, then averaged)
%   'fold'              : trial/non             (per cell, then averaged)
%   'delta'             : trial - non           (Hz; per cell, then averaged)
%
% errType:
%   'sd'  (default)  -> standard deviation of per-cell change
%   'sem'            -> sd/sqrt(n)
%
% Console table now includes paired t-test p for Trial vs Non-trial per rat + pooled.

p = inputParser;
p.addParameter('metric','percent',@(s) any(validatestring(s,{'percent','fold','delta'})));
p.addParameter('errType','sem',@(s) any(validatestring(s,{'sd','sem'})));
p.addParameter('barColor',[0.30 0.60 0.85]);
p.addParameter('epsDen',1e-12,@(x) isnumeric(x)&&isscalar(x)&&x>0); % protects %/fold from div-by-zero
p.parse(varargin{:});
metric   = p.Results.metric;
errType  = p.Results.errType;
barColor = p.Results.barColor;
epsDen   = p.Results.epsDen;

ratNames = {R.rat};
nR       = numel(R);

% ---------- collect per-rat per-cell rates ----------
perRat_non = cell(nR,1);
perRat_tri = cell(nR,1);
for rr = 1:nR
    tri = []; non = [];
    for d = 1:numel(R(rr).perDay)
        S = R(rr).perDay{d};
        if isempty(S), continue; end
        tested = isfinite(S.pVal);
        tri = [tri; S.meanTrialRate(tested)];
        non = [non; S.meanNonTrialRate(tested)];
    end
    perRat_tri{rr} = tri(:);
    perRat_non{rr} = non(:);
end
% pooled
allTri = vertcat(perRat_tri{:});
allNon = vertcat(perRat_non{:});

% ---------- per-cell paired change -> mean & error ----------
vals    = nan(nR,1);
errs    = nan(nR,1);
nCells  = zeros(nR,1);
meanNon = nan(nR,1);
meanTri = nan(nR,1);
p_pair  = nan(nR,1);  % paired t-test p for Trial vs Non-trial rates

for rr = 1:nR
    t = perRat_tri{rr}; n = perRat_non{rr};
    if isempty(t) || isempty(n), continue; end
    K = min(numel(t),numel(n));   % paired (tested cells)
    t = t(1:K); n = n(1:K);

    % per-cell change for bars
    switch metric
        case 'percent', d = 100*((t - n) ./ max(n,epsDen));
        case 'fold',    d =  (t ./ max(n,epsDen));
        otherwise,      d =  (t - n); % Hz
    end
    d = d(isfinite(d));
    if isempty(d), continue; end

    vals(rr) = mean(d,'omitnan');
    s        = std(d,'omitnan');
    errs(rr) = strcmp(errType,'sem') * (s/sqrt(numel(d))) + strcmp(errType,'sd') * s;
    nCells(rr) = numel(d);

    meanTri(rr) = mean(t,'omitnan');
    meanNon(rr) = mean(n,'omitnan');

    % paired t-test on Trial vs Non-trial rates for info table
    [~,p_pair(rr)] = ttest(t, n);  % per-rat paired test in Hz
end

% pooled “All”
Kall = min(numel(allTri), numel(allNon));
allTri_use = allTri(1:Kall);
allNon_use = allNon(1:Kall);
switch metric
    case 'percent', dAll = 100*((allTri_use - allNon_use) ./ max(allNon_use,epsDen));
    case 'fold',    dAll =  (allTri_use ./ max(allNon_use,epsDen));
    otherwise,      dAll =  (allTri_use - allNon_use);
end
dAll   = dAll(isfinite(dAll));
valAll = mean(dAll,'omitnan');
errAll = (strcmp(errType,'sem') * (std(dAll,'omitnan')/sqrt(numel(dAll)))) + ...
         (strcmp(errType,'sd')  *  std(dAll,'omitnan'));

% paired t-test on pooled Trial vs Non-trial rates
[~,p_pair_all] = ttest(allTri_use, allNon_use);
meanTri_all = mean(allTri_use,'omitnan');
meanNon_all = mean(allNon_use,'omitnan');

vals_plot = [vals; valAll];
errs_plot = [errs; errAll];
labels    = [ratNames, {'All'}];

% ---------- plot ----------
figure('Color','w','Position',[460 260 780 420]); hold on
b = bar(vals_plot, 'FaceColor', barColor, 'EdgeColor','none'); %#ok<NASGU>
errorbar(1:numel(vals_plot), vals_plot, errs_plot, errs_plot, '.k', 'LineWidth',1.2, 'CapSize',10);
yline(0,'k--','HandleVisibility','off');

xticks(1:numel(vals_plot)); xticklabels(labels); xtickangle(15);
ylabel(ylabel_for(metric));
title(sprintf('Change in firing rate (Non \\rightarrow Trial) — %s  (%s bars)', ...
      titleTag, upper(errType)));
box on

% value labels above/below bars
yl = ylim; pad = 0.04*(yl(2)-yl(1) + eps);
for i = 1:numel(vals_plot)
    txt = pretty_val(vals_plot(i), metric);
    if vals_plot(i) >= 0
        va = 'bottom'; ytxt = vals_plot(i) + errs_plot(i) + pad;
    else
        va = 'top';    ytxt = vals_plot(i) - errs_plot(i) - pad;
    end
    text(i, ytxt, txt, 'HorizontalAlignment','center', 'VerticalAlignment', va, 'FontSize',9);
end

% widen y-lims to fit labels
ymin = min(vals_plot - errs_plot); ymax = max(vals_plot + errs_plot);
rng  = ymax - ymin; if rng<=0, rng = 1; end
ylim([ymin - 0.10*rng, ymax + 0.12*rng]);

% ---------- console printout (matches bars) ----------
unitLabel = struct('percent','%','fold','x','delta','Hz');
fprintf('\n=== Non -> Trial change per rat (%s) with %s bars — values match the bars ===\n', metric, upper(errType));
fprintf('%-8s  nCells  BarMetric(%s)   %s_low   %s_high   Non_Hz   Trial_Hz   Delta_Hz   p_pair\n', ...
        'Rat', unitLabel.(metric), upper(errType), upper(errType));

for rr = 1:nR
    if isnan(vals(rr)), continue; end
    ciLo = vals(rr) - errs(rr);
    ciHi = vals(rr) + errs(rr);
    if isnan(p_pair(rr))
        ptxt = '   n/a ';
    else
        ptxt = sprintf('%.4g', p_pair(rr));
    end
    fprintf('%-8s  %6d  %10s   %8s  %8s   %7.4f   %8.4f   %8.4f   %7s\n', ...
        ratNames{rr}, nCells(rr), ...
        pretty_val(vals(rr), metric), ...
        pretty_val(ciLo, metric), pretty_val(ciHi, metric), ...
        meanNon(rr), meanTri(rr), (meanTri(rr)-meanNon(rr)), ptxt);
end

% pooled "All"
ciLo = valAll - errAll; ciHi = valAll + errAll;
if isnan(p_pair_all)
    ptxt_all = '   n/a ';
else
    ptxt_all = sprintf('%.4g', p_pair_all);
end
fprintf('%-8s  %6d  %10s   %8s  %8s   %7.4f   %8.4f   %8.4f   %7s\n\n', ...
    'All', numel(dAll), ...
    pretty_val(valAll, metric), ...
    pretty_val(ciLo,  metric), pretty_val(ciHi, metric), ...
    meanNon_all, meanTri_all, (meanTri_all - meanNon_all), ptxt_all);
end

% ---------- helpers ----------
function lab = ylabel_for(metric)
switch metric
    case 'percent', lab = 'Change from non-trial (%)';
    case 'fold',    lab = 'Fold change (Trial / Non-trial)';
    otherwise,      lab = '\Delta rate (Hz)   Trial - Non-trial';
end
end

function s = pretty_val(v, metric)
switch metric
    case 'percent', s = sprintf('%.1f%%', v);
    case 'fold',    s = sprintf('%.2f×',  v);
    otherwise,      s = sprintf('%.3f',   v);
end
end

function intervals = maskToIntervals(t, mask, dt)
% Convert a boolean mask on sample times t to [start end] intervals.
mask = mask(:); t = t(:);
dm = diff([false; mask; false]);
i1 = find(dm==1);
i2 = find(dm==-1)-1;
intervals = [t(i1)  t(i2)+dt];   % right edge extended by ~dt
end


function n = countInIntervals(spikes, intervals)
% Count spikes that fall inside any of the [start end] intervals.
if isempty(spikes) || isempty(intervals), n = 0; return; end
spikes = spikes(:);
inside = false(size(spikes));
for k = 1:size(intervals,1)
    inside = inside | (spikes>=intervals(k,1) & spikes<intervals(k,2));
end
n = sum(inside);
end

function drawViolin(x0, data, halfWidth)
% Simple single violin centered at x0 (vertical axis = Hz)
if nargin<3, halfWidth = 0.22; end
data = data(:); data = data(isfinite(data));
if isempty(data), return; end
N = 512;
yTop = max(data); yTop = max(yTop, eps) * 1.02;
y = linspace(0, yTop, N);
[f, y] = local_kde_reflect(data, y, []);   % uses your helper above

y = y(:);
f = [0; f(:); 0];
y = [0; y; yTop];

mx = max(f);  w = (f / max(mx,eps)) * halfWidth;

patch([x0 - w; flipud(x0 + w)], [y; flipud(y)], [0.70 0.80 0.95], ...
      'EdgeColor','none','FaceAlpha',0.70);

% mean ± SE ticks so the picture matches the numbers
mu = mean(data,'omitnan'); se = std(data,'omitnan')/sqrt(numel(data));
plot([x0-0.12 x0+0.12],[mu mu],'k-','LineWidth',1.2);
plot([x0 x0],[mu-se mu+se],'k-','LineWidth',1.2);
end

function Pop = accel_population_summary_split_POP(R, varargin)
% Population-level (across cells) speed-matched summary.
% Uses per-day S.pop.* fields (bin-wise population rates).
%
% Prints:
%   === Mean POPULATION firing rates (per rat; speed-matched bins) ===
% in the same spirit as your per-cell summary, but with nBins instead of nCells.

p = inputParser;
p.addParameter('nPerm', 500);
p.addParameter('titleTag', '');
p.addParameter('alphaStop', .1);
p.parse(varargin{:});
nPerm     = p.Results.nPerm;
titleTag  = p.Results.titleTag;
alphaStop = p.Results.alphaStop;

nR       = numel(R);
ratNames = {R.rat};

perRat_binDeltas  = cell(nR,1);
perRat_binTrialHz = cell(nR,1);
perRat_binNonHz   = cell(nR,1);
perRat_meanDelta   = nan(nR,1);
perRat_meanTrialHz = nan(nR,1);
perRat_meanNonHz   = nan(nR,1);
perRat_nBins       = zeros(nR,1);

for rr=1:nR
    dDel = []; dTri = []; dNon = [];
    for d=1:numel(R(rr).perDay)
        S = R(rr).perDay{d};
        if isempty(S) || ~isfield(S,'pop') || isempty(S.pop)
            continue;
        end

        % POP fields are now Hz per cell (we converted in accelBinMatched)
        dDel = [dDel; S.pop.deltaPerBin(:)];        %#ok<AGROW>
        dTri = [dTri; S.pop.trialRatePerBin(:)];    %#ok<AGROW>
        dNon = [dNon; S.pop.nonRatePerBin(:)];      %#ok<AGROW>
    end
    dDel = dDel(isfinite(dDel));
    dTri = dTri(isfinite(dTri));
    dNon = dNon(isfinite(dNon));

    perRat_binDeltas{rr}  = dDel;   % Hz per cell
    perRat_binTrialHz{rr} = dTri;   % Hz per cell
    perRat_binNonHz{rr}   = dNon;   % Hz per cell

    perRat_nBins(rr)       = numel(dDel);
    perRat_meanDelta(rr)   = mean(dDel,'omitnan');
    perRat_meanTrialHz(rr) = mean(dTri,'omitnan');
    perRat_meanNonHz(rr)   = mean(dNon,'omitnan');
end

delta_all = vertcat(perRat_binDeltas{:});
trial_all = vertcat(perRat_binTrialHz{:});
non_all   = vertcat(perRat_binNonHz{:});

% ---- EARLY-EXIT permutation per rat (bins as samples) ----
T_obs_rats  = nan(nR,1);
p_perm_rats = nan(nR,1);
T_perm_rats = cell(nR,1);

for rr=1:nR
    x = perRat_binDeltas{rr};
    x = x(isfinite(x));
    if isempty(x)
        continue;
    end
    T_obs_rats(rr) = mean(x,'omitnan');

    Tperm = zeros(nPerm,1);
    for b=1:nPerm
        flips   = (rand(size(x))<0.5)*2 - 1;
        Tperm(b)= mean(x .* flips,'omitnan');
    end
    p_perm_rats(rr) = mean(abs(Tperm) >= abs(T_obs_rats(rr)));
    T_perm_rats{rr} = Tperm;

    if p_perm_rats(rr) > alphaStop
        Pop = struct();
        Pop.early_exit   = true;
        Pop.reason       = sprintf('Rat %s exceeded alphaStop (p_{perm}=%.3g > %.3g).\n\n', ...
                                   ratNames{rr}, p_perm_rats(rr), alphaStop);
        Pop.offending_rat= ratNames{rr};
        Pop.p_perm_rat   = p_perm_rats(rr);
        Pop.nPerm        = nPerm;
        Pop.titleTag     = titleTag;
        Pop.partial      = struct('T_obs_rats', T_obs_rats, ...
                                  'p_perm_rats', p_perm_rats, ...
                                  'T_perm_rats', {T_perm_rats});
        return;
    end
end

% ---- pooled permutation on bin-level deltas ----
T_obs_all  = mean(delta_all,'omitnan');
T_perm_all = zeros(nPerm,1);
for b=1:nPerm
    flips        = (rand(size(delta_all))<0.5)*2 - 1;
    T_perm_all(b)= mean(delta_all .* flips,'omitnan');
end
p_perm_all = mean(abs(T_perm_all) >= abs(T_obs_all));

% ---- console print: POPULATION mean firing rates per rat ----
fprintf('\n=== Mean POPULATION firing rates (per rat; speed-matched bins) ===\n');
for rr=1:nR
    fprintf('%s: trial=%.4f Hz nontrial=%.4f Hz Δ=%.4f Hz nBins=%d p_perm=%.4g\n', ...
        ratNames{rr}, ...
        perRat_meanTrialHz(rr), perRat_meanNonHz(rr), perRat_meanDelta(rr), ...
        perRat_nBins(rr), p_perm_rats(rr));
end
fprintf('Pooled: trial=%.4f Hz nontrial=%.4f Hz Δ=%.4f Hz nBins=%d p_perm=%.4g\n', ...
    mean(trial_all,'omitnan'), mean(non_all,'omitnan'), mean(delta_all,'omitnan'), ...
    numel(delta_all), p_perm_all);

% ---- simple pooled permutation figure ----
fig = figure('Color','w','Position',[200 260 520 420]); hold on;
histogram(T_perm_all,'Normalization','percentage','EdgeColor','none');
xline(T_obs_all,'r-','LineWidth',2);
yL = ylim;
text(T_obs_all,0.9*yL(2),sprintf('obs=%.3f Hz\np_{perm}=%.4g',T_obs_all,p_perm_all), ...
    'Color','r','FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top');
xlabel('Mean \Delta rate under sign-flip null (POP bins)');
ylabel('Density');
title(sprintf('POPULATION permutation (sign flip) %s', titleTag));
box on;

% ---- pack outputs ----
Pop = struct();
Pop.early_exit = false;

Pop.perRat = struct( ...
    'name',           ratNames(:), ...
    'nBins',          num2cell(perRat_nBins(:)), ...
    'binDelta',       perRat_binDeltas(:), ...
    'binTrialHz',     perRat_binTrialHz(:), ...
    'binNontrialHz',  perRat_binNonHz(:), ...
    'meanDelta',      num2cell(perRat_meanDelta(:)), ...
    'meanTrialHz',    num2cell(perRat_meanTrialHz(:)), ...
    'meanNontrialHz', num2cell(perRat_meanNonHz(:)), ...
    'perm_Tobs',      num2cell(T_obs_rats(:)), ...
    'perm_p',         num2cell(p_perm_rats(:)) );

Pop.pooled = struct( ...
    'binDelta',      delta_all, ...
    'binTrialHz',    trial_all, ...
    'binNontrialHz', non_all, ...
    'meanDelta',     mean(delta_all,'omitnan'), ...
    'meanTrialHz',   mean(trial_all,'omitnan'), ...
    'meanNontrialHz',mean(non_all,'omitnan'), ...
    'perm_Tobs',     T_obs_all, ...
    'perm_p',        p_perm_all, ...
    'T_perm_all',    T_perm_all );

Pop.figure = fig;
end

function plot_rate_change_bars_POP(R, titleTag, varargin)
% POPULATION-level bar chart of change Non-trial -> Trial per rat (+ pooled).
% Samples = speed-matched POPULATION bins (S.pop.trialRatePerBin / nonRatePerBin).

p = inputParser;
p.addParameter('metric','percent',@(s) any(validatestring(s,{'percent','fold','delta'})));
p.addParameter('errType','sd',@(s) any(validatestring(s,{'sd','sem'})));
p.addParameter('barColor',[0.30 0.60 0.85]);
p.addParameter('epsDen',1e-12,@(x) isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
metric   = p.Results.metric;
errType  = p.Results.errType;
barColor = p.Results.barColor;
epsDen   = p.Results.epsDen;

ratNames = {R.rat};
nR       = numel(R);

perRat_non = cell(nR,1);
perRat_tri = cell(nR,1);
for rr = 1:nR
    tri = []; non = [];
    for d = 1:numel(R(rr).perDay)
        S = R(rr).perDay{d};
        if isempty(S) || ~isfield(S,'pop') || isempty(S.pop)
            continue;
        end

        % POP fields are now Hz per cell; just append
        tri = [tri; S.pop.trialRatePerBin(:)];
        non = [non; S.pop.nonRatePerBin(:)];
    end
    perRat_tri{rr} = tri(:);   % Hz per cell
    perRat_non{rr} = non(:);   % Hz per cell
end

allTri = vertcat(perRat_tri{:});
allNon = vertcat(perRat_non{:});

vals       = nan(nR,1);
errs       = nan(nR,1);
nBinsRat   = zeros(nR,1);
meanNon    = nan(nR,1);
meanTri    = nan(nR,1);
meanDelta  = nan(nR,1);    % Hz
p_ttest    = nan(nR,1);    % NEW: per-rat POP t-test p-values

for rr = 1:nR
    t = perRat_tri{rr};
    n = perRat_non{rr};
    if isempty(t) || isempty(n), continue; end

    K = min(numel(t),numel(n));
    t = t(1:K);
    n = n(1:K);

    % paired t-test on bin-wise POP rates for this rat
    if K >= 2
        [~, p_ttest(rr)] = ttest(t, n);
    else
        p_ttest(rr) = NaN;
    end

    deltaBins = t - n;
    deltaBins = deltaBins(isfinite(deltaBins));
    if isempty(deltaBins), continue; end

    nBinsRat(rr)  = numel(deltaBins);
    meanTri(rr)   = mean(t,'omitnan');
    meanNon(rr)   = mean(n,'omitnan');
    meanDelta(rr) = mean(deltaBins,'omitnan');  % Hz

    sHz = std(deltaBins,'omitnan');
    if strcmpi(errType,'sem')
        errHz = sHz / sqrt(numel(deltaBins));
    else
        errHz = sHz;
    end

    denom = max(meanNon(rr), epsDen);
    switch metric
        case 'percent'
            vals(rr) = 100 * (meanDelta(rr) / denom);
            errs(rr) = 100 * (errHz       / denom);
        case 'fold'
            vals(rr) = meanTri(rr) / denom;
            errs(rr) = errHz       / denom;
        otherwise
            vals(rr) = meanDelta(rr);
            errs(rr) = errHz;
    end
end

% pooled across rats
K_all      = min(numel(allTri), numel(allNon));
allTri_use = allTri(1:K_all);
allNon_use = allNon(1:K_all);

delta_all   = allTri_use - allNon_use;
delta_all   = delta_all(isfinite(delta_all));
meanTri_all = mean(allTri_use,'omitnan');
meanNon_all = mean(allNon_use,'omitnan');
meanDelta_all = mean(delta_all,'omitnan');

sHz_all = std(delta_all,'omitnan');
if strcmpi(errType,'sem')
    errHz_all = sHz_all / sqrt(numel(delta_all));
else
    errHz_all = sHz_all;
end

denomAll = max(meanNon_all, epsDen);
switch metric
    case 'percent'
        valAll = 100 * (meanDelta_all / denomAll);
        errAll = 100 * (errHz_all     / denomAll);
    case 'fold'
        valAll = meanTri_all / denomAll;
        errAll = errHz_all   / denomAll;
    otherwise
        valAll = meanDelta_all;
        errAll = errHz_all;
end

% pooled POP t-test
if K_all >= 2
    [~, p_ttest_all] = ttest(allTri_use, allNon_use);
else
    p_ttest_all = NaN;
end

vals_plot = [vals; valAll];
errs_plot = [errs; errAll];
labels    = [ratNames, {'All'}];

figure('Color','w','Position',[460 260 780 420]); hold on
bar(vals_plot, 'FaceColor', barColor, 'EdgeColor','none');
errorbar(1:numel(vals_plot), vals_plot, errs_plot, errs_plot, '.k', ...
    'LineWidth',1.2, 'CapSize',10);
yline(0,'k--','HandleVisibility','off');

xticks(1:numel(vals_plot)); xticklabels(labels); xtickangle(15);
ylabel(ylabel_for_POP(metric));
title(sprintf('POPULATION change Non \\rightarrow Trial — %s  (%s bars)', ...
      titleTag, upper(errType)));
box on

ymin = min(vals_plot - errs_plot); ymax = max(vals_plot + errs_plot);
rng  = ymax - ymin; if rng<=0, rng = 1; end
ylim([ymin - 0.10*rng, ymax + 0.12*rng]);

unitLabel = struct('percent','%','fold','x','delta','Hz');
fprintf('\n=== POPULATION Non -> Trial change per rat (%s) with %s bars ===\n', ...
    metric, upper(errType));
fprintf('%-8s  nBins  BarMetric(%s)   %s_low   %s_high   Non(Hz)   Trial(Hz)   Delta(Hz)   p_ttest\n', ...
        'Rat', unitLabel.(metric), upper(errType), upper(errType));

for rr = 1:nR
    if isnan(vals(rr)), continue; end
    ciLo = vals(rr) - errs(rr);
    ciHi = vals(rr) + errs(rr);
    fprintf('%-8s  %6d  %10s   %8s  %8s   %7.4f   %8.4f   %8.4f   %9.3g\n', ...
        ratNames{rr}, nBinsRat(rr), ...
        pretty_val_POP(vals(rr), metric), ...
        pretty_val_POP(ciLo, metric), pretty_val_POP(ciHi, metric), ...
        meanNon(rr), meanTri(rr), (meanTri(rr)-meanNon(rr)), p_ttest(rr));
end

ciLo = valAll - errAll; ciHi = valAll + errAll;
fprintf('%-8s  %6d  %10s   %8s  %8s   %7.4f   %8.4f   %8.4f   %9.3g\n', ...
    'All', numel(delta_all), ...
    pretty_val_POP(valAll, metric), ...
    pretty_val_POP(ciLo,  metric), pretty_val_POP(ciHi, metric), ...
    meanNon_all, meanTri_all, (meanTri_all - meanNon_all), p_ttest_all);
end

% ---- helpers ----
function lab = ylabel_for_POP(metric)
switch metric
    case 'percent', lab = 'POPULATION change from non-trial (%)';
    case 'fold',    lab = 'POPULATION fold change (Trial / Non-trial)';
    otherwise,      lab = 'POPULATION \Delta rate (Hz)   Trial - Non-trial';
end
end

function s = pretty_val_POP(v, metric)
switch metric
    case 'percent', s = sprintf('%.1f%%', v);
    case 'fold',    s = sprintf('%.2f×',  v);
    otherwise,      s = sprintf('%.3f',   v);
end
end

function plot_paired_means_slopegraph_POP(R, titleTag)
% plot_paired_means_slopegraph_POP
% Non-trial vs Trial POPULATION means per rat, using speed/accel-matched
% POP bins stored in S.pop.trialRatePerBin / S.pop.nonRatePerBin.
% RATES ARE CONVERTED TO Hz PER CELL using S.pop.nCells.

if nargin < 2
    titleTag = '';
end

ratNames = {R.rat};
nR       = numel(R);

meanNon_pop = nan(nR,1);   % mean POP non-trial rate per cell (Hz)
meanTri_pop = nan(nR,1);   % mean POP trial rate per cell (Hz)

for rr = 1:nR
    triBins = [];
    nonBins = [];

    for d = 1:numel(R(rr).perDay)
        S = R(rr).perDay{d};
        if isempty(S) || ~isfield(S,'pop') || isempty(S.pop)
            continue;
        end

        % raw POP rates (sum across cells for this day)
        tri = S.pop.trialRatePerBin(:);
        non = S.pop.nonRatePerBin(:);

        % convert to per-cell Hz using stored nCells if available
        nC = 1;
        if isfield(S.pop,'nCells') && ~isempty(S.pop.nCells) && S.pop.nCells > 0
            nC = S.pop.nCells;
        end

        triBins = [triBins; tri(:) ./ nC]; %#ok<AGROW>
        nonBins = [nonBins; non(:) ./ nC]; %#ok<AGROW>
    end

    meanTri_pop(rr) = mean(triBins,'omitnan');  % Hz per cell
    meanNon_pop(rr) = mean(nonBins,'omitnan');  % Hz per cell
end

% use only rats with finite POP means
valid = isfinite(meanNon_pop) & isfinite(meanTri_pop);
if ~any(valid)
    warning('plot_paired_means_slopegraph_POP: no finite POP means; skipping plot.');
    return;
end

figure('Color','w','Position',[320 280 520 420]); hold on

% slope lines Non -> Trial
for rr = find(valid)'
    plot([1 2], [meanNon_pop(rr) meanTri_pop(rr)], '-', ...
        'Color',[0.6 0.6 0.6], 'LineWidth',1.6);
end

% points at Non and Trial
scatter(ones(sum(valid),1), meanNon_pop(valid), 40, [0.35 0.65 0.95], 'filled');
scatter(2*ones(sum(valid),1), meanTri_pop(valid), 40, [0.95 0.45 0.45], 'filled');

xlim([0.8 2.2]);
xticks([1 2]);
xticklabels({'Non-trial','Trial'});

valsY = [meanNon_pop(valid); meanTri_pop(valid)];
ymin = min(valsY);
ymax = max(valsY);
pad  = 0.1 * max(ymax - ymin, eps);
ylim([ymin - pad, ymax + pad]);

yline(0,'k--');
box on
ylabel('Mean POPULATION rate per cell (Hz)');
title(sprintf('Non vs Trial (per-rat POPULATION means) %s', titleTag));

% rat labels on Non-trial side (only for valid rats)
vIdx = find(valid);
for ii = 1:numel(vIdx)
    rr = vIdx(ii);
    text(0.98, meanNon_pop(rr), ratNames{rr}, ...
        'HorizontalAlignment','right', 'Color',[0.4 0.4 0.4]);
end
end

function speedBin_dayLevel_heatmap(R, varargin)
p = inputParser;
p.addParameter('titleTag',''); p.parse(varargin{:});
titleTag = p.Results.titleTag;

ratNames = {R.rat};
nR = numel(R);
Pmat = nan(nR,3); Delta = nan(nR,3); N = zeros(nR,3);

for rr = 1:nR
  for d = 1:3
    S = R(rr).perDay{d};
    if isempty(S) || ~isfield(S,'dayLevel'), continue; end
    Pmat(rr,d) = S.dayLevel.p;
    Delta(rr,d)= S.dayLevel.delta;
    N(rr,d)    = S.dayLevel.nCells;
  end
end

Z = -log10(Pmat);
figure('Color','w','Position',[220 220 740 420]);
imagesc(Z); colorbar; caxis([0 5]);
colormap(parula);
xticks(1:3); xticklabels({'day−2','day−1','day 0'});
yticks(1:nR); yticklabels(ratNames);
title(sprintf('Day-level paired t across cells: -log10(p) %s', titleTag));
for rr=1:nR
  for d=1:3
    if ~isnan(Z(rr,d))
      txt = sprintf('n=%d \\Delta=%.3f', N(rr,d), Delta(rr,d));
      text(d, rr, txt, 'Color','w', 'HorizontalAlignment','center','FontSize',8);
    end
  end
end
end

function speedBin_dayLevel_bars(R, varargin)
p = inputParser; p.addParameter('titleTag',''); p.parse(varargin{:});
titleTag = p.Results.titleTag;

ratNames = {R.rat};
nR = numel(R);
D = nan(nR,3); SE = nan(nR,3);

for rr = 1:nR
  for d = 1:3
    S = R(rr).perDay{d};
    if isempty(S) || ~isfield(S,'dayLevel'), continue; end
    tested = isfinite(S.meanTrialRate) & isfinite(S.meanNonTrialRate);
    diffs  = S.meanTrialRate(tested) - S.meanNonTrialRate(tested);
    D(rr,d)= mean(diffs,'omitnan');
    SE(rr,d)= std(diffs,'omitnan')/sqrt(sum(isfinite(diffs)));
  end
end

figure('Color','w','Position',[220 260 1100 420]); hold on
ng = 3; bw = 0.22; off = [-bw, 0, bw];
for d = 1:3
  x = (1:nR) + off(d);
  bar(x, D(:,d), bw, 'EdgeColor','none');
  errorbar(x, D(:,d), SE(:,d), '.k', 'LineWidth',1.2, 'CapSize',10);
end
yline(0,'k--');
xticks(1:nR); xticklabels(ratNames); xtickangle(15);
legend({'day−2','day−1','day 0'},'Location','best');
ylabel('\Delta rate (Hz)  Trial − Non (across cells)');
title(sprintf('Day-level paired differences (means ± SE) %s', titleTag));
box on
end
