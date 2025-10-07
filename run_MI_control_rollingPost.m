function run_MI_control_rollingPost(WEdges)
% run_MI_control_rollingPost
% CS-aligned rolling removal:
%   Test   : remove spikes ONLY in the current window.
%   Control: remove same # spikes ONLY from OUTSIDE the full span.
% Optional tighter controls (speed / space / both), plus spike-count diagnostics.
% Group statistics are PAIRED t-tests of MI_rand vs MI_win (no "vs 0" anywhere).

% ---- analysis windows (allow negatives) ----
if nargin < 1 || isempty(WEdges)
    starts = -4:.5:16;                 % window start times (s)

    %starts = -2:.5:14;                 % window start times (s)
    %widths = ones(numel(starts),1)*2;  % 2-s wide windows
    widths = ones(numel(starts),1)*2;  % 2-s wide windows
    WEdges = [starts(:), starts(:)+widths(:)];
end

% ---- CONTROL MATCHING OPTIONS ----
ControlMatch = 'speed';   % 'none' | 'speed' | 'space' | 'speedspace'
NSpeedBins   = 200;       % used when ControlMatch includes 'speed'
SpaceBinSize = [];        % [] => use MI 'dim'; else numeric bin size (same units as pos)

% ---- config ----
ratNames  = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
velthresh = 4;    % cm/s
dim       = 2.5;
nIter     = 5;
alpha     = 0.05;

fprintf('\n===== Running MI_control_rollingPost (CS-aligned; control outside span; match=%s) =====\n', ControlMatch);

winCenters  = [];
curvesDelta = [];   % W × N_animals  (Δ = rand - win)
curvesRand  = [];   % W × N_animals
curvesWin   = [];   % W × N_animals
names       = {};

set(0,'DefaultFigureVisible','on');

for ii = 1:numel(ratNames)
    ratVar = ratNames{ii};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratVar))
        warning('Variable %s not found in base workspace. Skipping.', ratVar);
        continue;
    end
    rat = evalin('base', ratVar);

    % last 3 days up to An
    dateList = autoDateList(rat);
    idx = find(strcmp(dateList, rat.An));
    if isempty(idx) || idx < 3
        warning('%s missing enough days before An. Skipping...', ratVar);
        continue;
    end
    daysToUse = dateList(idx-2:idx);

    % needed fields
    spikes   = filterFieldsByDay(rat.Ca_peaks, daysToUse);
    pos      = filterFieldsByDay(rat.pos,       daysToUse);
    ts       = filterFieldsByDay(rat.Ca_ts,     daysToUse); %#ok<NASGU>
    cs_times = filterFieldsByDay(rat.CS_times,  daysToUse); % ACTUAL CS onsets

    % run rolling CS-aligned analysis with control-matching options
    out = MI_control_rollingPost(spikes, pos, velthresh, dim, ts, cs_times, ...
                                 nIter, WEdges, ControlMatch, NSpeedBins, SpaceBinSize);

    % spike-count parity (optional)
    print_spike_count_parity(out, ratVar);

    % aggregate per animal: Δ, Rand, Win curves
    [centersW, curveDeltaW, curveRandW, curveWinW, delta_allPoints, dayMeans] = ...
        aggregate_animal(out, rat, daysToUse); %#ok<ASGLU>  % dayMeans kept for possible later use

    if isempty(curveDeltaW)
        warning('%s: no usable data (check MI_win/MI_rand). Skipping.', ratVar);
        continue;
    end

    % harmonize x-axis
    if isempty(winCenters)
        winCenters = centersW(:);
    elseif numel(centersW) ~= numel(winCenters) || any(centersW(:) ~= winCenters(:))
        curveDeltaW = interp1(centersW(:), curveDeltaW(:), winCenters(:), 'linear', 'extrap');
        curveRandW  = interp1(centersW(:), curveRandW(:),  winCenters(:), 'linear', 'extrap');
        curveWinW   = interp1(centersW(:), curveWinW(:),   winCenters(:), 'linear', 'extrap');
    end

    curvesDelta(:, end+1) = curveDeltaW(:); %#ok<AGROW>
    curvesRand(:,  end+1) = curveRandW(:);  %#ok<AGROW>
    curvesWin(:,   end+1) = curveWinW(:);   %#ok<AGROW>
    names{end+1}          = ratVar;         %#ok<AGROW>

    % ---- per-animal within-animal per-window tests (paired across cells) ----
    [p_cellpair, q_cellpair, nWin, nSig_unc, nSig_fdr] = per_animal_window_tests_paired(out, rat, daysToUse, alpha); %#ok<ASGLU>
    muA  = mean(curveDeltaW,'omitnan');
    seA  = std(curveDeltaW,0,'omitnan')/sqrt(sum(~isnan(curveDeltaW)));
    fprintf('%s: mean Δ curve = %.3g ± %.3g (SEM across windows)\n', ratVar, muA, seA);
end

% ---------- Plot ----------
if isempty(curvesDelta)
    warning('No animals produced curves. Nothing to plot.');
    return;
end

set(0,'DefaultFigureVisible', 'on');
figure('Color','w'); hold on;
C = lines(size(curvesDelta,2));
for j = 1:size(curvesDelta,2)
    plot(winCenters, curvesDelta(:,j), '-', 'Color', [C(j,:) 0.75], 'LineWidth', 1.5);
end

% group mean ± SEM across animals (for Δ)
m  = nanmean(curvesDelta, 2);
se = nanstd(curvesDelta, 0, 2) ./ sqrt(size(curvesDelta,2));
xx = [winCenters; flipud(winCenters)];
yy = [m-se; flipud(m+se)];
fill(xx, yy, [0 0 0], 'FaceAlpha', 0.10, 'EdgeColor', 'none');
plot(winCenters, m, 'k-', 'LineWidth', 2.5);

% ----- Group per-window *paired* tests: MI_rand vs MI_win across animals -----
alpha = 0.05;
pvals = nan(numel(winCenters),1);
for w = 1:numel(winCenters)
    x = curvesRand(w,:); y = curvesWin(w,:);
    mask = ~isnan(x) & ~isnan(y);
    if nnz(mask) >= 2
        [~,pvals(w)] = ttest(x(mask), y(mask), 'Alpha', alpha);  % paired by animal
    end
end
qvals = bh_fdr(pvals);
sigIdx = find(qvals < alpha);
ymin = min(m-se);
for s = sigIdx(:)'
    plot(winCenters(s), ymin - 0.1*range(m), 'kv', 'MarkerFaceColor','k', 'MarkerSize',8);
end

xlabel('Window center relative to CS (s)');
ylabel('\Delta MI  (MI_{rand} - MI_{win})');
title(sprintf('Rolling control (CS-aligned; control outside span; match=%s)', ControlMatch));
grid on;
legend([names, {'Mean±SEM'}], 'Interpreter','none', 'Location','best');

outname = sprintf('MI_rollingPost_curves_CSaligned_match-%s.png', ControlMatch);
exportgraphics(gcf, outname, 'Resolution', 300);
fprintf('Saved figure: %s\n', outname);

% ---------- Group per-window table (paired across animals) ----------
fprintf('\n========== Group per-window results (paired t-test: MI_{rand} vs MI_{win} across animals) ==========\n');
fprintf('%10s  %12s  %12s  %12s  %8s\n', 'Center(s)', 'Mean Δ', 'SEM Δ', 'p_pair', 'q');
for w = 1:numel(winCenters)
    fprintf('%10.2f  %12.4g  %12.4g  %12.3g  %8.3f\n', ...
        winCenters(w), m(w), se(w), pvals(w), qvals(w));
end
fprintf('Significant windows (BH-FDR q<%.2f): %s\n', alpha, mat2str(winCenters(sigIdx),3));
fprintf('====================================================================================================\n');

end


% ------------------------------- helpers -------------------------------

function [centers, curveDelta, curveRand, curveWin, delta_allPoints, dayMeans] = aggregate_animal(out_struct, rat, daysToUse)
% Returns:
%  centers      : W×1 midpoints (from the first day encountered)
%  curveDelta   : W×1 mean over days×cells of (MI_rand − MI_win)
%  curveRand    : W×1 mean over days×cells of MI_rand
%  curveWin     : W×1 mean over days×cells of MI_win
%  delta_allPoints: vector of all pointwise Δ across windows×days×cells
%  dayMeans     : 1×Ndays mean Δ per day (mean over windows×cells)

centers = []; curveDelta = []; curveRand = []; curveWin = [];
delta_allPoints = []; dayMeans = [];

if isempty(out_struct) || ~isstruct(out_struct), return; end
fns = fieldnames(out_struct);
dayFns = fns(startsWith(fns,'MI_'));
if isempty(dayFns), return; end

Wref = [];               % number of windows (to align)
sumDelta = []; cntDelta = [];   % for Δ curve
sumRand  = []; cntRand  = [];   % for MI_rand curve
sumWin   = []; cntWin   = [];   % for MI_win curve

for di = 1:numel(dayFns)
    day = out_struct.(dayFns{di});
    if ~isfield(day,'winEdges') || ~isfield(day,'MI_rand') || ~isfield(day,'MI_win')
        continue;
    end

    W = size(day.winEdges,1);
    if isempty(centers)
        centers  = nanmean(day.winEdges,2);
        Wref     = W;
        sumDelta = zeros(W,1); cntDelta = zeros(W,1);
        sumRand  = zeros(W,1); cntRand  = zeros(W,1);
        sumWin   = zeros(W,1); cntWin   = zeros(W,1);
    end

    % ratemask (optional)
    maskField = sprintf('ratemask_%s', dayFns{di}(4:end));
    maskVec = true(1, size(day.MI_rand,2));
    if isfield(rat,'ratemask') && isfield(rat.ratemask, maskField)
        mv = rat.ratemask.(maskField)(:);
        if numel(mv) == size(day.MI_rand,2)
            maskVec = (mv == 1);
        end
    end

    MIrand = day.MI_rand(:, maskVec);
    MIwin  = day.MI_win(:,  maskVec);
    if isempty(MIrand) || isempty(MIwin), continue; end

    Delta = MIrand - MIwin;

    % day means across cells (per window)
    dayMeanDelta = nanmean(Delta, 2);
    dayMeanRand  = nanmean(MIrand, 2);
    dayMeanWin   = nanmean(MIwin,  2);

    sumDelta = sumDelta + dayMeanDelta;
    cntDelta = cntDelta + ~all(isnan(Delta),2);

    sumRand  = sumRand  + dayMeanRand;
    cntRand  = cntRand  + ~all(isnan(MIrand),2);

    sumWin   = sumWin   + dayMeanWin;
    cntWin   = cntWin   + ~all(isnan(MIwin),2);

    % collect ALL points (for all-points test A)
    delta_allPoints = [delta_allPoints; Delta(:)]; %#ok<AGROW>

    % day mean (for day-means test B)
    dayMeans = [dayMeans, mean(Delta,'all','omitnan')]; %#ok<AGROW>
end

if any(cntDelta>0)
    curveDelta = sumDelta ./ max(1,cntDelta);
    curveRand  = sumRand  ./ max(1,cntRand);
    curveWin   = sumWin   ./ max(1,cntWin);
else
    curveDelta = []; curveRand = []; curveWin = [];
end
end


function [p_unc, q_fdr, nW, nSig_unc, nSig_fdr] = per_animal_window_tests(out_struct, rat, daysToUse, alpha)
p_unc = []; q_fdr = []; nW = 0; nSig_unc = 0; nSig_fdr = 0;
if isempty(out_struct) || ~isstruct(out_struct), return; end
fns = fieldnames(out_struct);
dayFns = fns(startsWith(fns,'MI_'));
if isempty(dayFns), return; end

day0 = out_struct.(dayFns{1});
if ~isfield(day0,'winEdges'), return; end
centers = mean(day0.winEdges,2);
W = numel(centers);

Delta_byW = cell(W,1);
for di = 1:numel(dayFns)
    day = out_struct.(dayFns{di});
    if ~isfield(day,'winEdges') || ~isfield(day,'MI_rand') || ~isfield(day,'MI_win'), continue; end

    MIrand = day.MI_rand; MI_win = day.MI_win;
    if isempty(MIrand) || isempty(MI_win), continue; end

    maskField = sprintf('ratemask_%s', dayFns{di}(4:end));
    maskVec = true(1, size(MIrand,2));
    if isfield(rat,'ratemask') && isfield(rat.ratemask, maskField)
        mv = rat.ratemask.(maskField)(:);
        if numel(mv) == size(MIrand,2), maskVec = (mv==1); end
    end
    MIrand = MIrand(:,maskVec); MI_win = MI_win(:,maskVec);
    Delta  = MIrand - MI_win;

    for w = 1:W
        Delta_byW{w} = [Delta_byW{w}, Delta(w,:)]; %#ok<AGROW>
    end
end

p_unc = nan(W,1);
for w = 1:W
    x = Delta_byW{w}; x = x(~isnan(x));
    if numel(x) >= 2
        [~,p_unc(w)] = ttest(x, 0, 'Alpha', alpha);
    end
end
q_fdr = bh_fdr(p_unc); nW = W;
nSig_unc = sum(p_unc < alpha, 'omitnan');
nSig_fdr = sum(q_fdr < alpha, 'omitnan');
end

function q = bh_fdr(p)
q = nan(size(p));
ps = p(:); valid = ~isnan(ps); pv = ps(valid);
[sv, order] = sort(pv); m = numel(sv);
if m==0, q = reshape(q, size(p)); return; end
qv = sv .* m ./ (1:m)'; for i=m-1:-1:1, qv(i)=min(qv(i),qv(i+1)); end
q_full = nan(size(ps)); idx=zeros(size(order)); idx(order)=1:numel(order);
q_full(valid)=qv(idx); q = reshape(q_full, size(p));
end

function print_spike_count_parity(out_struct, label)
% Parity of kept spikes for TEST vs CONTROL, per window.
% Only counts cells where CONTROL succeeded (Iter_OK > 0).

if nargin < 2, label = ''; end
if isempty(out_struct) || ~isstruct(out_struct), return; end
fns = fieldnames(out_struct);
dayFns = fns(startsWith(fns,'MI_')); if isempty(dayFns), return; end

day0 = out_struct.(dayFns{1});
if ~isfield(day0,'winEdges') || ~isfield(day0,'counts'), return; end
centers = nanmean(day0.winEdges,2);  W = numel(centers);

K_keep_win  = [];  % W x cells (concat days)
K_keep_ctrl = [];
K_iter_ok   = [];  % W x cells

for di = 1:numel(dayFns)
    day = out_struct.(dayFns{di});
    if ~isfield(day,'counts'), continue; end
    K_keep_win  = [K_keep_win,  day.counts.N_keep_win];   %#ok<AGROW>
    K_keep_ctrl = [K_keep_ctrl, day.counts.N_keep_ctrl];  %#ok<AGROW>
    K_iter_ok   = [K_iter_ok,   day.counts.Iter_OK];      %#ok<AGROW>
end

fprintf('\n[spike parity %s]\n', label);
fprintf('%9s  %12s  %12s  %9s  %8s\n', ...
    'Center(s)', 'med keep(win)', 'med keep(ctrl)', '% parity', 'n OK');

for w = 1:W
    kw = K_keep_win(w, :);
    kc = K_keep_ctrl(w, :);
    ok = K_iter_ok(w, :) > 0;                         % control succeeded
    valid = ok & isfinite(kw) & isfinite(kc);
    if any(valid)
        pct_equal = 100 * mean(kw(valid) == kc(valid));
        n_ok = sum(valid);
        med_win  = median(kw(valid), 'omitnan');
        med_ctrl = median(kc(valid), 'omitnan');
    else
        pct_equal = NaN; n_ok = 0; med_win = NaN; med_ctrl = NaN;
    end
    fprintf('%9.2f  %12.1f  %12.1f  %8.1f%%  %8d\n', ...
        centers(w), med_win, med_ctrl, pct_equal, n_ok);
end
end

function [p_pair, q_fdr, nW, nSig_unc, nSig_fdr] = per_animal_window_tests_paired(out_struct, rat, daysToUse, alpha)
% For printing: within-animal per-window *paired* tests across cells.
% Pools cells across the animal's analyzed days for each window;
% paired t-test of MI_rand vs MI_win at each window.

p_pair = [];
q_fdr  = [];
nW = 0; nSig_unc = 0; nSig_fdr = 0;

if isempty(out_struct) || ~isstruct(out_struct), return; end
fns = fieldnames(out_struct);
dayFns = fns(startsWith(fns,'MI_'));
if isempty(dayFns), return; end

% find reference W and centers
day0 = out_struct.(dayFns{1});
if ~isfield(day0,'winEdges'), return; end
centers = mean(day0.winEdges,2);
W = numel(centers);

% build, for each window, paired lists (rand, win) across all days×cells
R_byW = cell(W,1);
W_byW = cell(W,1);
for di = 1:numel(dayFns)
    day = out_struct.(dayFns{di});
    if ~isfield(day,'winEdges') || ~isfield(day,'MI_rand') || ~isfield(day,'MI_win'), continue; end
    MIrand = day.MI_rand; MIwin = day.MI_win;
    if isempty(MIrand) || isempty(MIwin), continue; end

    % ratemask
    maskField = sprintf('ratemask_%s', dayFns{di}(4:end));
    maskVec = true(1, size(MIrand,2));
    if isfield(rat,'ratemask') && isfield(rat.ratemask, maskField)
        mv = rat.ratemask.(maskField)(:);
        if numel(mv) == size(MIrand,2), maskVec = (mv==1); end
    end
    MIrand = MIrand(:,maskVec); MIwin = MIwin(:,maskVec);

    % (optionally align to reference if window centers differ)
    centers_d = mean(day.winEdges,2);
    if numel(centers_d) ~= W || any(centers_d(:) ~= centers(:))
        % interpolate both matrices row-wise to reference centers
        MIrand = interp1(centers_d(:), MIrand, centers(:), 'linear', 'extrap');
        MIwin  = interp1(centers_d(:), MIwin,  centers(:), 'linear', 'extrap');
    end

    for w = 1:W
        R_byW{w} = [R_byW{w}, MIrand(w,:)]; %#ok<AGROW>
        W_byW{w} = [W_byW{w}, MIwin(w,:)];  %#ok<AGROW>
    end
end

p_pair = nan(W,1);
for w = 1:W
    x = R_byW{w}; y = W_byW{w};
    mask = isfinite(x) & isfinite(y);
    if nnz(mask) >= 2
        [~,p_pair(w)] = ttest(x(mask), y(mask), 'Alpha', alpha); % paired across cells
    end
end
q_fdr = bh_fdr(p_pair);
nW = W;
nSig_unc = sum(p_pair < alpha, 'omitnan');
nSig_fdr = sum(q_fdr < alpha, 'omitnan');
end
