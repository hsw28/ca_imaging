function run_MI_control_rollingPost(WEdges)
% run_MI_control_rollingPost
% CS-aligned rolling removal:
%   Test  : remove spikes ONLY in the current window.
%   Control: remove same # spikes ONLY from OUTSIDE the full span.
% Optional tighter controls (speed / space / both), plus spike-count diagnostics.


% ---- analysis windows (allow negatives) ----
if nargin < 2 || isempty(WEdges)
    starts = -4:.5:5;                 % window start times (s)
    widths = ones(numel(starts),1)*2;   % 2-s wide windows
    WEdges = [starts(:), starts(:)+widths(:)];
end

% ---- CONTROL MATCHING OPTIONS ----
ControlMatch = 'speed';        % 'none' | 'speed' | 'space' | 'speedspace'
NSpeedBins   = 50;  %10           % used when ControlMatch includes 'speed'
SpaceBinSize = [];            % [] => use MI 'dim'; else numeric bin size (same units as pos)

% ---- config ----
ratNames  = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
velthresh = 4;    % cm/s
dim       = 2.5;
nIter     = 5;
alpha     = 0.05;

fprintf('\n===== Running MI_control_rollingPost (CS-aligned; control outside span; match=%s) =====\n', ControlMatch);

winCenters   = [];
curves       = [];
names        = {};
allDelta     = [];
allDayMeans  = [];

set(0,'DefaultFigureVisible','on');

for ii = 1:numel(ratNames)
    ratVar = ratNames{ii}
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

    % ---- spike-count parity print (per animal) ----
    print_spike_count_parity(out, ratVar);

    % aggregate per animal
    [centersW, curveW, delta_allPoints, dayMeans] = aggregate_animal(out, rat, daysToUse);
    if isempty(curveW)
        warning('%s: no usable data (check MI_win/MI_rand). Skipping.', ratVar);
        continue;
    end

    if isempty(winCenters)
        winCenters = centersW(:);
    elseif numel(centersW) ~= numel(winCenters) || any(centersW(:) ~= winCenters(:))
        curveW = interp1(centersW(:), curveW(:), winCenters(:), 'linear', 'extrap');
    end

    curves(:, end+1) = curveW(:); %#ok<AGROW>
    names{end+1}     = ratVar;    %#ok<AGROW>

    allDelta     = [allDelta; delta_allPoints(:)]; %#ok<AGROW>
    allDayMeans  = [allDayMeans; dayMeans(:)];     %#ok<AGROW>

    % per-animal ALL-POINTS Δ vs 0
    [~,p_pair,~,stats_pair] = ttest(delta_allPoints, 0, 'Alpha', alpha);
    muA  = mean(curveW,'omitnan');
    seA  = std(curveW,0,'omitnan')/sqrt(sum(~isnan(curveW)));

    % per-window within-animal across-cells tests (print summary counts)
    [p_unc, q_fdr, nWin, nSig_unc, nSig_fdr] = per_animal_window_tests(out, rat, daysToUse, alpha);

    fprintf('%s: mean curve Δ = %.3f ± %.3f; all-points: t(%d)=%.3f, p=%.3g\n', ...
            ratVar, muA, seA, stats_pair.df, stats_pair.tstat, p_pair);
    fprintf('   per-window across-cells (uncorr/FDR): %d/%d and %d/%d; min p=%.3g\n', ...
            nSig_unc, nWin, nSig_fdr, nWin, min(p_unc,[],'omitnan'));
end

% ---------- Plot ----------
if isempty(curves)
    warning('No animals produced curves. Nothing to plot.');
    return;
end

set(0,'DefaultFigureVisible', 'on');
figure('Color','w'); hold on;
C = lines(size(curves,2));
for j = 1:size(curves,2)
    plot(winCenters, curves(:,j), '-', 'Color', [C(j,:) 0.75], 'LineWidth', 1.5);
end

% group mean ± SEM across animals
m  = nanmean(curves, 2);
se = nanstd(curves, 0, 2) ./ sqrt(size(curves,2));
xx = [winCenters; flipud(winCenters)];
yy = [m-se; flipud(m+se)];
fill(xx, yy, [0 0 0], 'FaceAlpha', 0.10, 'EdgeColor', 'none');
plot(winCenters, m, 'k-', 'LineWidth', 2.5);

% per-window group t-tests across animals' curves (vs 0) + BH-FDR
alpha = 0.05;
pvals = nan(size(m));
for w = 1:numel(winCenters)
    [~,pvals(w)] = ttest(curves(w,:), 0, 'Alpha', alpha);
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

% ---------- Printed tables ----------
fprintf('\n================ Group per-window results (ttest across animals vs 0) ================\n');
fprintf('%10s  %12s  %12s  %12s  %8s\n', 'Center(s)', 'Mean', 'SEM', 'p', 'q');
for w = 1:numel(winCenters)
    fprintf('%10.2f  %12.4g  %12.4g  %12.3g  %8.3f\n', ...
        winCenters(w), m(w), se(w), pvals(w), qvals(w));
end
fprintf('Significant windows (BH-FDR q<%.2f): %s\n', alpha, mat2str(winCenters(sigIdx),3));
fprintf('======================================================================================\n');

% ---------- Overall group tests ----------
% (A) ALL POINTS
x_all = allDelta(:); x_all = x_all(~isnan(x_all));
[~,p_all,~,stats_all] = ttest(x_all, 0, 'Alpha', alpha);
mu_all = mean(x_all,'omitnan');
se_all = std(x_all,0,'omitnan') / sqrt(numel(x_all));
fprintf('\n================ Overall group test — ALL POINTS =================\n');
fprintf('Mean Δ over all points: %g ± %g (SEM), N=%d\n', mu_all, se_all, numel(x_all));
fprintf('t(%d) = %.3f, p = %.3g\n', stats_all.df, stats_all.tstat, p_all);
fprintf('==================================================================\n');

% (B) DAY MEANS
x_days = allDayMeans(:); x_days = x_days(~isnan(x_days));
if ~isempty(x_days)
    [~,p_days,~,stats_days] = ttest(x_days, 0, 'Alpha', alpha);
    mu_days = mean(x_days,'omitnan');
    se_days = std(x_days,0,'omitnan') / sqrt(numel(x_days));
    fprintf('\n================ Overall group test — DAY MEANS ===================\n');
    fprintf('Mean Δ of animal×day means: %g ± %g (SEM), Ndays=%d\n', mu_days, se_days, numel(x_days));
    fprintf('t(%d) = %.3f, p = %.3g\n', stats_days.df, stats_days.tstat, p_days);
    fprintf('==================================================================\n');
else
    fprintf('\n[Info] No day means available (check MI_win/MI_rand presence per day).\n');
end
end


% ------------------------------- helpers -------------------------------

function [centers, curve, delta_allPoints, dayMeans] = aggregate_animal(out_struct, rat, daysToUse)
centers = []; curve = []; delta_allPoints = []; dayMeans = [];
if isempty(out_struct) || ~isstruct(out_struct), return; end
fns = fieldnames(out_struct);
dayFns = fns(startsWith(fns,'MI_'));
if isempty(dayFns), return; end

Wref = []; sumW = []; cntW = [];
for di = 1:numel(dayFns)
    day = out_struct.(dayFns{di});
    if ~isfield(day,'winEdges') || ~isfield(day,'MI_rand') || ~isfield(day,'MI_win'), continue; end
    W = size(day.winEdges,1);
    if isempty(centers)
        centers = nanmean(day.winEdges,2);
        Wref = W; sumW = zeros(W,1); cntW = zeros(W,1);
    end

    maskField = sprintf('ratemask_%s', dayFns{di}(4:end));
    maskVec = true(1, size(day.MI_rand,2));
    if isfield(rat,'ratemask') && isfield(rat.ratemask, maskField)
        mv = rat.ratemask.(maskField)(:);
        if numel(mv) == size(day.MI_rand,2), maskVec = (mv==1); end
    end

    MI_rand = day.MI_rand(:, maskVec);
    MI_win  = day.MI_win(:,  maskVec);
    if isempty(MI_rand) || isempty(MI_win), continue; end

    Delta  = MI_rand - MI_win;

    sumW = sumW + nanmean(Delta, 2);
    cntW = cntW + ~all(isnan(Delta),2);

    delta_allPoints = [delta_allPoints; Delta(:)]; %#ok<AGROW>
    dayMeans = [dayMeans, mean(Delta,'all','omitnan')]; %#ok<AGROW>
end

if any(cntW>0)
    curve = sumW ./ max(1,cntW);
else
    curve = [];
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
