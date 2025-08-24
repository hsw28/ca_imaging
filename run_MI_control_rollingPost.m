function run_MI_control_rollingPost(csusNum, WEdges)
% run_MI_control_rollingPost
% Calls MI_control_rollingPost per animal, aggregates curves across the last
% 3 days, plots one curve/animal (mean over days×cells of MI_rand−MI_win),
% overlays group mean±SEM, and prints the requested statistics:
%  - Per-animal paired t-test using *all points* (windows×days×cells):
%       ttest(MI_rand(:), MI_win_expanded(:))
%  - Group tests:
%       (A) All points: one-sample ttest of Δ = (MI_rand − MI_win) vs 0
%       (B) Day means: one-sample ttest across animal×day means of Δ vs 0
%
% Uses ratemask if available: rat.ratemask.ratemask_YYYY_MM_DD == 1.

%WEdges = [0 1; .5 1.5; 1 2; 1.5 2.5; 2 3; 2.5 3.5; 3 4; 3.5 4.5; 4 5];
%  WEdges = [0 2; 1 3; 2 4];

%WEdges = [0 1.5; .75 2.25; 1.5 3; 2.25 3.75];

%WEdges = [0 1; .5 1.5; 1 2; 1.5 2.5; 2 3; 2.5 3.5];

%WEdges = [-1 1; 0 2; 1 3; 2 4; 3 5; 4 6];

%WEdges = [0 2; 1 3; 2 4; 3 5; 4 6; 5 7; 8 10; 9 11; 10 12; 11 13];

%WEdges = [0 3; 1 4; 2 5; 3 6; 4 7; 5 8; 6 9; 7 10; 8 11; 9 12];
%WEdges = [0 3; 1 4; 2 5; 3 6; 4 7];


%WEdges = [0 2; 1 3; 2 4; 3 5; 4 6; 5 7; 8 10; 9 11; 10 12; 11 13];

starts = 0:.5:11;       % window start times
widths = ones(numel(starts),1)*2; % here always width 2
WEdges = [starts(:), starts(:)+widths(:)];

if nargin < 1, csusNum = 15; end

% ---- config ----

ratNames  = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
velthresh = 4;    % cm/s
dim       = 2.5;
nIter     = 100;
alpha     = 0.05;

fprintf('\n===== Running MI_control_rollingPost (csus%d) =====\n', csusNum);

winCenters = [];           % shared x-axis (filled after first animal)
curves     = [];           % W × N_animals (each col = animal curve)
names      = {};           % animal names included
allDelta   = [];           % collects ALL Δ across animals (for test A)
allDayMeans = [];          % collects animal×day means (for test B)

set(0,'DefaultFigureVisible','on');

for ii = 1:numel(ratNames)
    ratVar = ratNames{ii};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratVar))
        warning('Variable %s not found in base workspace. Skipping.', ratVar);
        continue;
    end
    rat = evalin('base', ratVar);

    % figure out last 3 days up to An
    dateList = autoDateList(rat);
    idx = find(strcmp(dateList, rat.An));
    if isempty(idx) || idx < 3
        warning('%s missing enough days before An. Skipping...', ratVar);
        continue;
    end
    daysToUse = dateList(idx-2:idx);

    % Pull fields needed to run the analysis
    spikes = filterFieldsByDay(rat.Ca_peaks, daysToUse);
    pos    = filterFieldsByDay(rat.pos,       daysToUse);
    ts     = filterFieldsByDay(rat.Ca_ts,     daysToUse);
    csus   = filterFieldsByDay(rat.(sprintf('csus%d', csusNum)), daysToUse);

    % We don't need actual MI to run the control here; pass an empty struct
    out = MI_control_rollingPost(spikes, pos, velthresh, dim, ts, csus, nIter, WEdges);

    fprintf('aggregating animal')
    % Aggregate this animal across days:
    [centersW, curveW, delta_allPoints, dayMeans] = aggregate_animal(out, rat, daysToUse);

    if isempty(curveW)
        warning('%s: no usable data (check MI_win/MI_rand fields). Skipping.', ratVar);
        continue;
    end

    % harmonize x-axis across animals
    if isempty(winCenters)
        winCenters = centersW(:);
    elseif numel(centersW) ~= numel(winCenters) || any(centersW(:) ~= winCenters(:))
        % interpolate onto shared axis
        curveW = interp1(centersW(:), curveW(:), winCenters(:), 'linear', 'extrap');
        % also align per-point deltas to shared axis if needed (we only need stats, already scalarized)
    end

    % store
    curves(:, end+1) = curveW(:); %#ok<AGROW>
    names{end+1}     = ratVar;    %#ok<AGROW>

    % collect for group stats
    allDelta   = [allDelta; delta_allPoints(:)]; %#ok<AGROW>
    allDayMeans = [allDayMeans; dayMeans(:)];    %#ok<AGROW>

    % -------- Per-animal stats (paired t-test MI_rand vs MI_win) --------
    % Note: test uses ALL points (windows×days×cells). This replicates each cell's MI_win
    % across windows; that's what you asked for.
    [p_pair,~,stats_pair] = ttest(delta_allPoints, 0, 'Alpha', alpha);  % equivalent to paired rand-vs-base
    muA  = mean(curveW,'omitnan');
    seA  = std(curveW,0,'omitnan')/sqrt(sum(~isnan(curveW)));

    % Also report per-window within-animal tests across cells (uncorr & FDR), useful to print
    [p_unc, q_fdr, nWin, nSig_unc, nSig_fdr] = per_animal_window_tests(out, rat, daysToUse, alpha);

    fprintf('%s: mean curve Δ = %.3f ± %.3f; paired all-points: t(%d)=%.3f, p=%.3g\n', ...
            ratVar, muA, seA, stats_pair, stats_pair, p_pair);
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


xlabel('Post-trial window center (s)');
ylabel('\Delta MI  (MI_{rand} - MI_{win})');
title(sprintf('Rolling-post control (csus%d): per-animal curves + mean±SEM', csusNum));
grid on;
legend([names, {'Mean±SEM'}], 'Interpreter','none', 'Location','best');

outname = sprintf('MI_rollingPost_curves_csus%d.png', csusNum);
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
% (A) ALL POINTS (windows×days×cells×animals)
x_all = allDelta(:); x_all = x_all(~isnan(x_all));
[~,p_all,~,stats_all] = ttest(x_all, 0, 'Alpha', alpha);
mu_all = mean(x_all,'omitnan');
se_all = std(x_all,0,'omitnan') / sqrt(numel(x_all));
fprintf('\n================ Overall group test — ALL POINTS =================\n');
fprintf('Mean Δ over all points: %g ± %g (SEM), N=%d\n', mu_all, se_all, numel(x_all));
fprintf('t(%d) = %.3f, p = %.3g\n', stats_all.df, stats_all.tstat, p_all);
fprintf('==================================================================\n');

% (B) DAY MEANS (animal×day means of Δ)
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


function [centers, curve, delta_allPoints, dayMeans] = aggregate_animal(out_struct, rat, daysToUse)
% Returns:
%  centers        : W×1 midpoints (from the first day encountered)
%  curve          : W×1 mean over days×cells of (MI_rand − MI_win)
%  delta_allPoints: vector of all pointwise Δ across windows×days×cells
%  dayMeans       : 1×Ndays mean Δ per day (mean over windows×cells)

centers = [];
curve   = [];
delta_allPoints = [];
dayMeans = [];

if isempty(out_struct) || ~isstruct(out_struct), return; end
fns = fieldnames(out_struct);
dayFns = fns(startsWith(fns,'MI_'));
if isempty(dayFns), return; end

Wref = [];     % number of windows (to align)
sumW = [];     % W×1 running sum for curve
cntW = [];     % W×1 running count for curve


for di = 1:numel(dayFns)
    day = out_struct.(dayFns{di});
    if ~isfield(day,'winEdges') || ~isfield(day,'MI_rand') || ~isfield(day,'MI_win')
        continue;
    end
    W = size(day.winEdges,1);
    if isempty(centers)
        centers = nanmean(day.winEdges,2);
        Wref = W;
        sumW = zeros(W,1);
        cntW = zeros(W,1);
    else
        if W ~= Wref
            % align this day's windows to the reference centers
            centers_d = nanmean(day.winEdges,2);
        else
            centers_d = centers;
        end
    end

    maskField = sprintf('ratemask_%s', dayFns{di}(4:end));
    maskVec = true(1, size(day.MI_rand,2));
    if isfield(rat,'ratemask') && isfield(rat.ratemask, maskField)
        mv = rat.ratemask.(maskField)(:);
        if numel(mv) == size(day.MI_rand,2)
            maskVec = (mv == 1);
        end
      end

    % Δ per window×cell for this day
    % Δ per window×cell for this day
    MI_rand = day.MI_rand(:, maskVec);   % W×Nc
    MI_win  = day.MI_win(:,  maskVec);   % W×Nc   <-- use ALL rows, not (1, :)

    if isempty(MI_rand) || isempty(MI_win), continue; end

    Delta  = MI_rand - MI_win;           % sign per your spec (control – removal)


    % update animal curve accumulators
    sumW = sumW + nanmean(Delta, 2);         % mean over cells for this day
    cntW = cntW + ~all(isnan(Delta),2);

    % collect ALL points (for all-points test A)
    delta_allPoints = [delta_allPoints; Delta(:)]; %#ok<AGROW>

    % day mean (for day-means test B)
    dayMeans = [dayMeans, mean(Delta,'all','omitnan')]; %#ok<AGROW>
end

if any(cntW>0)
    curve = sumW ./ max(1,cntW);             % mean over days of (mean over cells)
else
    curve = [];
end
end

function [p_unc, q_fdr, nW, nSig_unc, nSig_fdr] = per_animal_window_tests(out_struct, rat, daysToUse, alpha)
% For printing: within-animal per-window tests across cells (uncorr & FDR).
% Pools cells across the animal's analyzed days for each window; one-sample ttest vs 0.

p_unc = [];
q_fdr = [];
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

% build, for each window, the list of Δ across all (day×cells) for this animal
Delta_byW = cell(W,1);
for di = 1:numel(dayFns)
    day = out_struct.(dayFns{di});
    if ~isfield(day,'winEdges') || ~isfield(day,'MI_rand') || ~isfield(day,'MI_win'), continue; end
    centers_d = mean(day.winEdges,2);
    MIrand = day.MI_rand;
    MI_win = day.MI_win;
    if isempty(MIrand) || isempty(MI_win), continue; end

    % ratemask
    maskField = sprintf('ratemask_%s', dayFns{di}(4:end));
    maskVec = true(1, size(MIrand,2));
    if isfield(rat,'ratemask') && isfield(rat.ratemask, maskField)
        mv = rat.ratemask.(maskField)(:);
        if numel(mv) == size(MIrand,2), maskVec = (mv==1); end
    end
    MIrand = day.MI_rand(:,maskVec);   % W×Nc
    MI_win = day.MI_win(:,maskVec);    % W×Nc
    Delta  = MIrand - MI_win;          % W×Nc


    if numel(centers_d) ~= W || any(centers_d(:) ~= centers(:))
        % align rows to reference centers
        Delta = interp1(centers_d(:), Delta, centers(:), 'linear', 'extrap');
    end

    for w = 1:W
        Delta_byW{w} = [Delta_byW{w}, Delta(w,:)]; %#ok<AGROW>
    end
end

p_unc = nan(W,1);
for w = 1:W
    x = Delta_byW{w};
    x = x(~isnan(x));
    if numel(x) >= 2
        [~,p_unc(w)] = ttest(x, 0, 'Alpha', alpha);
    end
end
q_fdr = bh_fdr(p_unc);
nW = W;
nSig_unc = sum(p_unc < alpha, 'omitnan');
nSig_fdr = sum(q_fdr < alpha, 'omitnan');
end

function q = bh_fdr(p)
% Benjamini–Hochberg FDR; returns q-values shaped like p
q = nan(size(p));
ps = p(:);
valid = ~isnan(ps);
pv = ps(valid);
[sv, order] = sort(pv);
m = numel(sv);
if m==0, q = reshape(q, size(p)); return; end
qv = sv .* m ./ (1:m)';
for i = m-1:-1:1
    qv(i) = min(qv(i), qv(i+1));
end
q_full = nan(size(ps));
q_full(valid) = qv(invperm(order));
q = reshape(q_full, size(p));
end

function idx = invperm(order)
idx = zeros(size(order));
idx(order) = 1:numel(order);
end
