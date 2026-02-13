function OUT = run_task_space_stability_2x2(ratNames, varargin)
% RUN_TASK_SPACE_STABILITY_2X2
% Spatial stability analysis on a 2×2 grid, task vs non-task, robust to NaNs.
%
% Adds:
%   1) Flags to run only single-window and/or split-half metrics.
%   2) Plotting: ONLY 3 light day-lines per rat (first/last/middle or sampled),
%      plus 1 bold per-rat mean line.
%   3) Printed paired t-tests across rats (n=5) and across days (n~15 pooled),
%      done in z-space (Fisher z), which is the right space for averaging/testing.
%   4) Split-half uses multiple random splits (SplitHalfReps) -> yes, different
%      selections each rep via randperm(nW).

% ---------------------- Defaults (edit here if needed) -------------------
GridRC        = [2 2];
WinSec        = 2.0;
TraceWin      = [0 2];
BufferPost    = 0;
UseSpeedMask  = true;
VelThresh     = 4;        % cm/s
%CellNorm      = 'demean'; % 'none' | 'demean' | 'zscore' (applied to MAPS)  | 'meanrate'
CellNorm = 'meanrate';

MinWinFrames      = 8;
MinFramesPerBin   = 3;
SplitHalfReps     = 50;
RNGSeedBase       = 12345;

MinNCells_PV = 10;
MinNBins_SC  = 2;
MinStd       = 1e-12;

% NEW: run flags
RunSingle    = false;   % metrics 1a/1b
RunSplitHalf = true;   % metrics 2a/2b

% NEW: plotting controls
PlotMaxLightDays = 3;  % per rat
PlotLightDayMode = 'spread'; % 'spread' (first/mid/last) | 'random'

% NEW: stats controls
DoStats = true;
% ------------------------------------------------------------------------

% Parse varargin
p = inputParser;
addParameter(p,'RunSingle',RunSingle,@(x)islogical(x)&&isscalar(x));
addParameter(p,'RunSplitHalf',RunSplitHalf,@(x)islogical(x)&&isscalar(x));
addParameter(p,'SplitHalfReps',SplitHalfReps,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'PlotMaxLightDays',PlotMaxLightDays,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'PlotLightDayMode',PlotLightDayMode,@(s)ischar(s)||isstring(s));
addParameter(p,'DoStats',DoStats,@(x)islogical(x)&&isscalar(x));

parse(p,varargin{:});
RunSingle = p.Results.RunSingle;
RunSplitHalf = p.Results.RunSplitHalf;
SplitHalfReps = p.Results.SplitHalfReps;
PlotMaxLightDays = p.Results.PlotMaxLightDays;
PlotLightDayMode = char(p.Results.PlotLightDayMode);
DoStats = p.Results.DoStats;

ratNames = cellstr(ratNames);
nRats = numel(ratNames);

% Decide which metrics exist
% Always keep 4 slots for compatibility, but skip computation if flag off.
metricNames = {'m1a_PV_single','m1b_SC_single','m2a_PV_splithalf','m2b_SC_splithalf'};
useMetric = false(1,4);
if RunSingle
    useMetric(1) = true; % 1a
    useMetric(2) = true; % 1b
end
if RunSplitHalf
    useMetric(3) = true; % 2a
    useMetric(4) = true; % 2b
end

OUT = struct();
OUT.meta = struct('GridRC',GridRC,'WinSec',WinSec,'TraceWin',TraceWin,'BufferPost',BufferPost, ...
    'UseSpeedMask',UseSpeedMask,'VelThresh',VelThresh,'CellNorm',CellNorm, ...
    'MinWinFrames',MinWinFrames,'MinFramesPerBin',MinFramesPerBin,'SplitHalfReps',SplitHalfReps, ...
    'MinNCells_PV',MinNCells_PV,'MinNBins_SC',MinNBins_SC,'MinStd',MinStd, ...
    'RunSingle',RunSingle,'RunSplitHalf',RunSplitHalf, ...
    'PlotMaxLightDays',PlotMaxLightDays,'PlotLightDayMode',PlotLightDayMode, ...
    'DoStats',DoStats);

OUT.rats = struct( ...
    'animal', {}, 'days', {}, ...
    'z_task_byDay', {}, 'z_non_byDay', {}, ...
    'r_task_byDay', {}, 'r_non_byDay', {}, ...
    'z_task_mean', {}, 'z_non_mean', {}, ...
    'r_task_mean', {}, 'r_non_mean', {} );

% =============================== Main loop ===============================
for ii = 1:nRats
    ratVar = ratNames{ii};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratVar))
        warning('Variable %s not found in base workspace. Skipping.', ratVar);
        continue;
    end
    rat = evalin('base', ratVar);

    dateList = autoDateList_fallback(rat);
    if isempty(dateList)
        warning('[%s] No days found. Skipping.', ratVar);
        continue;
    end

    % ---- Select last 3 days up to An (inclusive) ----
    if isfield(rat,'An') && ~isempty(rat.An)
        idx = find(strcmp(dateList, rat.An), 1);
        if ~isempty(idx) && idx >= 3
            daysToUse = dateList(idx-2:idx);
        else
            daysToUse = dateList(max(1,numel(dateList)-2):end);
        end
    else
        daysToUse = dateList(max(1,numel(dateList)-2):end);
    end

    spikes_raw = filterFieldsByDay_fallback(rat.Ca_peaks, daysToUse);
    ts_raw     = filterFieldsByDay_fallback(rat.Ca_ts,    daysToUse);
    pos_raw    = filterFieldsByDay_fallback(rat.pos,      daysToUse);
    cs_raw     = filterFieldsByDay_fallback(rat.CS_times, daysToUse);

    [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, daysToUse);

    D = numel(ts);
    K = GridRC(1)*GridRC(2);

    % ---- Apply ratemask per day (keepCells where ratemask==1) ----
    for d = 1:D
        dlabel = daysToUse{d};
        keepCells = [];

        if isfield(rat,'ratemask')
            f = sprintf('ratemask_%s', dlabel);
            if isfield(rat.ratemask, f)
                rm = rat.ratemask.(f);
                keepCells = (rm == 1);
            end
        end

        if ~isempty(keepCells)
            % spikes{d} can be cell-of-cells or numeric matrix
            if iscell(spikes{d})
                % cell array: one entry per cell
                spikes{d} = spikes{d}(keepCells);
            elseif isnumeric(spikes{d})
                % numeric: rows are cells
                spikes{d} = spikes{d}(keepCells, :);
            else
                warning('[%s][%s] spikes day type not supported for ratemask filtering: %s', ...
                    ratVar, dlabel, class(spikes{d}));
            end
        end
    end

    z_task = nan(D, 4);
    z_non  = nan(D, 4);

    for d = 1:D
        t = ts{d}(:);
        if numel(t) < 10, continue; end

        edges = build_grid_edges_single_day(pos{d}, GridRC);
        [x, y] = interp_pos(pos{d}, t);

        S_counts = spikes_to_matrix(spikes{d}, t);

        if UseSpeedMask
            v = speed_cm_per_s(pos{d});
            v_i = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
            speed_ok = (v_i >= VelThresh);
        else
            speed_ok = true(size(t));
        end

        is_taskbuf = false(size(t));
        csd = cs{d}(:);
        if ~isempty(csd)
            for j = 1:numel(csd)
                t0 = csd(j) + TraceWin(1);
                t1 = csd(j) + TraceWin(1) + WinSec + BufferPost;
                is_taskbuf = is_taskbuf | (t >= t0 & t < t1);
            end
        end

        control_ok = ~is_taskbuf & speed_ok;

        taskWins = build_task_windows(csd, TraceWin(1), WinSec);

        [taskMaps, okTask] = window_maps_for_day( ...
            t, x, y, S_counts, taskWins, edges, ...
            speed_ok, UseSpeedMask, MinWinFrames, MinFramesPerBin, true(size(t)) );
        if ~okTask, continue; end
        nTaskW = size(taskMaps, 3);

        nonWins = sample_control_windows_strict(t, control_ok, WinSec, nTaskW, MinWinFrames, RNGSeedBase + 1000*ii + d);

        [nonMaps, okNon] = window_maps_for_day( ...
            t, x, y, S_counts, nonWins, edges, ...
            speed_ok, false, MinWinFrames, MinFramesPerBin, control_ok );
        if ~okNon, continue; end

        nNonW = size(nonMaps, 3);
        nW = min(nTaskW, nNonW);
        taskMaps = taskMaps(:,:,1:nW);
        nonMaps  = nonMaps(:,:,1:nW);

        [taskMaps, nonMaps] = normalize_maps_pair(taskMaps, nonMaps, CellNorm);

        % ---- metric 1a: PV single-window ----
        if useMetric(1)
            z_task(d,1) = metric_PV_single(taskMaps, K, MinNCells_PV, MinStd);
            z_non(d,1)  = metric_PV_single(nonMaps,  K, MinNCells_PV, MinStd);
        end

        % ---- metric 1b: single-cell single-window ----
        if useMetric(2)
            z_task(d,2) = metric_SC_single(taskMaps, K, MinNBins_SC, MinStd);
            z_non(d,2)  = metric_SC_single(nonMaps,  K, MinNBins_SC, MinStd);
        end

        % ---- metric 2a: PV split-half ----
        if useMetric(3)
            z_task(d,3) = metric_PV_splithalf(taskMaps, K, SplitHalfReps, RNGSeedBase + 2000*ii + d, MinNCells_PV, MinStd);
            z_non(d,3)  = metric_PV_splithalf(nonMaps,  K, SplitHalfReps, RNGSeedBase + 3000*ii + d, MinNCells_PV, MinStd);
        end

        % ---- metric 2b: single-cell split-half ----
        if useMetric(4)
            z_task(d,4) = metric_SC_splithalf(taskMaps, K, SplitHalfReps, RNGSeedBase + 4000*ii + d, MinNBins_SC, MinStd);
            z_non(d,4)  = metric_SC_splithalf(nonMaps,  K, SplitHalfReps, RNGSeedBase + 5000*ii + d, MinNBins_SC, MinStd);
        end
    end

    Rr = struct();
    Rr.animal = ratVar;
    Rr.days   = daysToUse;

    Rr.z_task_byDay = z_task;
    Rr.z_non_byDay  = z_non;

    Rr.r_task_byDay = tanh(z_task);
    Rr.r_non_byDay  = tanh(z_non);

    Rr.z_task_mean = mean(z_task, 1, 'omitnan');
    Rr.z_non_mean  = mean(z_non,  1, 'omitnan');
    Rr.r_task_mean = tanh(Rr.z_task_mean);
    Rr.r_non_mean  = tanh(Rr.z_non_mean);

    OUT.rats(end+1) = Rr; %#ok<AGROW>
end

% ================================ Plots =================================
plot_all_metrics(OUT, metricNames, useMetric, PlotMaxLightDays, PlotLightDayMode);

% ================================ Stats =================================
if DoStats
    print_stats(OUT, metricNames, useMetric);
end

end

% =========================================================================
%                                PLOTTING
% =========================================================================
function plot_all_metrics(OUT, metricNames, useMetric, maxLightDays, lightMode)
rats = OUT.rats;
if isempty(rats)
    warning('No rats to plot.');
    return;
end

mList = find(useMetric);
if isempty(mList)
    warning('No metrics selected to plot.');
    return;
end

nA = numel(rats);
cmap = lines(nA);
make_light = @(c,frac) (1-frac)*c + frac*[1 1 1];

figure('Color','w','Position',[120 120 1200 900]);

% choose subplot layout based on how many metrics are enabled
nM = numel(mList);
switch nM
    case 1, nR=1; nC=1;
    case 2, nR=1; nC=2;
    case 3, nR=2; nC=2;
    otherwise, nR=2; nC=2;
end

for idxM = 1:nM
    m = mList(idxM);
    subplot(nR,nC,idxM); hold on;

    % per-rat means (in r) computed via z-mean
    r_rat_task = nan(nA,1);
    r_rat_non  = nan(nA,1);

    % light day indices per rat
    lightIdx = cell(nA,1);

    for i = 1:nA
        zt = rats(i).z_task_byDay(:,m);
        zn = rats(i).z_non_byDay(:,m);

        % rat mean in z then tanh
        zt_mu = mean(zt, 'omitnan');
        zn_mu = mean(zn, 'omitnan');
        r_rat_task(i) = tanh(zt_mu);
        r_rat_non(i)  = tanh(zn_mu);

        % pick <= maxLightDays day indices with finite paired points
        mGood = isfinite(zt) & isfinite(zn);
        goodDays = find(mGood);
        if isempty(goodDays) || maxLightDays<=0
            lightIdx{i} = [];
        else
            if numel(goodDays) <= maxLightDays
                lightIdx{i} = goodDays(:);
            else
                switch lower(lightMode)
                    case 'random'
                        rng(9000+i+m);
                        sel = goodDays(randperm(numel(goodDays), maxLightDays));
                        lightIdx{i} = sort(sel(:));
                    otherwise % 'spread'
                        % spread across session: first/middle/last of goodDays
                        if maxLightDays == 1
                            lightIdx{i} = goodDays(round(end/2));
                        elseif maxLightDays == 2
                            lightIdx{i} = goodDays([1 end]);
                        else
                            mid = goodDays(round(numel(goodDays)/2));
                            lightIdx{i} = unique([goodDays(1); mid; goodDays(end)]);
                            % if still > maxLightDays (rare with duplicates), trim
                            if numel(lightIdx{i}) > maxLightDays
                                lightIdx{i} = lightIdx{i}(1:maxLightDays);
                            end
                        end
                end
            end
        end
    end

    keepA = isfinite(r_rat_task) & isfinite(r_rat_non);
    if ~any(keepA)
        title(metricNames{m}, 'Interpreter','none');
        axis off;
        continue;
    end

    % group bars: mean across animals in z then tanh (using rat-level means)
    z_group_task = mean(atanh_clip(r_rat_task(keepA)), 'omitnan');
    z_group_non  = mean(atanh_clip(r_rat_non(keepA)),  'omitnan');
    bar(1, tanh(z_group_task), 0.6, 'FaceColor',[0.35 0.70 1.00], 'EdgeColor','k');
    bar(2, tanh(z_group_non),  0.6, 'FaceColor',[0.85 0.55 0.25], 'EdgeColor','k');

    % 3 light day-lines per rat
    for i = 1:nA
        idxDays = lightIdx{i};
        if isempty(idxDays), continue; end
        c_light = make_light(cmap(i,:), 0.60);

        zt = rats(i).z_task_byDay(:,m);
        zn = rats(i).z_non_byDay(:,m);

        for dd = idxDays(:)'
            rt = tanh(zt(dd));
            rn = tanh(zn(dd));
            if isfinite(rt) && isfinite(rn)
                plot([1 2], [rt rn], '-', 'Color', c_light, 'LineWidth', 1.2);
                plot(1, rt, 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
                plot(2, rn, 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
            end
        end
    end

    % bold per-animal mean line
    for i = 1:nA
        if isfinite(r_rat_task(i)) && isfinite(r_rat_non(i))
            plot([1 2], [r_rat_task(i) r_rat_non(i)], '-', 'Color', cmap(i,:), 'LineWidth', 2.8);
            plot(1, r_rat_task(i), 'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 7);
            plot(2, r_rat_non(i),  'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 7);
        end
    end

    xlim([0.5 2.5]); xticks([1 2]);
    xticklabels({'task','non-task'});
    ylabel('Mean correlation (r)');
    title(metricNames{m}, 'Interpreter','none');
    yline(0,'k:');
    grid on; box on;
end
end

function z = atanh_clip(r)
r = max(min(r, 0.999999), -0.999999);
z = atanh(r);
end

% =========================================================================
%                                   STATS
% =========================================================================
function print_stats(OUT, metricNames, useMetric)
rats = OUT.rats;
if isempty(rats)
    fprintf('\n[run_task_space_stability_2x2] No rats; no stats.\n');
    return
end

mList = find(useMetric);
if isempty(mList)
    fprintf('\n[run_task_space_stability_2x2] No metrics selected; no stats.\n');
    return
end

fprintf('\n=== Task vs Non-task stability stats (paired t-tests; z-space) ===\n');

for m = mList
    % ----- rat-level paired test (n=rats) -----
    zT = nan(numel(rats),1);
    zN = nan(numel(rats),1);
    for i = 1:numel(rats)
        zt = rats(i).z_task_byDay(:,m);
        zn = rats(i).z_non_byDay(:,m);
        % rat mean in z (already)
        zT(i) = mean(zt,'omitnan');
        zN(i) = mean(zn,'omitnan');
    end
    keepR = isfinite(zT) & isfinite(zN);
    [pR, tR, dfR] = paired_t_from_vectors(zT(keepR), zN(keepR));

    % ----- day-level pooled paired test (n=days across all rats) -----
    zTd = []; zNd = [];
    for i = 1:numel(rats)
        zt = rats(i).z_task_byDay(:,m);
        zn = rats(i).z_non_byDay(:,m);
        keepD = isfinite(zt) & isfinite(zn);
        zTd = [zTd; zt(keepD)]; %#ok<AGROW>
        zNd = [zNd; zn(keepD)]; %#ok<AGROW>
    end
    [pD, tD, dfD] = paired_t_from_vectors(zTd, zNd);

    % Report means in BOTH z and r for readability
    mu_zT = mean(zT(keepR),'omitnan'); mu_zN = mean(zN(keepR),'omitnan');
    mu_rT = tanh(mu_zT);              mu_rN = tanh(mu_zN);

    fprintf('\n[%s]\n', metricNames{m});
    fprintf('  Rats  (n=%d): mean r task=%.3f  non=%.3f | mean z task=%.3f non=%.3f | paired t(%d)=%.3f  p=%.3g\n', ...
        nnz(keepR), mu_rT, mu_rN, mu_zT, mu_zN, dfR, tR, pR);
    fprintf('  Days  (n=%d): paired t(%d)=%.3f  p=%.3g   (pooled day pairs; z-space)\n', ...
        numel(zTd), dfD, tD, pD);
end
fprintf('\n');
end

function [p, tstat, df] = paired_t_from_vectors(a, b)
a = a(:); b = b(:);
m = isfinite(a) & isfinite(b);
a = a(m); b = b(m);
if numel(a) < 2
    p = NaN; tstat = NaN; df = NaN; return
end
d = b - a; % positive => non-task > task
n = numel(d);
mu = mean(d,'omitnan');
sd = std(d,0,'omitnan');
if ~isfinite(sd) || sd==0
    p = NaN; tstat = NaN; df = n-1; return
end
tstat = mu / (sd/sqrt(n));
df = n-1;
p = 2 * (1 - tcdf(abs(tstat), df));
end

% =========================================================================
%                             METRICS (z space)
% =========================================================================
function zmu = metric_PV_single(maps, K, minNcells, minStd)
[nc, ~, nW] = size(maps);
if nc < minNcells || nW < 2
    zmu = NaN; return;
end
z_list = [];
for k = 1:K
    X = squeeze(maps(:,k,:)); % Nc×nW
    if isempty(X), continue; end
    for i = 1:(nW-1)
        v1 = X(:,i);
        for j = (i+1):nW
            v2 = X(:,j);
            r = safe_corr_pairwise(v1, v2, minNcells, minStd);
            if isfinite(r), z_list(end+1,1) = atanh_clip(r); end %#ok<AGROW>
        end
    end
end
if isempty(z_list), zmu = NaN; else, zmu = mean(z_list,'omitnan'); end
end

function zmu = metric_SC_single(maps, K, minNbins, minStd)
[nc, ~, nW] = size(maps);
if K < minNbins || nc < 2 || nW < 2
    zmu = NaN; return;
end
z_cell = nan(nc,1);
for c = 1:nc
    X = squeeze(maps(c,:,:)); % K×nW
    z_list = [];
    for i = 1:(nW-1)
        v1 = X(:,i);
        for j = (i+1):nW
            v2 = X(:,j);
            r = safe_corr_pairwise(v1, v2, minNbins, minStd);
            if isfinite(r), z_list(end+1,1) = atanh_clip(r); end %#ok<AGROW>
        end
    end
    if ~isempty(z_list), z_cell(c) = mean(z_list,'omitnan'); end
end
if all(~isfinite(z_cell)), zmu = NaN; else, zmu = mean(z_cell,'omitnan'); end
end

function zmu = metric_PV_splithalf(maps, K, nRep, seed, minNcells, minStd)
[nc, ~, nW] = size(maps);
if nc < minNcells || nW < 4
    zmu = NaN; return;
end
rng(seed);

z_rep = nan(nRep,1);
for rr = 1:nRep
    idx = randperm(nW);               % <-- different random split each rep
    nA = floor(nW/2);
    A = idx(1:nA);
    B = idx(nA+1:end);
    if numel(A) < 2 || numel(B) < 2, continue; end

    z_k = nan(K,1);
    for k = 1:K
        XA = squeeze(maps(:,k,A));
        XB = squeeze(maps(:,k,B));
        vA = mean(XA, 2, 'omitnan');
        vB = mean(XB, 2, 'omitnan');
        r = safe_corr_pairwise(vA, vB, minNcells, minStd);
        if isfinite(r), z_k(k) = atanh_clip(r); end
    end
    if any(isfinite(z_k)), z_rep(rr) = mean(z_k,'omitnan'); end
end
if all(~isfinite(z_rep)), zmu = NaN; else, zmu = mean(z_rep,'omitnan'); end
end

function zmu = metric_SC_splithalf(maps, K, nRep, seed, minNbins, minStd)
[nc, ~, nW] = size(maps);
if K < minNbins || nc < 2 || nW < 4
    zmu = NaN; return;
end
rng(seed);

z_rep = nan(nRep,1);
for rr = 1:nRep
    idx = randperm(nW);               % <-- different random split each rep
    nA = floor(nW/2);
    A = idx(1:nA);
    B = idx(nA+1:end);
    if numel(A) < 2 || numel(B) < 2, continue; end

    z_cell = nan(nc,1);
    for c = 1:nc
        vA = mean(squeeze(maps(c,:,A)), 2, 'omitnan');
        vB = mean(squeeze(maps(c,:,B)), 2, 'omitnan');
        r  = safe_corr_pairwise(vA, vB, minNbins, minStd);
        if isfinite(r), z_cell(c) = atanh_clip(r); end
    end
    if any(isfinite(z_cell)), z_rep(rr) = mean(z_cell,'omitnan'); end
end
if all(~isfinite(z_rep)), zmu = NaN; else, zmu = mean(z_rep,'omitnan'); end
end

% =========================================================================
%                         WINDOW BUILDING + MAPS
% =========================================================================
function wins = build_task_windows(csd, tOffset, WinSec)
if isempty(csd)
    wins = zeros(0,2);
    return
end
csd = csd(:);
wins = [csd + tOffset, csd + tOffset + WinSec];
end

function wins = sample_control_windows_strict(t, okFrame, WinSec, nNeed, MinWinFrames, seed)
rng(seed);
t = t(:);
if nNeed < 1
    wins = zeros(0,2);
    return
end
okIdx = find(okFrame);
if numel(okIdx) < MinWinFrames
    wins = zeros(0,2);
    return
end
tEnd = t(end);

wins_out = nan(nNeed,2);
nFound = 0;
maxAttempts = max(2000, 50*nNeed);

for a = 1:maxAttempts
    s = okIdx(randi(numel(okIdx), 1, 1));
    t0 = t(s);
    t1 = t0 + WinSec;
    if t1 > tEnd, continue; end

    idx = find(t >= t0 & t < t1 & okFrame);
    if numel(idx) < MinWinFrames, continue; end

    nFound = nFound + 1;
    wins_out(nFound,:) = [t0 t1];
    if nFound >= nNeed, break; end
end

if nFound < 1
    wins = zeros(0,2);
else
    wins = wins_out(1:nFound,:);
    [~, ord] = sort(wins(:,1));
    wins = wins(ord,:);
end
end

function [maps, ok] = window_maps_for_day(t, x, y, S_counts, wins, edges, speed_ok, ...
    applySpeedMask, MinWinFrames, MinFramesPerBin, frame_ok)
ok = false;
maps = [];

if nargin < 11 || isempty(frame_ok)
    frame_ok = true(size(t));
end
if isempty(wins), return; end

t = t(:);
Nc = size(S_counts,1);
K  = (numel(edges.y)-1)*(numel(edges.x)-1);

dt = median(diff(t));
if ~isfinite(dt) || dt <= 0
    dt = max(eps, mean(diff(t), 'omitnan'));
end

goodMaps = cell(size(wins,1),1);
keepW = false(size(wins,1),1);

for w = 1:size(wins,1)
    t0 = wins(w,1); t1 = wins(w,2);

    idx = find(t >= t0 & t < t1 & frame_ok); % enforce allowed frames
    if isempty(idx), continue; end

    if applySpeedMask
        idx = idx(speed_ok(idx));
    end
    if numel(idx) < MinWinFrames, continue; end

    [~, kbin] = pos2bin(x(idx), y(idx), edges);

    M = nan(Nc, K);
    for k = 1:K
        hit = find(kbin == k);
        if numel(hit) < MinFramesPerBin, continue; end
        fr = sum(S_counts(:, idx(hit)), 2, 'omitnan') ./ (numel(hit)*dt);
        M(:,k) = fr;
    end

    if any(isfinite(M(:)))
        goodMaps{w} = M;
        keepW(w) = true;
    end
end

if ~any(keepW), return; end

Ms = goodMaps(keepW);
nW = numel(Ms);
maps = nan(Nc, K, nW);
for i = 1:nW
    maps(:,:,i) = Ms{i};
end
ok = true;
end

% =========================================================================
%                          MAP NORMALIZATION (Hz-space)
% =========================================================================
function [A, B] = normalize_maps_pair(A, B, mode)
mode = lower(string(mode));
if mode == "none"
    return
end

[nc, ~, ~] = size(A);
Ap = reshape(A, nc, []);
Bp = reshape(B, nc, []);
P  = [Ap, Bp];

mu = mean(P, 2, 'omitnan');

switch mode
    case "demean"
        A = A - reshape(mu, [nc 1 1]);
        B = B - reshape(mu, [nc 1 1]);

    case "zscore"
        sd = std(P, 0, 2, 'omitnan');
        sd(~isfinite(sd) | sd==0) = 1;
        A = (A - reshape(mu, [nc 1 1])) ./ reshape(sd, [nc 1 1]);
        B = (B - reshape(mu, [nc 1 1])) ./ reshape(sd, [nc 1 1]);
    case "meanrate"
      A = A ./ reshape(mu, [nc 1 1]);
      B = B ./ reshape(mu, [nc 1 1]);


    otherwise
        error('CellNorm must be ''none'',''demean'',''zscore''.');
end
end

% =========================================================================
%                               CORRELATION
% =========================================================================
function r = safe_corr_pairwise(a, b, minN, minStd)
a = a(:); b = b(:);
n = min(numel(a), numel(b));
a = a(1:n); b = b(1:n);

m = isfinite(a) & isfinite(b);
if nnz(m) < minN
    r = NaN; return;
end

aa = a(m); bb = b(m);
sa = std(aa); sb = std(bb);
if ~isfinite(sa) || ~isfinite(sb) || sa < minStd || sb < minStd
    r = NaN; return;
end

r = corr(aa, bb, 'type','Pearson');
end

% =========================================================================
%                               UTILITIES
% =========================================================================
function [edges] = build_grid_edges_single_day(posd, GridRC)
allx = posd.x(:); ally = posd.y(:);
allx = allx(isfinite(allx)); ally = ally(isfinite(ally));
if isempty(allx) || isempty(ally)
    edges.x = linspace(0, 1, GridRC(2)+1);
    edges.y = linspace(0, 1, GridRC(1)+1);
else
    edges.x = linspace(min(allx), max(allx), GridRC(2)+1);
    edges.y = linspace(min(ally), max(ally), GridRC(1)+1);
end
end

function [x_i, y_i] = interp_pos(posd, t)
tt = double(posd.t(:)); xx = double(posd.x(:)); yy = double(posd.y(:));
t  = double(t(:));
[ttu, ia] = unique(tt, 'stable');
xxu = xx(ia); yyu = yy(ia);
x_i = interp1(ttu, xxu, t, 'linear','extrap');
y_i = interp1(ttu, yyu, t, 'linear','extrap');
end

function v = speed_cm_per_s(posd)
t = double(posd.t(:)); x = double(posd.x(:)); y = double(posd.y(:));
n = min([numel(t), numel(x), numel(y)]);
t = t(1:n); x = x(1:n); y = y(1:n);
dt = diff(t); dt(end+1,1) = median(dt(dt>0),'omitnan');
dx = [diff(x); 0]; dy = [diff(y); 0];
v = hypot(dx,dy) ./ max(dt, eps);
end

function [rc_idx, k] = pos2bin(x, y, edges)
cx = discretize(x, edges.x);
cy = discretize(y, edges.y);
GridR = numel(edges.y)-1;
GridC = numel(edges.x)-1;
bad = isnan(cx) | isnan(cy) | cx<1 | cy<1 | cx>GridC | cy>GridR;
cx(bad) = NaN; cy(bad) = NaN;
rc_idx = [cy, cx];
k = nan(size(x));
m = isfinite(cx) & isfinite(cy);
if any(m)
    k(m) = sub2ind([GridR, GridC], cy(m), cx(m));
end
end

function S = spikes_to_matrix(daySpikes, t)
daySpikes = unwrap_spike_container(daySpikes);
if iscell(daySpikes)
    Nc = numel(daySpikes);
elseif isnumeric(daySpikes)
    Nc = size(daySpikes,1);
else
    error('Unsupported daySpikes type: %s', class(daySpikes));
end
t = double(t(:));
if numel(t) < 2
    S = zeros(Nc, numel(t), 'single');
    return;
end
dt = median(diff(t));
if ~isfinite(dt) || dt <= 0
    dt = max(eps, mean(diff(t),'omitnan'));
end
edges = [t - dt/2; t(end)+dt/2];
S = zeros(Nc, numel(t), 'single');
for c = 1:Nc
    st = extract_cell_spikes(daySpikes, c);
    if isempty(st), continue; end
    S(c,:) = histcounts(st, edges);
end
end

function obj = unwrap_spike_container(S)
if isstruct(S)
    f = fieldnames(S);
    if isempty(f), obj = []; return; end
    pick = [];
    for j = 1:numel(f)
        name = lower(f{j});
        if contains(name,'peak') || contains(name,'spike') || contains(name,'ca_peaks')
            pick = j; break;
        end
    end
    if isempty(pick), pick = 1; end
    obj = S.(f{pick});
else
    obj = S;
end
end

function st = extract_cell_spikes(container, c)
if iscell(container)
    st = container{c};
elseif isnumeric(container)
    if c > size(container,1), st = []; return; end
    st = container(c,:).';
else
    st = [];
end
st = double(st(:));
st = st(isfinite(st) & st > 0);
end

% =========================================================================
%                    DAY SELECTION + STRUCT STANDARDIZATION
% =========================================================================
function dateList = autoDateList_fallback(rat)
if exist('autoDateList','file') == 2
    try
        dateList = autoDateList(rat);
        return
    catch
    end
end
dateList = {};
if isfield(rat,'Ca_peaks') && isstruct(rat.Ca_peaks)
    dateList = fieldnames(rat.Ca_peaks);
elseif isfield(rat,'Ca_ts') && isstruct(rat.Ca_ts)
    dateList = fieldnames(rat.Ca_ts);
end
dateList = dateList(:);
end

function S = filterFieldsByDay_fallback(Sin, ~)
if exist('filterFieldsByDay','file') == 2
    try
        S = Sin;
        return
    catch
    end
end
S = Sin;
end

function [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, dayKeys)
spikes = to_daycells(spikes_raw, dayKeys);
ts_in  = to_daycells(ts_raw,    dayKeys);
pos_in = to_daycells(pos_raw,   dayKeys);
cs_in  = to_daycells(cs_raw,    dayKeys);

D = numel(dayKeys);
ts  = cell(1,D);
pos = cell(1,D);
cs  = cell(1,D);

for d = 1:D
    ts{d}  = coerce_ts_day(ts_in{d});
    pos{d} = coerce_pos_day(pos_in{d}, ts{d});
    cs{d}  = coerce_cs_day(cs_in{d});
end
end

function C = to_daycells(X, keys)
if iscell(X)
    C = X(1:min(numel(X), numel(keys)));
    if numel(C) < numel(keys), C(end+1:numel(keys)) = {[]}; end
    return
end
if isstruct(X)
    if numel(X) == 1
        fn = fieldnames(X);
        C = cell(1, numel(keys));
        for i = 1:numel(keys)
            k = keys{i};
            idx = find(strcmp(fn, k), 1);
            if isempty(idx), idx = find(contains(fn, k), 1); end
            if isempty(idx), C{i} = []; else, C{i} = X.(fn{idx}); end
        end
        return
    else
        C = arrayfun(@(j) X(j), 1:min(numel(X), numel(keys)), 'uni', 0);
        if numel(C) < numel(keys), C(end+1:numel(keys)) = {[]}; end
        return
    end
end
C = repmat({X}, 1, numel(keys));
end

function t = coerce_ts_day(td)
if iscell(td), if isempty(td), t=[]; return; else, td=td{1}; end, end
if isempty(td), t=[]; return; end

if isnumeric(td) && isvector(td)
    t = double(td(:));
    if ~isempty(t) && max(t,[],'omitnan') > 1e4, t = t/1000; end
    return
end

if isnumeric(td) && ismatrix(td) && ~isscalar(td)
    M = double(td);
    if isempty(M) || size(M,2) < 1
        t = []; return
    end
    col = min(2,size(M,2));
    if col < 1, t = []; return; end
    t = M(:, col);
    if ~isempty(t) && max(t,[],'omitnan') > 1e4, t = t/1000; end
    t = t(:);
    return
end

if istable(td)
    vn = lower(string(td.Properties.VariableNames));
    cname = vn(find(ismember(vn, ["time_ms","time","ts","t","timestamp","frame_ts","ca_ts"]),1));
    if isempty(cname), error('ts table: no time-like column found'); end
    t = double(td{:,find(vn==cname,1)});
    if contains(cname,"ms") || max(t,[],'omitnan') > 1e4, t = t/1000; end
    t = t(:);
    return
end

if isstruct(td)
    pref = ["time_ms","time","ts","t","timestamp","timestamps","frame_ts","ca_ts"];
    for j = 1:numel(pref)
        if isfield(td, pref{j})
            v = td.(pref{j});
            if isnumeric(v) && isvector(v) && numel(v) > 1
                t = double(v(:));
                if contains(pref(j),"ms") || max(t,[],'omitnan') > 1e4, t = t/1000; end
                return
            end
        end
    end
    error('ts struct: no numeric time vector found.');
end

error('Unsupported ts type: %s', class(td));
end

function P = coerce_pos_day(pd, tday)
col = @(v) v(:);
if iscell(pd), if isempty(pd), P=struct('t',[],'x',[],'y',[]); return; else, pd=pd{1}; end, end
if isempty(pd)
    if ~isempty(tday)
        P = struct('t',col(tday),'x',nan(numel(tday),1),'y',nan(numel(tday),1));
    else
        P = struct('t',[],'x',[],'y',[]);
    end
    return
end
if istable(pd)
    vn = lower(string(pd.Properties.VariableNames));
    tname = pick_name(vn, ["t","time","ts"]);
    xname = pick_name(vn, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"]);
    yname = pick_name(vn, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"]);
    t = []; if strlength(tname)>0, t = pd{:, find(vn==tname,1)}; end
    x = []; if strlength(xname)>0, x = pd{:, find(vn==xname,1)}; end
    y = []; if strlength(yname)>0, y = pd{:, find(vn==yname,1)}; end
    [t,x,y] = finalize_txy(t,x,y,tday);
    P = struct('t',col(t),'x',col(x),'y',col(y));
    return
end
if isnumeric(pd) && ismatrix(pd) && ~isscalar(pd)
    [t,x,y] = coerce_from_numeric(pd, tday);
    P = struct('t',col(t),'x',col(x),'y',col(y));
    return
end
if isstruct(pd)
    f = lower(string(fieldnames(pd)));
    t = get_field_if_exists(pd, pick_name(f, ["t","time","ts"]));
    x = get_field_if_exists(pd, pick_name(f, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"]));
    y = get_field_if_exists(pd, pick_name(f, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"]));
    [t,x,y] = finalize_txy(t,x,y,tday);
    P = struct('t',col(t),'x',col(x),'y',col(y));
    return
end
error('pos day: unsupported type %s', class(pd));
end

function name = pick_name(names, options)
name = "";
for k = 1:numel(options)
    idx = find(names==options(k),1);
    if ~isempty(idx), name = names(idx); return; end
end
end

function val = get_field_if_exists(S, name)
if strlength(name)==0 || ~isstruct(S), val=[]; return; end
fn = fieldnames(S);
idx = find(strcmpi(fn, char(name)),1);
if isempty(idx), val=[]; else, val = S.(fn{idx}); end
end

function [t,x,y] = coerce_from_numeric(M, tday)
[nr,nc] = size(M);
if nc==3
    t = M(:,1); x = M(:,2); y = M(:,3);
elseif nc==2
    x = M(:,1); y = M(:,2);
    if ~isempty(tday) && numel(tday)==nr, t = tday; else, t = (1:nr)'; end
else
    error('numeric matrix must be [n×3] or [n×2].');
end
[t,x,y] = finalize_txy(t,x,y,tday);
end

function [t,x,y] = finalize_txy(t,x,y,tday)
x = x(:); y = y(:);
if isempty(t)
    if ~isempty(tday) && numel(tday)==numel(x), t = tday(:); else, t = (1:numel(x))'; end
else
    t = t(:);
end
n = min([numel(t), numel(x), numel(y)]);
t = double(t(1:n)); x = double(x(1:n)); y = double(y(1:n));
end

function cs_vec = coerce_cs_day(csd)
if iscell(csd)
    if isempty(csd), cs_vec = []; return; end
    csd = csd{1};
end
if isempty(csd), cs_vec = []; return; end

if isnumeric(csd) && isvector(csd)
    cs = double(csd(:));
    cs = cs(isfinite(cs));
    if ~isempty(cs) && max(cs,[],'omitnan') > 1e4, cs = cs/1000; end
    cs_vec = cs(:);
    return
end

if istable(csd)
    vn = lower(string(csd.Properties.VariableNames));
    cname = pick_name(vn, ["cs_ms","cs_time","cstime","cs","onset","onsets","cs_onset", ...
                           "cs_time_ms","trial_cs","cue_onset","time","ts"]);
    if strlength(cname)==0, cs_vec = []; return; end
    tcol = double(csd{:, find(vn==cname,1)});
    tcol = tcol(isfinite(tcol));
    if ~isempty(tcol) && (contains(cname,"ms") || max(tcol,[],'omitnan') > 1e4), tcol = tcol/1000; end
    cs_vec = tcol(:);
    return
end

if isnumeric(csd) && ismatrix(csd) && ~isscalar(csd)
    M = double(csd);
    if isempty(M) || size(M,2) < 1 || size(M,1) < 1, cs_vec = []; return; end
    if size(M,2) == 1, cs = M(:,1); else, cs = M(:,2); end
    cs = cs(:);
    cs = cs(isfinite(cs));
    if ~isempty(cs) && max(cs,[],'omitnan') > 1e4, cs = cs/1000; end
    cs_vec = cs(:);
    return
end

cs_vec = [];
end
