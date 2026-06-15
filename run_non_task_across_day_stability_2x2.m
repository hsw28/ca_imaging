function OUT = run_non_task_across_day_stability_2x2(ratNames, varargin)
% RUN_NON_TASK_DAY_STABILITY_2X2
% Non-task spatial stability analysis on a 2x2 grid, robust to NaNs.
%
% Compares non-task activity maps across days:
%   1) An vs An-1
%   2) An vs An+2
%
% Metrics are analogous to run_task_space_stability_2x2, but the two
% conditions are day-pairs rather than task vs non-task within a day.
% Correlations are averaged/tested in Fisher z-space, then displayed as r.

% ---------------------- Defaults (edit here if needed) -------------------
GridRC        = [4 4];
WinSec        = 1;
TraceWin      = [0 10/7.5];
BufferPost    = 0;
UseSpeedMask  = true;
VelThresh     = 4;        % cm/s
CellNorm      = 'none'; % 'none' | 'demean' | 'zscore' | 'meanrate'

ControlWinsPerDay = 50;
MinWinFrames      = 8;
MinFramesPerBin   = 3;
SplitHalfReps     = 50;
RNGSeedBase       = 12345;

MinNCells_PV = 10;
MinNBins_SC  = 2;
MinStd       = 1e-12;

RunSingle    = true;   % metrics 1a/1b
RunSplitHalf = false;    % metrics 2a/2b

PlotMaxLightRats = inf;
DoStats = true;
% ------------------------------------------------------------------------

p = inputParser;
addParameter(p,'RunSingle',RunSingle,@(x)islogical(x)&&isscalar(x));
addParameter(p,'RunSplitHalf',RunSplitHalf,@(x)islogical(x)&&isscalar(x));
addParameter(p,'SplitHalfReps',SplitHalfReps,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'ControlWinsPerDay',ControlWinsPerDay,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'PlotMaxLightRats',PlotMaxLightRats,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'DoStats',DoStats,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});

RunSingle = p.Results.RunSingle;
RunSplitHalf = p.Results.RunSplitHalf;
SplitHalfReps = p.Results.SplitHalfReps;
ControlWinsPerDay = p.Results.ControlWinsPerDay;
PlotMaxLightRats = p.Results.PlotMaxLightRats;
DoStats = p.Results.DoStats;

ratNames = cellstr(ratNames);
nRats = numel(ratNames);

metricNames = {'m1a_PV_crossday_single','m1b_SC_crossday_single', ...
    'm2a_PV_crossday_splithalf','m2b_SC_crossday_splithalf'};
useMetric = false(1,4);
if RunSingle
    useMetric(1:2) = true;
end
if RunSplitHalf
    useMetric(3:4) = true;
end

pairNames = {'An_vs_AnMinus1','An_vs_AnPlus2'};

OUT = struct();
OUT.meta = struct('GridRC',GridRC,'WinSec',WinSec,'TraceWin',TraceWin,'BufferPost',BufferPost, ...
    'UseSpeedMask',UseSpeedMask,'VelThresh',VelThresh,'CellNorm',CellNorm, ...
    'ControlWinsPerDay',ControlWinsPerDay,'MinWinFrames',MinWinFrames, ...
    'MinFramesPerBin',MinFramesPerBin,'SplitHalfReps',SplitHalfReps, ...
    'MinNCells_PV',MinNCells_PV,'MinNBins_SC',MinNBins_SC,'MinStd',MinStd, ...
    'RunSingle',RunSingle,'RunSplitHalf',RunSplitHalf,'DoStats',DoStats);

OUT.pairNames = pairNames;
OUT.metricNames = metricNames;
OUT.rats = struct( ...
    'animal', {}, 'An', {}, 'pairs', {}, ...
    'z_byPair', {}, 'r_byPair', {}, ...
    'z_pair_mean', {}, 'r_pair_mean', {} );

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

    if isfield(rat,'An') && ~isempty(rat.An)
        idxAn = find(strcmp(dateList, rat.An), 1);
    else
        idxAn = [];
    end
    if isempty(idxAn)
        warning('[%s] Could not find rat.An in dateList. Skipping.', ratVar);
        continue;
    end

    pairDayIdx = { [idxAn, idxAn-1], [idxAn, idxAn+2] };
    pairDays = cell(1,2);
    validPair = false(1,2);
    for pidx = 1:2
        ids = pairDayIdx{pidx};
        if all(ids >= 1) && all(ids <= numel(dateList))
            pairDays{pidx} = dateList(ids);
            validPair(pidx) = true;
        else
            pairDays{pidx} = {};
        end
    end

    daysToUse = unique([pairDays{validPair}], 'stable');
    if isempty(daysToUse)
        warning('[%s] No requested day pairs exist around An=%s.', ratVar, rat.An);
        continue;
    end

    spikes_raw = filterFieldsByDay_fallback(rat.Ca_peaks, daysToUse);
    ts_raw     = filterFieldsByDay_fallback(rat.Ca_ts,    daysToUse);
    pos_raw    = filterFieldsByDay_fallback(rat.pos,      daysToUse);
    cs_raw     = filterFieldsByDay_fallback(rat.CS_times, daysToUse);
    [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, daysToUse);

    K = GridRC(1)*GridRC(2);
    z_pair = nan(2, 4);

    for pidx = 1:2
        if ~validPair(pidx), continue; end

        dAnLabel = pairDays{pidx}{1};
        dOtherLabel = pairDays{pidx}{2};
        iAn = find(strcmp(daysToUse, dAnLabel), 1);
        iOther = find(strcmp(daysToUse, dOtherLabel), 1);
        if isempty(iAn) || isempty(iOther), continue; end

        keepCells = pair_keep_cells(rat, dAnLabel, dOtherLabel, spikes{iAn}, spikes{iOther});
        spAn = filter_spikes_cells(spikes{iAn}, keepCells);
        spOther = filter_spikes_cells(spikes{iOther}, keepCells);

        mapsAn = non_task_maps_for_day(ts{iAn}, pos{iAn}, cs{iAn}, spAn, GridRC, WinSec, ...
            TraceWin, BufferPost, UseSpeedMask, VelThresh, ControlWinsPerDay, ...
            MinWinFrames, MinFramesPerBin, RNGSeedBase + 1000*ii + 10*pidx + 1);

        mapsOther = non_task_maps_for_day(ts{iOther}, pos{iOther}, cs{iOther}, spOther, GridRC, WinSec, ...
            TraceWin, BufferPost, UseSpeedMask, VelThresh, ControlWinsPerDay, ...
            MinWinFrames, MinFramesPerBin, RNGSeedBase + 1000*ii + 10*pidx + 2);

        if isempty(mapsAn) || isempty(mapsOther), continue; end
        [mapsAn, mapsOther] = normalize_maps_pair(mapsAn, mapsOther, CellNorm);

        if useMetric(1)
            z_pair(pidx,1) = metric_PV_cross_single(mapsAn, mapsOther, K, MinNCells_PV, MinStd);
        end
        if useMetric(2)
            z_pair(pidx,2) = metric_SC_cross_single(mapsAn, mapsOther, K, MinNBins_SC, MinStd);
        end
        if useMetric(3)
            z_pair(pidx,3) = metric_PV_cross_splithalf(mapsAn, mapsOther, K, SplitHalfReps, ...
                RNGSeedBase + 2000*ii + pidx, MinNCells_PV, MinStd);
        end
        if useMetric(4)
            z_pair(pidx,4) = metric_SC_cross_splithalf(mapsAn, mapsOther, K, SplitHalfReps, ...
                RNGSeedBase + 3000*ii + pidx, MinNBins_SC, MinStd);
        end
    end

    Rr = struct();
    Rr.animal = ratVar;
    Rr.An = rat.An;
    Rr.pairs = pairDays;
    Rr.z_byPair = z_pair;
    Rr.r_byPair = tanh(z_pair);
    Rr.z_pair_mean = mean(z_pair, 1, 'omitnan');
    Rr.r_pair_mean = tanh(Rr.z_pair_mean);

    OUT.rats(end+1) = Rr; %#ok<AGROW>
end

plot_all_metrics(OUT, metricNames, useMetric, PlotMaxLightRats);

if DoStats
    print_stats(OUT, metricNames, useMetric);
end

for i = 1:numel(OUT.rats)
    fprintf('%s  A-A: %.3f   A-B: %.3f\n', ...
        OUT.rats(i).animal, ...
        tanh(OUT.rats(i).z_byPair(1,2)), ...
        tanh(OUT.rats(i).z_byPair(2,2)));
end

end

% =========================================================================
%                                PLOTTING
% =========================================================================
function plot_all_metrics(OUT, metricNames, useMetric, maxLightRats)
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
nM = numel(mList);
switch nM
    case 1, nR=1; nC=1;
    case 2, nR=1; nC=2;
    otherwise, nR=2; nC=2;
end

for idxM = 1:nM
    m = mList(idxM);
    subplot(nR,nC,idxM); hold on;

    zA = nan(nA,1);
    zB = nan(nA,1);
    for i = 1:nA
        zA(i) = rats(i).z_byPair(1,m); % An vs An-1
        zB(i) = rats(i).z_byPair(2,m); % An vs An+2
    end
    keepA = isfinite(zA);
    keepB = isfinite(zB);

    if any(keepA)
        bar(1, tanh(mean(zA(keepA),'omitnan')), 0.6, 'FaceColor',[0.35 0.70 1.00], 'EdgeColor','k');
    end
    if any(keepB)
        bar(2, tanh(mean(zB(keepB),'omitnan')), 0.6, 'FaceColor',[0.85 0.55 0.25], 'EdgeColor','k');
    end

    nLight = min(nA, floor(maxLightRats));
    for i = 1:nLight
        c = cmap(i,:);
        cLight = make_light(c, 0.50);
        if isfinite(zA(i))
            plot(1, tanh(zA(i)), 'o', 'MarkerFaceColor', cLight, 'MarkerEdgeColor', cLight, 'MarkerSize', 5);
        end
        if isfinite(zB(i))
            plot(2, tanh(zB(i)), 'o', 'MarkerFaceColor', cLight, 'MarkerEdgeColor', cLight, 'MarkerSize', 5);
        end
        if isfinite(zA(i)) && isfinite(zB(i))
            plot([1 2], tanh([zA(i) zB(i)]), '-', 'Color', c, 'LineWidth', 2.5);
            plot([1 2], tanh([zA(i) zB(i)]), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor','k', ...
                'LineWidth',0.5, 'MarkerSize', 7);
        end
    end

    xlim([0.5 2.5]); xticks([1 2]);
    xticklabels({'An vs An-1','An vs An+2'});
    ylabel('Cross-day non-task stability (r)');
    title(metricNames{m}, 'Interpreter','none');
    yline(0,'k:');
    grid on; box on;
end
end

% =========================================================================
%                                   STATS
% =========================================================================
function print_stats(OUT, metricNames, useMetric)
rats = OUT.rats;
if isempty(rats)
    fprintf('\n[run_non_task_day_stability_2x2] No rats; no stats.\n');
    return
end

mList = find(useMetric);
fprintf('\n=== Non-task cross-day stability stats (paired t-tests; z-space) ===\n');
fprintf('Comparison: An-vs-An-1 stability against An-vs-An+2 stability.\n');

for m = mList
    zA = nan(numel(rats),1);
    zB = nan(numel(rats),1);
    for i = 1:numel(rats)
        zA(i) = rats(i).z_byPair(1,m);
        zB(i) = rats(i).z_byPair(2,m);
    end
    keep = isfinite(zA) & isfinite(zB);
    [pR, tR, dfR] = paired_t_from_vectors(zA(keep), zB(keep));

    mu_zA = mean(zA(keep),'omitnan'); mu_zB = mean(zB(keep),'omitnan');
    mu_rA = tanh(mu_zA);              mu_rB = tanh(mu_zB);

    fprintf('\n[%s]\n', metricNames{m});
    fprintf('  Rats (n=%d): mean r An~An-1=%.3f  An~An+2=%.3f | mean z %.3f vs %.3f | paired t(%d)=%.3f  p=%.3g\n', ...
        nnz(keep), mu_rA, mu_rB, mu_zA, mu_zB, dfR, tR, pR);
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
d = b - a; % positive => An~An+2 > An~An-1
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
function zmu = metric_PV_cross_single(A, B, K, minNcells, minStd)
[nc, ~, nWA] = size(A);
[ncB, ~, nWB] = size(B);
if nc ~= ncB || nc < minNcells || nWA < 1 || nWB < 1
    zmu = NaN; return;
end
z_list = [];
for k = 1:K
    XA = squeeze(A(:,k,:));
    XB = squeeze(B(:,k,:));
    for i = 1:nWA
        v1 = XA(:,i);
        for j = 1:nWB
            v2 = XB(:,j);
            r = safe_corr_pairwise(v1, v2, minNcells, minStd);
            if isfinite(r), z_list(end+1,1) = atanh_clip(r); end %#ok<AGROW>
        end
    end
end
if isempty(z_list), zmu = NaN; else, zmu = mean(z_list,'omitnan'); end
end

function zmu = metric_SC_cross_single(A, B, K, minNbins, minStd)
[nc, ~, nWA] = size(A);
[ncB, ~, nWB] = size(B);
if nc ~= ncB || K < minNbins || nc < 2 || nWA < 1 || nWB < 1
    zmu = NaN; return;
end
z_cell = nan(nc,1);
for c = 1:nc
    XA = squeeze(A(c,:,:));
    XB = squeeze(B(c,:,:));
    z_list = [];
    for i = 1:nWA
        v1 = XA(:,i);
        for j = 1:nWB
            v2 = XB(:,j);
            r = safe_corr_pairwise(v1, v2, minNbins, minStd);
            if isfinite(r), z_list(end+1,1) = atanh_clip(r); end %#ok<AGROW>
        end
    end
    if ~isempty(z_list), z_cell(c) = mean(z_list,'omitnan'); end
end
if all(~isfinite(z_cell)), zmu = NaN; else, zmu = mean(z_cell,'omitnan'); end
end

function zmu = metric_PV_cross_splithalf(A, B, K, nRep, seed, minNcells, minStd)
[nc, ~, nWA] = size(A);
[ncB, ~, nWB] = size(B);
if nc ~= ncB || nc < minNcells || nWA < 4 || nWB < 4
    zmu = NaN; return;
end
rng(seed);
z_rep = nan(nRep,1);
for rr = 1:nRep
    idxA = randperm(nWA);
    idxB = randperm(nWB);
    useA = idxA(1:floor(nWA/2));
    useB = idxB(1:floor(nWB/2));
    z_k = nan(K,1);
    for k = 1:K
        vA = mean(squeeze(A(:,k,useA)), 2, 'omitnan');
        vB = mean(squeeze(B(:,k,useB)), 2, 'omitnan');
        r = safe_corr_pairwise(vA, vB, minNcells, minStd);
        if isfinite(r), z_k(k) = atanh_clip(r); end
    end
    if any(isfinite(z_k)), z_rep(rr) = mean(z_k,'omitnan'); end
end
if all(~isfinite(z_rep)), zmu = NaN; else, zmu = mean(z_rep,'omitnan'); end
end

function zmu = metric_SC_cross_splithalf(A, B, K, nRep, seed, minNbins, minStd)
[nc, ~, nWA] = size(A);
[ncB, ~, nWB] = size(B);
if nc ~= ncB || K < minNbins || nc < 2 || nWA < 4 || nWB < 4
    zmu = NaN; return;
end
rng(seed);
z_rep = nan(nRep,1);
for rr = 1:nRep
    idxA = randperm(nWA);
    idxB = randperm(nWB);
    useA = idxA(1:floor(nWA/2));
    useB = idxB(1:floor(nWB/2));
    z_cell = nan(nc,1);
    for c = 1:nc
        vA = mean(squeeze(A(c,:,useA)), 2, 'omitnan');
        vB = mean(squeeze(B(c,:,useB)), 2, 'omitnan');
        r = safe_corr_pairwise(vA, vB, minNbins, minStd);
        if isfinite(r), z_cell(c) = atanh_clip(r); end
    end
    if any(isfinite(z_cell)), z_rep(rr) = mean(z_cell,'omitnan'); end
end
if all(~isfinite(z_rep)), zmu = NaN; else, zmu = mean(z_rep,'omitnan'); end
end

% =========================================================================
%                         NON-TASK WINDOW MAPS
% =========================================================================
function maps = non_task_maps_for_day(t, posd, csd, daySpikes, GridRC, WinSec, TraceWin, ...
    BufferPost, UseSpeedMask, VelThresh, nNeed, MinWinFrames, MinFramesPerBin, seed)
maps = [];
t = t(:);
if numel(t) < 10, return; end

edges = build_grid_edges_single_day(posd, GridRC);
[x, y] = interp_pos(posd, t);
S_counts = spikes_to_matrix(daySpikes, t);

if UseSpeedMask
    v = speed_cm_per_s(posd);
    v_i = interp1(posd.t(:), v(:), t, 'linear','extrap');
    speed_ok = (v_i >= VelThresh);
else
    speed_ok = true(size(t));
end

is_taskbuf = false(size(t));
csd = csd(:);
if ~isempty(csd)
    for j = 1:numel(csd)
        t0 = csd(j) + TraceWin(1);
        t1 = csd(j) + TraceWin(1) + WinSec + BufferPost;
        is_taskbuf = is_taskbuf | (t >= t0 & t < t1);
    end
end
control_ok = ~is_taskbuf & speed_ok;
nonWins = sample_control_windows_strict(t, control_ok, WinSec, nNeed, MinWinFrames, seed);

[maps, ok] = window_maps_for_day(t, x, y, S_counts, nonWins, edges, speed_ok, ...
    false, MinWinFrames, MinFramesPerBin, control_ok);
if ~ok, maps = []; end
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
    idx = find(t >= t0 & t < t1 & frame_ok);
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
%                          MAP NORMALIZATION
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
        mu(~isfinite(mu) | mu==0) = NaN;
        A = A ./ reshape(mu, [nc 1 1]);
        B = B ./ reshape(mu, [nc 1 1]);
    otherwise
        error('CellNorm must be ''none'', ''demean'', ''zscore'', or ''meanrate''.');
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

function z = atanh_clip(r)
r = max(min(r, 0.999999), -0.999999);
z = atanh(r);
end

% =========================================================================
%                               UTILITIES
% =========================================================================
function keepCells = pair_keep_cells(rat, dayA, dayB, spA, spB)
nA = spike_cell_count(spA);
nB = spike_cell_count(spB);
n = min(nA, nB);
keepCells = true(n,1);

if isfield(rat,'ratemask')
    fA = sprintf('ratemask_%s', dayA);
    fB = sprintf('ratemask_%s', dayB);
    if isfield(rat.ratemask, fA)
        rm = rat.ratemask.(fA);
        keepCells = keepCells & ratemask_to_keep(rm, n);
    end
    if isfield(rat.ratemask, fB)
        rm = rat.ratemask.(fB);
        keepCells = keepCells & ratemask_to_keep(rm, n);
    end
end
end

function keep = ratemask_to_keep(rm, n)
rm = logical(rm(:));
keep = false(n,1);
nUse = min(numel(rm), n);
keep(1:nUse) = rm(1:nUse);
end

function n = spike_cell_count(daySpikes)
daySpikes = unwrap_spike_container(daySpikes);
if iscell(daySpikes)
    n = numel(daySpikes);
elseif isnumeric(daySpikes)
    n = size(daySpikes,1);
else
    n = 0;
end
end

function out = filter_spikes_cells(daySpikes, keepCells)
out = unwrap_spike_container(daySpikes);
if isempty(keepCells), return; end
n = numel(keepCells);
if iscell(out)
    out = out(1:min(numel(out), n));
    keepCells = keepCells(1:numel(out));
    out = out(keepCells);
elseif isnumeric(out)
    out = out(1:min(size(out,1), n), :);
    keepCells = keepCells(1:size(out,1));
    out = out(keepCells,:);
end
end

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
    error('numeric matrix must be [n x 3] or [n x 2].');
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
