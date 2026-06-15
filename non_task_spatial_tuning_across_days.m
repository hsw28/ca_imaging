function R = non_task_spatial_tuning_across_days(ratNames, varargin)
% NON_TASK_SPATIAL_TUNING_ACROSS_DAYS
%
% Biological question:
%   Do cells preserve their spatial tuning across days/environments?
%
% Canonical workflow:
%   1) For each day, use all non-task data to build one spatial rate map per cell.
%   2) Smooth occupancy and spike-count maps, then divide to get rate maps.
%   3) For each cell, correlate the whole spatial map across days.
%   4) Compare An->An-1 versus An->An+2 across rats.
%
% Main output:
%   R(i).cell.z_An_AnMinus1 / r_An_AnMinus1 : per-cell cross-day stability
%   R(i).cell.z_An_AnPlus2  / r_An_AnPlus2  : per-cell cross-day stability
%   R(i).rat.z_An_AnMinus1 / r_An_AnMinus1 : rat mean in Fisher z/r space
%   R(i).rat.z_An_AnPlus2  / r_An_AnPlus2  : rat mean in Fisher z/r space

% ------------------------------- Options ---------------------------------
p = inputParser;
addParameter(p,'GridRC',[10 10], @(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'NumBins',[], @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=4));
addParameter(p,'TraceWin',[0 2], @(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'BufferPost',0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'UseSpeedMask',true, @(x)islogical(x)&&isscalar(x));
addParameter(p,'VelThresh',4, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'SmoothSigmaBins',1.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'SmoothKernelRadius',[], @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=1));
addParameter(p,'MinOccSecPerBin',0.05, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'MinValidBins',8, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'MinCellsPerRat',5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'DoPlots',true, @(x)islogical(x)&&isscalar(x));
addParameter(p,'DoStats',true, @(x)islogical(x)&&isscalar(x));
addParameter(p,'Verbose',true, @(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
opt = p.Results;

if ~isempty(opt.NumBins)
    opt.GridRC = best_factors_local(opt.NumBins);
end
opt.GridRC = double(opt.GridRC(:).');

ratNames = cellstr(ratNames);
R = struct([]);

for ii = 1:numel(ratNames)
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
    if ~isfield(rat,'An') || isempty(rat.An)
        warning('[%s] rat.An not found. Skipping.', ratVar);
        continue;
    end

    idxAn = find(strcmp(dateList, rat.An), 1);
    if isempty(idxAn)
        warning('[%s] rat.An=%s not found in dateList. Skipping.', ratVar, char(rat.An));
        continue;
    end

    pairIdx = [idxAn-1, idxAn, idxAn+2];
    pairOk = pairIdx >= 1 & pairIdx <= numel(dateList);
    daysNeeded = unique(dateList(pairIdx(pairOk)), 'stable');
    if isempty(daysNeeded)
        warning('[%s] No usable days around An=%s.', ratVar, char(rat.An));
        continue;
    end

    spikes_raw = filterFieldsByDay_fallback(rat.Ca_peaks, daysNeeded);
    ts_raw     = filterFieldsByDay_fallback(rat.Ca_ts,    daysNeeded);
    pos_raw    = filterFieldsByDay_fallback(rat.pos,      daysNeeded);
    cs_raw     = filterFieldsByDay_fallback(rat.CS_times, daysNeeded);
    [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, daysNeeded);

    mapsByDay = cell(numel(daysNeeded),1);
    validByDay = cell(numel(daysNeeded),1);
    cellKeepByDay = cell(numel(daysNeeded),1);
    occByDay = cell(numel(daysNeeded),1);
    edgesByDay = cell(numel(daysNeeded),1);

    for d = 1:numel(daysNeeded)
        dlabel = daysNeeded{d};
        cellKeepByDay{d} = ratemask_for_day(rat, dlabel, spikes{d});

        [mapsByDay{d}, validByDay{d}, occByDay{d}, edgesByDay{d}] = build_non_task_rate_maps_day( ...
            spikes{d}, ts{d}, pos{d}, cs{d}, opt);
    end

    idxLocalAn = find(strcmp(daysNeeded, dateList{idxAn}), 1);
    idxLocalMinus1 = [];
    idxLocalPlus2 = [];
    if idxAn-1 >= 1
        idxLocalMinus1 = find(strcmp(daysNeeded, dateList{idxAn-1}), 1);
    end
    if idxAn+2 <= numel(dateList)
        idxLocalPlus2 = find(strcmp(daysNeeded, dateList{idxAn+2}), 1);
    end

    [zMinus, rMinus, nBinsMinus] = correlate_day_cell_maps(mapsByDay, validByDay, cellKeepByDay, idxLocalAn, idxLocalMinus1, opt.MinValidBins);
    [zPlus,  rPlus,  nBinsPlus]  = correlate_day_cell_maps(mapsByDay, validByDay, cellKeepByDay, idxLocalAn, idxLocalPlus2,  opt.MinValidBins);

    R(ii).animal = ratVar; %#ok<AGROW>
    R(ii).meta.options = opt;
    R(ii).meta.dateList = dateList;
    R(ii).meta.An = rat.An;
    R(ii).meta.daysNeeded = daysNeeded;
    R(ii).meta.day_An = dateList{idxAn};
    R(ii).meta.day_AnMinus1 = day_or_empty(dateList, idxAn-1);
    R(ii).meta.day_AnPlus2 = day_or_empty(dateList, idxAn+2);
    R(ii).maps.rate_byDay = mapsByDay;
    R(ii).maps.valid_byDay = validByDay;
    R(ii).maps.cellKeep_byDay = cellKeepByDay;
    R(ii).maps.occupancy_byDay = occByDay;
    R(ii).maps.edges_byDay = edgesByDay;
    R(ii).cell.z_An_AnMinus1 = zMinus;
    R(ii).cell.r_An_AnMinus1 = rMinus;
    R(ii).cell.nBins_An_AnMinus1 = nBinsMinus;
    R(ii).cell.z_An_AnPlus2 = zPlus;
    R(ii).cell.r_An_AnPlus2 = rPlus;
    R(ii).cell.nBins_An_AnPlus2 = nBinsPlus;
    R(ii).rat.z_An_AnMinus1 = mean(zMinus(isfinite(zMinus)), 'omitnan');
    R(ii).rat.z_An_AnPlus2 = mean(zPlus(isfinite(zPlus)), 'omitnan');
    R(ii).rat.r_An_AnMinus1 = tanh(R(ii).rat.z_An_AnMinus1);
    R(ii).rat.r_An_AnPlus2 = tanh(R(ii).rat.z_An_AnPlus2);
    R(ii).rat.nCells_An_AnMinus1 = nnz(isfinite(zMinus));
    R(ii).rat.nCells_An_AnPlus2 = nnz(isfinite(zPlus));

    if opt.Verbose
        fprintf('[%s] %s~%s: r=%.3f nCells=%d | %s~%s: r=%.3f nCells=%d\n', ...
            ratVar, R(ii).meta.day_An, R(ii).meta.day_AnMinus1, R(ii).rat.r_An_AnMinus1, R(ii).rat.nCells_An_AnMinus1, ...
            R(ii).meta.day_An, R(ii).meta.day_AnPlus2, R(ii).rat.r_An_AnPlus2, R(ii).rat.nCells_An_AnPlus2);
    end
end

if opt.DoStats
    print_group_stats(R, opt);
end
if opt.DoPlots
    plot_group_summary(R);
end
end

% =========================================================================
%                             CORE ANALYSIS
% =========================================================================
function [rateMaps, validMaps, occSmooth, edges] = build_non_task_rate_maps_day(spikes, t, posd, csd, opt)
t = double(t(:));
Nc = spike_cell_count(spikes);
R = opt.GridRC(1);
C = opt.GridRC(2);
K = R*C;

rateMaps = nan(Nc, K);
validMaps = false(Nc, K);
occSmooth = nan(R, C);
edges = build_grid_edges_single_day(posd, opt.GridRC);

if numel(t) < 3 || Nc < 1
    return;
end

[x, y] = interp_pos(posd, t);
S = spikes_to_matrix(spikes, t);

dt = median(diff(t));
if ~isfinite(dt) || dt <= 0
    dt = max(eps, mean(diff(t),'omitnan'));
end

use = isfinite(t) & isfinite(x) & isfinite(y);

csd = csd(:);
if ~isempty(csd)
    is_taskbuf = false(size(t));
    for j = 1:numel(csd)
        t0 = csd(j) + opt.TraceWin(1);
        t1 = csd(j) + opt.TraceWin(2) + opt.BufferPost;
        is_taskbuf = is_taskbuf | (t >= t0 & t < t1);
    end
    use = use & ~is_taskbuf;
end

if opt.UseSpeedMask
    v = speed_cm_per_s(posd);
    v_i = interp1(posd.t(:), v(:), t, 'linear','extrap');
    use = use & (v_i >= opt.VelThresh);
end

if nnz(use) < 2
    return;
end

[~, kbin] = pos2bin(x(use), y(use), edges);
Suse = S(:, use);
kbin = kbin(:);

occ = zeros(R, C);
spikeCount = zeros(Nc, R, C);

for k = 1:K
    idx = find(kbin == k);
    if isempty(idx), continue; end
    [rr, cc] = ind2sub([R C], k);
    occ(rr, cc) = numel(idx) * dt;
    spikeCount(:, rr, cc) = sum(Suse(:, idx), 2, 'omitnan');
end

occSmooth = smooth2_nanaware(occ, opt.SmoothSigmaBins, opt.SmoothKernelRadius);
occValid = occSmooth >= opt.MinOccSecPerBin;

for c = 1:Nc
    sc = squeeze(spikeCount(c,:,:));
    scSmooth = smooth2_nanaware(sc, opt.SmoothSigmaBins, opt.SmoothKernelRadius);
    rm = scSmooth ./ occSmooth;
    rm(~occValid) = NaN;
    rateMaps(c,:) = reshape(rm, 1, []);
    validMaps(c,:) = isfinite(rateMaps(c,:));
end
end

function [z, r, nBins] = correlate_day_cell_maps(mapsByDay, validByDay, cellKeepByDay, idxA, idxB, minValidBins)
z = [];
r = [];
nBins = [];
if isempty(idxA) || isempty(idxB) || idxA < 1 || idxB < 1 || idxA > numel(mapsByDay) || idxB > numel(mapsByDay)
    return;
end
A = mapsByDay{idxA};
B = mapsByDay{idxB};
VA = validByDay{idxA};
VB = validByDay{idxB};
KA = cellKeepByDay{idxA};
KB = cellKeepByDay{idxB};
if isempty(A) || isempty(B)
    return;
end

Nc = min(size(A,1), size(B,1));
z = nan(Nc,1);
r = nan(Nc,1);
nBins = zeros(Nc,1);

for c = 1:Nc
    if c > numel(KA) || c > numel(KB) || ~KA(c) || ~KB(c)
        continue;
    end
    m = VA(c,:) & VB(c,:) & isfinite(A(c,:)) & isfinite(B(c,:));
    nBins(c) = nnz(m);
    if nBins(c) < minValidBins
        continue;
    end
    a = A(c,m).';
    b = B(c,m).';
    if std(a) <= 1e-12 || std(b) <= 1e-12
        continue;
    end
    rr = corr(a, b, 'type','Pearson');
    if isfinite(rr)
        r(c) = rr;
        z(c) = atanh_clip(rr);
    end
end
end

% =========================================================================
%                                  STATS
% =========================================================================
function print_group_stats(R, opt)
if isempty(R)
    fprintf('\n[canonical_non_task_spatial_tuning_across_days] No rats.\n');
    return;
end

zMinus = nan(numel(R),1);
zPlus = nan(numel(R),1);
nMinus = nan(numel(R),1);
nPlus = nan(numel(R),1);
for i = 1:numel(R)
    zMinus(i) = R(i).rat.z_An_AnMinus1;
    zPlus(i) = R(i).rat.z_An_AnPlus2;
    nMinus(i) = R(i).rat.nCells_An_AnMinus1;
    nPlus(i) = R(i).rat.nCells_An_AnPlus2;
end

keep = isfinite(zMinus) & isfinite(zPlus) & nMinus >= opt.MinCellsPerRat & nPlus >= opt.MinCellsPerRat;
[p, tstat, df] = paired_t_from_vectors(zMinus(keep), zPlus(keep));

fprintf('\n=== Canonical non-task spatial tuning across days ===\n');
fprintf('Per rat: mean Fisher z of per-cell full-map correlations.\n');
fprintf('Comparison: An~An-1 versus An~An+2.\n');
fprintf('Rats (n=%d): mean r An~An-1=%.3f  An~An+2=%.3f | mean z %.3f vs %.3f | paired t(%d)=%.3f  p=%.3g\n\n', ...
    nnz(keep), tanh(mean(zMinus(keep),'omitnan')), tanh(mean(zPlus(keep),'omitnan')), ...
    mean(zMinus(keep),'omitnan'), mean(zPlus(keep),'omitnan'), df, tstat, p);
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
%                                  PLOT
% =========================================================================
function plot_group_summary(R)
if isempty(R), return; end

nA = numel(R);
zMinus = nan(nA,1);
zPlus = nan(nA,1);
for i = 1:nA
    zMinus(i) = R(i).rat.z_An_AnMinus1;
    zPlus(i) = R(i).rat.z_An_AnPlus2;
end

cmap = lines(nA);
figure('Color','w','Position',[180 180 720 560]); hold on;

keepMinus = isfinite(zMinus);
keepPlus = isfinite(zPlus);
if any(keepMinus)
    bar(1, tanh(mean(zMinus(keepMinus),'omitnan')), 0.6, 'FaceColor',[0.35 0.70 1.00], 'EdgeColor','k');
end
if any(keepPlus)
    bar(2, tanh(mean(zPlus(keepPlus),'omitnan')), 0.6, 'FaceColor',[0.85 0.55 0.25], 'EdgeColor','k');
end

for i = 1:nA
    if isfinite(zMinus(i)) && isfinite(zPlus(i))
        plot([1 2], tanh([zMinus(i) zPlus(i)]), '-', 'Color', cmap(i,:), 'LineWidth', 2.5);
        plot([1 2], tanh([zMinus(i) zPlus(i)]), 'o', 'MarkerFaceColor', cmap(i,:), ...
            'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 7);
    elseif isfinite(zMinus(i))
        plot(1, tanh(zMinus(i)), 'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'MarkerSize', 7);
    elseif isfinite(zPlus(i))
        plot(2, tanh(zPlus(i)), 'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'MarkerSize', 7);
    end
end

xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'An vs An-1','An vs An+2'});
ylabel('Mean per-cell spatial map correlation (r)');
title('Canonical non-task spatial tuning stability');
yline(0,'k:');
grid on; box on;
end

% =========================================================================
%                             SMOOTHING HELPERS
% =========================================================================
function Y = smooth2_nanaware(X, sigma, radius)
X = double(X);
if sigma <= 0
    Y = X;
    return;
end
if isempty(radius)
    radius = max(1, ceil(3*sigma));
end

g = gaussian_kernel_1d(sigma, radius);
valid = isfinite(X);
X0 = X;
X0(~valid) = 0;

num = conv2(conv2(X0, g, 'same'), g.', 'same');
den = conv2(conv2(double(valid), g, 'same'), g.', 'same');
Y = num ./ den;
Y(den <= 0) = NaN;
end

function g = gaussian_kernel_1d(sigma, radius)
x = -radius:radius;
g = exp(-(x.^2) ./ (2*sigma.^2));
g = g ./ sum(g);
end

% =========================================================================
%                             DATA UTILITIES
% =========================================================================
function keepCells = ratemask_for_day(rat, dlabel, daySpikes)
n = spike_cell_count(daySpikes);
keepCells = true(n,1);
if isfield(rat,'ratemask')
    f = sprintf('ratemask_%s', dlabel);
    if isfield(rat.ratemask, f)
        rm = logical(rat.ratemask.(f));
        rm = rm(:);
        keepCells = false(n,1);
        nUse = min(n, numel(rm));
        keepCells(1:nUse) = rm(1:nUse);
    end
end
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

function z = atanh_clip(r)
r = max(min(r, 0.999999), -0.999999);
z = atanh(r);
end

function d = day_or_empty(dateList, idx)
if idx >= 1 && idx <= numel(dateList)
    d = dateList{idx};
else
    d = '';
end
end

function rc = best_factors_local(N)
validateattributes(N, {'numeric'}, {'scalar','integer','positive','finite'});
r = floor(sqrt(double(N)));
while r > 1 && mod(N, r) ~= 0
    r = r - 1;
end
c = N / r;
rc = [r, c];
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

function S = filterFieldsByDay_fallback(Sin, daysToUse)
if exist('filterFieldsByDay','file') == 2
    try
        S = filterFieldsByDay(Sin, daysToUse);
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
