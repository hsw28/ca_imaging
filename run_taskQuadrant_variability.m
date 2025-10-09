function R = run_taskQuadrant_variability(ratNames, varargin)
% run_taskQuadrant_variability
% Goal: During the TRACE window, split task frames into even/odd samples
% and compare task PV similarity "within quadrant" vs "across quadrants".
%
% Key idea: If across-quadrant similarity >= within-quadrant similarity,
% spatial context is disrupting (or mixing with) task coding.
%
% Usage:
%   R = run_taskQuadrant_variability({'rat0222','rat0307',...}, ...
%         'SubGridRC',[2 2], 'TraceWin',[0 2], 'UseSpeedMask',true, ...
%         'MinTraceEachHalf',1, 'MinTraceSamples',6, 'BinWeight','traceCount');
%
% Outputs per rat:
%   R(i).C.within.r / .z      : within-quadrant mean similarity (r and Fisher-z)
%   R(i).C.across.r / .z      : across-quadrant mean similarity
%   R(i).C.delta.r / .z       : within - across
%   R(i).C.p_perm             : strict permutation p-value
%   R(i).C.per_day            : day-wise within/across (z), valid quadrants, counts
%   R(i).meta                 : options and grid/ROI info
%
% Notes:
% - Uses only TRACE samples (t in CS+[TraceWin]).
% - Even/odd split is by frame index *within* the trace window per trial.
% - Quadrants come from the day’s trace ROI (percentile box with small margin).

% -------- Options --------
p = inputParser;
addParameter(p,'DaysMode','last3toAn');   % same idea as your runner
addParameter(p,'TraceWin',[0 2]);         % seconds from CS onset
addParameter(p,'SubGridRC',[2 2]);        % rows×cols (2×2 => quadrants)
addParameter(p,'ROIPrc',[5 95]);         % percentile ROI of trace samples
addParameter(p,'ROIMarginFrac',0.05);     % small padding on ROI box
addParameter(p,'VelThresh',4);            % cm/s (used only if UseSpeedMask)
addParameter(p,'UseSpeedMask',true);
addParameter(p,'CellNorm','demean');      % 'none' | 'demean' | 'zscore'
addParameter(p,'MinTraceEachHalf',1);     % >= this many even and odd samples
addParameter(p,'MinTraceSamples',3);      % total even+odd samples per quadrant
addParameter(p,'BinWeight','traceCount'); % 'none'|'traceCount'
addParameter(p,'NPerm',500);
addParameter(p,'DoPlots',true);
parse(p,varargin{:});
opt = p.Results;

ratNames = cellstr(ratNames);
R = struct([]);

for ii = 1:numel(ratNames)
    ratVar = ratNames{ii};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratVar))
        warning('Variable %s not found in base workspace. Skipping.', ratVar);
        continue;
    end
    rat = evalin('base', ratVar);

    % ----- which days -----
    dateList = autoDateList(rat);  % from your helper (Ca_traces fieldnames)
    switch lower(opt.DaysMode)
        case {'all','last3byca'}               % NEW: robust against An mismatch
            daysByCa = dateList;
            if strcmpi(opt.DaysMode,'last3byca') && numel(daysByCa) >= 3
                daysToUse = daysByCa(end-2:end);
            else
                daysToUse = daysByCa;
            end
        case 'last3toan'
            if isfield(rat,'An')
                idx = find(strcmp(dateList, rat.An), 1);
                if ~isempty(idx) && idx >= 3
                    daysToUse = dateList(idx-2:idx);
                else
                    % Fallback: still use dateList (not "dayXX") so labels match
                    daysToUse = dateList(max(1,end-2):end);
                end
            else
                daysToUse = dateList(max(1,end-2):end);
            end
        otherwise
            daysToUse = dateList;  % safe default
    end
    fprintf('[%s] daysToUse: %s\n', ratVar, strjoin(daysToUse, ', '));


    % ----- pull & standardize -----
    spikes_raw = filterFieldsByDay(rat.Ca_peaks, daysToUse);
    ts_raw     = filterFieldsByDay(rat.Ca_ts,    daysToUse);
    pos_raw    = filterFieldsByDay(rat.pos,      daysToUse);
    cs_raw     = filterFieldsByDay(rat.CS_times, daysToUse);
    [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, daysToUse);


    % ----- compute within/across quadrant similarity per day -----
    [z_within_day, z_across_day, per_day_dbg, grid_info] = ...
        quadrant_task_similarity_per_day(spikes, ts, pos, cs, opt);

    has_any = isfinite(z_within_day) | isfinite(z_across_day);
    zW = mean(z_within_day(has_any), 'omitnan');   % observed within (Fisher z)
    zA = mean(z_across_day(has_any), 'omitnan');   % observed across (Fisher z)
    delta_z = zW - zA;

    % strict permutation (permute odd-quadrant labels)
    [pR, pL, pTwo, delta_perm] = quadrant_strict_permutation(spikes, ts, pos, cs, opt, zW, zA);

    C.within.z = zW;      C.within.r = tanh(zW);
    C.across.z = zA;      C.across.r = tanh(zA);
    C.delta.z  = delta_z; C.delta.r  = tanh(delta_z);
    C.p_perm_right = pR;          % H1: within ≥ across
    C.p_perm_left  = pL;          % H1: within ≤ across  (disruption)
    C.p_perm_two   = pTwo;        % two-sided
    C.perm_delta   = delta_perm;  % optional: store null deltas for later QC
    C.per_day      = per_day_dbg;


    R(ii).animal = ratVar;
    R(ii).C = C;
    R(ii).meta.options = opt;
    R(ii).meta.days    = daysToUse;
    R(ii).meta.grid    = grid_info;
end

if opt.DoPlots, plot_quadrant_results(R); end
end

% ========== Core day-wise computation ==========
function [zW, zA, dbg_all, grid_info] = quadrant_task_similarity_per_day(spikes, ts, pos, cs, o)
D = numel(ts);
rows = o.SubGridRC(1); cols = o.SubGridRC(2); K = rows*cols;

zW = nan(D,1);   % within-quadrant Fisher z per day
zA = nan(D,1);   % across-quadrant Fisher z per day
dbg_all = repmat(struct('validK',[], 'n_even',[], 'n_odd',[], 'n_tot',[], ...
                        'weights',[], 'pairs',[], 'zW_i',[], 'zA_i',[], ...
                        'edges',[], 'ROI',[]), 1, D);

grid_info.SubGridRC = o.SubGridRC;
grid_info.ROIPrc = o.ROIPrc;
grid_info.ROIMarginFrac = o.ROIMarginFrac;

for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);
    S = spikes_to_matrix(spikes{d}, t);
    S = normalize_cells(S, o.CellNorm);

    % ----- per-day trace ROI and grid edges from trace samples -----
    K = o.SubGridRC(1)*o.SubGridRC(2);
    [edges, ~, ~, ROI] = build_grid_edges_from_trialROI_day( ...
        pos{d}, ts{d}, cs{d}, o.TraceWin, K, o.ROIPrc, o.ROIMarginFrac, o.SubGridRC);    grid_info.edges{d} = edges; grid_info.ROI{d} = ROI;

    % ----- optional speed mask -----
    if o.UseSpeedMask
        v = speed_cm_per_s(pos{d});
        v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
        mask_speed = (v >= o.VelThresh);
    else
        mask_speed = true(size(t));
    end

    PV_even = nan(size(S,1), K);
    PV_odd  = nan(size(S,1), K);
    n_even  = zeros(K,1);
    n_odd   = zeros(K,1);
    has_k   = false(K,1);

    csd = cs{d};
    for tr = 1:numel(csd)
        tidx = t >= csd(tr)+o.TraceWin(1) & t <= csd(tr)+o.TraceWin(2);
        tidx = tidx & mask_speed;
        idxs = find(tidx);
        if numel(idxs) < 2, continue; end

        idx_even = idxs(2:2:end);
        idx_odd  = idxs(1:2:end);
        [~, b_even] = pos2bin(x(idx_even), y(idx_even), edges);
        [~, b_odd ] = pos2bin(x(idx_odd ), y(idx_odd ), edges);

        for k = 1:K
            if any(b_even==k)
                ie = idx_even(b_even==k);
                PV_even(:,k) = nanmean([PV_even(:,k), mean(S(:,ie),2)],2);
                n_even(k) = n_even(k) + numel(ie);
                has_k(k) = true;
            end
            if any(b_odd==k)
                io = idx_odd(b_odd==k);
                PV_odd(:,k)  = nanmean([PV_odd(:,k),  mean(S(:,io),2)],2);
                n_odd(k) = n_odd(k) + numel(io);
                has_k(k) = true;
            end
        end
    end

    n_tot = n_even + n_odd;
    gate = (n_even >= o.MinTraceEachHalf) & (n_odd >= o.MinTraceEachHalf) & (n_tot >= o.MinTraceSamples);
    validK = find(has_k & gate);

    if isempty(validK)
        fprintf('[Quadrant d=%d] valid=[] | median n_tot=NaN | across pairs=0\n', d);
        continue;
    end

    % weights (default: trace counts)
    switch lower(o.BinWeight)
        case 'none'
            w = ones(K,1);
        case 'tracecount'
            w = n_tot;
        otherwise
            error('Unknown BinWeight: %s', o.BinWeight);
    end
    w(~isfinite(w))=0; w_valid = w(validK); if ~any(w_valid>0), w_valid(:)=1; end

    % WITHIN (per quadrant): corr(PV_even, PV_odd)
    zWi = nan(numel(validK),1);
    for i = 1:numel(validK)
        k = validK(i);
        zWi(i) = atanh(bound_r(safe_corr(PV_even(:,k), PV_odd(:,k))));
    end
    zW(d) = nansum(w_valid .* zWi) / max(1, nansum(w_valid));

    % ACROSS (between different quadrants): corr(mean(PV_even,PV_odd)_k1, _k2)
    if numel(validK) >= 2
        pairs = nchoosek(validK(:)',2);
        zAi = nan(size(pairs,1),1); wAi = nan(size(pairs,1),1);
        for i = 1:size(pairs,1)
            k1 = pairs(i,1); k2 = pairs(i,2);
            v1 = mean([PV_even(:,k1), PV_odd(:,k1)],2,'omitnan');
            v2 = mean([PV_even(:,k2), PV_odd(:,k2)],2,'omitnan');
            zAi(i) = atanh(bound_r(safe_corr(v1, v2)));
            wAi(i) = sqrt(w(k1)*w(k2));
        end
        wAi(~isfinite(wAi))=0; if ~any(wAi>0), wAi(:)=1; end
        zA(d) = nansum(wAi .* zAi) / max(1, nansum(wAi));
    else
        zA(d) = NaN;
        pairs = [];
    end

    % debug/record
    dbg_all(d).validK   = validK;
    dbg_all(d).n_even   = n_even;
    dbg_all(d).n_odd    = n_odd;
    dbg_all(d).n_tot    = n_tot;
    dbg_all(d).weights  = w;
    dbg_all(d).zW_i     = zWi;
    dbg_all(d).zA_i     = [];
    if ~isempty(pairs), dbg_all(d).pairs = pairs; dbg_all(d).zA_i = zAi; end
    dbg_all(d).edges    = edges;
    dbg_all(d).ROI      = ROI;

    pairs_count = max(numel(validK)*(numel(validK)-1)/2, 0);
    fprintf('[Quadrant d=%d] valid=%s | median n_tot=%g | across pairs=%d\n', ...
        d, mat2str(validK(:)'), median(n_tot(validK)), pairs_count);
end
end

% ========== Strict permutation (shuffle quadrant pairing of odd) ==========
function [p_right,p_left,p_two,delta_perm] = quadrant_strict_permutation(spikes, ts, pos, cs, o, zW_obs, zA_obs)
D = numel(ts);
rows = o.SubGridRC(1); cols = o.SubGridRC(2); K = rows*cols;

zW_perm = nan(o.NPerm,1);
if o.NPerm<=0, p_right=NaN; p_left=NaN; p_two=NaN; delta_perm=[]; return; end

for pidx = 1:o.NPerm
    zW_day = nan(D,1);
    for d = 1:D
        % Recompute PV_even/PV_odd for this day, then re-pair with a permutation
        t = ts{d}(:);
        [x, y] = interp_pos(pos{d}, t);
        S = normalize_cells(spikes_to_matrix(spikes{d}, t), o.CellNorm);

        [edges, ~] = build_traceROI_edges_day(pos{d}, ts{d}, cs{d}, ...
                                              o.TraceWin, o.ROIPrc, o.ROIMarginFrac, o.SubGridRC);

        if o.UseSpeedMask
            v = speed_cm_per_s(pos{d});
            v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
            mask_speed = (v >= o.VelThresh);
        else
            mask_speed = true(size(t));
        end

        PV_even = nan(size(S,1), K);
        PV_odd  = nan(size(S,1), K);
        n_even = zeros(K,1); n_odd = zeros(K,1); has_k=false(K,1);

        csd = cs{d};
        for tr = 1:numel(csd)
            tidx = t >= csd(tr)+o.TraceWin(1) & t <= csd(tr)+o.TraceWin(2);
            tidx = tidx & mask_speed;
            idxs = find(tidx); if numel(idxs) < 2, continue; end
            idx_even = idxs(2:2:end); idx_odd = idxs(1:2:end);
            [~, b_even] = pos2bin(x(idx_even), y(idx_even), edges);
            [~, b_odd ] = pos2bin(x(idx_odd ), y(idx_odd ), edges);

            for k = 1:K
                if any(b_even==k)
                    ie = idx_even(b_even==k);
                    PV_even(:,k) = nanmean([PV_even(:,k), mean(S(:,ie),2)],2);
                    n_even(k) = n_even(k) + numel(ie); has_k(k)=true;
                end
                if any(b_odd==k)
                    io = idx_odd(b_odd==k);
                    PV_odd(:,k)  = nanmean([PV_odd(:,k),  mean(S(:,io),2)],2);
                    n_odd(k) = n_odd(k) + numel(io);  has_k(k)=true;
                end
            end
        end

        n_tot = n_even + n_odd;
        gate = (n_even >= o.MinTraceEachHalf) & (n_odd >= o.MinTraceEachHalf) & (n_tot >= o.MinTraceSamples);
        validK = find(has_k & gate); if isempty(validK), continue; end

        % weights (trace counts)
        w = n_tot; w(~isfinite(w))=0;

        % permute odd-quadrant labels, keep within-quadrant weighting
        perm_valid = validK(randperm(numel(validK)));
        zP = nan(numel(validK),1); wP = nan(numel(validK),1);
        for j = 1:numel(validK)
            k1 = validK(j); k2 = perm_valid(j);
            zP(j) = atanh(max(min(safe_corr(PV_even(:,k1), PV_odd(:,k2)),0.99999),-0.99999));
            wP(j) = sqrt(w(k1)*w(k2));
        end
        wP(~isfinite(wP))=0; if ~any(wP>0), wP(:)=1; end
        zW_day(d) = nansum(wP .* zP) / max(1, nansum(wP));
    end
    zW_perm(pidx) = mean(zW_day(isfinite(zW_day)),'omitnan');
end

zW_perm   = zW_perm(isfinite(zW_perm));
delta_obs = zW_obs - zA_obs;          % << observed within minus across
delta_perm = zW_perm - zA_obs;        % << permuted within minus fixed across

p_right = mean(delta_perm >= delta_obs)           % H1: within ≥ across
p_left  = mean(delta_perm <= delta_obs);          % H1: within ≤ across  (disruption)
p_two   = 2*min(p_right, p_left);

p = p_two;
end

% ========== Plotting ==========
function plot_quadrant_results(R)
if isempty(R), return; end
within = arrayfun(@(s) s.C.within.r, R);
across = arrayfun(@(s) s.C.across.r, R);
mask = isfinite(within) & isfinite(across);

figure('Color','w','Position',[140 140 540 420]); hold on
idx = find(mask);
for ii = 1:numel(idx)
    k = idx(ii);
    plot([1 2], [within(k) across(k)], '-o', 'LineWidth', 1.2); hold on
end
bar(1, mean(within(mask),'omitnan'), 0.6, 'EdgeColor','k');
bar(2, mean(across(mask),'omitnan'), 0.6, 'EdgeColor','k');
xticks([1 2]); xticklabels({'within quadrant','across quadrants'});
ylabel('Mean PV correlation (r)'); grid on; box on
title(sprintf('Quadrant task stability (n=%d rats)', nnz(mask)));
% ---- Paired t-test on Fisher-z (within vs across) ----
wz = arrayfun(@(s) s.C.within.z, R);
az = arrayfun(@(s) s.C.across.z, R);
maskZ = isfinite(wz) & isfinite(az);

if any(maskZ)
    [~, p, ~, stats] = ttest(wz(maskZ), az(maskZ));  % paired t on z
    dz = mean(wz(maskZ)-az(maskZ)) / std(wz(maskZ)-az(maskZ));  % Cohen's dz

    % annotate above the bars
    ax = gca;
    yl = ylim(ax);
    y  = yl(2) + 0.05*range(yl);
    line(ax, [1 2], [y y], 'Color','k','LineWidth',1.2);
    text(1.5, y + 0.02*range(yl), ...
        sprintf('paired t (z): t(%d)=%.2f, p=%.3g', stats.df, stats.tstat, p), ...
        'HorizontalAlignment','center');
    % Also print permutation p's per rat:
    arrayfun(@(s) fprintf('%s perm p_right=%.3g  p_left=%.3g  p_two=%.3g\n', ...
        s.animal, s.C.p_perm_right, s.C.p_perm_left, s.C.p_perm_two), R);
    ylim(ax, [yl(1) y + 0.10*range(yl)]);

    % (optional) also show a nonparametric check on r
    wr = arrayfun(@(s) s.C.within.r, R);
    ar = arrayfun(@(s) s.C.across.r, R);
    try
        [p_sr,~,sr_stats] = signrank(wr(maskZ), ar(maskZ));
        fprintf('Wilcoxon signed-rank on r: W=%d, p=%.3g\n', sr_stats.signedrank, p_sr);
    catch
        % signrank may not exist in older toolboxes
    end
end

end

% =================== Utilities (self-contained) ===================
function days = auto_or_all_days(rat)
if isfield(rat,'dates'), days = rat.dates(:)';
elseif isfield(rat,'pos')
    if iscell(rat.pos), days = arrayfun(@(k) sprintf('day%02d',k), 1:numel(rat.pos), 'uni',0);
    elseif isstruct(rat.pos), days = arrayfun(@(k) sprintf('day%02d',k), 1:numel(rat.pos), 'uni',0);
    else, days = {'day01'};
    end
else, days = {'day01'};
end
end

function [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, dayKeys)
% Build 1×D cell arrays in the SAME order as dayKeys (e.g., {'2023_05_04',...})

spikes = to_daycells(spikes_raw, dayKeys);
ts_in  = to_daycells(ts_raw,   dayKeys);
pos_in = to_daycells(pos_raw,  dayKeys);
cs_in  = to_daycells(cs_raw,   dayKeys);

D = numel(dayKeys);
ts  = cell(1,D);
pos = cell(1,D);
cs  = cell(1,D);

for d = 1:D
    ts{d}  = coerce_ts_day(ts_in{d});
    pos{d} = coerce_pos_day(pos_in{d}, ts{d});
    cs{d}  = coerce_cs_day(cs_in{d}, ts{d});
end
end

function C = to_daycells(X, keys)
% Return cell array aligned to keys. Works for scalar structs with date fields.

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
            if isempty(idx)
                % relaxed fallback: match by substring if exact field not found
                idx = find(contains(fn, k), 1);
            end
            if isempty(idx), C{i} = []; else, C{i} = X.(fn{idx}); end
        end
        return
    else
        % array-of-structs fallback by order
        C = arrayfun(@(j) X(j), 1:min(numel(X), numel(keys)), 'uni', 0);
        if numel(C) < numel(keys), C(end+1:numel(keys)) = {[]}; end
        return
    end
end

% last resort: replicate scalar non-struct across days
C = repmat({X}, 1, numel(keys));
end

function C = to_cell(x)
if iscell(x), C=x; return; end
if isstruct(x) && numel(x)>1, C=arrayfun(@(k) x(k), 1:numel(x), 'uni',0); return; end
if isstruct(x) && numel(x)==1, C={x}; return; end
if istable(x) || isnumeric(x) || isstring(x) || ischar(x) || isempty(x), C={x}; return; end
error('Unsupported container type: %s', class(x));
end

function t = coerce_ts_day(td)
if iscell(td), if isempty(td), t=[]; return; else, td=td{1}; end, end
if isnumeric(td)&&isvector(td), t=double(td(:)); if max(t,[],'omitnan')>1e4, t=t/1000; end; return, end
if isnumeric(td)&&ismatrix(td)&&~isscalar(td)
    M=double(td); if size(M,2)>=2, t=M(:,2); else, t=M(:,1); end
    if max(t,[],'omitnan')>1e4, t=t/1000; end; t=t(:); return
end
if istable(td)
    vn=lower(string(td.Properties.VariableNames));
    cname=pick_name(vn, ["time_ms","time","ts","t","timestamp","frame_ts","ca_ts"]);
    if strlength(cname)==0, error('ts table: no time-like column found'); end
    t=double(td{:,find(vn==cname,1)}); if contains(cname,"ms")||max(t,[],'omitnan')>1e4, t=t/1000; end; t=t(:); return
end
if isstruct(td)
    cand={td}; fns=fieldnames(td);
    for k=1:numel(fns), if isstruct(td.(fns{k})), cand{end+1}=td.(fns{k}); end, end
    for c=1:numel(cand)
        C=cand{c}; fn=fieldnames(C);
        for j=1:numel(fn)
            v=C.(fn{j});
            if isnumeric(v)&&ismatrix(v)&&size(v,2)>=2&&size(v,1)>1
                t=double(v(:,2)); if max(t,[],'omitnan')>1e4, t=t/1000; end; t=t(:); return
            end
        end
    end
end
error('Unsupported ts type: %s', class(td));
end

function P = coerce_pos_day(pd, tday)
    if nargin<2, tday = []; end
    col = @(v) v(:);

    % unwrap single-cell containers
    if iscell(pd)
        if isempty(pd), P = struct('t',[],'x',[],'y',[]); return; end
        if numel(pd)==1, pd = pd{1}; else, pd = pd{1}; end
    end

    % empty fallback
    if isempty(pd)
        if ~isempty(tday)
            P = struct('t',col(tday),'x',nan(numel(tday),1),'y',nan(numel(tday),1));
        else
            P = struct('t',[],'x',[],'y',[]);
        end
        return
    end

    % table with any reasonable column names
    if istable(pd)
        vn = lower(string(pd.Properties.VariableNames));
        tname = pick_name(vn, ["t","time","ts"]);
        xname = pick_name(vn, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"]);
        yname = pick_name(vn, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"]);
        t=[]; if ~isempty(tname), t = pd{:, find(vn==tname,1)}; end
        x=[]; if ~isempty(xname), x = pd{:, find(vn==xname,1)}; end
        y=[]; if ~isempty(yname), y = pd{:, find(vn==yname,1)}; end
        [t,x,y] = finalize_txy(t,x,y,tday);
        P = struct('t',col(t),'x',col(x),'y',col(y));
        return
    end

    % numeric [n×3] (t,x,y) or [n×2] (x,y)
    if isnumeric(pd) && ismatrix(pd) && ~isscalar(pd)
        [t,x,y] = coerce_from_numeric(pd, tday);
        P = struct('t',col(t),'x',col(x),'y',col(y));
        return
    end

    % structs (possibly nested): try matrices first, then named fields
    if isstruct(pd)
        candidates = {pd};
        fns = fieldnames(pd);
        for k=1:numel(fns)
            if isstruct(pd.(fns{k})), candidates{end+1} = pd.(fns{k}); end
        end

        % (1) any numeric [n×2]/[n×3] inside
        for c=1:numel(candidates)
            C = candidates{c}; fn = fieldnames(C);
            for j=1:numel(fn)
                val = C.(fn{j});
                if isnumeric(val) && ismatrix(val) && ~isscalar(val)
                    try
                        [t,x,y] = coerce_from_numeric(val, tday);
                        P = struct('t',col(t),'x',col(x),'y',col(y));
                        return
                    catch
                    end
                end
            end
        end

        % (2) named fields (lots of aliases)
        t=[]; x=[]; y=[];
        for c=1:numel(candidates)
            C = candidates{c};
            f = lower(string(fieldnames(C)));
            getf = @(nm) get_field_if_exists(C, nm);
            if isempty(t), t = getf(pick_name(f, ["t","time","ts"])); end
            if isempty(x), x = getf(pick_name(f, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"])); end
            if isempty(y), y = getf(pick_name(f, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"])); end
        end
        [t,x,y] = finalize_txy(t,x,y,tday);
        P = struct('t',col(t),'x',col(x),'y',col(y));
        return
    end

    error('pos day: unsupported type %s', class(pd));
end

function [t,x,y] = coerce_from_numeric(M, tday)
    [nr,nc] = size(M);
    if nc==3
        t = M(:,1); x = M(:,2); y = M(:,3);
    elseif nc==2
        x = M(:,1); y = M(:,2); t = [];
    else
        error('numeric matrix must be [n×3] or [n×2], got [n×%d].', nc);
    end
    [t,x,y] = finalize_txy(t,x,y,tday);
end

function name = pick_name(names, options)
name=""; for k=1:numel(options), idx=find(names==options(k),1); if ~isempty(idx), name=names(idx); return; end, end
end

function val = get_field_if_exists(S, name)
if strlength(name)==0 || ~isstruct(S), val=[]; return; end
fn=fieldnames(S); idx=find(strcmpi(fn, char(name)),1); if isempty(idx), val=[]; else, val=S.(fn{idx}); end
end

function [t,x,y] = finalize_txy(t,x,y,tday)
if isempty(x)||isempty(y), error('pos day: missing x or y after coercion.'); end
x=x(:); y=y(:);
if isempty(t)
    if ~isempty(tday)&&numel(tday)==numel(x), t=tday(:); else, t=(1:numel(x))'; end
else, t=t(:);
end
n=min([numel(t),numel(x),numel(y)]); t=double(t(1:n)); x=double(x(1:n)); y=double(y(1:n));
end

function cs_vec = coerce_cs_day(csd, tday)
if iscell(csd), if isempty(csd), cs_vec=[]; return; else, csd=csd{1}; end, end
if isnumeric(csd) && isvector(csd), cs=double(csd(:)); if ~isempty(cs)&&max(cs,[],'omitnan')>1e4, cs=cs/1000; end, cs_vec=cs; return, end
if isnumeric(csd) && ismatrix(csd) && ~isscalar(csd)
    M=double(csd); if size(M,2)>=2, cs=M(:,2); else, cs=M(:,1); end
    if ~isempty(cs)&&max(cs,[],'omitnan')>1e4, cs=cs/1000; end, cs_vec=cs(:); return
end
if istable(csd)
    vn=lower(string(csd.Properties.VariableNames));
    cname=pick_name(vn, ["cs_ms","cs_time","cstime","cs","onset","onsets","cs_onset","cs_time_ms","trial_cs","cue_onset","time","ts"]);
    if strlength(cname)==0
        numMask=varfun(@isnumeric, csd, "OutputFormat","uniform"); idx=find(numMask,1);
        if isempty(idx), cs_vec=[]; return; end, tcol=double(csd{:,idx});
    else, tcol=double(csd{:, find(vn==cname,1)});
    end
    if contains(cname,"ms") || max(tcol,[],'omitnan')>1e4, tcol=tcol/1000; end
    cs_vec=tcol(:); return
end
if isstruct(csd)
    cand={csd}; fns=fieldnames(csd);
    for k=1:numel(fns), if isstruct(csd.(fns{k})), cand{end+1}=csd.(fns{k}); end, end
    for c=1:numel(cand)
        C=cand{c}; fn=fieldnames(C);
        for j=1:numel(fn)
            v=C.(fn{j});
            if isnumeric(v) && ismatrix(v) && size(v,1)>0
                if size(v,2)>=2, cs=double(v(:,2)); else, cs=double(v(:,1)); end
                if max(cs,[],'omitnan')>1e4, cs=cs/1000; end
                cs_vec=cs(:); return
            end
        end
    end
    cs_vec=[];
end
end

function [edges, ROI] = build_traceROI_edges_day(posd, tsd, csd, TraceWin, prc, marginFrac, SubGridRC)
% collect trace positions for this day
t = tsd(:); [x,y] = interp_pos(posd, t);
X=[]; Y=[];
for j=1:numel(csd)
    m = t >= csd(j)+TraceWin(1) & t <= csd(j)+TraceWin(2);
    if any(m), X=[X; x(m)]; Y=[Y; y(m)]; end
end
X = X(isfinite(X)); Y = Y(isfinite(Y));
if isempty(X) || isempty(Y)
    xmin = min(posd.x,[],'omitnan'); xmax = max(posd.x,[],'omitnan');
    ymin = min(posd.y,[],'omitnan'); ymax = max(posd.y,[],'omitnan');
else
    xmin = prctile(X, prc(1)); xmax = prctile(X, prc(2));
    ymin = prctile(Y, prc(1)); ymax = prctile(Y, prc(2));
end
dx = xmax-xmin; dy = ymax-ymin;
xmin = xmin - marginFrac*dx; xmax = xmax + marginFrac*dx;
ymin = ymin - marginFrac*dy; ymax = ymax + marginFrac*dy;

rows = SubGridRC(1); cols = SubGridRC(2);
edges.x = linspace(xmin, xmax, cols+1);
edges.y = linspace(ymin, ymax, rows+1);
ROI = struct('xmin',xmin,'xmax',xmax,'ymin',ymin,'ymax',ymax,'prc',prc,'marginFrac',marginFrac);
end

function [rc_idx, k] = pos2bin(x, y, edges)
[~, cx] = histc(x, edges.x); [~, cy] = histc(y, edges.y);
cx(cx<1 | cx>=numel(edges.x)) = NaN;
cy(cy<1 | cy>=numel(edges.y)) = NaN;
rc_idx = [cy, cx];
GridR = numel(edges.y)-1; GridC = numel(edges.x)-1;
k = nan(size(x));
m = isfinite(cx) & isfinite(cy);
k(m) = sub2ind([GridR, GridC], cy(m), cx(m));
end

function [x_i, y_i] = interp_pos(posd, t)
tt = double(posd.t(:)); xx = double(posd.x(:)); yy = double(posd.y(:));
t = double(t(:));
[ttu, ia] = unique(tt, 'stable'); xxu = xx(ia); yyu = yy(ia);
x_i = interp1(ttu, xxu, t, 'linear','extrap');
y_i = interp1(ttu, yyu, t, 'linear','extrap');
end

function S = spikes_to_matrix(daySpikes, t)
daySpikes = unwrap_spike_container(daySpikes);
if iscell(daySpikes), Nc=numel(daySpikes);
elseif isnumeric(daySpikes), Nc=size(daySpikes,1);
else, error('Unsupported daySpikes type: %s', class(daySpikes));
end
t=double(t(:));
if numel(t)<2, S=zeros(Nc, numel(t), 'single'); return; end
dt=median(diff(t)); if ~isfinite(dt)||dt<=0, dt=max(eps, mean(diff(t),'omitnan')); end
edges=[t - dt/2; t(end)+dt/2];
S=zeros(Nc, numel(t), 'single');
for c=1:Nc
    st = extract_cell_spikes(daySpikes, c);
    if isempty(st), continue; end
    S(c,:) = histcounts(st, edges);
end
end

function obj = unwrap_spike_container(S)
if isstruct(S)
    f=fieldnames(S); if isempty(f), obj=[]; return; end
    pick=[];
    for j=1:numel(f)
        name=lower(f{j});
        if contains(name,'peak') || contains(name,'spike') || contains(name,'ca_peaks'), pick=j; break, end
    end
    if isempty(pick), pick=1; end
    obj = S.(f{pick});
else
    obj=S;
end
end

function st = extract_cell_spikes(container, c)
if iscell(container), st=container{c};
elseif isnumeric(container), if c>size(container,1), st=[]; return; end, st=container(c,:).';
else, st=[]; end
st=double(st(:)); st=st(isfinite(st) & st>0);
end

function r = safe_corr(a,b)
if isempty(a)||isempty(b)||all(~isfinite(a))||all(~isfinite(b)), r=NaN; return; end
m = isfinite(a) & isfinite(b);
if nnz(m) < 3, r = NaN; return; end
r = corr(a(m), b(m), 'type','Pearson');
end

function x = bound_r(x), x = max(min(x,0.99999), -0.99999); end

function S = normalize_cells(S, mode)
switch lower(mode)
    case 'none'
    case 'demean'
        mu = mean(S,2,'omitnan'); S = S - mu;
    case 'zscore'
        mu = mean(S,2,'omitnan'); sd = std(S,0,2,'omitnan'); sd(sd==0|~isfinite(sd))=1;
        S = (S - mu) ./ sd;
    otherwise
        error('CellNorm must be ''zscore'', ''demean'', or ''none''.');
end
end

function v = speed_cm_per_s(posd)
t=double(posd.t(:)); x=double(posd.x(:)); y=double(posd.y(:));
n=min([numel(t),numel(x),numel(y)]); TXY=[t(1:n),x(1:n),y(1:n)];
try
    V = ca_velocity(TXY); v_times = double(V(2,:)).'; v_vals = double(V(1,:)).';
    [v_times_u, ia] = unique(v_times, 'stable'); v_vals_u = v_vals(ia);
    v = interp1(v_times_u, v_vals_u, t(1:n), 'linear', 'extrap');
    if any(~isfinite(v)), v = fillmissing(v,'nearest'); end
catch
    dt=diff(t(1:n)); dt(end+1,1)=median(dt(dt>0),'omitnan');
    dx=[diff(x(1:n));0]; dy=[diff(y(1:n));0]; v=hypot(dx,dy)./max(dt,eps);
end
end

function [edges, K, GridRC, ROI] = build_grid_edges_from_trialROI_day(posd, tsd, csd, TraceWin, NumBins, prc, marginFrac, GridRC_hint)
% Day-specific: build ROI from this day's trace positions and tile into NumBins.
if nargin<8 || isempty(GridRC_hint)
    GridRC = best_factors(NumBins);
else
    GridRC = GridRC_hint;
end
rows = GridRC(1); cols = GridRC(2);
K = rows*cols;

% collect positions inside the trace window for this day
t = tsd(:);
[x, y] = interp_pos(posd, t);
X = []; Y = [];
for j = 1:numel(csd)
    m = t >= csd(j)+TraceWin(1) & t <= csd(j)+TraceWin(2);
    if any(m), X = [X; x(m)]; Y = [Y; y(m)]; end %#ok<AGROW>
end
X = X(isfinite(X)); Y = Y(isfinite(Y));

if isempty(X) || isempty(Y)
    % fallback to this day's arena extents
    xmin = min(posd.x,[],'omitnan'); xmax = max(posd.x,[],'omitnan');
    ymin = min(posd.y,[],'omitnan'); ymax = max(posd.y,[],'omitnan');
else
    xmin = prctile(X, prc(1)); xmax = prctile(X, prc(2));
    ymin = prctile(Y, prc(1)); ymax = prctile(Y, prc(2));
end
dx = xmax - xmin; dy = ymax - ymin;
xmin = xmin - marginFrac*dx; xmax = xmax + marginFrac*dx;
ymin = ymin - marginFrac*dy; ymax = ymax + marginFrac*dy;

edges.x = linspace(xmin, xmax, cols+1);
edges.y = linspace(ymin, ymax, rows+1);
ROI = struct('xmin',xmin,'xmax',xmax,'ymin',ymin,'ymax',ymax,'prc',prc,'marginFrac',marginFrac);
end

function rc = best_factors(N)
% Return integer [rows cols] such that rows*cols == N and rows≈sqrt(N)
r = floor(sqrt(N));
while r > 1 && mod(N,r) ~= 0
    r = r - 1;
end
c = N / r;
rc = [r, c];
end
