function R = run_taskSpacePV_interference_TEMP(ratNames, varargin)
% RUN_TASKSPACEPV_INTERFERENCE
% End-to-end pipeline for testing Task↔Space population-vector interference.
%
% Implements three pieces:
%   (A) Task→Space interference time course:
%       r(t) = corr( PV_task(frame at t from CS), PV_control(same bin) )
%   (B) Space→Task intrusion:
%       within-bin (odd vs even) vs across-bin PV similarity during TRACE.
%   (C) Space stability with vs without task:
%       CTRL↔TASK same-bin similarity vs CTRL split-half in same bins.
%
% CONTROL spatial PVs are computed from NON-task time, excluding
%   [CS, CS+BufferPost] per trial (and optionally speed-masked).
%
% Usage:
%   R = run_taskSpacePV_interference({'rat0222','rat0307'}, 'DoPlots',false);
%
% Required fields in each rat struct (in base workspace):
%   rat.Ca_peaks, rat.Ca_ts, rat.pos, rat.CS_times
% External helpers (not defined here):
%   autoDateList(rat)            -> cellstr of day keys
%   filterFieldsByDay(struct,ks) -> pick fields by those day keys
%
% Key options (name/value pairs):
%   'DaysMode'       : {'last3toAn'|'all'|'last3byca'}         (default 'last3toAn')
%   'TraceWin'       : [t0 t1] seconds from CS (default [0 4])
%   'PrePostWin'     : [t0 t1] seconds from CS for A (default [-4 6])
%   'PreTestWin'     : [t0 t1] per-trial baseline for A (default [-2 0])
%   'BufferPost'     : seconds added after TraceWin for CTRL exclusion
%   'binSize'        : seconds per PV time bin (default 1/7.5)
%   'CellNorm'       : {'demean'|'zscore'|'none'} for PVs
%
%   Grid / ROI construction for space bins:
%   'NumBins'        : integer K (e.g., 6, 9, 12...). If empty, uses 'GridRC'
%   'GridRC'         : [rows cols] fallback when NumBins is empty
%   'UseTrialROI'    : bool. Use trace-period ROI to limit grid (default false)
%   'ROIByDay'       : bool. Build ROI per day (default true)
%   'ROIPrc'         : [lo hi] percentiles for ROI bounds
%   'ROIMarginFrac'  : fractional padding added to ROI box
%
%   CTRL occupancy & kinematics:
%   'MinOcc'         : minimum seconds per CTRL bin to keep
%   'VelThresh'      : cm/s threshold for speed masking
%   'UseSpeedMask'   : bool. Apply speed filter to samples
%   'CtrlSplitMode'  : 'interleaved' (placeholder; we use even/odd)
%
%   Intrusion (B) tuning:
%   'NPerm'              : group-level permutations
%   'MinTraceSamples'    : min total frames/bin in TRACE for B
%   'MinTraceEachHalf'   : min frames/bin in each half (odd/even)
%   'BinWeight'          : {'none'|'traceCount'|'ctrlReliability'|'traceCount*ctrlReliability'}
%   'AcrossMode'         : {'cross-halves'|'mean-of-halves'}
%
%   With/Without (C) tuning:
%   'NPerm_C'        : permutations for (C)
%   'NullMode_C'     : {'frame-redistribute'|'ctrl-derange-bins'}
%   'MinTraceSamples_C' : min TASK frames/bin for (C)
%   'MinCtrlOcc_C'      : require >= this many CTRL seconds (0 to skip)
%
% Output R (per animal):
%   R(i).A.t, .A.mean_r, .A.n, .A.ctrl_reliability, .A.pretrace (trialwise)
%   R(i).B.*  (within/across intrusion with permutation p-values)
%   R(i).C.*  (with vs without task, + distribution-level stats)
%   R(i).meta.* (options/days/grid/reliability)
%
% Recommended: keep default options first, then adjust thresholds once plots look sane.

% ------------------------------- Options ---------------------------------
p = inputParser;
addParameter(p,'DaysMode','last3toAn');
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'PrePostWin',[-4 6]);
addParameter(p,'PreTestWin',[-2 0]);
addParameter(p,'BufferPost',0);
addParameter(p,'binSize',1/7.5);
addParameter(p,'AcrossMode','cross-halves'); % 'cross-halves'|'mean-of-halves'
addParameter(p,'GridRC',[3 2]);              % fallback rows×cols if NumBins=[]
addParameter(p,'NumBins',3);                 % if given, auto factor into rows×cols
addParameter(p,'UseTrialROI',false);
addParameter(p,'ROIPrc',[5 95]);
addParameter(p,'ROIMarginFrac',0.05);
addParameter(p,'ROIByDay',true);
addParameter(p,'MinOcc',1/7.5);
addParameter(p,'VelThresh',4);
addParameter(p,'DoPlots',true);
addParameter(p,'NPerm',49);
addParameter(p,'UseSpeedMask',true);
addParameter(p,'CtrlSplitMode','interleaved');
addParameter(p,'CellNorm','demean');  % 'zscore'|'demean'|'none'
addParameter(p,'MinTraceSamples',3);
addParameter(p,'MinTraceEachHalf',1);
addParameter(p,'BinWeight','traceCount'); % 'none'|'traceCount'|'ctrlReliability'|'traceCount*ctrlReliability'
% Space-with-vs-without-task (C):
addParameter(p,'NPerm_C',49);
addParameter(p,'NullMode_C','frame-redistribute'); % 'frame-redistribute'|'ctrl-derange-bins'
addParameter(p,'MinTraceSamples_C',3);
addParameter(p,'MinCtrlOcc_C',0);
parse(p,varargin{:});
opt = p.Results;

ratNames = cellstr(ratNames);
R = struct([]);

% ========================== MAIN ANIMAL LOOP =============================
for ii = 1:numel(ratNames)
    ratVar = ratNames{ii};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratVar))
        warning('Variable %s not found in base workspace. Skipping.', ratVar);
        continue;
    end
    rat = evalin('base', ratVar);

    % ----- Select days (labels come from Ca_* fieldnames via autoDateList) -----
    dateList = autoDateList(rat); % external helper
    switch lower(opt.DaysMode)
        case {'all','last3byca'} % robust to An mismatch
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
                    daysToUse = dateList(max(1,end-2):end);
                end
            else
                daysToUse = dateList(max(1,end-2):end);
            end
        otherwise
            daysToUse = dateList;
    end
    fprintf('[%s] daysToUse: %s\n', ratVar, strjoin(daysToUse, ', '));

    % ----- Pull & standardize per-day containers (cells aligned to days) -----
    spikes_raw = filterFieldsByDay(rat.Ca_peaks, daysToUse); % external helper
    ts_raw     = filterFieldsByDay(rat.Ca_ts,    daysToUse);
    pos_raw    = filterFieldsByDay(rat.pos,      daysToUse);
    cs_raw     = filterFieldsByDay(rat.CS_times, daysToUse);
    [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, daysToUse);

    % ----- Build CONTROL spatial PVs (exclude Task+[0 BufferPost]) -----
    ctrl = compute_control_spatialPV(spikes, ts, pos, cs, ...
        'TraceWin',        opt.TraceWin, ...
        'BufferPost',      opt.BufferPost, ...
        'binSize',         opt.binSize, ...
        'GridRC',          opt.GridRC, ...
        'NumBins',         opt.NumBins, ...
        'UseTrialROI',     opt.UseTrialROI, ...
        'ROIByDay',        opt.ROIByDay, ...
        'ROIPrc',          opt.ROIPrc, ...
        'ROIMarginFrac',   opt.ROIMarginFrac, ...
        'MinOcc',          opt.MinOcc, ...
        'VelThresh',       opt.VelThresh, ...
        'UseSpeedMask',    opt.UseSpeedMask, ...
        'SplitMode',       opt.CtrlSplitMode, ...
        'CellNorm',        opt.CellNorm, ...
        'Label',           ratVar);

    % ----- (A) Task→Space interference time-course + per-trial paired test -----
    A = compute_task_to_space_timecourse(spikes, ts, pos, cs, ctrl, ...
        'PrePostWin',   opt.PrePostWin, ...
        'binSize',      opt.binSize, ...
        'VelThresh',    opt.VelThresh, ...
        'UseSpeedMask', opt.UseSpeedMask, ...
        'CellNorm',     opt.CellNorm, ...
        'TraceWin',     opt.TraceWin, ...
        'PreTestWin',   opt.PreTestWin, ...
        'Label',        ratVar);

    % Per-trial, bin-matched version (optional but very interpretable)
    A.binMatched = compute_task_to_space_trialwise_binMatched( ...
        spikes, ts, pos, cs, ctrl, ...
        'PreTestWin',   opt.PreTestWin, ...
        'TraceWin',     opt.TraceWin, ...
        'VelThresh',    opt.VelThresh, ...
        'UseSpeedMask', opt.UseSpeedMask, ...
        'CellNorm',     opt.CellNorm, ...
        'MinFramesPerBin', 2, ...
        'Label',        ratVar);

    if isfield(A,'pretrace') && ~isempty(A.pretrace.ttest)
        tt = A.pretrace.ttest;
        fprintf('[%s] Pre (%.1f–%.1fs) vs Trace (%.1f–%.1fs): paired t(%d)=%.2f, p=%.3g, nTrials=%d\n', ...
            ratVar, A.pretrace.win_pre(1), A.pretrace.win_pre(2), ...
            A.pretrace.win_trace(1), A.pretrace.win_trace(2), ...
            ifelse(isnan(tt.df), -1, tt.df), ifelse(isnan(tt.t), NaN, tt.t), ...
            ifelse(isnan(tt.p), NaN, tt.p), A.pretrace.n_pairs);
    end

    % ----- (B) Space→Task intrusion (TRACE within-bin vs across-bins) -----
    B = compute_space_to_task_intrusion(spikes, ts, pos, cs, ctrl, ...
        'TraceWin',           opt.TraceWin, ...
        'binSize',            opt.binSize, ...
        'NPerm',              opt.NPerm, ...
        'VelThresh',          opt.VelThresh, ...
        'UseSpeedMask',       opt.UseSpeedMask, ...
        'CellNorm',           opt.CellNorm, ...
        'MinTraceSamples',    opt.MinTraceSamples, ...
        'MinTraceEachHalf',   opt.MinTraceEachHalf, ...
        'BinWeight',          opt.BinWeight, ...
        'AcrossMode',         opt.AcrossMode, ...
        'Label',              ratVar);

        opt.TraceWin
        opt.VelThresh
        opt.UseSpeedMask
         opt.CellNorm
         opt.MinTraceSamples_C
         opt.MinCtrlOcc_C
         opt.NPerm_C
         opt.NullMode_C
    % ----- (C) Space stability WITH vs WITHOUT task (CTRL↔TASK vs CTRL split-half) -----
    C = compute_space_with_vs_without_task(spikes, ts, pos, cs, ctrl, ...
        'TraceWin',          opt.TraceWin, ...
        'VelThresh',         opt.VelThresh, ...
        'UseSpeedMask',      opt.UseSpeedMask, ...
        'CellNorm',          opt.CellNorm, ...
        'MinTraceSamples',   opt.MinTraceSamples_C, ...
        'MinCtrlOcc',        opt.MinCtrlOcc_C, ...
        'NPerm',             opt.NPerm_C, ...
        'NullMode',          opt.NullMode_C, ...
        'Label',             ratVar);

    % Add higher-powered per-bin stats (paired within-bin, Fisher-z)
    C = compute_space_with_without_binwise_stats(C, 'NPerm', 500);
    R(ii).C = C;

    % ----- Pack outputs per rat -----
    R(ii).animal = ratVar;
    R(ii).A = A;
    R(ii).B = B;
    R(ii).meta.options = opt;
    R(ii).meta.days = daysToUse;
    R(ii).meta.ctrl_reliability = ctrl.reliability;
    R(ii).meta.grid = ctrl.grid;
end

% ================================ PLOTS ==================================
plot_taskSpacePV_results(R); % time-course + intrusion per rat + group summary
try
    plot_spaceSimilarity_withWithout_group(R); % with vs without task (group)
catch ME
    warning('with/without task group plot failed: %s', ME.message);
end

% Optional distribution-level compare of WITHIN vs ACROSS during TRACE
stats = cdf_within_vs_across_trace(R); %#ok<NASGU>

end

% =========================================================================
% ========================== HELPER FUNCTIONS =============================
% =========================================================================

% ---------- (C) distribution-level stats: WITH vs WITHOUT task ----------
function C = compute_space_with_without_binwise_stats(C, varargin)
% Add per-day paired (within-bin) Fisher-z tests for WITH (CTRL↔TASK)
% vs WITHOUT (CTRL split-half) conditions on the SAME bins.
% Also includes:
%   - pooled paired t-test across all bins (animal-level)
%   - cluster-aware permutation p-values across days (preserve day clusters)
%
% Expects C to come from compute_space_with_vs_without_task, with:
%   C.byDay.withTask.r{d}    -> [Kx1] r per bin (NaN invalid)
%   C.byDay.withoutTask.r{d} -> [Kx1] r per bin (NaN invalid)
%
% Options:
%   'NPerm'      : # of sign-flip permutations across days (default 1000)
%   'ReportKS'   : supplemental two-sample KS on raw r (default true)
%   'ReportCliff': supplemental Cliff's delta on raw r (default true)

p = inputParser;
addParameter(p,'NPerm',1000);
addParameter(p,'ReportKS',true);
addParameter(p,'ReportCliff',true);
parse(p,varargin{:});
opt = p.Results;

if ~isfield(C,'byDay') || ~isfield(C.byDay,'withTask') || ~isfield(C.byDay,'withoutTask')
    warning('C.byDay.withTask/withoutTask not found. Nothing to do.');
    return
end

D = numel(C.byDay.withTask.r);
day = struct('nPairs',0,'t',NaN,'df',NaN,'p',NaN,'dz',NaN);

Z_with_all = []; Z_without_all = [];  % pooled across days (Fisher-z)
day_sizes = [];                       % bin counts per day (for info)

for d = 1:D
    rw = C.byDay.withTask.r{d};
    ru = C.byDay.withoutTask.r{d};
    if isempty(rw) || isempty(ru), continue; end

    % Keep only bins present/valid in BOTH
    m = isfinite(rw) & isfinite(ru);
    if ~any(m), continue; end
    rw = rw(m); ru = ru(m);

    % Fisher z transform (stable for correlations)
    zw = atanh(max(min(rw,0.99999),-0.99999));
    zu = atanh(max(min(ru,0.99999),-0.99999));
    delta = zw - zu;

    % Day-level paired test (WITH vs WITHOUT)
    if numel(delta) >= 2
        [~, p, ~, st] = ttest(zw, zu);
        dz = mean(delta) / std(delta);  % Cohen's dz for paired
        day(d) = struct('nPairs',numel(delta),'t',st.tstat,'df',st.df,'p',p,'dz',dz);
    else
        day(d) = struct('nPairs',numel(delta),'t',NaN,'df',NaN,'p',NaN,'dz',NaN);
    end

    % Accumulate for pooled paired test
    Z_with_all    = [Z_with_all; zw(:)]; %#ok<AGROW>
    Z_without_all = [Z_without_all; zu(:)]; %#ok<AGROW>
    day_sizes     = [day_sizes; numel(delta)]; %#ok<AGROW>

    % Optional supplemental raw-r descriptors
    if opt.ReportKS || opt.ReportCliff
        C.bin.byDay.stats(d).raw.n_with = numel(rw);
        C.bin.byDay.stats(d).raw.n_without = numel(ru);
        if opt.ReportKS
            [~, pks] = kstest2(rw, ru);
            C.bin.byDay.stats(d).raw.KS_p = pks;
        end
        if opt.ReportCliff
            C.bin.byDay.stats(d).raw.CliffDelta = cliff_delta(rw, ru);
        end
    end
end

% Pooled paired test across all bins (animal-level)
mask = isfinite(Z_with_all) & isfinite(Z_without_all);
if nnz(mask) >= 2
    [~, p_all, ~, st_all] = ttest(Z_with_all(mask), Z_without_all(mask));
    dz_all = mean(Z_with_all(mask)-Z_without_all(mask)) / std(Z_with_all(mask)-Z_without_all(mask));
else
    p_all = NaN; st_all = struct('tstat',NaN,'df',NaN); dz_all = NaN;
end

% Cluster-aware permutation across days: sign-flip each day’s mean Δz
NPerm = opt.NPerm;
perm_stat = nan(NPerm,1);

% Build per-day mean delta z
Dtot = numel(C.byDay.withTask.r);
day_mu = nan(Dtot,1);
for d = 1:Dtot
    rw = C.byDay.withTask.r{d}; ru = C.byDay.withoutTask.r{d};
    if isempty(rw) || isempty(ru), continue; end
    m = isfinite(rw) & isfinite(ru);
    if ~any(m), continue; end
    zw = atanh(max(min(rw(m),0.99999),-0.99999));
    zu = atanh(max(min(ru(m),0.99999),-0.99999));
    day_mu(d) = mean(zw-zu);
end
day_mu = day_mu(isfinite(day_mu));
obs_cluster_mean = mean(day_mu);

for b = 1:NPerm
    sgn = (rand(size(day_mu))>0.5)*2-1;  % random ± per day
    perm_stat(b) = mean(sgn .* day_mu);
end
p_right = mean(perm_stat >= obs_cluster_mean);
p_left  = mean(perm_stat <= obs_cluster_mean);
p_two   = 2*min(p_right, p_left);

% Pack results
C.distStats = struct();
C.distStats.perDay = day;                                  % array of per-day paired tests
C.distStats.pooled = struct('t',st_all.tstat,'df',st_all.df,'p',p_all,'dz',dz_all, ...
                            'nBins', nnz(mask));
C.distStats.cluster_perm = struct('NPerm',NPerm, ...
    'obs_meanDelta_z', obs_cluster_mean, ...
    'p_right', p_right, 'p_left', p_left, 'p_two', p_two);

% Tiny local Cliff’s delta helper (raw r)
function d = cliff_delta(x,y)
    x = x(:); y = y(:);
    x = x(isfinite(x)); y = y(isfinite(y));
    if isempty(x) || isempty(y), d = NaN; return; end
    nx = numel(x); ny = numel(y);
    gt = 0; lt = 0;
    for i=1:nx
        gt = gt + sum(x(i) > y);
        lt = lt + sum(x(i) < y);
    end
    d = (gt - lt) / (nx*ny);
end
end

% ---------- Time coercion utility ----------
function t = coerce_ts_day(td)
% Standardize a “time for day” container into a numeric vector (seconds).
% Accepts: vector, Nx2/3 numeric, table with time-like column, or struct.
% Auto-converts ms to s if values look like milliseconds.

if iscell(td), if isempty(td), t=[]; return; else, td=td{1}; end, end
if isnumeric(td) && isvector(td)
    t = double(td(:)); if max(t,[],'omitnan')>1e4, t=t/1000; end; return
end
if isnumeric(td) && ismatrix(td) && ~isscalar(td)
    M=double(td);
    if size(M,2)>=2
        t=M(:,2); if max(t,[],'omitnan')>1e4, t=t/1000; end; t=t(:); return
    else
        error('ts numeric matrix must have >=2 columns.');
    end
end
if istable(td)
    vn=lower(string(td.Properties.VariableNames));
    cname=pick_name(vn, ["time_ms","time","ts","t","timestamp","frame_ts","ca_ts"]);
    if strlength(cname)==0, error('ts table: no time-like column found'); end
    t=double(td{:,find(vn==cname,1)}); if contains(cname,"ms")||max(t,[],'omitnan')>1e4, t=t/1000; end; t=t(:); return
end
if isstruct(td)
    % heuristic dig through common shapes
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
    pref=["time_ms","time","ts","t","timestamp","timestamps","frame_ts","ca_ts"];
    for c=1:numel(cand)
        C=cand{c};
        for j=1:numel(pref)
            if isfield(C, pref{j})
                v=C.(pref{j});
                if isnumeric(v)&&isvector(v)&&numel(v)>1
                    t=double(v(:)); if contains(pref(j),"ms")||max(t,[],'omitnan')>1e4, t=t/1000; end; return
                end
            end
        end
    end
    error('ts struct: no numeric time vector or [n×2/3] matrix found.');
end
error('Unsupported ts type: %s', class(td));
end

% ---------- CTRL spatial PVs (non-task epochs) ----------
function CTRL = compute_control_spatialPV(spikes, ts, pos, cs, varargin)
% Build the spatial grid (arena-wide or trace-ROI), compute per-bin control
% PVs from non-task frames, and split-half reliability within control.
%
% Key behaviors:
%   - Excludes [CS .. CS+TraceWin(2)+BufferPost] for every trial.
%   - Optional speed masking and min occupancy gate (seconds/bin).
%   - Grid can be per-day (ROIByDay) or pooled; ROI can be trace-period ROI.

p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'BufferPost',4);
addParameter(p,'binSize',1/7.5);
addParameter(p,'GridRC',[2 3]);
addParameter(p,'NumBins',[]);
addParameter(p,'UseTrialROI',true);
addParameter(p,'ROIByDay',true);
addParameter(p,'ROIPrc',[5 95]);
addParameter(p,'ROIMarginFrac',0.05);
addParameter(p,'MinOcc',0.5);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'SplitMode','interleaved');
addParameter(p,'Label','');
addParameter(p,'CellNorm','demean');
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);

% ---- Build grid edges (per-day or pooled) ----
edgesByDay = cell(1,D);
ROIinfo    = cell(1,D);
if ~isempty(o.NumBins)
    rc_eff = best_factors(o.NumBins);
else
    rc_eff = o.GridRC;
end
K = rc_eff(1)*rc_eff(2);

if o.UseTrialROI && o.ROIByDay
    for d = 1:D
        [edgesByDay{d}, ~, ~, ROIinfo{d}] = build_grid_edges_from_trialROI_day( ...
            pos{d}, ts{d}, cs{d}, o.TraceWin, max(K,1), o.ROIPrc, o.ROIMarginFrac, rc_eff);
    end
elseif o.UseTrialROI && ~o.ROIByDay && ~isempty(o.NumBins)
    [edges1, ~, ~, ROI1] = build_grid_edges_from_trialROI( ...
        pos, ts, cs, o.TraceWin, o.NumBins, o.ROIPrc, o.ROIMarginFrac);
    edgesByDay(:) = {edges1};
    ROIinfo(:)    = {ROI1};
else
    [edges1, ~] = build_grid_edges(pos, rc_eff);
    edgesByDay(:) = {edges1};
end

% ---- CONTROL PVs per day ----
PV_by  = cell(1,D);
occ_by = cell(1,D);
rel_by = cell(1,D);

for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);

    % mask out Task + BufferPost
    is_taskbuf = false(size(t));
    csd = cs{d};
    if ~isempty(csd)
        for j = 1:numel(csd)
            is_taskbuf = is_taskbuf | (t >= csd(j)+o.TraceWin(1) & t <= csd(j)+o.TraceWin(2)+o.BufferPost);
        end
    end
    use = ~is_taskbuf;

    % speed mask (optional)
    if o.UseSpeedMask
        v = speed_cm_per_s(pos{d});
        v_i = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
        use = use & (v_i >= o.VelThresh);
    end

    % spikes to [cells × frames] matrix aligned to t
    S = spikes_to_matrix(spikes{d}, t);
    S = normalize_cells(S, o.CellNorm);

    % keep only control frames
    S   = S(:, use);
    t_u = t(use);
    x_u = x(use);
    y_u = y(use);

    % bin assignments
    [~, b] = pos2bin(x_u, y_u, edgesByDay{d});

    % occupancy (true seconds per bin)
    occ = zeros(K,1);
    if numel(t_u) >= 2
        dt = [diff(t_u); median(diff(t_u),'omitnan')];
    else
        dt = 0;
    end
    good = isfinite(b);
    if any(good)
        occ = occ + accumarray(b(good), dt(good), [K 1], @sum, 0);
    end

    % PV per bin (mean over time samples)
    Nc_d = size(S,1);
    PV   = nan(Nc_d, K);
    for k = 1:K
        sel = (b == k);
        if any(sel)
            PV(:,k) = mean(S(:,sel), 2, 'omitnan');
        end
    end

    % split-half reliability within CTRL
    rel = nan(K,1);
    sel_all = find(isfinite(b));
    for k = 1:K
        selk = sel_all(b(sel_all)==k);
        if numel(selk) < 4, continue; end
        idx_even = selk(2:2:end);
        idx_odd  = selk(1:2:end);
        v1 = mean(S(:,idx_even), 2, 'omitnan');
        v2 = mean(S(:,idx_odd ), 2, 'omitnan');
        rel(k) = safe_corr(v1, v2);
    end

    % occupancy gate
    PV(:, occ < o.MinOcc) = NaN;

    PV_by{d}  = PV;
    occ_by{d} = occ;
    rel_by{d} = rel;

    kept = all(isfinite(PV),1);
    if o.UseTrialROI && o.ROIByDay && ~isempty(ROIinfo{d})
        Rinfo = ROIinfo{d};
        fprintf(['[%s][CONTROL d=%d] Grid=%dx%d (K=%d via trialROI x=[%.2f %.2f], y=[%.2f %.2f]) | ', ...
                 'MinOcc=%.2fs | kept bins=%d/%d | split-half median r=%.2f\n'], ...
            L, d, rc_eff(1), rc_eff(2), K, Rinfo.xmin, Rinfo.xmax, Rinfo.ymin, Rinfo.ymax, ...
            o.MinOcc, nnz(kept), K, median(rel,'omitnan'));
    else
        fprintf('[%s][CONTROL d=%d] Grid=%dx%d (K=%d) | MinOcc=%.2fs | kept bins=%d/%d | split-half median r=%.2f\n', ...
            L, d, rc_eff(1), rc_eff(2), K, o.MinOcc, nnz(kept), K, median(rel,'omitnan'));
    end
end

% ---- Pack CTRL ----
CTRL.PV    = PV_by;                 % {day} : [Nc × K]
CTRL.occ   = occ_by;                % {day} : [K × 1] seconds
CTRL.grid.edges   = edgesByDay;     % {day} edges
CTRL.grid.GridRC  = rc_eff;
CTRL.grid.K       = K;
if o.UseTrialROI && o.ROIByDay
    CTRL.grid.origin = 'trialROI_byDay';
elseif o.UseTrialROI
    CTRL.grid.origin = 'trialROI_pooled';
else
    CTRL.grid.origin = 'arena';
end

CTRL.reliability_byDay = rel_by;    % {day}[K×1]
CTRL.reliability.r = vertcat(rel_by{:});
CTRL.reliability.z = atanh(max(min(CTRL.reliability.r,0.99999),-0.99999));

CTRL.params.MinOcc        = o.MinOcc;
CTRL.params.NumBins       = o.NumBins;
CTRL.params.GridRC_eff    = rc_eff;
CTRL.params.UseTrialROI   = o.UseTrialROI;
CTRL.params.ROIByDay      = o.ROIByDay;
CTRL.params.ROIPrc        = o.ROIPrc;
CTRL.params.ROIMarginFrac = o.ROIMarginFrac;
end

% ---------- (B) Space→Task intrusion ----------
function B = compute_space_to_task_intrusion(spikes, ts, pos, cs, CTRL, varargin)
% During TRACE only:
%   WITHIN: split-half (odd vs even) PV similarity within same bin.
%   ACROSS: similarity across different bins (two modes: cross-halves or mean-of-halves).
% Returns per-day & group Fisher-z means and permutation p-values on Δ (within−across).

p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'binSize',1/7.5);
addParameter(p,'NPerm',49);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'CellNorm','demean');
addParameter(p,'MinTraceSamples',6);
addParameter(p,'MinTraceEachHalf',1);
addParameter(p,'BinWeight','traceCount');     % 'none'|'traceCount'|'ctrlReliability'|'traceCount*ctrlReliability'
addParameter(p,'NullMode','derange-pairing'); % we derange valid bins
addParameter(p,'RecomputeAcrossInPerm',true);
addParameter(p,'AcrossMode','cross-halves');  % 'cross-halves'|'mean-of-halves'
addParameter(p,'Label','');
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);
K = CTRL.grid.K;

% Initialize output containers
B = struct();
B.byDay.within.z  = nan(D,1);
B.byDay.across.z  = nan(D,1);
B.byDay.within.r  = nan(D,1);
B.byDay.across.r  = nan(D,1);
B.byDay.valid     = false(D,1);
B.byDay.within.rBins = cell(D,1);  % per-bin WITHIN r for each day

% caches for permutations / reuse
PV_even_cache = cell(1,D);
PV_odd_cache  = cell(1,D);
validK_cache  = cell(1,D);
w_cache       = cell(1,D);

% keep all across-pair Fisher-z (k1,k2,z) for each day
pairs_cache = cell(1,D);

for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);
    S = spikes_to_matrix(spikes{d}, t);
    S = normalize_cells(S, o.CellNorm);

    % Speed mask
    if o.UseSpeedMask
        v = speed_cm_per_s(pos{d});
        v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
        mask_speed = (v >= o.VelThresh);
    else
        mask_speed = true(size(t));
    end

    edges_d = CTRL.grid.edges{d};
    PVd     = CTRL.PV{d};
    if isempty(PVd), continue; end

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

        % Even/odd halves (frame-level split)
        idx_even = clamp_idx(idxs(2:2:end), numel(t));
        idx_odd  = clamp_idx(idxs(1:2:end), numel(t));

        % Bin assignment for each half
        [~, b_even] = pos2bin(x(idx_even), y(idx_even), edges_d);
        [~, b_odd ] = pos2bin(x(idx_odd ), y(idx_odd ), edges_d);

        % Accumulate PV halves by bin
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

    % Gating: require both halves and CTRL PV validity
    n_tot    = n_even + n_odd;
    gate_ok  = (n_even >= o.MinTraceEachHalf) & (n_odd >= o.MinTraceEachHalf) ...
             & (n_tot >= o.MinTraceSamples);
    ctrl_ok  = all(isfinite(PVd),1)';
    validK   = find(has_k & gate_ok & ctrl_ok);
    if isempty(validK)
        fprintf('[%s][Intrusion d=%d] no valid bins.\n', L, d);
        continue;
    end

    % Bin weights for day means
    switch lower(o.BinWeight)
        case 'none',                          w = ones(K,1);
        case 'tracecount',                    w = n_tot;
        case 'ctrlreliability',               r = CTRL.reliability_byDay{d}; if isempty(r), r = zeros(K,1); end; w = max(r(:),0);
        case {'tracecount*ctrlreliability','tracecount×ctrlreliability'}
            r = CTRL.reliability_byDay{d}; if isempty(r), r = zeros(K,1); end
            w = n_tot .* max(r(:),0);
        otherwise, error('Unknown BinWeight: %s', o.BinWeight);
    end
    w(~isfinite(w)) = 0;
    w_valid = w(validK); if ~any(w_valid>0), w_valid = ones(numel(validK),1); end

    % WITHIN per-bin r (odd vs even same bin)
    zW = nan(numel(validK),1);
    rW_bins = nan(K,1);
    for iK = 1:numel(validK)
        k = validK(iK);
        r = safe_corr(PV_even(:,k), PV_odd(:,k));
        zW(iK) = atanh(max(min(r,0.99999),-0.99999));
        rW_bins(k) = r;
    end
    z_within = nansum(w_valid .* zW) / max(1, nansum(w_valid));

    % ACROSS: all bin pairs (two modes)
    zA = []; wA = []; pairs_cache{d} = [];
    if numel(validK) >= 2
        pairs = nchoosek(validK(:)',2);
        zA = nan(size(pairs,1),1);
        wA = nan(size(pairs,1),1);
        switch lower(o.AcrossMode)
            case 'cross-halves'
                for i = 1:size(pairs,1)
                    k1 = pairs(i,1); k2 = pairs(i,2);
                    r1 = safe_corr(PV_even(:,k1), PV_odd(:,k2));
                    r2 = safe_corr(PV_odd(:,k1),  PV_even(:,k2));
                    r  = mean([r1 r2],'omitnan');
                    zA(i) = atanh(max(min(r,0.99999),-0.99999));
                    wA(i) = sqrt(w(k1)*w(k2));
                    pairs_cache{d}(end+1, :) = [k1, k2, zA(i)]; %#ok<AGROW>
                end
            case 'mean-of-halves'
                v_even = PV_even; v_odd = PV_odd;
                for i = 1:size(pairs,1)
                    k1 = pairs(i,1); k2 = pairs(i,2);
                    v1 = mean([v_even(:,k1), v_odd(:,k1)],2,'omitnan');
                    v2 = mean([v_even(:,k2), v_odd(:,k2)],2,'omitnan');
                    r  = safe_corr(v1, v2);
                    zA(i) = atanh(max(min(r,0.99999),-0.99999));
                    wA(i) = sqrt(w(k1)*w(k2));
                    pairs_cache{d}(end+1, :) = [k1, k2, zA(i)]; %#ok<AGROW>
                end
            otherwise
                error('Unknown AcrossMode: %s', o.AcrossMode);
        end
        wA(~isfinite(wA)) = 0; if ~any(wA>0), wA(:)=1; end
        z_across = nansum(wA .* zA) / max(1, nansum(wA));
    else
        z_across = NaN;
    end

    % Record day summaries
    B.byDay.within.z(d) = z_within;                 B.byDay.within.r(d) = tanh(z_within);
    B.byDay.across.z(d) = z_across;                 B.byDay.across.r(d) = tanh(z_across);
    B.byDay.valid(d)    = true;
    B.byDay.within.rBins{d} = rW_bins;

    % caches for permutations
    PV_even_cache{d} = PV_even; PV_odd_cache{d} = PV_odd;
    validK_cache{d}  = validK;  w_cache{d}      = w;

    pairs_count = max(numel(validK)*(numel(validK)-1)/2, 0);
    fprintf('[%s][Intrusion d=%d] valid bins=%s | (n_even+n_odd) median=%g | across pairs=%d\n', ...
        L, d, mat2str(validK(:)'), median(n_tot(validK)), pairs_count);
end

% Group means across days
has_any = B.byDay.valid & isfinite(B.byDay.within.z) & isfinite(B.byDay.across.z);
z_within_mean = mean(B.byDay.within.z(has_any), 'omitnan');
z_across_mean = mean(B.byDay.across.z(has_any), 'omitnan');
delta_z       = z_within_mean - z_across_mean;

B.within.z = z_within_mean;    B.within.r = tanh(z_within_mean);
B.across.z = z_across_mean;    B.across.r = tanh(z_across_mean);
B.delta.z  = delta_z;          B.delta.r  = tanh(delta_z);

% Permutation test on Δz: derange across-bin pairing within day
NPerm = o.NPerm; perm_delta = nan(NPerm,1);
if NPerm > 0 && any(has_any)
    for b = 1:NPerm
        z_with_d = nan(D,1);
        z_acrs_d = nan(D,1);
        for d = 1:D
            validK = validK_cache{d};
            if numel(validK) < 2 || isempty(PV_even_cache{d}) || isempty(PV_odd_cache{d})
                continue;
            end
            % keep WITHIN as observed (tests ACROSS structure)
            z_with_d(d) = B.byDay.within.z(d);
            % ACROSS: derange validK labels
            pi = rand_derangement_idx(numel(validK));
            map = validK(pi);
            w = w_cache{d};
            wA = []; zA = [];
            switch lower(o.AcrossMode)
                case 'cross-halves'
                    for iK = 1:numel(validK)
                        k1 = validK(iK); k2 = map(iK);
                        r1 = safe_corr(PV_even_cache{d}(:,k1), PV_odd_cache{d}(:,k2));
                        r2 = safe_corr(PV_odd_cache{d}(:,k1),  PV_even_cache{d}(:,k2));
                        rr = mean([r1 r2],'omitnan');
                        zA(end+1,1) = atanh(max(min(rr,0.99999),-0.99999)); %#ok<AGROW>
                        wA(end+1,1) = sqrt(w(k1)*w(k2));                    %#ok<AGROW>
                    end
                case 'mean-of-halves'
                    for iK = 1:numel(validK)
                        k1 = validK(iK); k2 = map(iK);
                        v1 = mean([PV_even_cache{d}(:,k1), PV_odd_cache{d}(:,k1)],2,'omitnan');
                        v2 = mean([PV_even_cache{d}(:,k2), PV_odd_cache{d}(:,k2)],2,'omitnan');
                        rr = safe_corr(v1, v2);
                        zA(end+1,1) = atanh(max(min(rr,0.99999),-0.99999)); %#ok<AGROW>
                        wA(end+1,1) = sqrt(w(k1)*w(k2));                    %#ok<AGROW>
                    end
            end
            if isempty(wA) || ~any(isfinite(wA))
                z_acrs_d(d) = NaN;
            else
                wA(~isfinite(wA)) = 0; if ~any(wA>0), wA(:)=1; end
                z_acrs_d(d) = nansum(wA .* zA) / max(1, nansum(wA));
            end
        end
        m_with = mean(z_with_d(isfinite(z_with_d)), 'omitnan');
        m_acrs = mean(z_acrs_d(isfinite(z_acrs_d)), 'omitnan');
        perm_delta(b) = m_with - m_acrs;
    end
    if any(isfinite(perm_delta))
        B.p_perm_right = mean(perm_delta >= delta_z);
        B.p_perm_left  = mean(perm_delta <= delta_z);
        B.p_perm_two   = 2*min(B.p_perm_right, B.p_perm_left);
    else
        B.p_perm_right = NaN; B.p_perm_left = NaN; B.p_perm_two = NaN;
    end
else
    B.p_perm_right = NaN; B.p_perm_left = NaN; B.p_perm_two = NaN;
end

% expose caches (useful for later tables/plots)
B.cache = struct('pairs',{pairs_cache}, 'validK',{validK_cache}, 'mode',o.AcrossMode);
end

% ---------- Tiny utility: bound, dedupe, and sort indices ----------
function idx = clamp_idx(idx, N)
if isempty(idx), return; end
idx = idx(:);
idx = idx(isfinite(idx));
idx = idx(idx>=1 & idx<=N);
idx = unique(round(idx)); % ensure integer, unique, ascending
end

% ---------- Split pool into chunk sizes (greedy) ----------
function chunks = split_by_counts_ptr(pool, counts)
chunks = cell(numel(counts),1);
ptr = 1;
for j = 1:numel(counts)
    take = max(0, counts(j));
    if take==0 || ptr > numel(pool)
        chunks{j} = [];
    else
        last = min(numel(pool), ptr+take-1);
        chunks{j} = pool(ptr:last);
        ptr = last + 1;
    end
end
end

% ---------- Random derangement of 1:n (no fixed points) ----------
function pi = rand_derangement_idx(n)
if n<=1
    pi = 1:n;  % degenerate (caller should avoid n<=1 cases)
    return
end
ok = false;
while ~ok
    pi = randperm(n);
    ok = all(pi ~= 1:n);
end
end

% ---------- Terse if-else ----------
function out = ifelse(cond, a, b), if cond, out=a; else, out=b; end, end

% ---------- NaN-safe Pearson correlation with minimal length guard ----------
function r = safe_corr(a,b)
if isempty(a) || isempty(b) || all(~isfinite(a)) || all(~isfinite(b)), r = NaN; return; end
m = isfinite(a) & isfinite(b);
if nnz(m) < 3, r = NaN; return; end
r = corr(a(m), b(m), 'type','Pearson');
end

% ---------- Per-cell normalization ----------
function S = normalize_cells(S, mode)
switch lower(mode)
    case 'none'
        % leave raw counts
    case 'demean'
        mu = mean(S,2,'omitnan');
        S = S - mu;
    case 'zscore'
        mu = mean(S,2,'omitnan');
        sd = std(S,0,2,'omitnan'); sd(sd==0|~isfinite(sd)) = 1;
        S = (S - mu) ./ sd;
    otherwise
        error('CellNorm must be ''zscore'', ''demean'', or ''none''.');
end
end

% ---------- Interpolate position to timebase t ----------
function [x_i, y_i] = interp_pos(posd, t)
tt = double(posd.t(:)); xx = double(posd.x(:)); yy = double(posd.y(:));
if ~isnumeric(t), t = coerce_ts_day(t); else, t = double(t(:)); end
[ttu, ia] = unique(tt, 'stable'); xxu = xx(ia); yyu = yy(ia);
x_i = interp1(ttu, xxu, t, 'linear','extrap');
y_i = interp1(ttu, yyu, t, 'linear','extrap');
end

% ---------- Speed (cm/s) with robust fallback ----------
function v = speed_cm_per_s(posd)
t  = double(posd.t(:)); x  = double(posd.x(:)); y  = double(posd.y(:));
n  = min([numel(t), numel(x), numel(y)]);
TXY = [t(1:n), x(1:n), y(1:n)];
try
    V = ca_velocity(TXY);            % if available, use your velocity helper
    v_times = double(V(2,:)).'; v_vals = double(V(1,:)).';
    [v_times_u, ia] = unique(v_times, 'stable'); v_vals_u = v_vals(ia);
    v = interp1(v_times_u, v_vals_u, t(1:n), 'linear', 'extrap');
    if any(~isfinite(v)), v = fillmissing(v,'nearest'); end
catch
    % finite-difference fallback
    dt = diff(t(1:n)); dt(end+1,1) = median(dt(dt>0),'omitnan');
    dx = [diff(x(1:n)); 0]; dy = [diff(y(1:n)); 0];
    v = hypot(dx,dy) ./ max(dt, eps);
end
end

% ---------- Build arena-wide grid (fallback when no ROI) ----------
function [edges, K] = build_grid_edges(pos, GridRC)
allx = []; ally = [];
for d = 1:numel(pos)
    allx = [allx; pos{d}.x(:)];
    ally = [ally; pos{d}.y(:)];
end
allx = allx(isfinite(allx));  ally = ally(isfinite(ally));
if isempty(allx) || isempty(ally)
    edges.x = linspace(0, 1, GridRC(2)+1);
    edges.y = linspace(0, 1, GridRC(1)+1);
else
    edges.x = linspace(min(allx), max(allx), GridRC(2)+1);
    edges.y = linspace(min(ally), max(ally), GridRC(1)+1);
end
K = GridRC(1)*GridRC(2);
end

% ---------- Map (x,y) to linear bin index ----------
function [rc_idx, k] = pos2bin(x, y, edges)
cx = discretize(x, edges.x);  % 1..numel(edges.x)-1, or NaN
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

% ---------- Plotting: Task→Space timecourse + Intrusion ----------
function plot_taskSpacePV_results(R)
if isempty(R), warning('Empty results.'); return; end

TraceWin = R(1).meta.options.TraceWin;
tvec     = R(1).A.t; nR = numel(R);

% Per-rat timecourse & intrusion bars
for i = 1:nR
    if ~isfield(R(i),'A') || isempty(R(i).A), continue; end
    A = R(i).A; B = R(i).B;
    animal = R(i).animal;
    ctrlMean = mean(R(i).A.ctrl_reliability.r,'omitnan');

    figure('Color','w','Position',[120 120 820 560]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % (A) time course
    nexttile(1,[1 2]); hold on
    y = A.mean_r;
    yl = [-0.2 0.35];
    if any(isfinite(y))
        yl=[min(y,[],'omitnan')-0.05, max(y,[],'omitnan')+0.05];
    end
    patch([TraceWin(1) TraceWin(2) TraceWin(2) TraceWin(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.95 0.95 0.95], 'EdgeColor','none','DisplayName','Trace window');
    plot(A.t, y, 'k-', 'LineWidth', 2, 'DisplayName','corr(PV_{task}(t), PV_{ctrl})');
    yline(ctrlMean,'--','Color',[.5 .5 .5],'LineWidth',1.5,'DisplayName','ctrl split-half');
    xline(0,':','Color',[.5 .5 .5],'LineWidth',1.0,'DisplayName','CS');
    xlabel('Time from CS (s)'); ylabel('Correlation r');
    title(sprintf('[%s] Task→Space interference', animal), 'Interpreter','none');
    grid on; ylim(yl); legend('Location','best');

    % (B) intrusion bars
    nexttile(3,[1 2]); hold on
    if ~isempty(B)
        bar(1, B.within.r, 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k');
        bar(2, B.across.r, 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k');
        xticks([1 2]); xticklabels({'within (trace)','across (trace)'}); ylabel('Mean PV correlation (r)');
        txt = sprintf('\\Delta=%.3f, p_{perm}=%.3g', B.delta.r, B.p_perm_right);
        title(sprintf('[%s] Space→Task intrusion  (%s)', animal, txt), 'Interpreter','none');
        yline(0,'k:');
        ylim([min([0, B.within.r, B.across.r])-0.05, max([B.within.r, B.across.r])+0.10]);
        box on
    else
        axis off; text(0.5,0.5,'No trace bins available','HorizontalAlignment','center');
    end
end

% Group timecourse (mean ± SEM across rats)
T = numel(tvec); M = nan(nR, T); ctrlAll = nan(nR,1);
for i=1:nR
    if isempty(R(i).A), continue; end
    if numel(R(i).A.t) == T && all(R(i).A.t(:)==tvec(:))
        M(i,:) = R(i).A.mean_r(:);
        ctrlAll(i) = mean(R(i).A.ctrl_reliability.r,'omitnan');
    end
end

% Simple across-rat pre vs trace check on Fisher-z means
preMask   = (tvec >= -2 & tvec < -0);
traceMask = (tvec >=  2 & tvec <= 4);
z_group = atanh(M);
pre_z   = mean(z_group(:,preMask),  2,'omitnan');
trace_z = mean(z_group(:,traceMask),2,'omitnan');
ok = isfinite(pre_z) & isfinite(trace_z);
if any(ok)
    [~, p_preTrace, ~, st] = ttest(trace_z(ok), pre_z(ok));     % paired
    dz = mean(trace_z(ok)-pre_z(ok)) / std(trace_z(ok)-pre_z(ok)); % Cohen's dz
    fprintf('Timecourse: trace vs pre (z): t(%d)=%.2f, p=%.3g, dz=%.2f\n', st.df, st.tstat, p_preTrace, dz);
end

okR = any(isfinite(M),2) & isfinite(ctrlAll);
if any(okR)
    mu = nanmean(M(okR,:),1); se = nanstd(M(okR,:),0,1)./sqrt(nnz(okR));
    figure('Color','w','Position',[140 140 860 420]); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % group timecourse
    nexttile(1); hold on
    yl = [min(mu-se,[],'omitnan')-0.05, max(mu+se,[],'omitnan')+0.05];
    patch([TraceWin(1) TraceWin(2) TraceWin(2) TraceWin(1)], [yl(1) yl(1) yl(2) yl(2)], [0.95 0.95 0.95], 'EdgeColor','none');
    xx = [tvec; flipud(tvec)]; yy = [mu(:)-se(:); flipud(mu(:)+se(:))];
    fill(xx, yy, [0 0 0], 'FaceAlpha',0.08, 'EdgeColor','none'); plot(tvec, mu, 'k-', 'LineWidth',2.5);
    yline(mean(ctrlAll(okR)),'--','Color',[.5 .5 .5],'LineWidth',1.5); xline(0,':','Color',[.5 .5 .5]);
    xlabel('Time from CS (s)'); ylabel('Mean correlation r'); title(sprintf('Task→Space (group mean ± SEM; n=%d)', nnz(okR)));
    grid on; ylim(yl);

    % (B) per-rat per-day within vs across overlay
    nexttile(2); hold on
    bar(1, 0, 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k');   % legend color
    bar(2, 0, 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k');
    xticks([1 2]); xticklabels({'within','across'}); ylabel('Correlation r'); box on
    title('Space→Task (per-day lines, per-rat means)');

    cmap = lines(nR);
    make_light = @(c,frac) (1 - frac)*c + frac*[1 1 1];

    all_w = []; all_a = [];
    for i = 1:nR
        if ~isfield(R(i),'B') || isempty(R(i).B) || ~isfield(R(i).B,'byDay'), continue; end
        w_i = R(i).B.byDay.within.r(:);
        a_i = R(i).B.byDay.across.r(:);
        v_i = R(i).B.byDay.valid(:);
        m   = v_i & isfinite(w_i) & isfinite(a_i);
        if ~any(m), continue; end

        c_base  = cmap(i,:);
        c_light = make_light(c_base, 0.60);

        % per-day lines
        for d = find(m).'
            plot([1 2], [w_i(d) a_i(d)], '-', 'Color', c_light, 'LineWidth', 1);
            plot(1, w_i(d), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
            plot(2, a_i(d), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
        end

        % per-rat means
        w_mu = mean(w_i(m), 'omitnan');
        a_mu = mean(a_i(m), 'omitnan');
        plot([1 2], [w_mu a_mu], '-', 'Color', c_base, 'LineWidth', 2.5);
        plot(1, w_mu, 'o', 'MarkerFaceColor', c_base, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 6);
        plot(2, a_mu, 'o', 'MarkerFaceColor', c_base, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 6);

        all_w(end+1) = w_mu; %#ok<AGROW>
        all_a(end+1) = a_mu; %#ok<AGROW>
    end

    % overlay group means
    if ~isempty(all_w), bar(1, mean(all_w,'omitnan'), 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k'); end
    if ~isempty(all_a), bar(2, mean(all_a,'omitnan'), 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k'); end

    yline(0,'k:'); ylim auto

    % Paired t across animal-days on Fisher-z
    wz_days = []; az_days = [];
    for i = 1:nR
        if ~isfield(R(i),'B') || isempty(R(i).B) || ~isfield(R(i).B,'byDay'), continue; end
        Zi_w = R(i).B.byDay.within.z(:);
        Zi_a = R(i).B.byDay.across.z(:);
        v    = R(i).B.byDay.valid(:);
        m    = v & isfinite(Zi_w) & isfinite(Zi_a);
        wz_days = [wz_days; Zi_w(m)]; %#ok<AGROW>
        az_days = [az_days; Zi_a(m)]; %#ok<AGROW>
    end
    n_pairs = numel(wz_days);
    if n_pairs >= 2
        [~, p, ~, stats] = ttest(wz_days, az_days);    % paired on z
        dz = mean(wz_days-az_days) / std(wz_days-az_days);
        ax = gca;
        yl = ylim(ax); y = yl(2) + 0.05*range(yl);
        line(ax, [1 2], [y y], 'Color','k','LineWidth',1.2);
        text(1.5, y + 0.02*range(yl), ...
            sprintf('paired t across days: n=%d, t(%d)=%.2f, p=%.3g, dz=%.2f', n_pairs, stats.df, stats.tstat, p, dz), ...
            'HorizontalAlignment','center');
        ylim(ax, [yl(1) y + 0.08*range(yl)]);
    end

    % CDFs per rat of CTRL↔TASK binwise r during TRACE (panel C)
    try
        plot_spaceTask_C_traceCDF_per_rat(R);
    catch ME
        warning('CDF plot failed: %s', ME.message);
    end
end
end

% ---------- (A) Task→Space timecourse + per-trial paired test ----------
function A = compute_task_to_space_timecourse(spikes, ts, pos, cs, CTRL, varargin)
p = inputParser;
addParameter(p,'PrePostWin',[-4 16]);
addParameter(p,'TraceWin',[0 2]);      % Trial TRACE window (for per-trial summary)
addParameter(p,'PreTestWin',[-2 0]);   % Trial PRE window (paired with TRACE)
addParameter(p,'binSize',1/7.5);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'CellNorm','demean');
addParameter(p,'Label','');
parse(p,varargin{:});
o = p.Results;

% Build relative time axis for pooled time course
t_vec = (o.PrePostWin(1):o.binSize:o.PrePostWin(2))';
T = numel(t_vec);
sum_z = zeros(T,1);  cnt = zeros(T,1);

% Masks for per-trial windows (relative to CS)
isPreMask   = (t_vec >= o.PreTestWin(1) & t_vec < o.PreTestWin(2));
isTraceMask = (t_vec >= o.TraceWin(1)   & t_vec < o.TraceWin(2));

% Per-trial accumulators (Fisher-z means)
z_pre_trials   = [];
z_trace_trials = [];

dbg = struct('n_attempt',0,'n_ok',0, ...
             'n_speedFail',0,'n_oob',0,'n_ctrlMissing',0,'n_corrNaN',0);

D = numel(ts);
for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);
    S = spikes_to_matrix(spikes{d}, t);
    S = normalize_cells(S, o.CellNorm);

    if o.UseSpeedMask
        v = speed_cm_per_s(pos{d});
        v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
        use_speed = (v >= o.VelThresh);
    else
        use_speed = true(size(t));
    end

    edges_d = CTRL.grid.edges{d};
    PVd     = CTRL.PV{d};
    if isempty(PVd), continue; end

    csd = cs{d};
    for tr = 1:numel(csd)
        t0 = csd(tr);

        % tmp per-trial collectors
        z_pre_tp   = [];
        z_trace_tp = [];

        % loop over relative bins
        for i = 1:T
            win_t = t0 + t_vec(i);
            [~, idx] = min(abs(t - win_t));
            dbg.n_attempt = dbg.n_attempt + 1;

            if ~use_speed(idx), dbg.n_speedFail = dbg.n_speedFail + 1; continue; end

            [~, b] = pos2bin(x(idx), y(idx), edges_d);
            if isnan(b) || b<1 || b>size(PVd,2)
                dbg.n_oob = dbg.n_oob + 1; continue;
            end
            ctrl_pv = PVd(:,b);
            if any(~isfinite(ctrl_pv))
                dbg.n_ctrlMissing = dbg.n_ctrlMissing + 1; continue;
            end

            r = safe_corr(S(:,idx), ctrl_pv);
            if ~isfinite(r)
                dbg.n_corrNaN = dbg.n_corrNaN + 1; continue;
            end

            z = atanh(max(min(r,0.99999),-0.99999));

            % pooled time course
            sum_z(i) = sum_z(i) + z;
            cnt(i)   = cnt(i) + 1;
            dbg.n_ok = dbg.n_ok + 1;

            % per-trial assignment
            if isPreMask(i)
                z_pre_tp(end+1) = z; %#ok<AGROW>
            elseif isTraceMask(i)
                z_trace_tp(end+1) = z; %#ok<AGROW>
            end
        end

        % finalize per-trial means if both windows had data
        if ~isempty(z_pre_tp) && ~isempty(z_trace_tp)
            z_pre_trials(end+1)   = mean(z_pre_tp,   'omitnan'); %#ok<AGROW>
            z_trace_trials(end+1) = mean(z_trace_tp, 'omitnan'); %#ok<AGROW>
        end
    end
end

% Pack pooled time course
mean_z = sum_z ./ max(1,cnt);
A.t = t_vec;
A.mean_z = mean_z;
A.mean_r = tanh(mean_z);
A.n = cnt;

% carry reliability (from CTRL)
A.ctrl_reliability.r = CTRL.reliability.r;
A.ctrl_reliability.z = CTRL.reliability.z;

% Per-rat paired t-test across trials (pre > trace on Fisher-z)
A.pretrace.win_pre   = o.PreTestWin;
A.pretrace.win_trace = o.TraceWin;
A.pretrace.z_pre     = z_pre_trials(:);
A.pretrace.z_trace   = z_trace_trials(:);
A.pretrace.r_pre     = tanh(z_pre_trials(:));
A.pretrace.r_trace   = tanh(z_trace_trials(:));
A.pretrace.n_pairs   = numel(z_pre_trials);

if A.pretrace.n_pairs >= 2
  [h,p,~,stats] = ttest(A.pretrace.z_pre, A.pretrace.z_trace, 'Tail','right');
  A.pretrace.delta_z       = A.pretrace.z_pre - A.pretrace.z_trace;        % per-trial Δz
  A.pretrace.mean_delta_z  = mean(A.pretrace.delta_z,'omitnan');
  A.pretrace.cohens_dz     = A.pretrace.mean_delta_z ./ std(A.pretrace.delta_z,'omitnan');
  A.pretrace.tail = 'right';
  A.pretrace.ttest = struct('h',h,'p',p,'t',stats.tstat,'df',stats.df);
  A.pretrace.frac_trace_below_pre = mean(A.pretrace.z_trace < A.pretrace.z_pre);
else
  A.pretrace.ttest = struct('h',NaN,'p',NaN,'t',NaN,'df',NaN);
end

A.debug = dbg;
A.debug.ctrl_valid_bins = arrayfun(@(d) find(all(isfinite(CTRL.PV{d}),1)), 1:D, 'uni',0);
A.debug.ctrl_occ_sec    = CTRL.occ;
A.debug.ctrl_MinOcc     = CTRL.params.MinOcc;
A.debug.time_binSize    = o.binSize;
end

% ---------- Build ROI grid from pooled trace positions (all days) ----------
function [edges, K, GridRC, ROI] = build_grid_edges_from_trialROI(pos, ts, cs, TraceWin, NumBins, prc, marginFrac)
% Collect trace-period positions across days, bound ROI by percentiles, pad
% the box, and tile it into NumBins bins using near-square [rows cols].

% collect trace samples
X = []; Y = [];
for d = 1:numel(pos)
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);
    csd = cs{d};
    for j = 1:numel(csd)
        m = t >= csd(j)+TraceWin(1) & t <= csd(j)+TraceWin(2);
        if any(m)
            X = [X; x(m)]; %#ok<AGROW>
            Y = [Y; y(m)]; %#ok<AGROW>
        end
    end
end
X = X(isfinite(X)); Y = Y(isfinite(Y));

% Fallback to arena extents if no trace samples
if isempty(X) || isempty(Y)
    allx = []; ally = [];
    for d = 1:numel(pos), allx = [allx; pos{d}.x(:)]; ally = [ally; pos{d}.y(:)]; end
    allx = allx(isfinite(allx)); ally = ally(isfinite(ally));
    if isempty(allx) || isempty(ally)
        xmin = 0; xmax = 1; ymin = 0; ymax = 1;
    else
        xmin = min(allx); xmax = max(allx);
        ymin = min(ally); ymax = max(ally);
    end
else
    lo = prc(1); hi = prc(2);
    xmin = prctile(X, lo); xmax = prctile(X, hi);
    ymin = prctile(Y, lo); ymax = prctile(Y, hi);
end

% pad ROI
dx = xmax - xmin; dy = ymax - ymin;
xmin = xmin - marginFrac*dx; xmax = xmax + marginFrac*dx;
ymin = ymin - marginFrac*dy; ymax = ymax + marginFrac*dy;

% tile ROI into near-square grid
GridRC = best_factors(NumBins);
rows = GridRC(1); cols = GridRC(2);
edges.x = linspace(xmin, xmax, cols+1);
edges.y = linspace(ymin, ymax, rows+1);
K = rows*cols;

ROI = struct('xmin',xmin,'xmax',xmax,'ymin',ymin,'ymax',ymax,'prc',prc,'marginFrac',marginFrac);
end

% ---------- Choose near-square [rows cols] factors ----------
function rc = best_factors(N)
r = floor(sqrt(N));
while r > 1 && mod(N,r) ~= 0
    r = r - 1;
end
c = N / r;
rc = [r, c];
end

% ---------- Align raw inputs to day cells ----------
function [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, dayKeys)
% Convert various shapes (structs/cells/tables) into 1×D cell arrays
% aligned to the provided dayKeys order.

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

% ---------- Generic to-day-cells converter ----------
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
            if isempty(idx), idx = find(contains(fn, k), 1); end % relaxed match
            if isempty(idx), C{i} = []; else, C{i} = X.(fn{idx}); end
        end
        return
    else
        C = arrayfun(@(j) X(j), 1:min(numel(X), numel(keys)), 'uni', 0);
        if numel(C) < numel(keys), C(end+1:numel(keys)) = {[]}; end
        return
    end
end
% scalar fallback: replicate across days
C = repmat({X}, 1, numel(keys));
end

% ---------- Day-specific trace ROI grid ----------
function [edges, K, GridRC, ROI] = build_grid_edges_from_trialROI_day(posd, tsd, csd, TraceWin, NumBins, prc, marginFrac, GridRC_hint)
if nargin<8 || isempty(GridRC_hint)
    GridRC = best_factors(NumBins);
else
    GridRC = GridRC_hint;
end
rows = GridRC(1); cols = GridRC(2);
K = rows*cols;

% trace samples for this day
t = tsd(:);
[x, y] = interp_pos(posd, t);
X = []; Y = [];
for j = 1:numel(csd)
    m = t >= csd(j)+TraceWin(1) & t <= csd(j)+TraceWin(2);
    if any(m), X = [X; x(m)]; Y = [Y; y(m)]; end %#ok<AGROW>
end

% ROI bounds
if isempty(X) || isempty(Y)
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

% ---------- (C) WITH vs WITHOUT task ----------
function C = compute_space_with_vs_without_task(spikes, ts, pos, cs, CTRL, varargin)
% OLD/WORKING variant, now with [OLDCHK] printouts to expose gating & pairing.

p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',true);
addParameter(p,'CellNorm','demean');
addParameter(p,'MinTraceSamples',3);
addParameter(p,'MinCtrlOcc',0);
addParameter(p,'NPerm',49);
addParameter(p,'NullMode','frame-redistribute');
addParameter(p,'Label','');
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts); K = CTRL.grid.K;

r_bin_byDay   = cell(1,D);
z_day         = nan(D,1);
n_valid_day   = zeros(D,1);
r_cell_byDay  = cell(1,D);

perm_cache = struct('validK',{cell(1,D)}, 'task_idx',{cell(1,D)}, ...
                    'task_cnt',{cell(1,D)}, 'PV_ctrl',{cell(1,D)});

fprintf('[OLDCHK] Params: TraceWin=[%g %g] MinTrace=%d MinCtrlOcc=%g UseSpeedMask=%d VelThresh=%g\n', ...
    o.TraceWin, o.MinTraceSamples, o.MinCtrlOcc, o.UseSpeedMask, o.VelThresh);

for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);
    S = spikes_to_matrix(spikes{d}, t);
    S = normalize_cells(S, o.CellNorm);

    if o.UseSpeedMask
        v = speed_cm_per_s(pos{d});
        v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
        speed_ok = (v >= o.VelThresh);
    else
        speed_ok = true(size(t));
    end

    edges_d = CTRL.grid.edges{d};
    PV_ctrl = CTRL.PV{d};
    if isempty(PV_ctrl), fprintf('[OLDCHK][d=%d] PV_ctrl empty\n',d); continue; end
    Nc = size(PV_ctrl,1);

    % --- TASK PV ---
    PV_task = nan(Nc, K); n_task = zeros(K,1); task_idx_by_bin = cell(K,1);
    csd = cs{d};
    for tr = 1:numel(csd)
        tidx = (t >= csd(tr)+o.TraceWin(1) & t <= csd(tr)+o.TraceWin(2)) & speed_ok;
        ii = find(tidx);
        if isempty(ii), continue; end
        [~, b] = pos2bin(x(ii), y(ii), edges_d);
        for k = 1:K
            sel = (b == k);
            if any(sel)
                ii_k = ii(sel);
                PV_task(:,k) = nanmean([PV_task(:,k), mean(S(:, ii_k), 2, 'omitnan')], 2);
                n_task(k)    = n_task(k) + numel(ii_k);
                task_idx_by_bin{k} = [task_idx_by_bin{k}; ii_k]; %#ok<AGROW>
            end
        end
    end

    ok_ctrl = all(isfinite(PV_ctrl),1)';
    if o.MinCtrlOcc > 0 && isfield(CTRL,'occ') && numel(CTRL.occ) >= d
        ok_ctrl = ok_ctrl & (CTRL.occ{d}(:) >= o.MinCtrlOcc);
    end
    ok_task = (n_task >= o.MinTraceSamples);
    validK  = find(ok_ctrl & ok_task);

    fprintf('[OLDCHK][d=%d] n_task>0 bins=%d | ok_ctrl=%d | ok_task=%d | validK=%d\n', ...
        d, nnz(n_task>0), nnz(ok_ctrl), nnz(ok_task), numel(validK));

    n_valid_day(d) = numel(validK);
    if isempty(validK), r_bin_byDay{d} = nan(K,1); continue; end

    r_bin = nan(K,1);
    for k = validK'
        r_bin(k) = safe_corr(PV_ctrl(:,k), PV_task(:,k));
    end
    r_bin_byDay{d} = r_bin;

    zz = atanh(max(min(r_bin(validK),0.99999),-0.99999));
    z_day(d) = mean(zz,'omitnan');

    r_cell = nan(Nc,1); M = validK(:)';
    if ~isempty(M)
        for c = 1:Nc
            r_cell(c) = safe_corr(PV_ctrl(c,M).', PV_task(c,M).');
        end
    end
    r_cell_byDay{d} = r_cell;

    perm_cache.validK{d}   = validK;
    perm_cache.PV_ctrl{d}  = PV_ctrl;
    perm_cache.task_idx{d} = task_idx_by_bin;
    perm_cache.task_cnt{d} = n_task;

    fprintf('[OLDCHK][d=%d] finite r_task_bins=%d (of K=%d)\n', d, nnz(isfinite(r_bin)), K);
end

C.bin.byDay.r        = r_bin_byDay;
C.bin.byDay.z_mean   = z_day;
C.bin.byDay.n_valid  = n_valid_day;
C.bin.group.z_mean   = mean(z_day(isfinite(z_day)),'omitnan');
C.bin.group.r_mean   = tanh(C.bin.group.z_mean);

with_byDay   = cell(1,D);
without_byDay= cell(1,D);
z_ctrl_day   = nan(D,1);

for d = 1:D
    r_task_bins = r_bin_byDay{d};
    if isempty(r_task_bins), fprintf('[OLDCHK][d=%d] r_task_bins empty\n',d); continue; end
    validK = find(isfinite(r_task_bins));
    if isempty(validK), fprintf('[OLDCHK][d=%d] validK from task is empty\n',d); continue; end

    with_byDay{d} = r_task_bins(validK);

    r_ctrl = [];
    if isfield(CTRL,'reliability_byDay') && numel(CTRL.reliability_byDay) >= d ...
            && ~isempty(CTRL.reliability_byDay{d})
        r_ctrl = CTRL.reliability_byDay{d}(:);
    end
    if isempty(r_ctrl), r_ctrl = nan(CTRL.grid.K,1); end
    without_byDay{d} = r_ctrl(validK);

    fprintf('[OLDCHK][d=%d] paired bins=%d | ctrl finite over validK=%d\n', ...
        d, numel(validK), nnz(isfinite(r_ctrl(validK))));

    z_ctrl_day(d) = mean(atanh(max(min(without_byDay{d},0.99999),-0.99999)),'omitnan');
end

z_with_day = C.bin.byDay.z_mean;
z_with_anim = mean(z_with_day(isfinite(z_with_day)),'omitnan');
z_ctrl_anim = mean(z_ctrl_day(isfinite(z_ctrl_day)),'omitnan');

C.byDay.withTask.r    = with_byDay;
C.byDay.withoutTask.r = without_byDay;
C.withTask.z          = z_with_anim;
C.withTask.r          = tanh(z_with_anim);
C.withoutTask.z       = z_ctrl_anim;
C.withoutTask.r       = tanh(z_ctrl_anim);

fprintf('[OLDCHK] paired bins per day: ');
fprintf('%d ', arrayfun(@(dd) numel(C.byDay.withTask.r{dd}), 1:numel(C.byDay.withTask.r)));
fprintf('\n');

% --- permutations (unchanged from your old) ---
Nperm = o.NPerm; z_perm = nan(Nperm,1);
if Nperm > 0 && any(isfinite(z_day))
    switch lower(o.NullMode)
        case 'frame-redistribute'
            for pidx = 1:Nperm
                z_day_p = nan(D,1);
                for d = 1:D
                    validK = perm_cache.validK{d}; if isempty(validK), continue; end
                    PV_ctrl = perm_cache.PV_ctrl{d};
                    taskIdx = perm_cache.task_idx{d};
                    cnt     = perm_cache.task_cnt{d};
                    pool = [];
                    for k = validK(:).', pool = [pool; clamp_idx(taskIdx{k}, numel(ts{d}))]; end %#ok<AGROW>
                    pool = unique(pool); if isempty(pool), continue; end
                    pool = pool(randperm(numel(pool)));
                    req = max(0, round(cnt(:)));
                    req(~ismember((1:K).', validK)) = 0;
                    if sum(req) > numel(pool)
                        deficit = sum(req) - numel(pool);
                        for j = numel(req):-1:1
                            take = min(req(j), deficit);
                            req(j) = req(j) - take; deficit = deficit - take;
                            if deficit<=0, break; end
                        end
                    end
                    chunks = split_by_counts_ptr(pool, req);
                    t_day = ts{d}(:);
                    S = spikes_to_matrix(spikes{d}, t_day);
                    S = normalize_cells(S, o.CellNorm);
                    PV_task_p = nan(size(S,1), K);
                    for j = 1:numel(validK)
                        k = validK(j);
                        seg = clamp_idx(chunks{k}, size(S,2));
                        if ~isempty(seg)
                            PV_task_p(:,k) = mean(S(:, seg), 2, 'omitnan');
                        end
                    end
                    r_bin_p = nan(K,1);
                    for k = validK(:).'
                        r_bin_p(k) = safe_corr(PV_ctrl(:,k), PV_task_p(:,k));
                    end
                    z_day_p(d) = mean(atanh(max(min(r_bin_p(validK),0.99999),-0.99999)),'omitnan');
                end
                z_perm(pidx) = mean(z_day_p(isfinite(z_day_p)),'omitnan');
            end
        case 'ctrl-derange-bins'
            for pidx = 1:Nperm
                z_day_p = nan(D,1);
                for d = 1:D
                    validK = perm_cache.validK{d}; if numel(validK)<2, continue; end
                    PV_ctrl = perm_cache.PV_ctrl{d};
                    iiByBin = perm_cache.task_idx{d};
                    t_day = ts{d}(:);
                    S = spikes_to_matrix(spikes{d}, t_day);
                    S = normalize_cells(S, o.CellNorm);
                    PV_task = nan(size(PV_ctrl,1), K);
                    for j = 1:numel(validK)
                        k = validK(j);
                        seg = clamp_idx(iiByBin{k}, size(S,2));
                        if ~isempty(seg), PV_task(:,k) = mean(S(:, seg), 2, 'omitnan'); end
                    end
                    pi = rand_derangement_idx(numel(validK));
                    map = validK(pi);
                    r_bin_p = nan(K,1);
                    for j = 1:numel(validK)
                        k = validK(j); kp = map(j);
                        r_bin_p(k) = safe_corr(PV_ctrl(:,kp), PV_task(:,k));
                    end
                    z_day_p(d) = mean(atanh(max(min(r_bin_p(validK),0.99999),-0.99999)),'omitnan');
                end
                z_perm(pidx) = mean(z_day_p(isfinite(z_day_p)),'omitnan');
            end
        otherwise
            error('Unknown NullMode for C: %s', o.NullMode);
    end
end

z_obs = C.bin.group.z_mean;
if ~isempty(z_perm) && any(isfinite(z_perm))
    p_right = mean(z_perm >= z_obs);
    p_left  = mean(z_perm <= z_obs);
    p_two   = 2*min(p_right, p_left);
else
    p_right = NaN; p_left = NaN; p_two = NaN;
end

C.perm.group.z_mean   = z_perm;
C.perm.p_right        = p_right;
C.perm.p_left         = p_left;
C.perm.p_two          = p_two;
C.perm.mode           = o.NullMode;
C.perm.NPerm          = o.NPerm;

C.tests = struct();
C.tests.dayZ = z_day(:);
C.params = o; C.label = L;
C.cell.byDay.r = r_cell_byDay;
C.cell.group.r_median = median(vertcat(r_cell_byDay{:}), 'omitnan');

fprintf('[%s][C] obs z=%.3f (r=%.3f), CTRL(wo-task) z=%.3f (r=%.3f), perm mean z=%.3f (p_right=%.3g)\n', ...
    L, z_obs, tanh(z_obs), mean(z_ctrl_day,'omitnan'), tanh(mean(z_ctrl_day,'omitnan')), ...
    mean(z_perm,'omitnan'), p_right);
end


% ---------- WITH/WITHOUT group plot ----------
function plot_spaceSimilarity_withWithout_group(R)
% Bars = group means; light lines = per-rat per-day, dark lines = per-rat mean.
if isempty(R), return; end
nR = numel(R);
cmap = lines(nR);
make_light = @(c,frac) (1-frac)*c + frac*[1 1 1];

with_all    = nan(1,nR);
without_all = nan(1,nR);
with_days_all    = [];
without_days_all = [];

figure('Color','w','Position',[160 160 860 520]); hold on

for i = 1:nR
    c_base  = cmap(i,:);
    c_light = make_light(c_base, 0.60);

    with_days_i = []; without_days_i = [];
    if isfield(R(i),'C') && isfield(R(i).C,'byDay') && ...
       isfield(R(i).C.byDay,'withTask') && isfield(R(i).C.byDay.withTask,'r') && ...
       isfield(R(i).C.byDay,'withoutTask') && isfield(R(i).C.byDay.withoutTask,'r')

        wcell = R(i).C.byDay.withTask.r;
        ucell = R(i).C.byDay.withoutTask.r;

        if iscell(wcell) && iscell(ucell)
            nd = min(numel(wcell), numel(ucell));
            for d = 1:nd
                wv = wcell{d}; uv = ucell{d};
                wv = wv(isfinite(wv)); uv = uv(isfinite(uv));
                if ~isempty(wv) && ~isempty(uv)
                    with_days_i(end+1)    = mean(wv,'omitnan'); %#ok<AGROW>
                    without_days_i(end+1) = mean(uv,'omitnan'); %#ok<AGROW>
                end
            end
        end
    end

    % per-day light lines
    nd_i = min(numel(with_days_i), numel(without_days_i));
    for d = 1:nd_i
        plot([1 2], [with_days_i(d) without_days_i(d)], '-', ...
             'Color', c_light, 'LineWidth', 1.0);
        plot(1, with_days_i(d),    'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
        plot(2, without_days_i(d), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
    end

    % per-rat mean (dark)
    if nd_i >= 1
        w_mu = mean(with_days_i,'omitnan');
        u_mu = mean(without_days_i,'omitnan');
        plot([1 2], [w_mu u_mu], '-', 'Color', c_base, 'LineWidth', 2.5);
        plot(1, w_mu, 'o', 'MarkerFaceColor', c_base, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 6);
        plot(2, u_mu, 'o', 'MarkerFaceColor', c_base, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 6);
    end

    % fallback per-rat scalars if present
    if isfield(R(i),'C') && isfield(R(i).C,'withTask') && isfield(R(i).C.withTask,'r')
        with_all(i) = R(i).C.withTask.r;
    end
    if isfield(R(i),'C') && isfield(R(i).C,'withoutTask') && isfield(R(i).C.withoutTask,'r')
        without_all(i) = R(i).C.withoutTask.r;
    end

    with_days_all    = [with_days_all, with_days_i]; %#ok<AGROW>
    without_days_all = [without_days_all, without_days_i]; %#ok<AGROW>
end

% Group bars
if all(isfinite(with_all)) && all(isfinite(without_all))
    bar(1, mean(with_all,'omitnan'),    0.6, 'FaceColor',[0.35 0.7 1.0], 'EdgeColor','k');
    bar(2, mean(without_all,'omitnan'), 0.6, 'FaceColor',[0.85 0.55 0.25], 'EdgeColor','k');
else
    bar(1, mean(with_days_all,'omitnan'),    0.6, 'FaceColor',[0.35 0.7 1.0], 'EdgeColor','k');
    bar(2, mean(without_days_all,'omitnan'), 0.6, 'FaceColor',[0.85 0.55 0.25], 'EdgeColor','k');
end

% Paired t across animal-days (Fisher-z)
mask = isfinite(with_days_all) & isfinite(without_days_all);
nDays = nnz(mask);
if nDays >= 2
    zw = atanh(max(min(with_days_all(mask),   0.99999), -0.99999));
    zu = atanh(max(min(without_days_all(mask),0.99999), -0.99999));
    [~, p, ~, st] = ttest(zw, zu);
    dz = mean(zw-zu)/std(zw-zu);
    yl = ylim; y = yl(2) + 0.06*range(yl);
    line([1 2], [y y], 'Color','k','LineWidth',1.2);
    text(1.5, y + 0.02*range(yl), sprintf('paired t across days (z): t(%d)=%.2f, p=%.3g, dz=%.2f, n=%d', st.df, st.tstat, p, dz, nDays), 'HorizontalAlignment','center');
    ylim([yl(1) y + 0.10*range(yl)]);
end

xticks([1 2]); xticklabels({'with task','without task'});
ylabel('Mean PV correlation (r)');
title('Space similarity: with vs without task (observed only)');
yline(0,'k:'); grid on; box on
end

% ---------- TRACE CDFs: WITHIN vs ACROSS ----------
function stats = cdf_within_vs_across_trace(R, Nperm)
if nargin<2, Nperm = 20000; end

A = numel(R);
within_by_rat = cell(A,1);
across_by_rat = cell(A,1);

% collect per-rat distributions across days
for i = 1:A
    if ~isfield(R(i),'B') || isempty(R(i).B), continue; end
    Bi = R(i).B;

    % WITHIN: r per bin across all days
    rW = [];
    if isfield(Bi,'byDay') && isfield(Bi.byDay,'within') && iscell(Bi.byDay.within.rBins)
        for d = 1:numel(Bi.byDay.within.rBins)
            v = Bi.byDay.within.rBins{d};
            if ~isempty(v), rW = [rW; v(:)]; end %#ok<AGROW>
        end
    end
    rW = rW(isfinite(rW));

    % ACROSS: per-pair z (convert to r)
    rA = [];
    if isfield(Bi,'cache') && isfield(Bi.cache,'pairs') && ~isempty(Bi.cache.pairs)
        for d = 1:numel(Bi.cache.pairs)
            P = Bi.cache.pairs{d};
            if isempty(P), continue; end
            z = P(:,3); z = z(isfinite(z));
            rA = [rA; tanh(z(:))]; %#ok<AGROW>
        end
    end
    rA = rA(isfinite(rA));

    within_by_rat{i} = rW;
    across_by_rat{i} = rA;
end

% drop rats missing either distribution
keep = cellfun(@(x) ~isempty(x), within_by_rat) & cellfun(@(x) ~isempty(x), across_by_rat);
within_by_rat = within_by_rat(keep);
across_by_rat = across_by_rat(keep);
Aeff = numel(within_by_rat);
if Aeff==0
    warning('No rats with both within & across data.'); stats = struct(); return;
end

% pooled vectors (for quick CDF & tests)
rW_all = vertcat(within_by_rat{:});
rA_all = vertcat(across_by_rat{:});

% CDF plot
figure('Color','w','Position',[220 220 720 540]); hold on
for i = 1:Aeff
    [fW,xW] = ecdf(within_by_rat{i}); plot(xW,fW,'-','Color',[0.7 0.85 1.0], 'LineWidth',1);
    [fA,xA] = ecdf(across_by_rat{i}); plot(xA,fA,'-','Color',[1.0 0.8 0.7], 'LineWidth',1);
end
[fW,xW] = ecdf(rW_all);  plot(xW,fW,'-','Color',[0.1 0.35 0.9], 'LineWidth',2.5);
[fA,xA] = ecdf(rA_all);  plot(xA,fA,'-','Color',[0.85 0.35 0.15], 'LineWidth',2.5);
xlabel('PV correlation r (TRACE)'); ylabel('F(r \leq x)');
title(sprintf('TRACE CDFs: WITHIN vs ACROSS  (rats=%d, nW=%d, nA=%d)', Aeff, numel(rW_all), numel(rA_all)));
legend({'rat within','rat across','pooled within','pooled across'}, 'Location','southeast'); legend boxoff
grid on; xlim([-1 1]); ylim([0 1]);

% quick pooled tests (anti-conservative)
try, [~, p_ks]   = kstest2(rW_all, rA_all);   catch, p_ks = NaN; end
try, p_wil = ranksum(rW_all, rA_all, 'tail','right'); catch, p_wil = NaN; end % within > across?

% Cliff's delta on r and Δz on Fisher z
cliff = cliffs_delta(rW_all, rA_all);
zW = atanh(bound_r(rW_all)); zA = atanh(bound_r(rA_all));
delta_z = mean(zW,'omitnan') - mean(zA,'omitnan');

% hierarchical permutation across rats: sign-flip per-rat Δz
rat_stat_obs = nan(Aeff,1);
for i = 1:Aeff
    zWi = atanh(bound_r(within_by_rat{i}));
    zAi = atanh(bound_r(across_by_rat{i}));
    rat_stat_obs(i) = mean(zWi,'omitnan') - mean(zAi,'omitnan');
end
obs = mean(rat_stat_obs,'omitnan');

perm = nan(Nperm,1);
for b = 1:Nperm
    s = (rand(Aeff,1)>0.5)*2 - 1;  % sign flip per rat
    perm(b) = mean(s .* rat_stat_obs, 'omitnan');
end
p_right = mean(perm >= obs);
p_two   = 2*min(p_right, 1-p_right);

% pack
stats = struct();
stats.pooled = struct('KS_p', p_ks, 'ranksum_right_p', p_wil, ...
                      'CliffsDelta_r', cliff, 'delta_z', delta_z, ...
                      'n_within', numel(rW_all), 'n_across', numel(rA_all));
stats.hier_perm = struct('obs_delta_z', obs, 'p_right', p_right, 'p_two', p_two, ...
                         'r_diff_approx', tanh(obs), 'rats', Aeff);

fprintf('\nTRACE CDF (WITHIN vs ACROSS):\n');
fprintf('  Pooled: KS p=%.3g, Wilcoxon p_right=%.3g, CliffΔ=%.3f, Δz=%.3f (r≈%.3f)\n', ...
        p_ks, p_wil, cliff, delta_z, tanh(delta_z));
fprintf('  Hierarchical permutation: obs=%.3f (r≈%.3f), p_right=%.4g, p_two=%.4g, rats=%d\n\n', ...
        obs, tanh(obs), p_right, p_two, Aeff);

% small inner helpers
function d = cliffs_delta(x,y)
    x = x(isfinite(x)); y = y(isfinite(y));
    if isempty(x)||isempty(y), d=NaN; return; end
    nx=numel(x); ny=numel(y);
    allv=[x(:); y(:)]; [ranks,~]=tiedrank(allv);
    rx=ranks(1:nx); U = nx*ny + nx*(nx+1)/2 - sum(rx);
    d = (2*U)/(nx*ny) - 1;
end
function r = bound_r(r), r = max(min(r,0.999999),-0.999999); end
end

% ---------- Task→Space per-trial (bin-matched) ----------
function A2 = compute_task_to_space_trialwise_binMatched(spikes, ts, pos, cs, CTRL, varargin)
% Per-trial, compare PRE vs TRACE correlations to CTRL within only the bins
% visited in BOTH windows (bin-matched), then do a paired t-test across trials.

p = inputParser;
addParameter(p,'PreTestWin',[-2 0]);
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',true);
addParameter(p,'CellNorm','demean');    % 'none'|'demean'|'zscore'
addParameter(p,'MinFramesPerBin',2);    % require >= this many in BOTH windows
addParameter(p,'Label','');
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);
K = CTRL.grid.K;

% containers
z_pre_trials   = [];
z_trace_trials = [];
n_bins_trial   = [];
day_of_trial   = [];

for d = 1:D
    t  = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);
    S  = spikes_to_matrix(spikes{d}, t);
    S  = normalize_cells(S, o.CellNorm);

    if o.UseSpeedMask
        v = speed_cm_per_s(pos{d});
        v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
        speed_ok = (v >= o.VelThresh);
    else
        speed_ok = true(size(t));
    end

    edges_d  = CTRL.grid.edges{d};
    PV_ctrl  = CTRL.PV{d};            % [Nc x K]
    if isempty(PV_ctrl), continue; end
    ctrl_ok  = all(isfinite(PV_ctrl),1)';  % gate to bins with a control map

    csd = cs{d}(:);
    for tr = 1:numel(csd)
        t0 = csd(tr);

        idx_pre   = find((t >= t0+o.PreTestWin(1) & t < t0+o.PreTestWin(2)) & speed_ok);
        idx_trace = find((t >= t0+o.TraceWin(1)   & t < t0+o.TraceWin(2))   & speed_ok);
        if numel(idx_pre) < o.MinFramesPerBin || numel(idx_trace) < o.MinFramesPerBin
            continue
        end

        % bin assign and counts
        [~, b_pre]   = pos2bin(x(idx_pre),   y(idx_pre),   edges_d);
        [~, b_trace] = pos2bin(x(idx_trace), y(idx_trace), edges_d);
        pre_cnt   = accumarray(max(1,b_pre(~isnan(b_pre))),  1, [K 1], @sum, 0);
        trace_cnt = accumarray(max(1,b_trace(~isnan(b_trace))),1, [K 1], @sum, 0);

        % matched bins: CTRL-valid and sufficient frames in BOTH windows
        matched = find(ctrl_ok & pre_cnt >= o.MinFramesPerBin & trace_cnt >= o.MinFramesPerBin);
        if isempty(matched), continue; end

        z_pre_k   = nan(numel(matched),1);
        z_trace_k = nan(numel(matched),1);

        for j = 1:numel(matched)
            k = matched(j);
            seg_pre   = idx_pre(b_pre   == k);
            seg_trace = idx_trace(b_trace == k);

            if numel(seg_pre)   >= o.MinFramesPerBin
                PV_pre   = mean(S(:, seg_pre),   2, 'omitnan');
                z_pre_k(j) = atanh(max(min(safe_corr(PV_pre, PV_ctrl(:,k)), 0.99999), -0.99999));
            end
            if numel(seg_trace) >= o.MinFramesPerBin
                PV_tr    = mean(S(:, seg_trace), 2, 'omitnan');
                z_trace_k(j) = atanh(max(min(safe_corr(PV_tr, PV_ctrl(:,k)), 0.99999), -0.99999));
            end
        end

        m = isfinite(z_pre_k) & isfinite(z_trace_k);
        if ~any(m), continue; end

        z_pre_trials(end+1,1)   = mean(z_pre_k(m),   'omitnan'); %#ok<AGROW>
        z_trace_trials(end+1,1) = mean(z_trace_k(m), 'omitnan'); %#ok<AGROW>
        n_bins_trial(end+1,1)   = nnz(m);                             %#ok<AGROW>
        day_of_trial(end+1,1)   = d;                                  %#ok<AGROW>
    end
end

% pack + stats
A2 = struct();
A2.trial.z_pre     = z_pre_trials;
A2.trial.z_trace   = z_trace_trials;
A2.trial.r_pre     = tanh(z_pre_trials);
A2.trial.r_trace   = tanh(z_trace_trials);
A2.trial.delta_z   = z_pre_trials - z_trace_trials;   % positive => disruption
A2.trial.n_bins    = n_bins_trial;
A2.trial.day_id    = day_of_trial;
A2.n_trials        = numel(z_pre_trials);
A2.label           = L;
A2.params          = o;

if A2.n_trials >= 2
    [~, p, ~, st] = ttest(z_pre_trials, z_trace_trials, 'Tail','right');
    A2.stats = struct( ...
        't', st.tstat, 'df', st.df, 'p_right', p, ...
        'mean_delta_z', mean(A2.trial.delta_z,'omitnan'), ...
        'cohens_dz', mean(A2.trial.delta_z,'omitnan')/std(A2.trial.delta_z,0,'omitnan'), ...
        'frac_trace_below_pre', mean(z_trace_trials < z_pre_trials) ...
    );
else
    A2.stats = struct('t',NaN,'df',NaN,'p_right',NaN,'mean_delta_z',NaN,'cohens_dz',NaN,'frac_trace_below_pre',NaN);
end

fprintf('[%s][binMatched] n_trials=%d | mean(Δz)=%.3f, dz=%.2f | t(%g)=%.2f, p_right=%.3g | median bins/trial=%g\n', ...
    L, A2.n_trials, A2.stats.mean_delta_z, A2.stats.cohens_dz, A2.stats.df, A2.stats.t, A2.stats.p_right, ...
    median(A2.trial.n_bins,'omitnan'));
end

% ---------- Position coercion ----------
function P = coerce_pos_day(pd, tday)
% Standardize a “pos for day” container into struct with fields t,x,y.
% Accepts: cell, table, numeric [n×2] or [n×3], or nested struct.

if nargin<2, tday = []; end
col = @(v) v(:);
if iscell(pd), if isempty(pd), P=struct('t',[],'x',[],'y',[]); return; end, if numel(pd)==1, pd=pd{1}; else, pd=pd{1}; end, end
if isempty(pd)
    if ~isempty(tday), P=struct('t',col(tday),'x',nan(numel(tday),1),'y',nan(numel(tday),1));
    else, P=struct('t',[],'x',[],'y',[]); end, return
end
if istable(pd)
    vn=lower(string(pd.Properties.VariableNames));
    tname = pick_name(vn, ["t","time","ts"]);
    xname = pick_name(vn, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"]);
    yname = pick_name(vn, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"]);
    t=[]; if ~isempty(tname), t = pd{:, find(vn==tname,1)}; end
    x=[]; if ~isempty(xname), x = pd{:, find(vn==xname,1)}; end
    y=[]; if ~isempty(yname), y = pd{:, find(vn==yname,1)}; end
    [t,x,y] = finalize_txy(t,x,y,tday);
    P = struct('t',col(t),'x',col(x),'y',col(y)); return
end
if isnumeric(pd) && ismatrix(pd) && ~isscalar(pd)
    [t,x,y] = coerce_from_numeric(pd, tday);
    P = struct('t',col(t),'x',col(x),'y',col(y)); return
end
if isstruct(pd)
    candidates={pd}; fns=fieldnames(pd);
    for k=1:numel(fns), if isstruct(pd.(fns{k})), candidates{end+1}=pd.(fns{k}); end, end
    for c=1:numel(candidates)
        C=candidates{c}; fn=fieldnames(C);
        for j=1:numel(fn)
            val=C.(fn{j});
            if isnumeric(val)&&ismatrix(val)&&~isscalar(val)
                try
                    [t,x,y]=coerce_from_numeric(val, tday);
                    P=struct('t',col(t),'x',col(x),'y',col(y)); return
                catch, end
            end
        end
    end
    t=[]; x=[]; y=[];
    for c=1:numel(candidates)
        C=candidates{c}; f=lower(string(fieldnames(C)));
        getf=@(nm) get_field_if_exists(C, nm);
        if isempty(t), t=getf(pick_name(f, ["t","time","ts"])); end
        if isempty(x), x=getf(pick_name(f, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"])); end
        if isempty(y), y=getf(pick_name(f, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"])); end
    end
    [t,x,y]=finalize_txy(t,x,y,tday);
    P=struct('t',col(t),'x',col(x),'y',col(y)); return
end
error('pos day: unsupported type %s', class(pd));
end

% ---------- Numeric pos coercion ----------
function [t,x,y] = coerce_from_numeric(M, tday)
[nr,nc]=size(M);
if nc==3, t=M(:,1); x=M(:,2); y=M(:,3);
elseif nc==2, x=M(:,1); y=M(:,2); if ~isempty(tday)&&numel(tday)==nr, t=tday; else, t=(1:nr)'; end
else, error('numeric matrix must be [n×3] or [n×2], got [n×%d].', nc);
end
[t,x,y] = finalize_txy(t,x,y,tday);
end

% ---------- Pick best variable name from options ----------
function name = pick_name(names, options)
name=""; for k=1:numel(options), idx=find(names==options(k),1); if ~isempty(idx), name=names(idx); return; end, end
end

% ---------- Safe struct field get ----------
function val = get_field_if_exists(S, name)
if strlength(name)==0 || ~isstruct(S), val=[]; return; end
fn=fieldnames(S); idx=find(strcmpi(fn, char(name)),1);
if isempty(idx), val=[]; else, val=S.(fn{idx}); end
end

% ---------- Final clean of t,x,y vectors ----------
function [t,x,y] = finalize_txy(t,x,y,tday)
if isempty(x)||isempty(y), error('pos day: missing x or y after coercion.'); end
x=x(:); y=y(:);
if isempty(t)
    if ~isempty(tday)&&numel(tday)==numel(x), t=tday(:); else, t=(1:numel(x))'; end
else, t=t(:);
end
n=min([numel(t),numel(x),numel(y)]); t=double(t(1:n)); x=double(x(1:n)); y=double(y(1:n));
end

% ---------- CS time coercion ----------
function cs_vec = coerce_cs_day(csd, tday)
if iscell(csd), if isempty(csd), cs_vec=[]; return; else, csd=csd{1}; end, end
if isnumeric(csd) && isvector(csd), cs=double(csd(:)); if ~isempty(cs) && max(cs,[],'omitnan')>1e4, cs=cs/1000; end, cs_vec=cs; return, end
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
    for c=1:numel(cand)
        C=cand{c}; fn=fieldnames(C);
        for j=1:numel(fn)
            v=C.(fn{j});
            if isnumeric(v)&&isvector(v)&&numel(v)>0
                cs=double(v(:)); if max(cs,[],'omitnan')>1e4, cs=cs/1000; end, cs_vec=cs; return
            end
        end
    end
    cs_vec=[]; return
end
cs_vec=[];
end

% ---------- Spikes -> [cells × frames] aligned to t ----------
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

% ---------- Choose a spike subfield if spikes given as struct ----------
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

% ---------- Extract one cell’s spike times from generic container ----------
function st = extract_cell_spikes(container, c)
if iscell(container), st=container{c};
elseif isnumeric(container), if c>size(container,1), st=[]; return; end, st=container(c,:).';
else, st=[]; end
st=double(st(:)); st=st(isfinite(st) & st>0);
end
