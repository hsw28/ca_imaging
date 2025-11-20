function R = run_taskSpacePV_interference_DEP(ratNames, varargin)
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

addParameter(p,'TemporalBinMode','frames');   % 'time' | 'frames'
addParameter(p,'FramesPerBin',1);          % integer; if set, overrides NumTimeBins
addParameter(p,'NumTimeBins',15);           % used when TemporalBinMode='time' or FramesPerBin=[]


addParameter(p,'VelThresh',4);
addParameter(p,'DoPlots',true);
addParameter(p,'NPerm',49);
addParameter(p,'UseSpeedMask',true);
addParameter(p,'CtrlSplitMode','interleaved');
addParameter(p,'CellNorm','demean');  % 'zscore'|'demean'|'none'
addParameter(p,'MinTraceSamples',0);
addParameter(p,'MinTraceEachHalf',0);
addParameter(p,'BinWeight','traceCount'); % 'none'|'traceCount'|'ctrlReliability'|'traceCount*ctrlReliability'

% Space-with-vs-without-task (C):
addParameter(p,'NPerm_C',49);
addParameter(p,'NullMode_C','frame-redistribute'); % 'frame-redistribute'|'ctrl-derange-bins'
addParameter(p,'MinTraceSamples_C',0);
addParameter(p,'MinCtrlOcc_C',0);

% --- Across-trials (AT) options ---
addParameter(p,'DoAcrossTrials',true);  % enable/disable AT analysis
addParameter(p,'NPerm_AT',49);         % permutations for AT (per-day, optional)
addParameter(p,'MinTrials_AT',5);       % minimum trials/day to include AT for that day

%not using
addParameter(p,'NumBins',3);                 % if given, auto factor into rows×cols
addParameter(p,'UseTrialROI',false);
addParameter(p,'ROIPrc',[5 95]);
addParameter(p,'ROIMarginFrac',0.05);
addParameter(p,'ROIByDay',true);
addParameter(p,'MinOcc',1/7.5);
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
    %    fprintf('[%s] Pre (%.1f–%.1fs) vs Trace (%.1f–%.1fs): paired t(%d)=%.2f, p=%.3g, nTrials=%d\n', ...
    %        ratVar, A.pretrace.win_pre(1), A.pretrace.win_pre(2), ...
    %        A.pretrace.win_trace(1), A.pretrace.win_trace(2), ...
    %        ifelse(isnan(tt.df), -1, tt.df), ifelse(isnan(tt.t), NaN, tt.t), ...
    %        ifelse(isnan(tt.p), NaN, tt.p), A.pretrace.n_pairs);
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
    C = compute_space_with_without_binwise_stats(C, 'NPerm', 49);
    R(ii).C = C;

    % ----- (AT) Across-trials even→even: within-place vs across-place -----
    if opt.DoAcrossTrials
      AT = compute_space_to_task_acrossTrials(spikes, ts, pos, cs, ctrl, ...
              'TraceWin',[0 2], 'NPerm',49, 'MinTrials',5, ...
              'CellNorm','demean', 'Label', ratVar);

        R(ii).AT = AT;
    else
        R(ii).AT = struct();
    end


    % ----- Pack outputs per rat -----
    R(ii).animal = ratVar;
    R(ii).A = A;
    R(ii).B = B;
    R(ii).AT = AT;   % already set in the block above if DoAcrossTrials=true

    R(ii).meta.options = opt;
    R(ii).meta.days = daysToUse;
    R(ii).meta.ctrl_reliability = ctrl.reliability;
    R(ii).meta.grid = ctrl.grid;
end

% ================================ PLOTS ==================================
%plot_taskSpacePV_results(R); % time-course + intrusion per rat + group summary

summarize_taskSpacePV_clean(R);
% Across-trials (AT) plots
plot_space_to_task_acrossTrials_perAnimal(R);
plot_space_to_task_acrossTrials_group(R);
cdf_acrossTrials_within_vs_across_AT(R);


%try
%    plot_spaceSimilarity_withWithout_group(R); % with vs without task (group)
%catch ME
%    warning('with/without task group plot failed: %s', ME.message);
%end

% Optional distribution-level compare of WITHIN vs ACROSS during TRACE
%stats = cdf_within_vs_across_trace(R); %#ok<NASGU>

% Across-trials (AT) per-animal and group summaries
%try
%    plot_space_to_task_acrossTrials_perAnimal(R);
%catch ME
%    warning('Across-trials per-animal plot failed: %s', ME.message);
%end

%try
%    plot_space_to_task_acrossTrials_group(R);
%catch ME
%    warning('Across-trials group plot failed: %s', ME.message);
%end

end



% =========================================================================
% ========================== HELPER FUNCTIONS =============================
% =========================================================================

% ---------- (C) distribution-level stats: WITH vs WITHOUT task ----------
function C = compute_space_with_without_binwise_stats(C, varargin)
% Add per-day paired (within-bin) Fisher-z tests for WITH (CTRL↔TASK)
% vs WITHOUT (CTRL split-half). Also pooled paired t and cluster-aware perms.

p = inputParser;
addParameter(p,'NPerm',49);
addParameter(p,'ReportKS',true);
addParameter(p,'ReportCliff',true);
parse(p,varargin{:});
opt = p.Results;

if ~isfield(C,'byDay') || ~isfield(C.byDay,'withTask') || ~isfield(C.byDay,'withoutTask')
    warning('C.byDay.withTask/withoutTask not found. Nothing to do.'); return
end

D = numel(C.byDay.withTask.r);
dayStats = repmat(struct('nPairs',0,'t',NaN,'df',NaN,'p',NaN,'dz',NaN), D, 1);
Z_with_all = []; Z_without_all = [];

for d = 1:D
    rw = C.byDay.withTask.r{d};  ru = C.byDay.withoutTask.r{d};
    if isempty(rw) || isempty(ru), continue; end
    m = isfinite(rw) & isfinite(ru);
    if ~any(m), continue; end

    zw = atanh(max(min(rw(m), 0.99999), -0.99999));
    zu = atanh(max(min(ru(m), 0.99999), -0.99999));
    delta = zw - zu;

    if numel(delta) >= 2
        [~, p, ~, st] = ttest(zw, zu);
        dz = mean(delta)/std(delta);
        dayStats(d) = struct('nPairs',numel(delta),'t',st.tstat,'df',st.df,'p',p,'dz',dz);
    else
        dayStats(d).nPairs = numel(delta);
    end

    Z_with_all    = [Z_with_all; zw(:)]; %#ok<AGROW>
    Z_without_all = [Z_without_all; zu(:)]; %#ok<AGROW>

    % optional raw-r descriptors
    if opt.ReportKS || opt.ReportCliff
        if ~isfield(C,'bin'), C.bin = struct(); end
        if ~isfield(C.bin,'byDay'), C.bin.byDay = struct(); end
        if ~isfield(C.bin.byDay,'stats'), C.bin.byDay.stats = repmat(struct(), D, 1); end

        C.bin.byDay.stats(d).raw = struct();
        C.bin.byDay.stats(d).raw.n_with = nnz(isfinite(rw));
        C.bin.byDay.stats(d).raw.n_without = nnz(isfinite(ru));
        if opt.ReportKS
            try [~, pks] = kstest2(rw(isfinite(rw)), ru(isfinite(ru))); catch, pks = NaN; end
            C.bin.byDay.stats(d).raw.KS_p = pks;
        end
        if opt.ReportCliff
            C.bin.byDay.stats(d).raw.CliffDelta = local_cliff_delta(rw, ru);
        end
    end
end

% pooled paired test (bins pooled across days)
mask = isfinite(Z_with_all) & isfinite(Z_without_all);
if nnz(mask) >= 2
    [~, p_all, ~, st_all] = ttest(Z_with_all(mask), Z_without_all(mask));
    dz_all = mean(Z_with_all(mask) - Z_without_all(mask)) ./ std(Z_with_all(mask) - Z_without_all(mask));
else
    p_all = NaN; st_all = struct('tstat',NaN,'df',NaN); dz_all = NaN;
end

% cluster-aware permutation: sign-flip per-day mean Δz
day_mu = nan(D,1);
for d = 1:D
    rw = C.byDay.withTask.r{d};  ru = C.byDay.withoutTask.r{d};
    if isempty(rw) || isempty(ru), continue; end
    m = isfinite(rw) & isfinite(ru);
    if ~any(m), continue; end
    zw = atanh(max(min(rw(m), 0.99999), -0.99999));
    zu = atanh(max(min(ru(m), 0.99999), -0.99999));
    day_mu(d) = mean(zw - zu, 'omitnan');
end
day_mu = day_mu(isfinite(day_mu));
obs_cluster_mean = mean(day_mu,'omitnan');

NPerm = opt.NPerm;
perm_stat = nan(NPerm,1);
for b = 1:NPerm
    s = (rand(size(day_mu))>0.5)*2 - 1;
    perm_stat(b) = mean(s .* day_mu);
end
p_right = (sum(perm_stat >= obs_cluster_mean) + 1) / (nnz(isfinite(perm_stat)) + 1);
p_left  = (sum(perm_stat <= obs_cluster_mean) + 1) / (nnz(isfinite(perm_stat)) + 1);
p_two   = 2*min(p_right, p_left);

% pack
if ~isfield(C,'distStats'), C.distStats = struct(); end
C.distStats.perDay  = dayStats;
C.distStats.pooled  = struct('t',st_all.tstat,'df',st_all.df,'p',p_all,'dz',dz_all, 'nBins', nnz(mask));
C.distStats.cluster_perm = struct('NPerm',NPerm,'obs_meanDelta_z',obs_cluster_mean, ...
                                  'p_right',p_right,'p_left',p_left,'p_two',p_two);

    function d = local_cliff_delta(x,y)
        x=x(:); y=y(:); x=x(isfinite(x)); y=y(isfinite(y));
        if isempty(x)||isempty(y), d=NaN; return, end
        nx=numel(x); ny=numel(y); gt=0; lt=0;
        for i=1:nx, gt=gt+sum(x(i)>y); lt=lt+sum(x(i)<y); end
        d=(gt-lt)/(nx*ny);
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
    %    fprintf(['[%s][CONTROL d=%d] Grid=%dx%d (K=%d via trialROI x=[%.2f %.2f], y=[%.2f %.2f]) | ', ...
    %             'MinOcc=%.2fs | kept bins=%d/%d | split-half median r=%.2f\n'], ...
    %        L, d, rc_eff(1), rc_eff(2), K, Rinfo.xmin, Rinfo.xmax, Rinfo.ymin, Rinfo.ymax, ...
    %        o.MinOcc, nnz(kept), K, median(rel,'omitnan'));
    else
    %    fprintf('[%s][CONTROL d=%d] Grid=%dx%d (K=%d) | MinOcc=%.2fs | kept bins=%d/%d | split-half median r=%.2f\n', ...
    %        L, d, rc_eff(1), rc_eff(2), K, o.MinOcc, nnz(kept), K, median(rel,'omitnan'));
    end
end

% ---- Pack CTRL ----
CTRL.PV    = PV_by;                 % {day} : [Nc × K]
CTRL.occ   = occ_by;                % {day} : [K × 1] seconds
CTRL.grid.edges   = edgesByDay;     % {day} edges
CTRL.grid.GridRC  = rc_eff;
CTRL.grid.K       = K;

% NEW: store per-day distributions for CDFs
AT_byDay_within_list_z  = cell(D,1);  % Fisher-z values for WITHIN pairs
AT_byDay_across_list_z  = cell(D,1);  % Fisher-z values for ACROSS pairs


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
%BEST_FACTORS Return near-square [rows cols] integer factors of N.
% Picks the largest divisor <= sqrt(N) for rows, and cols = N/rows.

validateattributes(N, {'numeric'}, {'scalar','integer','positive','finite'});
r = floor(sqrt(double(N)));
while r > 1 && mod(N, r) ~= 0
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

function C = compute_space_with_vs_without_task(spikes, ts, pos, cs, CTRL, varargin)
% CURRENT variant, now with [NEWCHK] printouts to trace differences & recompute path.

p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',true);
addParameter(p,'CellNorm','demean');     % 'zscore'|'demean'|'none'
addParameter(p,'MinTraceSamples',0);     % min TASK frames per bin
addParameter(p,'MinCtrlOcc',0);          % extra CTRL occupancy seconds to require
addParameter(p,'NPerm',49);
addParameter(p,'NullMode','frame-redistribute');
addParameter(p,'Label','');
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts); K = CTRL.grid.K;

AT_byDay_within_list_z  = cell(D,1);
AT_byDay_across_list_z  = cell(D,1);


r_bin_byDay   = cell(1,D);
z_day         = nan(D,1);
n_valid_day   = zeros(D,1);
r_cell_byDay  = cell(1,D);

perm_cache = struct('validK',{cell(1,D)}, 'task_idx',{cell(1,D)}, ...
                    'task_cnt',{cell(1,D)}, 'PV_ctrl',{cell(1,D)});

fprintf('[NEWCHK] Params: TraceWin=[%g %g] MinTrace=%d MinCtrlOcc=%g UseSpeedMask=%d VelThresh=%g NPerm=%d\n', ...
    o.TraceWin, o.MinTraceSamples, o.MinCtrlOcc, o.UseSpeedMask, o.VelThresh, o.NPerm);

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
    if isempty(PV_ctrl), fprintf('[NEWCHK][d=%d] PV_ctrl empty\n',d); continue; end
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

    fprintf('[NEWCHK][d=%d] n_task>0 bins=%d | ok_ctrl=%d | ok_task=%d | validK=%d\n', ...
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

    fprintf('[NEWCHK][d=%d] finite r_task_bins=%d (of K=%d)\n', d, nnz(isfinite(r_bin)), K);
end

C.bin.byDay.r        = r_bin_byDay;
C.bin.byDay.z_mean   = z_day;
C.bin.byDay.n_valid  = n_valid_day;
C.bin.group.z_mean   = mean(z_day(isfinite(z_day)),'omitnan');
C.bin.group.r_mean   = tanh(C.bin.group.z_mean);

with_byDay    = cell(1,D);
without_byDay = cell(1,D);
z_ctrl_day    = nan(D,1);

for d = 1:D
    r_task_bins = C.bin.byDay.r{d};
    if isempty(r_task_bins)
        fprintf('[NEWCHK][d=%d] no task bins (r_task_bins empty)\n', d);
        continue;
    end

    % precomp reliability
    if isfield(CTRL,'reliability_byDay') && numel(CTRL.reliability_byDay) >= d ...
            && ~isempty(CTRL.reliability_byDay{d})
        r_ctrl = CTRL.reliability_byDay{d}(:); src = 'precomp';
    else
        r_ctrl = nan(CTRL.grid.K,1); src = 'none';
    end

    mask = isfinite(r_task_bins) & isfinite(r_ctrl);

    % --- recompute CTRL split-half if no overlap
    if ~any(mask)
        t  = ts{d}(:);
        [x, y] = interp_pos(pos{d}, t);
        S  = spikes_to_matrix(spikes{d}, t);
        S  = normalize_cells(S, o.CellNorm);

        is_task = false(size(t));
        if ~isempty(cs{d})
            for j = 1:numel(cs{d})
                is_task = is_task | (t >= cs{d}(j)+o.TraceWin(1) & t <= cs{d}(j)+o.TraceWin(2));
            end
        end
        use = ~is_task;

        if o.UseSpeedMask
            v = speed_cm_per_s(pos{d});
            v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
            use = use & (v >= o.VelThresh);
        end

        [~, b_ctrl] = pos2bin(x(use), y(use), CTRL.grid.edges{d});
        S_ctrl = S(:, use);

        r_ctrl2 = nan(K,1);
        for k = 1:K
            idx = find(b_ctrl == k);
            if numel(idx) < 4, continue; end
            idx_odd  = idx(1:2:end);
            idx_even = idx(2:2:end);
            v1 = mean(S_ctrl(:, idx_odd),  2, 'omitnan');
            v2 = mean(S_ctrl(:, idx_even), 2, 'omitnan');
            r  = safe_corr(v1, v2);
            if ~isfinite(r)
                n = numel(idx); h = floor(n/2);
                v1 = mean(S_ctrl(:, idx(1:h)),    2, 'omitnan');
                v2 = mean(S_ctrl(:, idx(h+1:n)),  2, 'omitnan');
                r  = safe_corr(v1, v2);
            end
            r_ctrl2(k) = r;
        end

        r_ctrl = r_ctrl2; src = 'recomputed';
        mask   = isfinite(r_task_bins) & isfinite(r_ctrl);
    end

    n_task_fin  = nnz(isfinite(r_task_bins));
    n_ctrl_fin  = nnz(isfinite(r_ctrl));
    n_paired    = nnz(mask);
    fprintf('[NEWCHK][d=%d] task finite=%d | ctrl finite=%d (%s) | paired=%d\n', ...
        d, n_task_fin, n_ctrl_fin, src, n_paired);

    if ~any(mask), continue; end
    with_byDay{d}    = r_task_bins(mask);
    without_byDay{d} = r_ctrl(mask);

    z_ctrl_day(d) = mean(atanh(max(min(without_byDay{d},0.99999),-0.99999)), 'omitnan');
end

% extra quick checks (unchanged)
for d = 1:D
    r_task_bins = C.bin.byDay.r{d};
    fprintf('[NEWCHK][chk d=%d] validK=%d | task finite=%d\n', d, C.bin.byDay.n_valid(d), nnz(isfinite(r_task_bins)));
    if isfield(CTRL,'reliability_byDay') && numel(CTRL.reliability_byDay) >= d && ~isempty(CTRL.reliability_byDay{d})
        fprintf('                   ctrl finite(precomp)=%d\n', nnz(isfinite(CTRL.reliability_byDay{d})));
    else
        fprintf('                   ctrl finite(precomp)=0 (missing)\n');
    end
end

C.byDay.withTask.r    = with_byDay;
C.byDay.withoutTask.r = without_byDay;

z_with_day  = C.bin.byDay.z_mean;
z_with_anim = mean(z_with_day(isfinite(z_with_day)),'omitnan');
z_ctrl_anim = mean(z_ctrl_day(isfinite(z_ctrl_day)),'omitnan');

C.withTask.z    = z_with_anim;     C.withTask.r    = tanh(z_with_anim);
C.withoutTask.z = z_ctrl_anim;     C.withoutTask.r = tanh(z_ctrl_anim);

fprintf('[NEWCHK] paired bins per day: ');
fprintf('%d ', arrayfun(@(dd) numel(C.byDay.withTask.r{dd}), 1:numel(C.byDay.withTask.r)));
fprintf('\n');

% --- permutations (same logic; default NPerm=49 here) ---
Nperm = o.NPerm;
z_perm = nan(Nperm,1);
if Nperm > 0 && any(isfinite(z_day))
    switch lower(o.NullMode)
        case 'frame-redistribute'
            for pidx = 1:Nperm
                z_day_p = nan(D,1);
                for d = 1:D
                    validK = perm_cache.validK{d};
                    if isempty(validK), continue; end
                    PV_ctrl = perm_cache.PV_ctrl{d};
                    taskIdx = perm_cache.task_idx{d};
                    cnt     = perm_cache.task_cnt{d};

                    pool = [];
                    for k = validK(:).', pool = [pool; clamp_idx(taskIdx{k}, numel(ts{d}))]; end %#ok<AGROW>
                    pool = unique(pool);
                    if isempty(pool), continue; end
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
                    validK = perm_cache.validK{d};
                    if numel(validK) < 2, continue; end
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
fprintf('[NEWCHK] perm for trace disrupts space AGAINST NULL\n');
if ~isempty(z_perm) && any(isfinite(z_perm))
    p_right = mean(z_perm >= z_obs);
    p_left  = mean(z_perm <= z_obs);
    p_two   = 2*min(p_right, p_left);  % keep semicolon if you want silence
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

end



% ---------- helper: split-half CTRL from frames (even vs odd) ----------
function [r_ctrl, src] = local_ctrl_split_half_from_frames(daySpikes, t, posd, csd, edges_d, ...
        cellNorm, useSpeed, velThresh, excludeTrace, traceWin, K)

r_ctrl = [];
src = 'none';

t = t(:);
[x, y] = interp_pos(posd, t);
S = spikes_to_matrix(daySpikes, t);
S = normalize_cells(S, cellNorm);

if useSpeed
    v = speed_cm_per_s(posd);
    v = interp1(posd.t(:), v(:), t, 'linear','extrap');
    speed_ok = (v >= velThresh);
else
    speed_ok = true(size(t));
end

mask = speed_ok;
if excludeTrace && ~isempty(csd)
    inTrace = false(size(t));
    csd = csd(:);
    for j = 1:numel(csd)
        inTrace = inTrace | (t >= csd(j)+traceWin(1) & t <= csd(j)+traceWin(2));
    end
    mask = mask & ~inTrace;
end

frames = find(mask);
if numel(frames) < 4
    src = 'split-half:insufficient_frames';
    return
end

fe = frames(2:2:end);
fo = frames(1:2:end);
if isempty(fe) || isempty(fo)
    src = 'split-half:no_even_or_odd';
    return
end

% build PVs for even and odd halves
PV_e = nan(size(S,1), K);
PV_o = nan(size(S,1), K);

[~, bE] = pos2bin(x(fe), y(fe), edges_d);
[~, bO] = pos2bin(x(fo), y(fo), edges_d);

for k = 1:K
    selE = (bE == k);
    selO = (bO == k);
    if any(selE)
        PV_e(:,k) = mean(S(:, fe(selE)), 2, 'omitnan');
    end
    if any(selO)
        PV_o(:,k) = mean(S(:, fo(selO)), 2, 'omitnan');
    end
end

% split-half r per bin
r_ctrl = nan(K,1);
for k = 1:K
    r_ctrl(k) = safe_corr(PV_e(:,k), PV_o(:,k));
end

src = sprintf('split-half:%s', ternary(excludeTrace,'nonTRACE','allframes'));
end

% tiny utility
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


function plot_space_with_vs_without_task_group(R)
% Build per-animal and per-animal-day summaries from C.tests.dayZ and C.perm.group.z_mean.

% --- collect per-animal-day z and mapping ---
animalNames = {R.animal};
nA = numel(R);

dayZ_all = [];     % [Ndays_total x 1]
anim_all = [];     % indices mapping day -> animal
perm_means = [];   % per-animal perm mean (just for display pairing)

for i = 1:nA
    if ~isfield(R(i),'C') || isempty(R(i).C), continue; end
    z_d = R(i).C.tests.dayZ;
    z_d = z_d(isfinite(z_d));
    dayZ_all = [dayZ_all; z_d]; %#ok<AGROW>
    anim_all = [anim_all; i*ones(numel(z_d),1)]; %#ok<AGROW>
end

% per-animal observed means (across its days)
obs_animZ = accumarray(anim_all, dayZ_all, [nA 1], @mean, NaN);
obs_animR = tanh(obs_animZ);

% group permutation pool (already group-mean over days; we’ll display its mean)
perm_group_z = [];
for i = 1:nA
    if ~isfield(R(i),'C') || isempty(R(i).C) || isempty(R(i).C.perm.group.z_mean), continue; end
    perm_group_z = [perm_group_z; R(i).C.perm.group.z_mean(:)]; %#ok<AGROW>
end
perm_group_mean_z = mean(perm_group_z,'omitnan');
perm_group_mean_r = tanh(perm_group_mean_z);

% --- t-tests (vs 0) ---
% Across animals (n = #animals with data)
maskA = isfinite(obs_animZ);
p_anim = NaN; df_anim = NaN; t_anim = NaN;
if any(maskA)
    [~, p_anim, ~, st] = ttest(obs_animZ(maskA), 0);   % Fisher-z vs 0
    df_anim = st.df; t_anim = st.tstat;
end

% Across animal-days (n = total valid days)
maskD = isfinite(dayZ_all);
p_day = NaN; df_day = NaN; t_day = NaN;
if any(maskD)
    [~, p_day, ~, st2] = ttest(dayZ_all(maskD), 0);
    df_day = st2.df; t_day = st2.tstat;
end

% --- Plot (bars + lines) ---
figure('Color','w','Position',[160 160 860 480]); hold on

% Positions
xObs = 1; xNull = 2;

% Light lines for animal-days: connect Obs to PermMean (same for all; draw toward group perm)
% We plot observed (day) against group perm mean as a reference:
for i = 1:nA
    z_d = []; if isfield(R(i),'C') && ~isempty(R(i).C), z_d = R(i).C.tests.dayZ; end
    z_d = z_d(isfinite(z_d));
    for j = 1:numel(z_d)
        plot([xObs xNull], [tanh(z_d(j)), perm_group_mean_r], '-', 'Color',[0.6 0.6 0.9], 'LineWidth',0.8);
    end
end

% Dark lines for per-animal averages:
for i = 1:nA
    if ~isfinite(obs_animR(i)), continue; end
    plot([xObs xNull], [obs_animR(i), perm_group_mean_r], '-', 'Color',[0.1 0.2 0.8], 'LineWidth',2.0);
end

% Bars: group means
bar(xObs, mean(obs_animR(isfinite(obs_animR)),'omitnan'), 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k');
bar(xNull, perm_group_mean_r,                          0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k');

% Cosmetics
xlim([0.5 2.5]); xticks([xObs xNull]); xticklabels({'Observed (CTRL↔TASK)', 'Permuted mean'});
ylabel('Spatial PV similarity r'); yline(0,'k:');
title(sprintf('Space stability with vs without task  |  t_anim=%.2f (df=%g, p=%.3g);  t_day=%.2f (df=%g, p=%.3g)', ...
    t_anim, df_anim, p_anim, t_day, df_day, p_day));

box on; grid on;

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
if nargin<2, Nperm = 49; end

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

%fprintf('[%s][binMatched] n_trials=%d | mean(Δz)=%.3f, dz=%.2f | t(%g)=%.2f, p_right=%.3g | median bins/trial=%g\n', ...
%    L, A2.n_trials, A2.stats.mean_delta_z, A2.stats.cohens_dz, A2.stats.df, A2.stats.t, A2.stats.p_right, ...
%    median(A2.trial.n_bins,'omitnan'));
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







% ---------- (B) Space→Task intrusion ----------
% ---------- (B) Space→Task intrusion ----------
function B = compute_space_to_task_intrusion(spikes, ts, pos, cs, CTRL, varargin)
% Compute within-bin vs across-bin task PV similarity and permutation tests.
% Default "AcrossMode" is 'cross-halves' (symmetric odd–even across bins).

p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'binSize',1/7.5);
addParameter(p,'NPerm',49);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'CellNorm','demean');          % 'zscore'|'demean'|'none'
addParameter(p,'MinTraceSamples',6);
addParameter(p,'MinTraceEachHalf',1);
addParameter(p,'BinWeight','traceCount');     % 'none'|'traceCount'|'ctrlReliability'|'traceCount*ctrlReliability'
addParameter(p,'NullMode','derange-pairing'); % 'derange-pairing'|'frame-redistribute'
addParameter(p,'RecomputeAcrossInPerm',true);
addParameter(p,'AcrossMode','cross-halves');  % 'cross-halves'|'mean-of-halves'
addParameter(p,'Label','');
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);
K = CTRL.grid.K;

% ---- Container scaffolding ----
B = struct();
B.byDay = struct();
B.byDay.within = struct('r',[],'z',[],'rBins',[]);   % rBins will be {D×1} cells
B.byDay.across = struct('r',[],'z',[]);
B.byDay.valid  = false(D,1);
B.byDay.within.rBins = cell(D,1);                   % <- correct per-day cell init

B.perm = struct();
B.perm.byDay = struct('delta', {cell(1,D)});        % cell{d} -> [Nperm×1] Δz for that day

z_within_day = nan(D,1);
z_across_day = nan(D,1);
has_any      = false(D,1);

% caches for permutations
PV_even_cache = cell(1,D);
PV_odd_cache  = cell(1,D);
validK_cache  = cell(1,D);
w_cache       = cell(1,D);

% frame index caches (for frame-redistribute)
even_idx_cache = cell(1,D);
odd_idx_cache  = cell(1,D);
n_even_cache   = cell(1,D);
n_odd_cache    = cell(1,D);

for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);
    S = spikes_to_matrix(spikes{d}, t);
    S = normalize_cells(S, o.CellNorm);

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

    even_idx_by_bin = cell(K,1);
    odd_idx_by_bin  = cell(K,1);

    csd = cs{d};
    for tr = 1:numel(csd)
        tidx = t >= csd(tr)+o.TraceWin(1) & t <= csd(tr)+o.TraceWin(2);
        tidx = tidx & mask_speed;
        idxs = find(tidx);
        if numel(idxs) < 2, continue; end

        idx_even = clamp_idx(idxs(2:2:end), numel(t));
        idx_odd  = clamp_idx(idxs(1:2:end), numel(t));

        [~, b_even] = pos2bin(x(idx_even), y(idx_even), edges_d);
        [~, b_odd ] = pos2bin(x(idx_odd ), y(idx_odd ), edges_d);

        for k = 1:K
            if any(b_even==k)
                ie = idx_even(b_even==k);
                if ~isempty(ie)
                    PV_even(:,k) = nanmean([PV_even(:,k), mean(S(:,ie),2)],2);
                    n_even(k) = n_even(k) + numel(ie);
                    has_k(k) = true;
                    even_idx_by_bin{k} = [even_idx_by_bin{k}; ie]; %#ok<AGROW>
                end
            end
            if any(b_odd==k)
                io = idx_odd(b_odd==k);
                if ~isempty(io)
                    PV_odd(:,k)  = nanmean([PV_odd(:,k),  mean(S(:,io),2)],2);
                    n_odd(k)  = n_odd(k) + numel(io);
                    has_k(k) = true;
                    odd_idx_by_bin{k} = [odd_idx_by_bin{k}; io]; %#ok<AGROW>
                end
            end
        end
    end

    % gate bins by trace samples + control PV availability
    n_tot = n_even + n_odd;
    gate_counts = (n_even >= o.MinTraceEachHalf) & (n_odd >= o.MinTraceEachHalf) & (n_tot >= o.MinTraceSamples);
    validK = find(has_k & gate_counts & all(isfinite(PVd),1)');
    if isempty(validK)
        fprintf('[%s][Intrusion d=%d] no valid bins.\n', L, d);
        continue;
    end
    has_any(d) = true;

    % bin weights
    switch lower(o.BinWeight)
        case 'none',                               w = ones(K,1);
        case 'tracecount',                         w = n_tot;
        case 'ctrlreliability',                    r = CTRL.reliability_byDay{d}; if isempty(r), r = zeros(K,1); end, w = max(r(:),0);
        case {'tracecount*ctrlreliability','tracecount×ctrlreliability'}
                                                   r = CTRL.reliability_byDay{d}; if isempty(r), r = zeros(K,1); end, w = n_tot .* max(r(:),0);
        otherwise, error('Unknown BinWeight: %s', o.BinWeight);
    end
    w(~isfinite(w)) = 0;
    w_valid = w(validK); if ~any(w_valid>0), w_valid = ones(numel(validK),1); end

    % observed WITHIN (Fisher-z weighted mean)
    zW = nan(numel(validK),1);
    for iK = 1:numel(validK)
        k = validK(iK);
        zW(iK) = atanh(max(min(safe_corr(PV_even(:,k), PV_odd(:,k)),0.99999),-0.99999));
    end
    z_within_day(d) = nansum(w_valid .* zW) / max(1, nansum(w_valid));

    % store per-bin r for this day
    rW_bins = nan(K,1);
    for iK = 1:numel(validK)
        k = validK(iK);
        rW_bins(k) = safe_corr(PV_even(:,k), PV_odd(:,k));
    end
    B.byDay.within.rBins{d} = rW_bins;

    % observed ACROSS
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
                    zA(i) = mean([atanh(max(min(r1,0.99999),-0.99999)), atanh(max(min(r2,0.99999),-0.99999))], 'omitnan');
                    wA(i) = sqrt(w(k1)*w(k2));
                end
            case 'mean-of-halves'
                for i = 1:size(pairs,1)
                    k1 = pairs(i,1); k2 = pairs(i,2);
                    v1 = mean([PV_even(:,k1), PV_odd(:,k1)],2,'omitnan');
                    v2 = mean([PV_even(:,k2), PV_odd(:,k2)],2,'omitnan');
                    zA(i) = atanh(max(min(safe_corr(v1, v2),0.99999),-0.99999));
                    wA(i) = sqrt(w(k1)*w(k2));
                end
            otherwise
                error('Unknown AcrossMode: %s', o.AcrossMode);
        end
        wA(~isfinite(wA)) = 0; if ~any(wA>0), wA(:)=1; end
        z_across_day(d) = nansum(wA .* zA) / max(1, nansum(wA));
    else
        z_across_day(d) = NaN;
    end

    % caches
    PV_even_cache{d} = PV_even;
    PV_odd_cache{d}  = PV_odd;
    validK_cache{d}  = validK;
    w_cache{d}       = w;

    even_idx_cache{d} = even_idx_by_bin;
    odd_idx_cache{d}  = odd_idx_by_bin;
    n_even_cache{d}   = n_even;
    n_odd_cache{d}    = n_odd;
end

% mean across days (Fisher z)
m_valid        = isfinite(z_within_day) & (B.byDay.valid | true); %#ok<NASGU>
z_within_mean  = mean(z_within_day(isfinite(z_within_day)),'omitnan');
z_across_mean  = mean(z_across_day(isfinite(z_across_day)),'omitnan');
delta_z        = z_within_mean - z_across_mean;

% ------- permutations -------
Nperm = o.NPerm;
perm_delta        = nan(Nperm,1);
zW_perm_mean_all  = nan(Nperm,1);
zA_perm_mean_all  = nan(Nperm,1);

if Nperm>0 && any(isfinite(z_within_day))
    switch lower(o.NullMode)

    case 'derange-pairing'
        for pidx = 1:Nperm
            zW_perm_day = nan(D,1);
            zA_perm_day = nan(D,1);
            for d = 1:D
                validK = validK_cache{d};
                if isempty(validK), continue; end

                w = w_cache{d};
                w_valid = w(validK); if ~any(w_valid>0), w_valid(:)=1; end
                PV_even = PV_even_cache{d};
                PV_odd  = PV_odd_cache{d};

                % derange the mapping O<-E across bins
                Kvalid = numel(validK);
                if Kvalid < 1, continue; end
                pi_idx = rand_derangement_idx(Kvalid);
                map_to = validK(pi_idx);

                % WITHIN under null
                zP = nan(Kvalid,1); wP = nan(Kvalid,1);
                for j = 1:Kvalid
                    k1 = validK(j); k2 = map_to(j);
                    zP(j) = atanh(max(min(safe_corr(PV_even(:,k1), PV_odd(:,k2)),0.99999),-0.99999));
                    wP(j) = sqrt(w(k1)*w(k2));
                end
                wP(~isfinite(wP))=0; if ~any(wP>0), wP(:)=1; end
                zW_perm_day(d) = nansum(wP .* zP) / max(1,nansum(wP));

                % ACROSS under null (optional recompute)
                if o.RecomputeAcrossInPerm && Kvalid >= 2
                    pairs = nchoosek(validK(:)',2);
                    zA = nan(size(pairs,1),1);
                    wA = nan(size(pairs,1),1);
                    switch lower(o.AcrossMode)
                        case 'cross-halves'
                            for i = 1:size(pairs,1)
                                k1 = pairs(i,1); k2 = pairs(i,2);
                                k1o = map_to(validK==k1);
                                k2o = map_to(validK==k2);
                                r1 = safe_corr(PV_even(:,k1), PV_odd(:,k2o));
                                r2 = safe_corr(PV_odd(:,k1o),  PV_even(:,k2));
                                zA(i) = mean([atanh(max(min(r1,0.99999),-0.99999)), atanh(max(min(r2,0.99999),-0.99999))], 'omitnan');
                                wA(i) = sqrt(w(k1)*w(k2));
                            end
                        case 'mean-of-halves'
                            for i = 1:size(pairs,1)
                                k1 = pairs(i,1); k2 = pairs(i,2);
                                k1o = map_to(validK==k1);
                                k2o = map_to(validK==k2);
                                v1 = mean([PV_even(:,k1), PV_odd(:,k1o)],2,'omitnan');
                                v2 = mean([PV_even(:,k2), PV_odd(:,k2o)],2,'omitnan');
                                zA(i) = atanh(max(min(safe_corr(v1, v2),0.99999),-0.99999));
                                wA(i) = sqrt(w(k1)*w(k2));
                            end
                    end
                    wA(~isfinite(wA))=0; if ~any(wA>0), wA(:)=1; end
                    zA_perm_day(d) = nansum(wA .* zA) / max(1, nansum(wA));
                else
                    zA_perm_day(d) = z_across_day(d);
                end

                % ---- save per-day Δz for this permutation properly ----
                if ~isfield(B.perm,'byDay') || isempty(B.perm.byDay.delta)
                    B.perm.byDay.delta = cell(1,D);
                end
                if isfinite(zW_perm_day(d)) && isfinite(zA_perm_day(d))
                    B.perm.byDay.delta{d}(end+1,1) = zW_perm_day(d) - zA_perm_day(d); %#ok<AGROW>
                end
            end
            zW_perm_mean_all(pidx) = mean(zW_perm_day(isfinite(zW_perm_day)),'omitnan');
            zA_perm_mean_all(pidx) = mean(zA_perm_day(isfinite(zA_perm_day)),'omitnan');
            perm_delta(pidx)       = zW_perm_mean_all(pidx) - zA_perm_mean_all(pidx);
        end

    case 'frame-redistribute'
        for pidx = 1:Nperm
            zW_perm_day = nan(D,1);
            zA_perm_day = nan(D,1);

            for d = 1:D
                validK = validK_cache{d};
                if isempty(validK), continue; end
                w = w_cache{d};
                w_valid = w(validK); if ~any(w_valid>0), w_valid(:)=1; end

                even_idx_by_bin = even_idx_cache{d};
                odd_idx_by_bin  = odd_idx_cache{d};
                ne = n_even_cache{d};
                no = n_odd_cache{d};

                E_pool = vertcat(even_idx_by_bin{validK});
                O_pool = vertcat(odd_idx_by_bin{validK});
                E_pool = unique(clamp_idx(E_pool, numel(ts{d})));
                O_pool = unique(clamp_idx(O_pool, numel(ts{d})));
                if isempty(E_pool) || isempty(O_pool), continue; end

                % randomize/split to match original counts
                E_perm = E_pool(randperm(numel(E_pool)));
                O_perm = O_pool(randperm(numel(O_pool)));
                ce = max(0, round(ne(validK))); co = max(0, round(no(validK)));

                if sum(ce) > numel(E_perm)
                    deficit = sum(ce) - numel(E_perm);
                    for j = numel(ce):-1:1
                        take = min(ce(j), deficit);
                        ce(j) = ce(j) - take; deficit = deficit - take;
                        if deficit <= 0, break; end
                    end
                end
                if sum(co) > numel(O_perm)
                    deficit = sum(co) - numel(O_perm);
                    for j = numel(co):-1:1
                        take = min(co(j), deficit);
                        co(j) = co(j) - take; deficit = deficit - take;
                        if deficit <= 0, break; end
                    end
                end

                E_chunks = split_by_counts_ptr(E_perm, ce);
                O_chunks = split_by_counts_ptr(O_perm, co);

                % rebuild permuted PVs
                t_day = ts{d}(:);
                Sday  = spikes_to_matrix(spikes{d}, t_day);
                Sday  = normalize_cells(Sday, o.CellNorm);
                PV_even_p = nan(size(Sday,1), K);
                PV_odd_p  = nan(size(Sday,1), K);
                for j = 1:numel(validK)
                    k = validK(j);
                    segE = clamp_idx(E_chunks{j}, size(Sday,2));
                    segO = clamp_idx(O_chunks{j}, size(Sday,2));
                    if ~isempty(segE), PV_even_p(:,k) = mean(Sday(:, segE), 2, 'omitnan'); end
                    if ~isempty(segO), PV_odd_p(:,k)  = mean(Sday(:, segO),  2, 'omitnan'); end
                end

                % WITHIN (perm)
                zW_vec = nan(numel(validK),1);
                for j = 1:numel(validK)
                    k = validK(j);
                    zW_vec(j) = atanh(max(min(safe_corr(PV_even_p(:,k), PV_odd_p(:,k)),0.99999),-0.99999));
                end
                zW_perm_day(d) = nansum(w_valid .* zW_vec) / max(1,nansum(w_valid));

                % ACROSS (perm) per AcrossMode
                if numel(validK) >= 2
                    pairs = nchoosek(validK(:)', 2);
                    zA_vec = nan(size(pairs,1),1);
                    wA     = nan(size(pairs,1),1);
                    switch lower(o.AcrossMode)
                        case 'cross-halves'
                            for iPair = 1:size(pairs,1)
                                k1 = pairs(iPair,1); k2 = pairs(iPair,2);
                                r1 = safe_corr(PV_even_p(:,k1), PV_odd_p(:,k2));
                                r2 = safe_corr(PV_odd_p(:,k1),  PV_even_p(:,k2));
                                zA_vec(iPair) = mean([atanh(max(min(r1,0.99999),-0.99999)), atanh(max(min(r2,0.99999),-0.99999))], 'omitnan');
                                wA(iPair)     = sqrt(w(k1)*w(k2));
                            end
                        case 'mean-of-halves'
                            for iPair = 1:size(pairs,1)
                                k1 = pairs(iPair,1); k2 = pairs(iPair,2);
                                v1 = mean([PV_even_p(:,k1), PV_odd_p(:,k1)],2,'omitnan');
                                v2 = mean([PV_even_p(:,k2), PV_odd_p(:,k2)],2,'omitnan');
                                zA_vec(iPair) = atanh(max(min(safe_corr(v1, v2),0.99999),-0.99999));
                                wA(iPair)     = sqrt(w(k1)*w(k2));
                            end
                    end
                    wA(~isfinite(wA))=0; if ~any(wA>0), wA(:)=1; end
                    zA_perm_day(d) = nansum(wA .* zA_vec) / max(1, nansum(wA));
                end

                % save per-day Δz for this permutation
                if isfinite(zW_perm_day(d)) && isfinite(zA_perm_day(d))
                    B.perm.byDay.delta{d}(end+1,1) = zW_perm_day(d) - zA_perm_day(d); %#ok<AGROW>
                end
            end

            zW_perm_mean_all(pidx) = mean(zW_perm_day(isfinite(zW_perm_day)),'omitnan');
            zA_perm_mean_all(pidx) = mean(zA_perm_day(isfinite(zA_perm_day)),'omitnan');
            perm_delta(pidx)       = zW_perm_mean_all(pidx) - zA_perm_mean_all(pidx);
        end

    otherwise
        error('Unknown NullMode: %s', o.NullMode);
    end
end

fprintf('Obs: within_z=%.3f, across_z=%.3f, Δ_z=%.3f\n', z_within_mean, z_across_mean, delta_z);
fprintf('Null (perm): mean(within_z)=%.3f (sd=%.3f), mean(across_z)=%.3f (sd=%.3f), mean(Δ_z)=%.3f (sd=%.3f)\n', ...
        mean(zW_perm_mean_all,'omitnan'), std(zW_perm_mean_all,'omitnan'), ...
        mean(zA_perm_mean_all,'omitnan'), std(zA_perm_mean_all,'omitnan'), ...
        mean(perm_delta,'omitnan'),       std(perm_delta,'omitnan'));

% p-values
perm_delta = perm_delta(isfinite(perm_delta));
delta_obs  = delta_z;
if isempty(perm_delta)
    p_right = NaN; p_left = NaN; p_two = NaN;
else
    p_right = mean(perm_delta >= delta_obs);
    p_left  = mean(perm_delta <= delta_obs);
    p_two   = 2*min(p_right, p_left);
end

% per-day exposure for across-day testing
B.byDay.within.z = z_within_day;
B.byDay.across.z = z_across_day;
B.byDay.within.r = tanh(z_within_day);
B.byDay.across.r = tanh(z_across_day);
B.byDay.valid    = isfinite(z_within_day) & (isfinite(z_within_day)|isfinite(z_across_day));

% pack
B.within.r = tanh(z_within_mean);  B.within.z = z_within_mean;
B.across.r = tanh(z_across_mean);  B.across.z = z_across_mean;
B.delta.r  = tanh(delta_z);        B.delta.z  = delta_z;

B.p_perm_right = p_right;
B.p_perm_left  = p_left;
B.p_perm_two   = p_two;
B.perm_delta   = perm_delta;

B.null = struct('mode', o.NullMode, 'recomputeAcross', logical(o.RecomputeAcrossInPerm), ...
                'NPerm', o.NPerm, 'AcrossMode', o.AcrossMode);

B.notes = struct('MinTraceSamples',o.MinTraceSamples, ...
                 'MinTraceEachHalf',o.MinTraceEachHalf, ...
                 'BinWeight',o.BinWeight);
end






function B = compute_space_to_task_acrossTrials(spikes, ts, pos, cs, CTRL, varargin)
% compute_space_to_task_acrossTrials  (no nested helpers)
% Uses file-scope helpers (interp_pos, spikes_to_matrix, normalize_cells,
% pos2bin, choose_half_vec, choose_frames, pair_corr_with_flags_fileScope).

p = inputParser;
addParameter(p,'Halves','even');
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'CellNorm','demean');
addParameter(p,'NPerm',49);
addParameter(p,'MinTrials',2);
addParameter(p,'Debug',true);
addParameter(p,'Label','');

% soft-guard / audit options
addParameter(p,'MinNCorr',10);
addParameter(p,'MinStdBin',1e-3);
addParameter(p,'TopRInspect',10);
addParameter(p,'DoBounds',true);
addParameter(p,'DoSpearman',true);
addParameter(p,'DoLOO',false);
addParameter(p,'DupMode','ts-only');  % 'ts-only' or 'hybrid'

parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);
K = CTRL.grid.K;

zW_day_all   = nan(D,1);   zA_day_all   = nan(D,1);
zW_day_filt  = nan(D,1);   zA_day_filt  = nan(D,1);
has_any      = false(D,1);

keep_counts = struct('W_all',zeros(D,1),'A_all',zeros(D,1), ...
                     'W_keep',zeros(D,1),'A_keep',zeros(D,1));
flag_counts = struct('lowstd',zeros(D,1),'dup',zeros(D,1),'bound',zeros(D,1));

B.audit = cell(D,1);

for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);
    S = spikes_to_matrix(spikes{d}, t);
    S = normalize_cells(S, o.CellNorm);

    if o.UseSpeedMask
        v = speed_cm_per_s(pos{d});
        v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
        mask_speed = (v >= o.VelThresh);
    else
        mask_speed = true(size(t));
    end

    edges_d = CTRL.grid.edges{d};
    csd     = cs{d};
    nTrials = numel(csd);
    if nTrials < o.MinTrials, continue; end

    PV_even_trials = cell(nTrials,1);
    PV_odd_trials  = cell(nTrials,1);
    SrcFrames_even = cell(nTrials, K);
    SrcFrames_odd  = cell(nTrials, K);

    % Build per-trial PVs and frame sources
    for tr = 1:nTrials
        idx = (t >= csd(tr)+o.TraceWin(1) & t <= csd(tr)+o.TraceWin(2)) & mask_speed;
        frames = find(idx);
        if numel(frames) < 2, continue; end

        idx_even = frames(2:2:end);
        idx_odd  = frames(1:2:end);

        % even
        PV_e = nan(size(S,1), K);
        if ~isempty(idx_even)
            [~, bE] = pos2bin(x(idx_even), y(idx_even), edges_d);
            for k = 1:K
                sel = (bE == k);
                if any(sel)
                    f = idx_even(sel);
                    PV_e(:,k) = mean(S(:, f), 2, 'omitnan');
                    SrcFrames_even{tr,k} = f;
                else
                    SrcFrames_even{tr,k} = [];
                end
            end
        end
        PV_even_trials{tr} = PV_e;

        % odd
        PV_o = nan(size(S,1), K);
        if ~isempty(idx_odd)
            [~, bO] = pos2bin(x(idx_odd), y(idx_odd), edges_d);
            for k = 1:K
                sel = (bO == k);
                if any(sel)
                    f = idx_odd(sel);
                    PV_o(:,k) = mean(S(:, f), 2, 'omitnan');
                    SrcFrames_odd{tr,k} = f;
                else
                    SrcFrames_odd{tr,k} = [];
                end
            end
        end
        PV_odd_trials{tr} = PV_o;
    end

    % Which trials are valid for the chosen halves?
    switch lower(o.Halves)
        case 'even', PV_trials = PV_even_trials;
        case 'odd',  PV_trials = PV_odd_trials;
        case 'both', PV_trials = []; % handled later
        otherwise, error('Halves must be even|odd|both.');
    end

    if ~strcmpi(o.Halves,'both')
        valid_trials = find(~cellfun(@(M) isempty(M) || all(~isfinite(M(:))), PV_trials));
    else
        vt_even = find(~cellfun(@(M) isempty(M) || all(~isfinite(M(:))), PV_even_trials));
        vt_odd  = find(~cellfun(@(M) isempty(M) || all(~isfinite(M(:))), PV_odd_trials));
        valid_trials = union(vt_even, vt_odd);
    end
    if numel(valid_trials) < 2, continue; end
    has_any(d) = true;

    hasBin_even = false(numel(valid_trials), K);
    hasBin_odd  = false(numel(valid_trials), K);
    for ii = 1:numel(valid_trials)
        tIdx = valid_trials(ii);
        if ~isempty(PV_even_trials{tIdx}), hasBin_even(ii,:) = any(isfinite(PV_even_trials{tIdx}), 1); end
        if ~isempty(PV_odd_trials{tIdx}),  hasBin_odd(ii,:)  = any(isfinite(PV_odd_trials{tIdx}), 1);  end
    end
    if strcmpi(o.Halves,'even'), hasBin = hasBin_even;
    elseif strcmpi(o.Halves,'odd'), hasBin = hasBin_odd;
    else, hasBin = hasBin_even | hasBin_odd;
    end

    if o.Debug
        bins_per_trial = sum(hasBin,2);
        nCells = NaN;
        for tTry = valid_trials(:).'
            Mi = PV_even_trials{tTry}; if isempty(Mi), Mi = PV_odd_trials{tTry}; end
            if ~isempty(Mi), nCells = size(Mi,1); break; end
        end
        fprintf('   [AT DEBUG d=%d %s] validTrials=%d | bins/trial: mean=%.2f (min=%d, max=%d) | K=%d | nCells=%s\n', ...
            d, o.Halves, numel(valid_trials), mean(bins_per_trial), min(bins_per_trial), max(bins_per_trial), ...
            K, mat2str(nCells));

        n_within_est = 0;
        for k = 1:K, n_within_est = n_within_est + nchoosek(sum(hasBin(:,k)),2); end
        n_across_est = 0;
        for i1 = 1:numel(valid_trials)-1
            for i2 = i1+1:numel(valid_trials)
                n_across_est = n_across_est + ...
                    (sum(hasBin(i1,:)) * sum(hasBin(i2,:)) - sum(hasBin(i1,:) & hasBin(i2,:)));
            end
        end
        fprintf('   [AT DEBUG d=%d] est WITHIN pairs=%d | est ACROSS mismatched pairs=%d\n', d, n_within_est, n_across_est);
    end

    % Accumulators
    ZW=[]; ZA=[]; FW=[]; FA=[]; NW=[]; NA=[]; W_info=[]; A_info=[];
    across_checked = 0;

    % Tally + realized counters
    drop_counts = struct('kept',0,'lowstd',0,'tsdup',0,'algdup_unrobust',0,'affine',0,'bound',0);
    realized_counts = struct('within',0,'across',0);

    % ---------- WITHIN ----------
    for k = 1:K
        trials_k = find(hasBin(:,k));
        if numel(trials_k) < 2, continue; end
        for i1 = 1:numel(trials_k)-1
            t1 = valid_trials(trials_k(i1));
            v1 = choose_half_vec(PV_even_trials, PV_odd_trials, t1, k, o.Halves);
            for i2 = i1+1:numel(trials_k)
                t2 = valid_trials(trials_k(i2));
                v2 = choose_half_vec(PV_even_trials, PV_odd_trials, t2, k, o.Halves);

                f1 = choose_frames(SrcFrames_even, SrcFrames_odd, t1, k, o.Halves);
                f2 = choose_frames(SrcFrames_even, SrcFrames_odd, t2, k, o.Halves);

                [r,nC,flags,ex] = pair_corr_with_flags_fileScope( ...
                    v1, v2, o.MinNCorr, o.MinStdBin, f1, f2, S, o.DoSpearman, o.DoLOO, o.DoBounds);
                realized_counts.within = realized_counts.within + 1;

                % Timestamp duplicate?
                ts1 = unique(f1(:)); ts2 = unique(f2(:));
                tsdup = ~isempty(ts1) && ~isempty(ts2) && numel(ts1)==numel(ts2) && all(ts1==ts2);

                % Flags/context
                is_lowstd   = flags(1) ~= 0;
                is_dup_alg  = flags(2) ~= 0;
                over_bound  = flags(3) ~= 0;

                has_loo   = isfield(ex,'r_loo')   && isfinite(ex.r_loo);
                loo_close = has_loo && (abs(ex.r_loo - r) <= 0.02);

                has_spear = isfield(ex,'r_spear') && isfinite(ex.r_spear);
                spear_hi  = has_spear && (ex.r_spear >= 0.98);
                robust    = loo_close || spear_hi;

                has_bound   = isfield(ex,'r_bound') && isfinite(ex.r_bound);
                r_bound_val = NaN; if has_bound, r_bound_val = ex.r_bound; end

                has_R2         = isfield(ex,'R2') && isfinite(ex.R2);
                affine_perfect = has_R2 && (ex.R2 > 0.9999);

                diffnorm = NaN;
                if isfield(ex,'diffnorm') && isfinite(ex.diffnorm), diffnorm = ex.diffnorm; end

                % --------- Decide keep (single policy) ---------
                keep = true;
                if is_lowstd
                    keep = false;
                elseif tsdup
                    keep = false;
                elseif affine_perfect && isfinite(diffnorm) && (diffnorm <= 1e-12)
                    keep = false;
                elseif over_bound && has_bound && (r > r_bound_val + 0.25) && ~robust
                    keep = false;
                elseif is_dup_alg && ~robust
                    keep = false;
                end

                % ---- flag tallies (independent of keep) ----
                flags_seen.lowstd = flags(1) ~= 0;
                flags_seen.dup    = flags(2) ~= 0;   % algorithmic/affine "dup" flag
                flags_seen.bound  = flags(3) ~= 0;

                % accumulate “was flagged” counts (for audit, not drop tallies)
                flag_counts.lowstd(d) = flag_counts.lowstd(d) + flags_seen.lowstd;
                flag_counts.dup(d)    = flag_counts.dup(d)    + flags_seen.dup;
                flag_counts.bound(d)  = flag_counts.bound(d)  + flags_seen.bound;

                % ---- final keep decision already computed above -> variable `keep` ----
                if keep
                    drop_counts.kept = drop_counts.kept + 1;
                else
                    % increment exactly one drop bucket based on the reason that made keep=false
                    if flags_seen.lowstd
                        drop_counts.lowstd = drop_counts.lowstd + 1;
                    elseif tsdup
                        drop_counts.tsdup = drop_counts.tsdup + 1;
                    elseif affine_perfect && isfinite(diffnorm) && (diffnorm <= 1e-12)
                        drop_counts.affine = drop_counts.affine + 1;
                    elseif flags_seen.bound && has_bound && (r > r_bound_val + 0.25) && ~robust
                        drop_counts.bound = drop_counts.bound + 1;
                    elseif (~robust) && is_dup_alg
                        drop_counts.algdup_unrobust = drop_counts.algdup_unrobust + 1;
                    else
                        % optional: a catch-all if you want to know when none matched
                        % drop_counts.other = getfieldwithdefault(drop_counts,'other',0) + 1;
                    end
                end

                if ~isfinite(r) || ~keep, continue; end
                ZW(end+1,1) = atanh(max(min(r,0.999999),-0.999999));
                NW(end+1,1) = nC;
                FW(end+1,:) = flags;
                W_info(end+1,:) = [t1 t2 k k r ex.r_spear ex.r_loo ex.r_bound];
            end
        end
    end

    % ---------- ACROSS ----------
    for i1 = 1:numel(valid_trials)-1
        t1 = valid_trials(i1);
        bins1 = find(hasBin(i1,:)); if isempty(bins1), continue; end
        for i2 = i1+1:numel(valid_trials)
            t2 = valid_trials(i2);
            bins2 = find(hasBin(i2,:)); if isempty(bins2), continue; end
            for k1 = bins1
                for k2 = bins2
                    if k2 == k1, continue; end
                    across_checked = across_checked + 1;

                    v1 = choose_half_vec(PV_even_trials, PV_odd_trials, t1, k1, o.Halves);
                    v2 = choose_half_vec(PV_even_trials, PV_odd_trials, t2, k2, o.Halves);
                    f1 = choose_frames(SrcFrames_even, SrcFrames_odd, t1, k1, o.Halves);
                    f2 = choose_frames(SrcFrames_even, SrcFrames_odd, t2, k2, o.Halves);

                    [r,nC,flags,ex] = pair_corr_with_flags_fileScope( ...
                        v1, v2, o.MinNCorr, o.MinStdBin, f1, f2, S, o.DoSpearman, o.DoLOO, o.DoBounds);
                    realized_counts.across = realized_counts.across + 1;

                    % Timestamp duplicate?
                    ts1 = unique(f1(:)); ts2 = unique(f2(:));
                    tsdup = ~isempty(ts1) && ~isempty(ts2) && numel(ts1)==numel(ts2) && all(ts1==ts2);

                    % Flags/context
                    is_lowstd   = flags(1) ~= 0;
                    is_dup_alg  = flags(2) ~= 0;
                    over_bound  = flags(3) ~= 0;

                    has_loo   = isfield(ex,'r_loo')   && isfinite(ex.r_loo);
                    loo_close = has_loo && (abs(ex.r_loo - r) <= 0.02);

                    has_spear = isfield(ex,'r_spear') && isfinite(ex.r_spear);
                    spear_hi  = has_spear && (ex.r_spear >= 0.98);
                    robust    = loo_close || spear_hi;

                    has_bound   = isfield(ex,'r_bound') && isfinite(ex.r_bound);
                    r_bound_val = NaN; if has_bound, r_bound_val = ex.r_bound; end

                    has_R2         = isfield(ex,'R2') && isfinite(ex.R2);
                    affine_perfect = has_R2 && (ex.R2 > 0.9999);

                    diffnorm = NaN;
                    if isfield(ex,'diffnorm') && isfinite(ex.diffnorm), diffnorm = ex.diffnorm; end

                    % --------- Decide keep (single policy) ---------
                    keep = true;
                    if is_lowstd
                        keep = false;
                    elseif tsdup
                        keep = false;
                    elseif affine_perfect && isfinite(diffnorm) && (diffnorm <= 1e-12)
                        keep = false;
                    elseif over_bound && has_bound && (r > r_bound_val + 0.25) && ~robust
                        keep = false;
                    elseif is_dup_alg && ~robust
                        keep = false;
                    end

                    % Tally drop reasons / kept
                    % ---- flag tallies (independent of keep) ----
                    flags_seen.lowstd = flags(1) ~= 0;
                    flags_seen.dup    = flags(2) ~= 0;   % algorithmic/affine "dup" flag
                    flags_seen.bound  = flags(3) ~= 0;

                    % accumulate “was flagged” counts (for audit, not drop tallies)
                    flag_counts.lowstd(d) = flag_counts.lowstd(d) + flags_seen.lowstd;
                    flag_counts.dup(d)    = flag_counts.dup(d)    + flags_seen.dup;
                    flag_counts.bound(d)  = flag_counts.bound(d)  + flags_seen.bound;

                    % ---- final keep decision already computed above -> variable `keep` ----
                    if keep
                        drop_counts.kept = drop_counts.kept + 1;
                    else
                        % increment exactly one drop bucket based on the reason that made keep=false
                        if flags_seen.lowstd
                            drop_counts.lowstd = drop_counts.lowstd + 1;
                        elseif tsdup
                            drop_counts.tsdup = drop_counts.tsdup + 1;
                        elseif affine_perfect && isfinite(diffnorm) && (diffnorm <= 1e-12)
                            drop_counts.affine = drop_counts.affine + 1;
                        elseif flags_seen.bound && has_bound && (r > r_bound_val + 0.25) && ~robust
                            drop_counts.bound = drop_counts.bound + 1;
                        elseif (~robust) && is_dup_alg
                            drop_counts.algdup_unrobust = drop_counts.algdup_unrobust + 1;
                        else
                            % optional: a catch-all if you want to know when none matched
                            % drop_counts.other = getfieldwithdefault(drop_counts,'other',0) + 1;
                        end
                    end

                    if ~isfinite(r) || ~keep, continue; end
                    ZA(end+1,1) = atanh(max(min(r,0.999999),-0.999999));
                    NA(end+1,1) = nC;
                    FA(end+1,:) = flags;
                    A_info(end+1,:) = [t1 t2 k1 k2 r ex.r_spear ex.r_loo ex.r_bound];
                end
            end
        end
    end

    % ==== AFTER finishing both WITHIN and ACROSS loops (per day d) ====
    % Tally + assert
    n_realized = realized_counts.within + realized_counts.across;
    n_dropped  = drop_counts.lowstd + drop_counts.tsdup + drop_counts.algdup_unrobust ...
                 + drop_counts.affine + drop_counts.bound;

    assert(drop_counts.kept + n_dropped == n_realized, ...
           'Tally mismatch: kept(%d)+dropped(%d) ~= realized(%d)', ...
           drop_counts.kept, n_dropped, n_realized);

    if o.Debug
        % ACROSS-pairs audit you already compute elsewhere
        fprintf('   [AT DEBUG d=%d] ACROSS pairs k1~=k2 verified (n=%d)\n', d, across_checked);

        % Show pre-keep (realized) vs post-keep (ZW/ZA) counts
        fprintf('   [AT DEBUG d=%d] realized (pre-keep) WITHIN n=%d | ACROSS n=%d\n', ...
            d, realized_counts.within, realized_counts.across);
        fprintf('   [AT DEBUG d=%d] kept (post-keep)     WITHIN n=%d | ACROSS n=%d\n', ...
            d, numel(ZW), numel(ZA));

        % Single unified kept/dropped line
        fprintf(['   [AT DEBUG d=%d] kept=%d | dropped=%d (lowstd=%d tsdup=%d ' ...
                 'algdup_unrobust=%d affine=%d bound=%d)\n'], ...
            d, drop_counts.kept, n_dropped, drop_counts.lowstd, drop_counts.tsdup, ...
            drop_counts.algdup_unrobust, drop_counts.affine, drop_counts.bound);
    end

    % ==== (then keep your existing overlap and r-summary prints that follow) ====

    % ---------- Day-level z, flags, and audit ----------
    % NEW: save the kept per-pair distributions for CDFs
    AT_byDay_within_list_z{d} = ZW(:);
    AT_byDay_across_list_z{d} = ZA(:);


    zW_day_all(d) = mean(ZW,'omitnan');
    zA_day_all(d) = mean(ZA,'omitnan');

    lowstd_W = (size(FW,2)>=1) & FW(:,1);  lowstd_A = (size(FA,2)>=1) & FA(:,1);
    dup_W    = (size(FW,2)>=2) & FW(:,2);  dup_A    = (size(FA,2)>=2) & FA(:,2);
    bnd_W    = (size(FW,2)>=3) & FW(:,3);  bnd_A    = (size(FA,2)>=3) & FA(:,3);

    flag_counts.lowstd(d) = nnz(lowstd_W) + nnz(lowstd_A);
    flag_counts.dup(d)    = nnz(dup_W)    + nnz(dup_A);
    flag_counts.bound(d)  = nnz(bnd_W)    + nnz(bnd_A);

    keepW = ~(lowstd_W | dup_W | bnd_W);
    keepA = ~(lowstd_A | dup_A | bnd_A);
    keep_counts.W_all(d)  = numel(ZW);
    keep_counts.A_all(d)  = numel(ZA);
    keep_counts.W_keep(d) = nnz(keepW);
    keep_counts.A_keep(d) = nnz(keepA);

    zW_day_filt(d) = mean(ZW(keepW),'omitnan');
    zA_day_filt(d) = mean(ZA(keepA),'omitnan');

    % High-r audit
    [~, idxW] = sort(tanh(ZW), 'descend'); idxW = idxW(1:min(o.TopRInspect, numel(idxW)));
    [~, idxA] = sort(tanh(ZA), 'descend'); idxA = idxA(1:min(o.TopRInspect, numel(idxA)));
    audit_lines = strings(0,1);
    fmt = '      [TOP %s] r=%.4f | keep=%d | flags=[lowstd:%d dup:%d bound:%d] | t1=%d t2=%d k1=%d k2=%d | Spearman=%.3f | LOO=%.3f | r_bound=%.3f';
    for q = 1:numel(idxW)
        i = idxW(q);
        audit_lines(end+1,1) = sprintf(fmt,'WITHIN',tanh(ZW(i)), keepW(i), FW(i,1),FW(i,2),FW(i,3), ...
            W_info(i,1),W_info(i,2),W_info(i,3),W_info(i,4), W_info(i,6), W_info(i,7), W_info(i,8));
    end
    for q = 1:numel(idxA)
        i = idxA(q);
        audit_lines(end+1,1) = sprintf(fmt,'ACROSS',tanh(ZA(i)), keepA(i), FA(i,1),FA(i,2),FA(i,3), ...
            A_info(i,1),A_info(i,2),A_info(i,3),A_info(i,4), A_info(i,6), A_info(i,7), A_info(i,8));
    end
    B.audit{d} = audit_lines;
    if o.Debug && ~isempty(audit_lines)
        fprintf('   [AT DEBUG d=%d] High-r audit (top %d W, top %d A):\n', d, numel(idxW), numel(idxA));
        for s = 1:numel(audit_lines), fprintf('%s\n', audit_lines(s)); end
    end
end

% --------- Per-rat summaries (ALL) ----------
valid_all  = has_any & isfinite(zW_day_all)  & isfinite(zA_day_all);
zW_mu_all  = mean(zW_day_all(valid_all),'omitnan');
zA_mu_all  = mean(zA_day_all(valid_all),'omitnan');
delta_all  = (zW_day_all(valid_all) - zA_day_all(valid_all));
delta_obs_all = mean(delta_all,'omitnan');

Np = max(49, 2*o.NPerm);
perm_all = nan(Np,1);
for b = 1:Np
    s = (rand(size(delta_all))>0.5)*2-1; perm_all(b) = mean(s .* delta_all);
end
p_two_all = (sum(abs(perm_all) >= abs(delta_obs_all)) + 1) / (Np + 1);

% --------- Per-rat summaries (FILTERED) ----------
valid_f = has_any & isfinite(zW_day_filt) & isfinite(zA_day_filt);
zW_mu_f = mean(zW_day_filt(valid_f),'omitnan');
zA_mu_f = mean(zA_day_filt(valid_f),'omitnan');
delta_f = (zW_day_filt(valid_f) - zA_day_filt(valid_f));
delta_obs_f = mean(delta_f,'omitnan');

perm_f = nan(Np,1);
for b = 1:Np
    s = (rand(size(delta_f))>0.5)*2-1; perm_f(b) = mean(s .* delta_f);
end
p_two_f = (sum(abs(perm_f) >= abs(delta_obs_f)) + 1) / (Np + 1);

% --------- Pack outputs ----------
B.byDay.valid       = has_any;
B.byDay.within.z    = zW_day_all;         B.byDay.across.z    = zA_day_all;
B.byDay.within.r    = tanh(zW_day_all);   B.byDay.across.r    = tanh(zA_day_all);
B.byDay.within_f.z  = zW_day_filt;        B.byDay.across_f.z  = zA_day_filt;
B.byDay.within_f.r  = tanh(zW_day_filt);  B.byDay.across_f.r  = tanh(zA_day_filt);

B.counts = struct('within_all', keep_counts.W_all, 'across_all', keep_counts.A_all, ...
                  'within_kept', keep_counts.W_keep, 'across_kept', keep_counts.A_keep);
B.flags  = struct('lowstd', flag_counts.lowstd, 'dup', flag_counts.dup, 'bound', flag_counts.bound);

B.within.z = zW_mu_all;  B.within.r = tanh(zW_mu_all);
B.across.z = zA_mu_all;  B.across.r = tanh(zA_mu_all);
B.delta.z  = delta_obs_all; B.delta.r = tanh(delta_obs_all);
B.perm_all = struct('delta_null', perm_all, 'delta_obs', delta_obs_all, 'p_two', p_two_all);

B.filtered.within.z = zW_mu_f;  B.filtered.within.r = tanh(zW_mu_f);
B.filtered.across.z = zA_mu_f;  B.filtered.across.r = tanh(zA_mu_f);
B.filtered.delta.z  = delta_obs_f; B.filtered.delta.r = tanh(delta_obs_f);
B.filtered.perm = struct('delta_null', perm_f, 'delta_obs', delta_obs_f, 'p_two', p_two_f);

% NEW: expose per-day lists (both z and r) for CDF plots
B.byDay.within.list_z = AT_byDay_within_list_z;
B.byDay.across.list_z = AT_byDay_across_list_z;
B.byDay.within.list_r = cellfun(@tanh, AT_byDay_within_list_z, 'uni', 0);
B.byDay.across.list_r = cellfun(@tanh, AT_byDay_across_list_z, 'uni', 0);

if any(valid_all)
    fprintf('[AT SUMMARY %s] days=%d | ALL Δz=%.3f (r≈%.3f), p_two=%.3g | FILTERED Δz=%.3f (r≈%.3f), p_two=%.3g\n', ...
        L, nnz(valid_all), delta_obs_all, tanh(delta_obs_all), p_two_all, ...
        delta_obs_f, tanh(delta_obs_f), p_two_f);
end
end

% ---------- tiny selectors (no state) ----------
function v = choose_half_vec(PVe, PVo, tIdx, k, halves)
    if strcmpi(halves,'even') || (strcmpi(halves,'both') && ~isempty(PVe{tIdx}) && any(isfinite(PVe{tIdx}(:))))
        v = PVe{tIdx}(:,k);
    else
        v = PVo{tIdx}(:,k);
    end
end
function f = choose_frames(Fe, Fo, tIdx, k, halves)
    if strcmpi(halves,'even')
        f = Fe{tIdx,k};
    else
        f = Fo{tIdx,k};
    end
end

function [r, nC, flags, extras] = pair_corr_with_flags_fileScope(v1, v2, minN, minStd, ...
                                            frames1, frames2, S_full, doSpearman, doLOO, doBounds)
% flags: [lowstd dup bound]
flags = [false false false];
extras = struct('r_spear',NaN,'r_loo',NaN,'r_bound',NaN,'a',NaN,'b',NaN,'R2',NaN,...
                'diffnorm',NaN,'maxabs',NaN);

% basic guards
s1 = std(v1,0,'omitnan');
s2 = std(v2,0,'omitnan');
if s1 < minStd || s2 < minStd
    flags(1) = true;
end

[r, nC] = safe_corr(v1, v2, minN, minStd);
if ~isfinite(r)
    return
end

% duplicate / affine check for very high r
if r > 0.98
    diffv = v1 - v2;
    % NaN-robust L2 norm and max|.| without 'omitnan'
    dm = diffv(isfinite(diffv));
    extras.diffnorm = sqrt(sum(dm.^2));                  % ||Δ||_2
    if isempty(dm), extras.maxabs = NaN; else, extras.maxabs = max(abs(dm)); end

    m = isfinite(v1) & isfinite(v2);
    if nnz(m) >= 3
        X = [v1(m) ones(nnz(m),1)];
        y =  v2(m);
        ab = X\y;
        yhat = X*ab;
        res = y - yhat;
        extras.a = ab(1);
        extras.b = ab(2);
        vy = var(y,0,'omitnan');
        if vy > 0
            extras.R2 = 1 - var(res,0,'omitnan')/vy;
        end
        % mark as duplicate/degenerate if virtually identical or perfectly affine
        if (~isnan(extras.maxabs) && extras.maxabs < 1e-9) || (isfinite(extras.R2) && extras.R2 > 0.9999)
            flags(2) = true;
        end
    end
end

% Spearman (optional)
if doSpearman
    m = isfinite(v1) & isfinite(v2);
    if nnz(m) >= minN
        try
            extras.r_spear = corr(v1(m), v2(m), 'Type','Spearman');
        catch
            % ignore
        end
    end
end

% Leave-one-out (optional)
if doLOO
    m = isfinite(v1) & isfinite(v2);
    idx = find(m);
    nn = numel(idx);
    if nn >= max(minN,12)
        rs = nan(nn,1);
        for jj = 1:nn
            mm = m;
            mm(idx(jj)) = false;
            try
                rs(jj) = corr(v1(mm), v2(mm), 'Type','Pearson','Rows','pairwise');
            catch
                rs(jj) = NaN;
            end
        end
        extras.r_loo = mean(rs,'omitnan');
    end
end

% reliability bound (optional)
if doBounds
    rel1 = split_half_reliability_fileScope(frames1, S_full);
    rel2 = split_half_reliability_fileScope(frames2, S_full);
    if ~isfinite(rel1), rel1 = NaN; end
    if ~isfinite(rel2), rel2 = NaN; end
    rel1 = max(rel1,0);
    rel2 = max(rel2,0);
    extras.r_bound = sqrt(rel1 * rel2);
    if isfinite(extras.r_bound) && r > extras.r_bound + 0.05
        flags(3) = true;
    end
end
end


function rb = split_half_reliability_fileScope(frames, S_full)
rb = NaN;
if isempty(frames) || numel(frames) < 4, return; end
f = unique(frames(:));
fe = f(2:2:end); fo = f(1:2:end);
if isempty(fe) || isempty(fo), return; end
v1 = mean(S_full(:, fe),2,'omitnan');
v2 = mean(S_full(:, fo),2,'omitnan');
rb = safe_corr(v1, v2, 5, 0);  % min guards: low cost
end





% ===================== local helpers (debug-safe) ======================
function [z, nOverlaps] = local_within_same_bin(PV_trials, valid_trials, hasBin, K, minN, minStd)
    if nargin < 5, minN = 10; end
    if nargin < 6, minStd = 1e-3; end
    z = []; nOverlaps = [];
    for k = 1:K
        trials_k = find(hasBin(:,k));   % trials that have data for bin k
        if numel(trials_k) < 2, continue; end
        for i1 = 1:numel(trials_k)-1
            t1 = valid_trials(trials_k(i1));
            v1 = PV_trials{t1}(:,k);
            for i2 = i1+1:numel(trials_k)
                t2 = valid_trials(trials_k(i2));
                v2 = PV_trials{t2}(:,k);
                [r, nC] = safe_corr(v1, v2, minN, minStd);
                if isfinite(r)
                    z(end+1,1) = atanh(max(min(r,0.99999),-0.99999)); %#ok<AGROW>
                    nOverlaps(end+1,1) = nC; %#ok<AGROW>
                end
            end
        end
    end
end


% signature change
function [z, nOverlaps, pairs_used] = local_across_diff_bins(PV_trials, valid_trials, hasBin, K, minN, minStd)
    if nargin < 5, minN = 10; end
    if nargin < 6, minStd = 1e-3; end
    z = []; nOverlaps = []; pairs_used = [];  % <-- new
    for i1 = 1:numel(valid_trials)-1
        t1 = valid_trials(i1);
        bins1 = find(hasBin(i1,:));
        if isempty(bins1), continue; end
        for i2 = i1+1:numel(valid_trials)
            t2 = valid_trials(i2);
            bins2 = find(hasBin(i2,:));
            if isempty(bins2), continue; end
            for k1 = bins1
                v1 = PV_trials{t1}(:,k1);
                for k2 = bins2
                    if k2 == k1, continue; end    % enforce different bins
                    v2 = PV_trials{t2}(:,k2);
                    [r, nC] = safe_corr(v1, v2, minN, minStd);
                    if isfinite(r)
                        z(end+1,1) = atanh(max(min(r,0.99999),-0.99999));
                        nOverlaps(end+1,1) = nC;
                        pairs_used(end+1,:) = [k1, k2];     %#ok<AGROW>
                    end
                end
            end
        end
    end
end


function PV_out = shuffle_bins(PV_in)
    % Shuffle columns (bins) independently within each trial matrix
    PV_out = PV_in;
    for i = 1:numel(PV_in)
        M = PV_in{i};
        if isempty(M), continue; end
        K = size(M,2);
        PV_out{i} = M(:, randperm(K));
    end
end


function plot_space_to_task_acrossTrials_perAnimal(R)
% Show within vs across (r) for AT, one figure per animal, plus day dots.

for i = 1:numel(R)
    if ~isfield(R(i),'AT') || isempty(R(i).AT), continue; end
    AT = R(i).AT;
    if ~isfield(AT,'within') || ~isfield(AT,'across'), continue; end
    figure('Color','w','Position',[180 180 640 420]); hold on
    bar(1, AT.within.r, 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k');
    bar(2, AT.across.r, 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k');
    xticks([1 2]); xticklabels({'within (same bin, across trials)','across (diff bins, across trials)'});
    ylabel('Mean PV correlation (r)');
    title(sprintf('[%s] Across-trials even→even (\\Delta=%.3f)', R(i).animal, AT.delta.r), 'Interpreter','none');
    yline(0,'k:'); box on; grid on

    % show per-day dots if present
    if isfield(AT,'byDay') && isfield(AT.byDay,'within') && isfield(AT.byDay,'across')
        zw = AT.byDay.within.z; za = AT.byDay.across.z;
        w = tanh(zw); a = tanh(za);
        m = isfinite(w) & isfinite(a);
        plot(ones(nnz(m),1), w(m), 'ko', 'MarkerFaceColor',[0.1 0.4 0.9], 'MarkerSize',5);
        plot(2*ones(nnz(m),1), a(m), 'ko', 'MarkerFaceColor',[0.9 0.4 0.1], 'MarkerSize',5);
        for j = find(m).'
            plot([1 2],[w(j) a(j)], '-', 'Color',[0 0 0 0.25]);
        end
    end
end
end

function plot_space_to_task_acrossTrials_group(R)
% Group bar of AT within vs across using per-animal means; show a paired t on Fisher-z.

W = []; A = []; % per-animal means (z)
for i = 1:numel(R)
    if ~isfield(R(i),'AT') || isempty(R(i).AT), continue; end
    AT = R(i).AT;

    if isfield(AT,'byDay') && isfield(AT.byDay,'within') && isfield(AT.byDay,'across') ...
            && isfield(AT.byDay.within,'z') && isfield(AT.byDay.across,'z')
        zw = AT.byDay.within.z;
        za = AT.byDay.across.z;
        m = isfinite(zw) & isfinite(za);
        if any(m)
            W(end+1) = mean(zw(m),'omitnan'); %#ok<AGROW>
            A(end+1) = mean(za(m),'omitnan'); %#ok<AGROW>
        end
    elseif isfield(AT,'within') && isfield(AT,'across') ...
            && isfield(AT.within,'z') && isfield(AT.across,'z')
        % fallback to animal-level scalars if per-day missing
        W(end+1) = AT.within.z;  %#ok<AGROW>
        A(end+1) = AT.across.z;  %#ok<AGROW>
    end
end

if isempty(W) || isempty(A)
    warning('Across-trials group: no data to plot.');
    return
end

[~, p, ~, st] = ttest(W, A);  % paired on Fisher-z
dz = mean(W-A)/std(W-A);      % Cohen’s dz
Wr = tanh(mean(W));
Ar = tanh(mean(A));

figure('Color','w','Position',[200 200 700 420]); hold on
bar(1, Wr, 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k');
bar(2, Ar, 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k');
xticks([1 2]); xticklabels({'within (across trials)','across (across trials)'});
ylabel('Mean PV correlation (r)');

% jittered paired lines (in r, for readability)
wr = tanh(W(:));
ar = tanh(A(:));
for i = 1:numel(wr)
    plot([1 2],[wr(i) ar(i)], '-', 'Color',[0 0 0 0.25]);
    plot(1, wr(i), 'ko', 'MarkerFaceColor',[0.1 0.4 0.9], 'MarkerSize',5);
    plot(2, ar(i), 'ko', 'MarkerFaceColor',[0.9 0.4 0.1], 'MarkerSize',5);
end

title(sprintf('Across-trials even→even | paired t (z): t(%d)=%.2f, p=%.3g, dz=%.2f, n=%d', ...
    st.df, st.tstat, p, dz, numel(W)));
yline(0,'k:'); box on; grid on
end

function summarize_taskSpacePV_clean(R)
% Prints standardized stats + makes bar plots for:
% (B) Task within vs across (odd↔even)
% (AT) Task within vs across (even↔even across trials)
% (C) WITH vs WITHOUT task

do_one('B: Task within vs across (odd↔even)', ...
    @get_B_dayZ, @get_B_dayPermNull, R);

do_one('AT: Task within vs across (even↔even across trials)', ...
    @get_AT_dayZ, @get_AT_dayPermNull, R);

do_one('C: WITH vs WITHOUT task', ...
    @get_C_dayZ, @get_C_dayPermNull, R);
end

% --------- per-condition extractors ---------
% --------- per-condition extractors (now return did) ---------
function [Zw, Za, aid, did] = get_B_dayZ(R)
Zw=[]; Za=[]; aid=[]; did=[];
for i=1:numel(R)
    if ~isfield(R(i),'B')||~isfield(R(i).B,'byDay'), continue; end
    zW = R(i).B.byDay.within.z(:);
    zA = R(i).B.byDay.across.z(:);
    v  = R(i).B.byDay.valid(:);
    m  = v & isfinite(zW) & isfinite(zA);
    Zw = [Zw; zW(m)]; Za=[Za; zA(m)];
    aid = [aid; i*ones(nnz(m),1)];
    did = [did; find(m)];  % day index within that animal
end
end

function dayNull = get_B_dayPermNull(R)
dayNull = cell(numel(R),1);
for i=1:numel(R)
    if isfield(R(i),'B') && isfield(R(i).B,'perm') && isfield(R(i).B.perm,'byDay')
        dayNull{i} = R(i).B.perm.byDay.delta;  % cell{day} -> [P x 1] Δz
    end
end
end

function [Zw, Za, aid, did] = get_AT_dayZ(R)
Zw=[]; Za=[]; aid=[]; did=[];
for i=1:numel(R)
    if ~isfield(R(i),'AT')||~isfield(R(i).AT,'byDay'), continue; end
    zW = R(i).AT.byDay.within.z(:);
    zA = R(i).AT.byDay.across.z(:);
    m  = isfinite(zW) & isfinite(zA);
    Zw = [Zw; zW(m)]; Za=[Za; zA(m)];
    aid = [aid; i*ones(nnz(m),1)];
    did = [did; find(m)];
end
end

function dayNull = get_AT_dayPermNull(R)
dayNull = cell(numel(R),1);
for i=1:numel(R)
    if isfield(R(i),'AT') && isfield(R(i).AT,'day')
        D = numel(R(i).AT.day);
        tmp = cell(D,1);
        for d=1:D
            if isfield(R(i).AT.day(d),'perm') && isfield(R(i).AT.day(d).perm,'delta_null')
                tmp{d} = R(i).AT.day(d).perm.delta_null(:);
            end
        end
        dayNull{i} = tmp;
    end
end
end

function [Zw, Za, aid, did] = get_C_dayZ(R)
Zw=[]; Za=[]; aid=[]; did=[];
for i=1:numel(R)
    if ~isfield(R(i),'C')||~isfield(R(i).C,'bin')||~isfield(R(i).C.bin,'byDay'), continue; end
    z_with = R(i).C.bin.byDay.z_mean(:);   % WITH
    if ~isfield(R(i).C,'byDay')||~isfield(R(i).C.byDay,'withoutTask'), continue; end
    ur = R(i).C.byDay.withoutTask.r;
    z_without = nan(numel(ur),1);
    for d=1:numel(ur)
        if ~isempty(ur{d})
            z_without(d) = mean(atanh(max(min(ur{d}(:),0.99999),-0.99999)),'omitnan');
        end
    end
    m = isfinite(z_with) & isfinite(z_without);
    Zw  = [Zw;  z_with(m)];
    Za  = [Za;  z_without(m)];
    aid = [aid; i*ones(nnz(m),1)];
    did = [did; find(m)];
end
end

function dayNull = get_C_dayPermNull(R)
dayNull = cell(numel(R),1);
for i=1:numel(R)
    if isfield(R(i),'C') && isfield(R(i).C,'distStats') && isfield(R(i).C.distStats,'perm')
        dayNull{i} = R(i).C.distStats.perm.byDay; % cell{day} -> [P x 1]
    end
end
end

% --------- core printer + plotter ----------
function do_one(titleStr, pullZ, pullNull, R)
[Zw_all, Za_all, aid] = pullZ(R);
animals = {R.animal};

% 1) Paired t across ALL DAYS (two-sided)
[~, p_all, ~, st] = ttest(Zw_all, Za_all);
dz = mean(Zw_all-Za_all)/std(Zw_all-Za_all);
fprintf('\n== %s ==\n', titleStr);
fprintf('All days paired t: nDays=%d, t(%d)=%.2f, p=%.3g, dz=%.2f\n', numel(Zw_all), st.df, st.tstat, p_all, dz);

% 2) Per-animal paired t (pooling that animal’s days)
for i=1:numel(animals)
    m = (aid==i);
    if nnz(m) >= 2
        [~, p, ~, st_i] = ttest(Zw_all(m), Za_all(m));
        dz_i = mean(Zw_all(m)-Za_all(m))/std(Zw_all(m)-Za_all(m));
        fprintf('  [%s] paired t across its days: t(%d)=%.2f, p=%.3g, dz=%.2f, nDays=%d\n', ...
            animals{i}, st_i.df, st_i.tstat, p, dz_i, nnz(m));
    end
end

% 3) Per-animal two-sided sign-flip (across days) on Δz
fprintf('  Sign-flip across days (two-sided):\n');
for i=1:numel(animals)
    m = (aid==i);
    if nnz(m) >= 2
        d = Zw_all(m) - Za_all(m);
        Np = 49;
        s = (rand(Np, numel(d))>0.5)*2-1;
        perm = mean(s.*d',2);
        obs  = mean(d);
        p_two = (sum(abs(perm) >= abs(obs)) + 1) / (Np + 1);
        fprintf('    [%s] p_two=%.4g (obs Δz=%.3f)\n', animals{i}, p_two, obs);
    end
end

% 4) Bin-shuffle per day → aggregate (two-sided):
fprintf('  Bin-shuffle per day → aggregate (two-sided):\n');

% Pull the z-values (Zw, Za) *and* the mapping (aid, did)
[Zw, Za, aid, did] = pullZ(R);                 % <-- keep Zw & Za
dayNull = pullNull(R);                          % cell{animal}{day} = [P x 1] Δz samples

for i = 1:numel({R.animal})
    keep = (aid == i);
    if nnz(keep) < 2 || isempty(dayNull) || numel(dayNull) < i || isempty(dayNull{i})
        continue;
    end

    % observed per-day Δz for this animal
    obs_day = Zw(keep) - Za(keep);             % Δz per day
    did_i   = did(keep);                       % day indices within this animal

    % retain only days that have a stored null
    good = false(size(did_i));
    for k = 1:numel(did_i)
        d0 = did_i(k);
        good(k) = d0 <= numel(dayNull{i}) && ~isempty(dayNull{i}{d0});
    end
    obs_day = obs_day(good);
    did_i   = did_i(good);
    if numel(did_i) < 2, continue; end

    % align across days by taking the first P samples common to all days
    P = min( cellfun(@numel, dayNull{i}(did_i)) );
    P = min(P, 490);
    if P < 5, continue; end

    permd = nan(P,1);
    for p = 1:P
        mu = nan(numel(did_i),1);
        for k = 1:numel(did_i)
            mu(k) = dayNull{i}{did_i(k)}(p);   % this is a Δz sample for that day
        end
        permd(p) = mean(mu, 'omitnan');        % group (per-animal) Δz for perm p
    end

    % observed group Δz (mean across that animal’s days)
    obs  = mean(obs_day, 'omitnan');

    % two-sided p with add-one correction
    p_two = (sum(abs(permd) >= abs(obs)) + 1) / (numel(permd) + 1);
    fprintf('    [%s] p_two=%.4g (obs Δz=%.3f)\n', R(i).animal, p_two, obs);
end



% ------- Plot bars with bold per-animal and light per-day lines (in r) -------
% ------- Plot bars with per-rat colors (dark=rat mean, light=per-day), in r -------
figure('Color','w','Position',[200 200 820 480]); hold on
bar(1, 0, 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k');   % legend color
bar(2, 0, 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k');
xticks([1 2]); xticklabels({'within','across'}); ylabel('Correlation r'); box on
title(titleStr);

% convert the per-day z to r and keep a day index per animal
[Zw, Za, aid, did] = pullZ(R);           % << note: we’ll update pullZ to also return did
rw = tanh(Zw);  ra = tanh(Za);

nA   = numel({R.animal});
cmap = lines(nA);
make_light = @(c,frac) (1 - frac)*c + frac*[1 1 1];

grp_w = []; grp_a = [];
for i = 1:nA
    m = (aid==i) & isfinite(rw) & isfinite(ra);
    if ~any(m), continue; end
    c_base  = cmap(i,:);
    c_light = make_light(c_base, 0.60);

    % per-day light lines/dots
    idx = find(m);
    for j = 1:numel(idx)
        plot([1 2], [rw(idx(j)) ra(idx(j))], '-', 'Color', c_light, 'LineWidth', 1);
        plot(1, rw(idx(j)), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
        plot(2, ra(idx(j)), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
    end

    % per-rat means (dark)
    w_mu = mean(rw(m),'omitnan');
    a_mu = mean(ra(m),'omitnan');
    plot([1 2], [w_mu a_mu], '-', 'Color', c_base, 'LineWidth', 2.5);
    plot(1, w_mu, 'o', 'MarkerFaceColor', c_base, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 6);
    plot(2, a_mu, 'o', 'MarkerFaceColor', c_base, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 6);

    grp_w(end+1) = w_mu; %#ok<AGROW>
    grp_a(end+1) = a_mu; %#ok<AGROW>
end

% overlay group means
if ~isempty(grp_w), bar(1, mean(grp_w,'omitnan'), 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k'); end
if ~isempty(grp_a), bar(2, mean(grp_a,'omitnan'), 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k'); end
yline(0,'k:'); ylim auto
end


% ---------- NaN-safe Pearson correlation with guards for low variance ----------
% ---------- NaN-safe Pearson correlation (flexible signature) ----------
% Usage patterns supported:
%   r           = safe_corr(a,b)
%   [r,n]       = safe_corr(a,b)
%   r           = safe_corr(a,b,minN,minStd)
%   [r,n]       = safe_corr(a,b,minN,minStd)
%
% - minN   : minimum overlapping finite samples required (default 3)
% - minStd : minimum std threshold; if either vector's std < minStd,
%            returns r=0 (instead of NaN from corr) (default 1e-12)
function varargout = safe_corr(a,b,varargin)

    % defaults
    if numel(varargin) >= 1 && ~isempty(varargin{1}), minN   = varargin{1}; else, minN   = 3;     end
    if numel(varargin) >= 2 && ~isempty(varargin{2}), minStd = varargin{2}; else, minStd = 1e-12; end

    % basic finite mask & count
    m = isfinite(a) & isfinite(b);
    nC = nnz(m);

    if nC < minN
        r = NaN;
    else
        a = a(m); b = b(m);

        sa = std(a); sb = std(b);
        if ~isfinite(sa) || ~isfinite(sb) || sa < minStd || sb < minStd
            % treat constant vectors as zero-similarity rather than NaN
            r = 0;
        else
            r = corr(a, b, 'type','Pearson');
        end
    end

    % outputs: 1st=r, optional 2nd=nC
    if nargout <= 1
        varargout = {r};
    else
        varargout = {r, nC};
    end
end

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

% ---------- tiny helpers used in B and C ----------
function ii = clamp_idx(ii, n)
% Clamp indices to [1..n] and drop NaNs/zeros/negatives.
    if isempty(ii), ii = []; return; end
    ii = ii(:);
    ii = ii(isfinite(ii));           % drop NaN/Inf
    ii = round(ii);                  % ensure integer
    ii(ii < 1) = 1;
    ii(ii > n) = n;
end

function chunks = split_by_counts_ptr(pool, counts)
% Split a vector of frame indices "pool" into cell chunks of given sizes.
% Extra safety: if pool is short, we truncate later counts gracefully.
    counts = counts(:);
    counts(~isfinite(counts)) = 0;
    counts = max(0, round(counts));
    chunks = cell(numel(counts),1);
    if isempty(pool) || ~any(counts), return; end
    p = 1;
    for i = 1:numel(counts)
        c = counts(i);
        if c <= 0 || p > numel(pool)
            chunks{i} = [];
        else
            q = min(p + c - 1, numel(pool));
            chunks{i} = pool(p:q);
            p = q + 1;
        end
    end
end

function pi = rand_derangement_idx(K)
% Return a derangement of 1:K (no fixed points). For K<=2, fall back safely.
    if K <= 1
        pi = 1:K;
        return
    elseif K == 2
        pi = [2 1];
        return
    end
    % Sattolo’s algorithm for a single-cycle permutation (guaranteed no fixed points)
    pi = 1:K;
    for i = K:-1:2
        j = randi([1 i-1]);   % j < i
        tmp = pi(i);
        pi(i) = pi(j);
        pi(j) = tmp;
    end
    % Rarely, Sattolo can still produce a fixed point if implemented incorrectly;
    % sanity check and repair (very low probability path).
    fixed = find(pi == 1:K);
    if ~isempty(fixed)
        % swap any fixed item with its neighbor
        for f = fixed(:).'
            s = mod(f, K) + 1;
            [pi(f), pi(s)] = deal(pi(s), pi(f));
        end
    end
end

function stats = cdf_acrossTrials_within_vs_across_AT(R, Nperm)
% CDF of across-trials (AT) WITHIN vs ACROSS distributions (r), with
% pooled KS, right-tailed Wilcoxon (within > across?), Cliff's delta,
% Δz, and a hierarchical sign-flip permutation across animals.
if nargin<2, Nperm = 49; end

A = numel(R);
within_by_rat = cell(A,1);
across_by_rat = cell(A,1);

% Collect per-rat distributions across days
for i = 1:A
    if ~isfield(R(i),'AT') || ~isfield(R(i).AT,'byDay'), continue; end
    ATi = R(i).AT.byDay;

    rW = [];
    if isfield(ATi,'within') && isfield(ATi.within,'list_r') && iscell(ATi.within.list_r)
        for d = 1:numel(ATi.within.list_r)
            v = ATi.within.list_r{d};
            if ~isempty(v), rW = [rW; v(:)]; end %#ok<AGROW>
        end
    end

    rA = [];
    if isfield(ATi,'across') && isfield(ATi.across,'list_r') && iscell(ATi.across.list_r)
        for d = 1:numel(ATi.across.list_r)
            v = ATi.across.list_r{d};
            if ~isempty(v), rA = [rA; v(:)]; end %#ok<AGROW>
        end
    end

    within_by_rat{i} = rW(isfinite(rW));
    across_by_rat{i} = rA(isfinite(rA));
end

% Keep rats that have both dists
keep = cellfun(@(x) ~isempty(x), within_by_rat) & cellfun(@(x) ~isempty(x), across_by_rat);
within_by_rat = within_by_rat(keep);
across_by_rat = across_by_rat(keep);
Aeff = numel(within_by_rat);

if Aeff==0
    warning('AT CDF: no animals with both WITHIN and ACROSS distributions.');
    stats = struct();
    return
end

% Pooled vectors (for CDF & simple tests)
rW_all = vertcat(within_by_rat{:});
rA_all = vertcat(across_by_rat{:});

% ---- Plot CDFs ----
figure('Color','w','Position',[220 220 720 540]); hold on
for i = 1:Aeff
    [fW,xW] = ecdf(within_by_rat{i}); plot(xW,fW,'-','Color',[0.7 0.85 1.0],'LineWidth',1);
    [fA,xA] = ecdf(across_by_rat{i}); plot(xA,fA,'-','Color',[1.0 0.8 0.7],'LineWidth',1);
end
[fW,xW] = ecdf(rW_all);  plot(xW,fW,'-','Color',[0.1 0.35 0.9],'LineWidth',2.5);
[fA,xA] = ecdf(rA_all);  plot(xA,fA,'-','Color',[0.85 0.35 0.15],'LineWidth',2.5);
xlabel('PV correlation r (AT: across trials)'); ylabel('F(r \leq x)');
title(sprintf('AT CDFs: WITHIN vs ACROSS  (rats=%d, nW=%d, nA=%d)', Aeff, numel(rW_all), numel(rA_all)));
legend({'rat within','rat across','pooled within','pooled across'}, 'Location','southeast'); legend boxoff
grid on; xlim([-1 1]); ylim([0 1]);

% ---- Simple pooled tests (anti-conservative) ----
try, [~, p_ks]   = kstest2(rW_all, rA_all);          catch, p_ks = NaN; end
try, p_wil = ranksum(rW_all, rA_all,'tail','right');  catch, p_wil = NaN; end % within > across?

% Cliff's delta on r and Δz on Fisher z
cliff = local_cliffs_delta(rW_all, rA_all);
zW = atanh(bound_r(rW_all)); zA = atanh(bound_r(rA_all));
delta_z = mean(zW,'omitnan') - mean(zA,'omitnan');

% ---- Hierarchical permutation across rats: sign-flip per-rat Δz ----
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

% ---- Pack ----
stats = struct();
stats.pooled = struct('KS_p', p_ks, 'ranksum_right_p', p_wil, ...
                      'CliffsDelta_r', cliff, 'delta_z', delta_z, ...
                      'n_within', numel(rW_all), 'n_across', numel(rA_all));
stats.hier_perm = struct('obs_delta_z', obs, 'p_right', p_right, 'p_two', p_two, ...
                         'r_diff_approx', tanh(obs), 'rats', Aeff);

fprintf('\nAT CDF (WITHIN vs ACROSS across trials):\n');
fprintf('  Pooled: KS p=%.3g, Wilcoxon p_right=%.3g, CliffΔ=%.3f, Δz=%.3f (r≈%.3f)\n', ...
        p_ks, p_wil, cliff, delta_z, tanh(delta_z));
fprintf('  Hierarchical permutation: obs=%.3f (r≈%.3f), p_right=%.4g, p_two=%.4g, rats=%d\n\n', ...
        obs, tanh(obs), p_right, p_two, Aeff);

% ---- small helpers ----
function d = local_cliffs_delta(x,y)
    x = x(isfinite(x)); y = y(isfinite(y));
    if isempty(x)||isempty(y), d=NaN; return; end
    nx=numel(x); ny=numel(y);
    allv=[x(:); y(:)]; [ranks,~]=tiedrank(allv);
    rx=ranks(1:nx); U = nx*ny + nx*(nx+1)/2 - sum(rx);
    d = (2*U)/(nx*ny) - 1;
end
function r = bound_r(r), r = max(min(r,0.999999),-0.999999); end
end
