function R = run_task_to_space_interference(ratNames, varargin)
% RUN_TASK_TO_SPACE_INTERFERENCE  (C) Space stability with vs without task
% Time-binned version so you can compare against the space→task analysis.
%
%   CTRL↔TASK similarity is computed per (temporal bin, spatial bin), then
%   averaged on the Fisher-z scale within each day. WITHOUT-task is the
%   CTRL split-half reliability in the same spatial bins.
%
% CONTROL spatial PVs are computed from NON-task frames, excluding
%   [CS .. CS+TraceWin(2)+BufferPost] for every trial. Optional speed mask.
%
% Usage:
%   R = run_task_to_space_interference({'rat0222','rat0307'}, 'DoPlots',false);
%
% Required fields per rat (in base workspace):
%   rat.Ca_peaks, rat.Ca_ts, rat.pos, rat.CS_times
%
% Key options (name/value pairs):
%   'DaysMode'        : {'last3toAn'|'all'|'last3byca'}  (default 'last3toAn')
%   'TraceWin'        : [t0 t1] seconds from CS          (default [0 2])
%   'BufferPost'      : seconds excluded after TraceWin  (default 0)
%   'CellNorm'        : {'demean'|'zscore'|'none'}       (default 'demean')
%   'GridRC'          : [rows cols] for spatial grid     (default [3 2])
%   'NumBins'         : integer, overrides GridRC        (default [])
%
%   Frame-aligned gates (sampling ~7.5 Hz):
%   'MinCtrlFrames'   : min CTRL frames per spatial bin          (default 2)
%   'MinTraceFrames'  : min TASK frames per (temporal, spatial)  (default 2)
%   'TimeBins'        : number of temporal bins in TraceWin      (default 15)
%   'FramesPerBin'    : override: fixed frames per temporal bin  (default [])
%
%   CTRL kinematics:
%   'VelThresh'       : cm/s threshold for speed mask    (default 4)
%   'UseSpeedMask'    : bool                             (default true)
%
%   Permutation (C):
%   'NPerm_C'         : permutations for group null      (default 49)
%   'NullMode_C'      : {'frame-redistribute'|'ctrl-derange-bins'} (default 'frame-redistribute')
%
% Output per animal:
%   R(i).C.* and R(i).meta.*
%
% Notes:
%   - Temporal binning is by frame index, not seconds. If FramesPerBin is
%     empty, frames inside TraceWin are split as evenly as possible into
%     TimeBins chunks per trial.
%   - Day-level z is the mean of per-(b,k) Fisher z values.

% ------------------------------- Options ---------------------------------
p = inputParser;
addParameter(p,'DaysMode','last3toAn');
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'BufferPost',0);
addParameter(p,'GridRC',[2 2]);
addParameter(p,'NumBins',[]);

addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',true);
addParameter(p,'CellNorm','demean');  % 'zscore'|'demean'|'none'

% frame-based gates
addParameter(p,'MinCtrlFrames',2);
addParameter(p,'MinTraceFrames',2);
addParameter(p,'TimeBins',15);
addParameter(p,'FramesPerBin',[]);

% (C) perms
addParameter(p,'NPerm_C',500);
addParameter(p,'NullMode_C','frame-redistribute');

addParameter(p,'DoPlots',true);
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

    % ----- Select days (by Ca_* fieldnames via autoDateList) -----
    dateList = autoDateList(rat); % external helper
    switch lower(opt.DaysMode)
        case {'all','last3byca'}
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

    % ----- Pull & standardize per-day containers -----
    spikes_raw = filterFieldsByDay(rat.Ca_peaks, daysToUse); % external
    ts_raw     = filterFieldsByDay(rat.Ca_ts,    daysToUse);
    pos_raw    = filterFieldsByDay(rat.pos,      daysToUse);
    cs_raw     = filterFieldsByDay(rat.CS_times, daysToUse);
    [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, daysToUse);

    % ----- CONTROL spatial PVs (exclude Task+[0 BufferPost]) -----
    ctrl = compute_control_spatialPV(spikes, ts, pos, cs, ...
        'TraceWin',        opt.TraceWin, ...
        'BufferPost',      opt.BufferPost, ...
        'GridRC',          opt.GridRC, ...
        'NumBins',         opt.NumBins, ...
        'MinCtrlFrames',   opt.MinCtrlFrames, ...
        'VelThresh',       opt.VelThresh, ...
        'UseSpeedMask',    opt.UseSpeedMask, ...
        'CellNorm',        opt.CellNorm, ...
        'Label',           ratVar);

    % ----- Space stability WITH vs WITHOUT task (time-binned) -----
    C = compute_space_with_vs_without_task(spikes, ts, pos, cs, ctrl, ...
        'TraceWin',        opt.TraceWin, ...
        'VelThresh',       opt.VelThresh, ...
        'UseSpeedMask',    opt.UseSpeedMask, ...
        'CellNorm',        opt.CellNorm, ...
        'MinTraceFrames',  opt.MinTraceFrames, ...
        'TimeBins',        opt.TimeBins, ...
        'FramesPerBin',    opt.FramesPerBin, ...
        'NPerm',           opt.NPerm_C, ...
        'NullMode',        opt.NullMode_C, ...
        'Label',           ratVar);

    % Per-day paired z stats
    C = compute_space_with_without_binwise_stats(C, 'NPerm', 49);
    R(ii).C = C;

    % ----- Pack meta -----
    R(ii).animal = ratVar;
    R(ii).meta.options = opt;
    R(ii).meta.days = daysToUse;
    R(ii).meta.ctrl_reliability = ctrl.reliability;
    R(ii).meta.grid = ctrl.grid;
end

% ================================ PLOTS ==================================
summarize_taskSpacePV_clean(R);
end

% =========================================================================
% ========================== HELPER FUNCTIONS =============================
% =========================================================================

% ---------- CTRL spatial PVs (non-task frames) ----------
function CTRL = compute_control_spatialPV(spikes, ts, pos, cs, varargin)
p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'BufferPost',0);
addParameter(p,'GridRC',[3 2]);
addParameter(p,'NumBins',[]);
addParameter(p,'MinCtrlFrames',2);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'Label','');
addParameter(p,'CellNorm','demean');
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);

% ---- Build grid edges (per-day or pooled) ----
if ~isempty(o.NumBins)
    rc_eff = best_factors(o.NumBins);
else
    rc_eff = o.GridRC;
end
K = rc_eff(1)*rc_eff(2);

edgesByDay = cell(1,D);
for d = 1:D
    [edgesByDay{d}] = build_grid_edges_single_day(pos{d}, rc_eff);
end

% ---- CONTROL PVs per day ----
PV_by  = cell(1,D);
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
    x_u = x(use);
    y_u = y(use);

    % bin assignments (spatial)
    [~, b] = pos2bin(x_u, y_u, edgesByDay{d});

    % PV per bin (mean over frames)
    Nc_d = size(S,1);
    PV   = nan(Nc_d, K);
    for k = 1:K
        idx = find(b == k);
        if numel(idx) >= o.MinCtrlFrames
            PV(:,k) = mean(S(:,idx), 2, 'omitnan');  % frame-aligned mean
        end
    end

    % split-half reliability within CTRL (even/odd frames inside each bin)
    rel = nan(K,1);
    for k = 1:K
        idx = find(b == k);
        if numel(idx) < max(4, 2*o.MinCtrlFrames), continue; end
        idx_even = idx(2:2:end);
        idx_odd  = idx(1:2:end);
        v1 = mean(S(:,idx_even), 2, 'omitnan');
        v2 = mean(S(:,idx_odd ), 2, 'omitnan');
        rel(k) = safe_corr(v1, v2);
    end

    PV_by{d}  = PV;
    rel_by{d} = rel;

    kept = all(isfinite(PV),1);
    fprintf('[%s][CTRL d=%d] kept bins %d/%d | split-half median r=%.2f\n', L, d, nnz(kept), K, median(rel,'omitnan'));
end

% ---- Pack CTRL ----
CTRL.PV    = PV_by;                 % {day} : [Nc × K]
CTRL.grid.edges   = edgesByDay;     % {day} edges
CTRL.grid.GridRC  = rc_eff;
CTRL.grid.K       = K;

CTRL.reliability_byDay = rel_by;    % {day}[K×1]
CTRL.reliability.r = vertcat(rel_by{:});
CTRL.reliability.z = atanh(max(min(CTRL.reliability.r,0.99999),-0.99999));

CTRL.params.MinCtrlFrames = o.MinCtrlFrames;
CTRL.params.NumBins       = o.NumBins;
CTRL.params.GridRC_eff    = rc_eff;
end

% ---------- WITH vs WITHOUT task (C) : TIME-BINNED ----------
function C = compute_space_with_vs_without_task(spikes, ts, pos, cs, CTRL, varargin)
p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',true);
addParameter(p,'CellNorm','demean');
addParameter(p,'MinTraceFrames',2);
addParameter(p,'TimeBins',15);
addParameter(p,'FramesPerBin',[]);
addParameter(p,'NPerm',49);
addParameter(p,'NullMode','frame-redistribute');
addParameter(p,'Label','');
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts); K = CTRL.grid.K;

r_bin_byDay   = cell(1,D);         % mean r over time-bins per spatial bin (for display)
z_day         = nan(D,1);          % day-level mean of all (b,k) Fisher z
n_valid_day   = zeros(D,1);
perZ_list_byDay = cell(1,D);       % store list of per-(b,k) z's for sanity checks

for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);
    S = spikes_to_matrix(spikes{d}, t);
    S = normalize_cells(S, o.CellNorm);

    % speed mask
    if o.UseSpeedMask
        v = speed_cm_per_s(pos{d});
        v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
        speed_ok = (v >= o.VelThresh);
    else
        speed_ok = true(size(t));
    end

    edges_d = CTRL.grid.edges{d};
    PV_ctrl = CTRL.PV{d};
    if isempty(PV_ctrl), continue; end
    Nc = size(PV_ctrl,1);

    % ---- Build per-(trial, temporal-bin) frame indices within TraceWin ----
    csd = cs{d}(:);
    B = o.TimeBins;
    listZ = [];                    % all per-(b,k) Fisher z for this day
    r_bin_accum = nan(K, B);       % r for each spatial bin and time bin

    for tr = 1:numel(csd)
        t0 = csd(tr) + o.TraceWin(1);
        t1 = csd(tr) + o.TraceWin(2);
        idx_all = find(t >= t0 & t < t1);
        if isempty(idx_all), continue; end

        % frame-aligned temporal splits
        if ~isempty(o.FramesPerBin)
            step = max(1, round(o.FramesPerBin));
            edges_idx = 1:step:(numel(idx_all)+1);
            if edges_idx(end) ~= numel(idx_all)+1
                edges_idx(end+1) = numel(idx_all)+1;
            end
        else
            edges_idx = round(linspace(1, numel(idx_all)+1, B+1));
            edges_idx(1) = 1; edges_idx(end) = numel(idx_all)+1;
        end
        nTB = numel(edges_idx)-1;   % actual number of temporal bins this trial

        for bti = 1:nTB
            if edges_idx(bti) >= edges_idx(bti+1), continue; end
            seg = idx_all(edges_idx(bti):edges_idx(bti+1)-1);
            if o.UseSpeedMask, seg = seg(speed_ok(seg)); end
            if numel(seg) < o.MinTraceFrames, continue; end

            % spatial bin per frame
            [~, b_space] = pos2bin(x(seg), y(seg), edges_d);

            % per-bin TASK PV for this temporal chunk
            for k = 1:K
                kk = (b_space == k);
                if nnz(kk) >= o.MinTraceFrames && all(isfinite(PV_ctrl(:,k)))
                    pv_task = mean(S(:, seg(kk)), 2, 'omitnan');
                    r = safe_corr(PV_ctrl(:,k), pv_task);
                    r_bin_accum(k, bti) = nanmean([r_bin_accum(k, bti), r]); % accumulate across trials if same (b,k)
                end
            end
        end
    end

    % per-day aggregation on Fisher-z over all valid (b,k)
    r_valid = r_bin_accum(isfinite(r_bin_accum));
    if isempty(r_valid)
        r_bin_byDay{d} = nan(K,1);
        n_valid_day(d) = 0;
        continue;
    end

    z_list = atanh(max(min(r_bin_accum(:),0.99999),-0.99999));
    z_list = z_list(isfinite(z_list));
    perZ_list_byDay{d} = z_list;
    z_day(d) = mean(z_list, 'omitnan');

    % also keep a simple spatial-bin mean r (over time bins) for display
    r_spatial = nan(K,1);
    for k = 1:K
        rv = r_bin_accum(k, :);
        rv = rv(isfinite(rv));
        if ~isempty(rv)
            r_spatial(k) = tanh(mean(atanh(max(min(rv,0.99999),-0.99999))));
        end
    end
    r_bin_byDay{d} = r_spatial;
    n_valid_day(d) = numel(z_list);
end

C.bin.byDay.r        = r_bin_byDay;          % mean r per spatial bin (over time)
C.bin.byDay.z_mean   = z_day;                % day-level Fisher-z mean over (b,k)
C.bin.byDay.n_valid  = n_valid_day;
C.bin.byDay.z_list   = perZ_list_byDay;      % stash raw per-(b,k) z
C.bin.group.z_mean   = mean(z_day(isfinite(z_day)),'omitnan');
C.bin.group.r_mean   = tanh(C.bin.group.z_mean);

% WITHOUT task = CTRL split-half reliability (precomputed if available)
with_byDay    = cell(1,D);
without_byDay = cell(1,D);
z_ctrl_day    = nan(D,1);

for d = 1:D
    r_task_bins = C.bin.byDay.r{d};
    if isempty(r_task_bins), continue; end

    if isfield(CTRL,'reliability_byDay') && numel(CTRL.reliability_byDay) >= d ...
            && ~isempty(CTRL.reliability_byDay{d})
        r_ctrl = CTRL.reliability_byDay{d}(:);
    else
        r_ctrl = nan(CTRL.grid.K,1);
    end

    mask = isfinite(r_task_bins) & isfinite(r_ctrl);
    if ~any(mask), continue; end
    with_byDay{d}    = r_task_bins(mask);
    without_byDay{d} = r_ctrl(mask);

    z_ctrl_day(d) = mean(atanh(max(min(without_byDay{d},0.99999),-0.99999)), 'omitnan');
end

% day-level delta (WITH minus WITHOUT) in Fisher-z
delta_day = z_day - z_ctrl_day;      % [D x 1], may have NaNs
C.bin.byDay.delta_z = delta_day;
C.bin.group.delta_z_mean = mean(delta_day(isfinite(delta_day)), 'omitnan');
C.bin.group.delta_r_mean = tanh(C.bin.group.delta_z_mean);



% --- PERMUTATION: spatial-bin derangement null for Δz (per day, then group) ---
P = o.NPerm;
perm_byDay = cell(D,1);
if P > 0
    for d = 1:D
        if ~isfinite(z_day(d)) || ~isfinite(z_ctrl_day(d)), continue; end

        % Collect the r values per spatial bin across time-bins that built z_day(d).
        % We saved per-day, per-spatial mean r in C.bin.byDay.r{d}; we also saved the raw
        % per-(b,k) z in C.bin.byDay.z_list{d}, but for a clean spatial derangement we
        % re-evaluate z_with from the r_bin_accum we just used.
        % To do that, recreate a matrix of r over (k, time-bin) for this day:
        % (We cached r_bin_accum as 'r_bin_byDay{d}' only after averaging over time-bins;
        % so we’ll rebuild a minimal proxy from S and CTRL now.)

        % Lightweight reconstitution: use stored CTRL PVs + day’s TraceWin frames and edges
        % to compute per-(k,b) TASK PVs again (same as above but without speed/occupancy checks)
        % Then derange CTRL’s columns and recompute z_with^perm.
        edges_d = CTRL.grid.edges{d};
        PV_ctrl = CTRL.PV{d};
        if isempty(PV_ctrl), continue; end
        Nc = size(PV_ctrl,1);

        % Build a list of (k,b) task vectors actually used earlier:
        csd = cs{d}(:);
        B = o.TimeBins;
        % Collect per-(k,b) TASK PV means in a cell; skip ones that would have been invalid
        TASK_kb = cell(CTRL.grid.K, B);
        t = ts{d}(:);
        [x, y] = interp_pos(pos{d}, t);
        S = spikes_to_matrix(spikes{d}, t); %#ok<NBRAK>
        S = normalize_cells(S, o.CellNorm);
        if o.UseSpeedMask
            v = speed_cm_per_s(pos{d});
            v = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
            speed_ok = (v >= o.VelThresh);
        else
            speed_ok = true(size(t));
        end

        for tr = 1:numel(csd)
            t0 = csd(tr) + o.TraceWin(1);
            t1 = csd(tr) + o.TraceWin(2);
            idx_all = find(t >= t0 & t < t1);
            if isempty(idx_all), continue; end
            edges_idx = round(linspace(1, numel(idx_all)+1, B+1));
            edges_idx([1 end]) = [1, numel(idx_all)+1];
            for bti = 1:(numel(edges_idx)-1)
                seg = idx_all(edges_idx(bti):edges_idx(bti+1)-1);
                if o.UseSpeedMask, seg = seg(speed_ok(seg)); end
                if numel(seg) < o.MinTraceFrames, continue; end
                [~, b_space] = pos2bin(x(seg), y(seg), edges_d);
                for k = 1:CTRL.grid.K
                    hit = (b_space == k);
                    if nnz(hit) >= o.MinTraceFrames && all(isfinite(PV_ctrl(:,k)))
                        vec = mean(S(:, seg(hit)), 2, 'omitnan');
                        % concatenate across trials (we'll average later in r / z)
                        if isempty(TASK_kb{k,bti})
                            TASK_kb{k,bti} = vec;
                        else
                            % average PVs across trials for same (k,b)
                            TASK_kb{k,bti} = nanmean([TASK_kb{k,bti}, vec], 2);
                        end
                    end
                end
            end
        end

        % observed z_with_day(d) was already computed; now build perm nulls by deranging CTRL bins
        validK = find(all(isfinite(PV_ctrl),1));
        if numel(validK) < 2, continue; end

        perm_d = nan(P,1);
        for pidx = 1:P
            % derangement of validK
            map = validK(randperm(numel(validK)));
            if any(map==validK), % enforce derangement
                for rep = 1:5
                    map = validK(randperm(numel(validK)));
                    if ~any(map==validK), break; end
                end
                if any(map==validK)  % last resort rotate
                    map = validK([2:end,1]);
                end
            end
            % recompute r over all available (k,b) using CTRL(:, map(k)) vs TASK_kb{k,b}
            r_list = [];
            for bti = 1:B
                for ii = 1:numel(validK)
                    k = validK(ii); kp = map(ii);
                    tv = TASK_kb{k,bti};
                    if ~isempty(tv)
                        r_list(end+1,1) = safe_corr(PV_ctrl(:,kp), tv); %#ok<AGROW>
                    end
                end
            end
            if isempty(r_list), continue; end
            zw_perm = atanh(max(min(r_list,0.99999),-0.99999));
            z_with_perm = mean(zw_perm,'omitnan');
            perm_d(pidx) = z_with_perm - z_ctrl_day(d);
        end
        perm_byDay{d} = perm_d;
    end

    % --- Δ-based permutation test (bin-derangement) --------------------------
    % observed animal-level statistic = mean Δz across this animal's valid days
    delta_obs = mean(delta_day(isfinite(delta_day)), 'omitnan');

    % build a permutation distribution of the same statistic from perm_byDay
    % (each entry perm_byDay{d} is a [P x 1] vector of Δz^perm for day d)
    % choose the common number of perms Pcommon across days
    valid = ~cellfun('isempty', perm_byDay);
    if any(valid)
        Pcommon = min(cellfun(@(v) sum(isfinite(v)), perm_byDay(valid)));
    else
        Pcommon = 0;
    end

    if Pcommon >= 5
        delta_perm_anim = nan(Pcommon,1);
        for p = 1:Pcommon
            dvals = nan(sum(valid),1);
            idx = 1;
            for d = 1:D
                if ~valid(d), continue; end
                if numel(perm_byDay{d}) >= p && isfinite(perm_byDay{d}(p))
                    dvals(idx) = perm_byDay{d}(p); idx = idx + 1;
                end
            end
            delta_perm_anim(p) = mean(dvals(isfinite(dvals)), 'omitnan');  % same stat as delta_obs
        end

        pr = mean(delta_perm_anim >= delta_obs);
        pl = mean(delta_perm_anim <= delta_obs);
        pt = 2*min(pr, pl);
    else
        delta_perm_anim = [];
        pr = NaN; pl = NaN; pt = NaN;
    end

    fprintf('p value against shuff')
    pt

    % store Δ-based null + p-values
    C.delta = struct();
    C.delta.obs_z     = delta_obs;
    C.delta.obs_r     = tanh(delta_obs);
    C.delta.perm_z    = delta_perm_anim;
    C.delta.p_right   = pr;
    C.delta.p_left    = pl;
    C.delta.p_two     = pt;
end

% --- Make WITH/WITHOUT daywise r available to downstream stats/plots ---
C.byDay.withTask.r    = with_byDay;      % cell{day} of r per spatial bin
C.byDay.withoutTask.r = without_byDay;   % cell{day} of r per spatial bin (CTRL split-half)
end

% ---------- (C) per-day paired z + pooled + sign-flip ----------
function C = compute_space_with_without_binwise_stats(C, varargin)
p = inputParser; addParameter(p,'NPerm',49); parse(p,varargin{:}); opt = p.Results;
if ~isfield(C,'byDay') || ~isfield(C.byDay,'withTask') || ~isfield(C.byDay,'withoutTask')
    warning('C.byDay.withTask/withoutTask not found.'); return
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
        [~, p, ~, st] = ttest(zw, zu); dz = mean(delta)/std(delta);
        dayStats(d) = struct('nPairs',numel(delta),'t',st.tstat,'df',st.df,'p',p,'dz',dz);
    else
        dayStats(d).nPairs = numel(delta);
    end
    Z_with_all    = [Z_with_all; zw(:)]; %#ok<AGROW>
    Z_without_all = [Z_without_all; zu(:)]; %#ok<AGROW>
end
mask = isfinite(Z_with_all) & isfinite(Z_without_all);
if nnz(mask) >= 2
    [~, p_all, ~, st_all] = ttest(Z_with_all(mask), Z_without_all(mask));
    dz_all = mean(Z_with_all(mask) - Z_without_all(mask)) ./ std(Z_with_all(mask) - Z_without_all(mask));
else
    p_all = NaN; st_all = struct('tstat',NaN,'df',NaN); dz_all = NaN;
end
% cluster-aware permutation: sign-flip on per-day mean Δz
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
day_mu = day_mu(isfinite(day_mu)); obs_cluster_mean = mean(day_mu,'omitnan');
NPerm = opt.NPerm; perm_stat = nan(NPerm,1);
for b = 1:NPerm
    s = (rand(size(day_mu))>0.5)*2 - 1; perm_stat(b) = mean(s .* day_mu);
end
p_right = (sum(perm_stat >= obs_cluster_mean) + 1) / (nnz(isfinite(perm_stat)) + 1);
p_left  = (sum(perm_stat <= obs_cluster_mean) + 1) / (nnz(isfinite(perm_stat)) + 1);
p_two   = 2*min(p_right, p_left);
if ~isfield(C,'distStats'), C.distStats = struct(); end
C.distStats.perDay  = dayStats;
C.distStats.pooled  = struct('t',st_all.tstat,'df',st_all.df,'p',p_all,'dz',dz_all, 'nBins', nnz(mask));
C.distStats.cluster_perm = struct('NPerm',NPerm,'obs_meanDelta_z',obs_cluster_mean, ...
                                  'p_right',p_right,'p_left',p_left,'p_two',p_two);
end

% ========================== UTILITIES ==========================
function [edges] = build_grid_edges_single_day(posd, GridRC)
allx = posd.x(:); ally = posd.y(:);
allx = allx(isfinite(allx));  ally = ally(isfinite(ally));
if isempty(allx) || isempty(ally)
    edges.x = linspace(0, 1, GridRC(2)+1);
    edges.y = linspace(0, 1, GridRC(1)+1);
else
    edges.x = linspace(min(allx), max(allx), GridRC(2)+1);
    edges.y = linspace(min(ally), max(ally), GridRC(1)+1);
end
end

function rc = best_factors(N)
validateattributes(N, {'numeric'}, {'scalar','integer','positive','finite'});
r = floor(sqrt(double(N)));
while r > 1 && mod(N, r) ~= 0, r = r - 1; end
c = N / r; rc = [r, c];
end

function [x_i, y_i] = interp_pos(posd, t)
tt = double(posd.t(:)); xx = double(posd.x(:)); yy = double(posd.y(:));
t = double(t(:));
[ttu, ia] = unique(tt, 'stable'); xxu = xx(ia); yyu = yy(ia);
x_i = interp1(ttu, xxu, t, 'linear','extrap');
y_i = interp1(ttu, yyu, t, 'linear','extrap');
end

function v = speed_cm_per_s(posd)
t  = double(posd.t(:)); x  = double(posd.x(:)); y  = double(posd.y(:));
n  = min([numel(t), numel(x), numel(y)]);
dt = diff(t(1:n)); dt(end+1,1) = median(dt(dt>0),'omitnan');
dx = [diff(x(1:n)); 0]; dy = [diff(y(1:n)); 0];
v = hypot(dx,dy) ./ max(dt, eps);
end

function [rc_idx, k] = pos2bin(x, y, edges)
cx = discretize(x, edges.x); cy = discretize(y, edges.y);
GridR = numel(edges.y)-1; GridC = numel(edges.x)-1;
bad = isnan(cx) | isnan(cy) | cx<1 | cy<1 | cx>GridC | cy>GridR;
cx(bad) = NaN; cy(bad) = NaN;
rc_idx = [cy, cx]; k = nan(size(x));
m = isfinite(cx) & isfinite(cy);
if any(m), k(m) = sub2ind([GridR, GridC], cy(m), cx(m)); end
end

function [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, dayKeys)
spikes = to_daycells(spikes_raw, dayKeys);
ts_in  = to_daycells(ts_raw,   dayKeys);
pos_in = to_daycells(pos_raw,  dayKeys);
cs_in  = to_daycells(cs_raw,   dayKeys);
D = numel(dayKeys);
ts  = cell(1,D); pos = cell(1,D); cs  = cell(1,D);
for d = 1:D
    ts{d}  = coerce_ts_day(ts_in{d});
    pos{d} = coerce_pos_day(pos_in{d}, ts{d});
    cs{d}  = coerce_cs_day(cs_in{d}, ts{d});
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
if isnumeric(td) && isvector(td)
    t = double(td(:)); if max(t,[],'omitnan')>1e4, t=t/1000; end; return
end
if isnumeric(td) && ismatrix(td) && ~isscalar(td)
    M=double(td); t = M(:, min(2,size(M,2))); if max(t,[],'omitnan')>1e4, t=t/1000; end; t=t(:); return
end
if istable(td)
    vn=lower(string(td.Properties.VariableNames));
    cname = vn( find(ismember(vn, ["time_ms","time","ts","t","timestamp","frame_ts","ca_ts"]),1) );
    if isempty(cname), error('ts table: no time-like column found'); end
    t=double(td{:,find(vn==cname,1)}); if contains(cname,"ms")||max(t,[],'omitnan')>1e4, t=t/1000; end; t=t(:); return
end
if isstruct(td)
    pref=["time_ms","time","ts","t","timestamp","timestamps","frame_ts","ca_ts"];
    for j=1:numel(pref)
        if isfield(td, pref{j})
            v=td.(pref{j}); if isnumeric(v)&&isvector(v)&&numel(v)>1
                t=double(v(:)); if contains(pref(j),"ms")||max(t,[],'omitnan')>1e4, t=t/1000; end; return
            end
        end
    end
    error('ts struct: no numeric time vector found.');
end
error('Unsupported ts type: %s', class(td));
end

function P = coerce_pos_day(pd, tday)
if nargin<2, tday = []; end
col = @(v) v(:);
if iscell(pd), if isempty(pd), P=struct('t',[],'x',[],'y',[]); return; else, pd=pd{1}; end, end
if isempty(pd)
    if ~isempty(tday), P=struct('t',col(tday),'x',nan(numel(tday),1),'y',nan(numel(tday),1));
    else, P=struct('t',[],'x',[],'y',[]); end, return
end
if istable(pd)
    vn=lower(string(pd.Properties.VariableNames));
    tname = pick_name(vn, ["t","time","ts"]);
    xname = pick_name(vn, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"]);
    yname = pick_name(vn, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"]);
    t = []; if ~isempty(tname), t = pd{:, find(vn==tname,1)}; end
    x = []; if ~isempty(xname), x = pd{:, find(vn==xname,1)}; end
    y = []; if ~isempty(yname), y = pd{:, find(vn==yname,1)}; end
    [t,x,y] = finalize_txy(t,x,y,tday);
    P = struct('t',col(t),'x',col(x),'y',col(y)); return
end
if isnumeric(pd) && ismatrix(pd) && ~isscalar(pd)
    [t,x,y] = coerce_from_numeric(pd, tday);
    P = struct('t',col(t),'x',col(x),'y',col(y)); return
end
if isstruct(pd)
    t=[]; x=[]; y=[]; f=lower(string(fieldnames(pd)));
    t = get_field_if_exists(pd, pick_name(f, ["t","time","ts"]));
    x = get_field_if_exists(pd, pick_name(f, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"]));
    y = get_field_if_exists(pd, pick_name(f, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"]));
    [t,x,y]=finalize_txy(t,x,y,tday);
    P=struct('t',col(t),'x',col(x),'y',col(y)); return
end
error('pos day: unsupported type %s', class(pd));
end

function name = pick_name(names, options)
name=""; for k=1:numel(options), idx=find(names==options(k),1); if ~isempty(idx), name=names(idx); return; end, end
end

function val = get_field_if_exists(S, name)
if strlength(name)==0 || ~isstruct(S), val=[]; return; end
fn=fieldnames(S); idx=find(strcmpi(fn, char(name)),1);
if isempty(idx), val=[]; else, val=S.(fn{idx}); end
end

function [t,x,y] = coerce_from_numeric(M, tday)
[nr,nc]=size(M);
if nc==3, t=M(:,1); x=M(:,2); y=M(:,3);
elseif nc==2, x=M(:,1); y=M(:,2); if ~isempty(tday)&&numel(tday)==nr, t=tday; else, t=(1:nr)'; end
else, error('numeric matrix must be [n×3] or [n×2].');
end
[t,x,y] = finalize_txy(t,x,y,tday);
end

function [t,x,y] = finalize_txy(t,x,y,tday)
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
if istable(csd)
    vn=lower(string(csd.Properties.VariableNames));
    cname = pick_name(vn, ["cs_ms","cs_time","cstime","cs","onset","onsets","cs_onset","cs_time_ms","trial_cs","cue_onset","time","ts"]);
    if strlength(cname)==0, cs_vec=[]; return; end
    tcol=double(csd{:, find(vn==cname,1)}); if contains(cname,"ms")||max(tcol,[],'omitnan')>1e4, tcol=tcol/1000; end
    cs_vec=tcol(:); return
end
if isnumeric(csd) && ismatrix(csd) && ~isscalar(csd)
    M=double(csd); cs = M(:, min(2,size(M,2))); if max(cs,[],'omitnan')>1e4, cs=cs/1000; end, cs_vec=cs(:); return
end
if isstruct(csd), cs_vec=[]; return, end
cs_vec=[];
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
    pick=[]; for j=1:numel(f)
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

% ---------- Corr + normalization ----------
function varargout = safe_corr(a,b,varargin)
a = a(:);
b = b(:);
if numel(a) ~= numel(b), n = min(numel(a), numel(b)); a = a(1:n); b = b(1:n); end
if numel(varargin) >= 1 && ~isempty(varargin{1}), minN   = varargin{1}; else, minN   = 3;     end
if numel(varargin) >= 2 && ~isempty(varargin{2}), minStd = varargin{2}; else, minStd = 1e-12; end
m = isfinite(a) & isfinite(b);
nC = nnz(m);
if nC < minN
    r = NaN;
else
    a = a(m); b = b(m);
    sa = std(a); sb = std(b);
    if ~isfinite(sa) || ~isfinite(sb) || sa < minStd || sb < minStd
        r = 0;
    else
        r = corr(a, b, 'type','Pearson');
    end
end
if nargout <= 1, varargout = {r};
else, varargout = {r, nC};
end
end

function S = normalize_cells(S, mode)
switch lower(mode)
    case 'none'
    case 'demean'
        mu = mean(S,2,'omitnan'); S = S - mu;
    case 'zscore'
        mu = mean(S,2,'omitnan'); sd = std(S,0,2,'omitnan'); sd(sd==0|~isfinite(sd)) = 1; S = (S - mu) ./ sd;
    otherwise
        error('CellNorm must be ''zscore'',''demean'',''none''.');
end
end

% ---------- Small helpers used above ----------
function summarize_taskSpacePV_clean(R)
% SUMMARIZE_TASKSPACEPV_CLEAN
% Group viz + stats for (C) WITH vs WITHOUT task.
% - All averaging/testing done in Fisher-z (atanh), bars plotted in r.
% - Each rat contributes one mean (across its valid days) -> group mean.
%
% Expected fields per animal (produced by your pipeline):
%   R(i).C.bin.byDay.z_mean     : per-day WITH-task mean Fisher-z
%   R(i).C.byDay.withoutTask.r  : cell{day} of r-values (CTRL split-half) per bin
%
% This function is self-contained and includes its helpers below.

titleStr = 'C: WITH vs WITHOUT task';

% ---------- harvest per-day z (WITH) and z (WITHOUT) ----------
[Zw_all, Za_all, aid, did] = get_C_dayZ(R);  % column vectors; aligned by (animal, day)

if isempty(Zw_all) || isempty(Za_all)
    warning('No WITH/WITHOUT day-level data found to summarize.');
    return
end

animals = {R.animal}';

% ---------- PRIMARY STATS: per-rat (subject-level) ----------
% Build per-rat mean z (WITH/WITHOUT), then paired t across rats.
nA = numel(animals);
rat_mean_with_z    = nan(nA,1);
rat_mean_without_z = nan(nA,1);

for i = 1:nA
    m = (aid == i);
    if any(m)
        rat_mean_with_z(i)    = mean(Zw_all(m), 'omitnan');
        rat_mean_without_z(i) = mean(Za_all(m), 'omitnan');
    end
end

keepA = isfinite(rat_mean_with_z) & isfinite(rat_mean_without_z);
Zw_rats = rat_mean_with_z(keepA);
Za_rats = rat_mean_without_z(keepA);

% Paired t across rats (Fisher-z)
[pA, tA, dfA, dzA] = paired_t_z(Zw_rats, Za_rats);

% ---------- SECONDARY (OPTIONAL): pooled per-day paired t in z ----------
keepD = isfinite(Zw_all) & isfinite(Za_all);
Zw_days = Zw_all(keepD);
Za_days = Za_all(keepD);
[pD, tD, dfD, dzD] = paired_t_z(Zw_days, Za_days);

% ---------- PLOT ----------
figure('Color','w','Position',[160 160 900 540]); hold on

% 1) Per-day light lines, grouped by rat
cmap = lines(nA);
make_light = @(c,frac) (1-frac)*c + frac*[1 1 1];

% Convert the specific day's WITH/WITHOUT z back to r for display
rw_day = tanh(Zw_all);
ru_day = tanh(Za_all);

for i = 1:nA
    c_base  = cmap(i,:);
    c_light = make_light(c_base, 0.60);

    m = (aid == i) & isfinite(rw_day) & isfinite(ru_day);
    did_i = find(m);
    for j = 1:numel(did_i)
        plot([1 2], [rw_day(did_i(j)) ru_day(did_i(j))], '-', ...
            'Color', c_light, 'LineWidth', 1.0);
        plot(1, rw_day(did_i(j)), 'o', 'MarkerFaceColor', c_light, ...
            'MarkerEdgeColor', c_light, 'MarkerSize', 4);
        plot(2, ru_day(did_i(j)), 'o', 'MarkerFaceColor', c_light, ...
            'MarkerEdgeColor', c_light, 'MarkerSize', 4);
    end
end

% 2) Per-rat dark means (in r, from tanh of per-rat mean z)
rat_with_r    = tanh(rat_mean_with_z);
rat_without_r = tanh(rat_mean_without_z);

for i = 1:nA
    if isfinite(rat_with_r(i)) && isfinite(rat_without_r(i))
        plot([1 2], [rat_with_r(i) rat_without_r(i)], '-', 'Color', cmap(i,:), 'LineWidth', 2.5);
        plot(1, rat_with_r(i),    'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
        plot(2, rat_without_r(i), 'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
    end
end

% 3) Group bars (mean across rats in z, back-transformed)
bar(1, tanh(mean(Zw_rats,'omitnan')), 0.6, 'FaceColor',[0.35 0.70 1.00], 'EdgeColor','k');
bar(2, tanh(mean(Za_rats,'omitnan')), 0.6, 'FaceColor',[0.85 0.55 0.25], 'EdgeColor','k');

% Cosmetics
xlim([0.5 2.5]); xticks([1 2]);
xticklabels({'with task','without task'});  % <- fixed label order
ylabel('Mean PV correlation (r)');
title(titleStr);
yline(0,'k:'); grid on; box on

% Stats annotation (primary: per-rat; secondary: per-day)
yl = ylim; y = yl(2) + 0.06*range(yl);
line([1 2], [y y], 'Color','k', 'LineWidth', 1.2);
txt1 = sprintf('per-rat paired t (z): t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d', dfA, tA, pA, dzA, numel(Zw_rats));
text(1.5, y + 0.02*range(yl), txt1, 'HorizontalAlignment','center');
ylim([yl(1) y + 0.10*range(yl)]);

% Also print to console a concise summary, including the secondary per-day test
fprintf('\n== %s ==\n', titleStr);
fprintf('Primary (per-rat):  paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d\n', dfA, tA, pA, dzA, numel(Zw_rats));
fprintf('Secondary (per-day): paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nDays=%d\n', dfD, tD, pD, dzD, numel(Zw_days));

end

% ========================= HELPERS =========================

function [Zw, Za, aid, did] = get_C_dayZ(R)
% Return aligned per-day Fisher-z values:
%   Zw : WITH-task (already z in R(i).C.bin.byDay.z_mean)
%   Za : WITHOUT-task (convert day’s r cell to z by averaging within day in z)
Zw=[]; Za=[]; aid=[]; did=[];

for i = 1:numel(R)
    has_with = isfield(R(i),'C') && isfield(R(i).C,'bin') && isfield(R(i).C.bin,'byDay') ...
               && isfield(R(i).C.bin.byDay,'z_mean');
    has_without = isfield(R(i),'C') && isfield(R(i).C,'byDay') ...
                  && isfield(R(i).C.byDay,'withoutTask') && isfield(R(i).C.byDay.withoutTask,'r');

    if ~has_with || ~has_without, continue; end

    z_with = R(i).C.bin.byDay.z_mean(:);        % z per day (WITH)
    ur     = R(i).C.byDay.withoutTask.r;        % cell{day} vectors of r (WITHOUT)

    % Convert WITHOUT r-vectors to per-day Fisher-z means
    z_without = nan(numel(ur),1);
    for d = 1:numel(ur)
        if ~isempty(ur{d})
            r_d = ur{d}(:);
            r_d = r_d(isfinite(r_d));
            if ~isempty(r_d)
                z_without(d) = mean(atanh(max(min(r_d,0.999999),-0.999999)), 'omitnan');
            end
        end
    end

    % Align days that have both WITH and WITHOUT
    m = isfinite(z_with) & isfinite(z_without);
    Zw  = [Zw;  z_with(m)];
    Za  = [Za;  z_without(m)];
    aid = [aid; i*ones(nnz(m),1)];
    % for completeness (not used here), the within-animal day indices:
    did_local = (1:numel(z_with)).';
    did = [did; did_local(m)];
end
end

function [p, t, df, dz] = paired_t_z(z1, z2)
% Paired t-test in z space + Cohen's dz for paired design.
% Returns two-sided p (MATLAB ttest default), t statistic, df, and dz.
if numel(z1) ~= numel(z2)
    error('paired_t_z: inputs must be same length.');
end
mask = isfinite(z1) & isfinite(z2);
z1 = z1(mask); z2 = z2(mask);

if numel(z1) < 2
    p = NaN; t = NaN; df = 0; dz = NaN; return
end

[~, p, ~, st] = ttest(z1, z2);   % two-sided
t  = st.tstat;
df = st.df;

d  = z1 - z2;
dz = mean(d,'omitnan') / std(d, 0, 'omitnan');  % Cohen's dz (paired)
end
