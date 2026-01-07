function R = run_task_to_space_interference(ratNames, varargin)
% RUN_TASK_TO_SPACE_INTERFERENCE  Space stability with vs without task
%
% TEMPLATE version:
%   - Build a CTRL spatial template from HALF of non-task frames (CTRL_A).
%   - Compare template to:
%       (1) other half of non-task frames (CTRL_B)  -> WITHOUT task
%       (2) TASK frames in TraceWin (time-binned)   -> WITH task
%   - All averaging/testing done in Fisher-z (atanh).
%
% Required per rat (in base workspace):
%   rat.Ca_peaks, rat.Ca_ts, rat.pos, rat.CS_times
%
% Options (name/value):
%   'DaysMode'        : {'last3toAn'|'all'|'last3byca'}  default 'last3toAn'
%   'TraceWin'        : [t0 t1] seconds from CS          default [0 2]
%   'BufferPost'      : seconds excluded after TraceWin  default 0
%   'CellNorm'        : {'demean'|'zscore'|'none'}       default 'none'
%   'GridRC'          : [rows cols] for spatial grid     default [2 2]
%   'NumBins'         : integer, overrides GridRC        default []
%   'MinCtrlFrames'   : min CTRL frames per spatial bin *per half* default 10
%   'MinTraceFrames'  : min TASK frames per (timebin,spatialbin)    default 2
%   'TimeBins'        : number of temporal bins in TraceWin          default 15
%   'FramesPerBin'    : fixed frames per temporal bin override       default []
%   'CtrlSplitMode'   : {'random'|'evenodd'}             default 'random'
%   'CtrlSplitSeed'   : [] or scalar                     default []
%   'VelThresh'       : cm/s speed threshold             default 4
%   'UseSpeedMask'    : bool                             default true
%
% Permutations:
%   'NPerm_C'         : group null perms                 default 500
%   'NullMode_C'      : (reserved; currently uses flip within-day pairs)
%   'DeltaPermType'   : {'derange'|'flip'|'both'}        default 'flip'
%
% Single-cell (corr across spatial bins per cell):
%   'DoSingleCell'    : bool                             default true
%   'NPerm_SC'        : perms                            default 500
%   'NullMode_SC'     : {'flip-cells'|'template-derange-bins'} default 'flip-cells'
%
% Output:
%   R(i).C.* and R(i).meta.*

% ------------------------------- Options ---------------------------------
p = inputParser;
addParameter(p,'DaysMode','last3toAn');
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'BufferPost',0);
addParameter(p,'GridRC',[2 2]);
addParameter(p,'NumBins',[]);

addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',true);
addParameter(p,'CellNorm','none');  % 'zscore'|'demean'|'none'

addParameter(p,'MinCtrlFrames',10);
addParameter(p,'MinTraceFrames',2);
addParameter(p,'TimeBins',15);
addParameter(p,'FramesPerBin',[]);

addParameter(p,'CtrlSplitMode','random', @(s) any(strcmpi(s,{'random','evenodd'})));
addParameter(p,'CtrlSplitSeed',[], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)));

addParameter(p,'NPerm_C',500);
addParameter(p,'NullMode_C','frame-redistribute'); %#ok<NASGU>
addParameter(p,'DeltaPermType','flip', @(s) any(strcmpi(s,{'derange','flip','both'})));

addParameter(p,'DoSingleCell',true, @(x) islogical(x)&&isscalar(x));
addParameter(p,'NPerm_SC',500);
addParameter(p,'NullMode_SC','flip-cells', @(s) any(strcmpi(s,{'flip-cells','template-derange-bins'})));

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

    % ----- Select days -----
    dateList = autoDateList_fallback(rat);
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
                    daysToUse = dateList(max(1,numel(dateList)-2):end);
                end
            else
                daysToUse = dateList(max(1,numel(dateList)-2):end);
            end
        otherwise
            daysToUse = dateList;
    end
    fprintf('[%s] daysToUse: %s\n', ratVar, strjoin(daysToUse, ', '));

    % ----- Pull & standardize per-day containers -----
    spikes_raw = filterFieldsByDay_fallback(rat.Ca_peaks, daysToUse);
    ts_raw     = filterFieldsByDay_fallback(rat.Ca_ts,    daysToUse);
    pos_raw    = filterFieldsByDay_fallback(rat.pos,      daysToUse);
    cs_raw     = filterFieldsByDay_fallback(rat.CS_times, daysToUse);
    [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, daysToUse);

    % ----- CONTROL spatial PV template from CTRL_A; WITHOUT = CTRL_B vs template -----
    ctrl = compute_control_spatialPV_template(spikes, ts, pos, cs, ...
        'TraceWin',        opt.TraceWin, ...
        'BufferPost',      opt.BufferPost, ...
        'GridRC',          opt.GridRC, ...
        'NumBins',         opt.NumBins, ...
        'MinCtrlFrames',   opt.MinCtrlFrames, ...
        'VelThresh',       opt.VelThresh, ...
        'UseSpeedMask',    opt.UseSpeedMask, ...
        'CellNorm',        opt.CellNorm, ...
        'CtrlSplitMode',   opt.CtrlSplitMode, ...
        'CtrlSplitSeed',   opt.CtrlSplitSeed, ...
        'Label',           ratVar);

    % ----- WITH task: TASK vs template -----
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
        'NPerm_SC',        opt.NPerm_SC, ...
        'NullMode_SC',     opt.NullMode_SC, ...
        'Label',           ratVar, ...
        'DeltaPermType',   opt.DeltaPermType);

    % Per-day paired z stats (quick)
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
if opt.DoPlots
    summarize_taskSpacePV_clean(R);
    plot_taskToSpace_deltaPermHist_perRat(R, 'PermType', opt.DeltaPermType);

    if opt.DoSingleCell
        summarize_taskSpacePV_singleCell_clean(R);
        plot_taskToSpace_SC_deltaPermHist_perRat(R);
    end
end

for i = 1:numel(R)
    fprintf('\n=== %s ===\n', R(i).animal);

    hasSC = isfield(R(i),'C') && isfield(R(i).C,'SC') && ~isempty(R(i).C.SC);
    fprintf('has C.SC: %d\n', hasSC);

    if hasSC && isfield(R(i).C.SC,'byDay') ...
            && isfield(R(i).C.SC.byDay,'withTask') && isfield(R(i).C.SC.byDay.withTask,'zmean')
        disp('withTask zmean:'); disp(R(i).C.SC.byDay.withTask.zmean(:)');
    else
        disp('withTask zmean: MISSING');
    end

    if hasSC && isfield(R(i).C.SC,'byDay') ...
            && isfield(R(i).C.SC.byDay,'withoutTask') && isfield(R(i).C.SC.byDay.withoutTask,'zmean')
        disp('withoutTask zmean:'); disp(R(i).C.SC.byDay.withoutTask.zmean(:)');
    else
        disp('withoutTask zmean: MISSING');
    end

    if hasSC && isfield(R(i).C.SC,'delta') && isfield(R(i).C.SC.delta,'perm_z') && ~isempty(R(i).C.SC.delta.perm_z)
        fprintf('SC.delta.perm_z length: %d\n', numel(R(i).C.SC.delta.perm_z));
    else
        disp('SC.delta.perm_z: MISSING/EMPTY');
    end
end

end

% =========================================================================
% ========================== HELPER FUNCTIONS =============================
% =========================================================================

% ---------- CTRL spatial PV TEMPLATE from half of non-task frames ----------
function CTRL = compute_control_spatialPV_template(spikes, ts, pos, cs, varargin)
p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'BufferPost',0);
addParameter(p,'GridRC',[3 2]);
addParameter(p,'NumBins',[]);
addParameter(p,'MinCtrlFrames',10);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'Label','');
addParameter(p,'CellNorm','demean');
addParameter(p,'CtrlSplitMode','random');
addParameter(p,'CtrlSplitSeed',[]);
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);

% ---- Build grid edges ----
if ~isempty(o.NumBins)
    rc_eff = best_factors(o.NumBins);
else
    rc_eff = o.GridRC;
end
K = rc_eff(1)*rc_eff(2);

edgesByDay = cell(1,D);
for d = 1:D
    edgesByDay{d} = build_grid_edges_single_day(pos{d}, rc_eff);
end

PV_template_by = cell(1,D);    % template from CTRL_A
PV_ctrlB_by    = cell(1,D);    % map from CTRL_B
rel_by         = cell(1,D);    % corr(template_k, ctrlB_k) per bin

if ~isempty(o.CtrlSplitSeed)
    rng(o.CtrlSplitSeed);
end

for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);

    % mask out Task + BufferPost => CONTROL frames
    is_taskbuf = false(size(t));
    csd = cs{d};
    if ~isempty(csd)
        for j = 1:numel(csd)
            is_taskbuf = is_taskbuf | (t >= csd(j)+o.TraceWin(1) & t <= csd(j)+o.TraceWin(2)+o.BufferPost);
        end
    end
    use = ~is_taskbuf;

    if o.UseSpeedMask
        v = speed_cm_per_s(pos{d});
        v_i = interp1(pos{d}.t(:), v(:), t, 'linear','extrap');
        use = use & (v_i >= o.VelThresh);
    end

    S = spikes_to_matrix(spikes{d}, t);
    S = normalize_cells(S, o.CellNorm);

    S   = S(:, use);
    x_u = x(use);
    y_u = y(use);

    [~, b] = pos2bin(x_u, y_u, edgesByDay{d});

    Nc_d = size(S,1);
    PVtempl = nan(Nc_d, K);
    PVb     = nan(Nc_d, K);
    rel     = nan(K,1);

    for k = 1:K
        idx = find(b == k);
        if numel(idx) < 2*o.MinCtrlFrames
            continue;
        end

        [idxA, idxB] = split_indices(idx, o.CtrlSplitMode);

        if numel(idxA) < o.MinCtrlFrames || numel(idxB) < o.MinCtrlFrames
            continue;
        end

        vA = mean(S(:, idxA), 2, 'omitnan');
        vB = mean(S(:, idxB), 2, 'omitnan');

        PVtempl(:,k) = vA;
        PVb(:,k)     = vB;
        rel(k)       = safe_corr(vA, vB);
    end

    PV_template_by{d} = PVtempl;
    PV_ctrlB_by{d}    = PVb;
    rel_by{d}         = rel;

    kept = all(isfinite(PVtempl),1) & all(isfinite(PVb),1);
    fprintf('[%s][CTRL-TEMPLATE d=%d] kept bins %d/%d | template↔ctrlB median r=%.2f\n', ...
        L, d, nnz(kept), K, median(rel,'omitnan'));
end

CTRL.PV_template = PV_template_by;
CTRL.PV_ctrlB    = PV_ctrlB_by;

CTRL.grid.edges  = edgesByDay;
CTRL.grid.GridRC = rc_eff;
CTRL.grid.K      = K;

CTRL.reliability_byDay = rel_by;
CTRL.reliability.r     = vertcat(rel_by{:});
CTRL.reliability.z     = atanh(max(min(CTRL.reliability.r,0.99999),-0.99999));

CTRL.params.MinCtrlFrames   = o.MinCtrlFrames;
CTRL.params.NumBins         = o.NumBins;
CTRL.params.GridRC_eff      = rc_eff;
CTRL.params.CtrlSplitMode   = o.CtrlSplitMode;
CTRL.params.CtrlSplitSeed   = o.CtrlSplitSeed;

CTRL.maps.template_byDay = PV_template_by;
CTRL.maps.ctrlB_byDay    = PV_ctrlB_by;
end

function [idxA, idxB] = split_indices(idx, mode)
idx = idx(:);
switch lower(mode)
    case 'evenodd'
        idxA = idx(1:2:end);
        idxB = idx(2:2:end);
    otherwise
        rp = idx(randperm(numel(idx)));
        nA = floor(numel(rp)/2);
        idxA = rp(1:nA);
        idxB = rp(nA+1:end);
end
end

% ---------- WITH vs WITHOUT task (C) ----------
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
addParameter(p,'NullMode','frame-redistribute'); %#ok<NASGU>
addParameter(p,'Label','');
addParameter(p,'NPerm_SC',500);
addParameter(p,'NullMode_SC','flip-cells', @(s) any(strcmpi(s,{'flip-cells','template-derange-bins'})));
addParameter(p,'DeltaPermType','flip', @(s) any(strcmpi(s,{'derange','flip','both'})));
parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);
K = CTRL.grid.K;

r_bin_byDay     = cell(1,D);
z_day           = nan(D,1);
n_valid_day     = zeros(D,1);
perZ_list_byDay = cell(1,D);

C = struct();

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

    edges_d  = CTRL.grid.edges{d};
    PV_templ = CTRL.PV_template{d};
    if isempty(PV_templ), continue; end

    csd = cs{d}(:);
    B = o.TimeBins;
    r_bin_accum = nan(K, B);

    for tr = 1:numel(csd)
        t0 = csd(tr) + o.TraceWin(1);
        t1 = csd(tr) + o.TraceWin(2);
        idx_all = find(t >= t0 & t < t1);
        if isempty(idx_all), continue; end

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
        nTB = numel(edges_idx)-1;

        for bti = 1:nTB
            if edges_idx(bti) >= edges_idx(bti+1), continue; end
            seg = idx_all(edges_idx(bti):edges_idx(bti+1)-1);
            if o.UseSpeedMask, seg = seg(speed_ok(seg)); end
            if numel(seg) < o.MinTraceFrames, continue; end

            [~, b_space] = pos2bin(x(seg), y(seg), edges_d);

            for k = 1:K
                hit = (b_space == k);
                if nnz(hit) >= o.MinTraceFrames && all(isfinite(PV_templ(:,k)))
                    pv_task = mean(S(:, seg(hit)), 2, 'omitnan');
                    r = safe_corr(PV_templ(:,k), pv_task);
                    r_bin_accum(k, bti) = nanmean([r_bin_accum(k, bti), r]);
                end
            end
        end
    end

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

    r_spatial = nan(K,1);
    for k = 1:K
        rv = r_bin_accum(k,:);
        rv = rv(isfinite(rv));
        if ~isempty(rv)
            r_spatial(k) = tanh(mean(atanh(max(min(rv,0.99999),-0.99999))));
        end
    end
    r_bin_byDay{d} = r_spatial;
    n_valid_day(d) = numel(z_list);

    % ---------------- SINGLE-CELL ----------------
    if isfield(CTRL,'maps') && numel(CTRL.maps.template_byDay) >= d && ~isempty(CTRL.maps.template_byDay{d}) && ...
       numel(CTRL.maps.ctrlB_byDay)    >= d && ~isempty(CTRL.maps.ctrlB_byDay{d})

        mapA = CTRL.maps.template_byDay{d};  % Nc×K
        mapB = CTRL.maps.ctrlB_byDay{d};     % Nc×K
        Nc = size(mapA,1);

        % WITHOUT: per-cell corr across bins
        z_sc_without = nan(Nc,1);
        for c = 1:Nc
            r = safe_corr(mapA(c,:).', mapB(c,:).', 3, 1e-12);
            if isfinite(r), z_sc_without(c) = atanh(max(min(r,0.99999),-0.99999)); end
        end

        % --- WITH task (FIXED): pool across TRIALS within each timebin to fill bins ---
        % Build pooled TASK maps per timebin: taskMap_tb(:,:,bti) is Nc×K
        B  = o.TimeBins;
        Nc = size(mapA,1);

        % Pool frame indices by (timebin, spatialbin) across ALL trials
        idxPool = cell(B, K);   % each entry holds frame indices into t/S

        csd = cs{d}(:);
        for tr = 1:numel(csd)
            t0 = csd(tr) + o.TraceWin(1);
            t1 = csd(tr) + o.TraceWin(2);
            idx_all = find(t >= t0 & t < t1);
            if isempty(idx_all), continue; end

            % same temporal split logic
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
            nTB = numel(edges_idx)-1;

            for bti = 1:nTB
                seg = idx_all(edges_idx(bti):edges_idx(bti+1)-1);
                if isempty(seg), continue; end

                if o.UseSpeedMask
                    seg = seg(speed_ok(seg));
                end
                if numel(seg) < o.MinTraceFrames
                    continue
                end

                [~, b_space] = pos2bin(x(seg), y(seg), edges_d);

                % add frames into pools by spatial bin
                for k = 1:K
                    hit = (b_space == k);
                    if nnz(hit) >= o.MinTraceFrames && all(isfinite(mapA(:,k)))
                        idxPool{bti,k} = [idxPool{bti,k}; seg(hit)]; %#ok<AGROW>
                    end
                end
            end
        end

        % Convert pooled indices -> pooled task maps
        taskMap_tb = nan(Nc, K, B);  % Nc×K×B
        for bti = 1:B
            for k = 1:K
                idxk = idxPool{bti,k};
                if numel(idxk) >= o.MinTraceFrames
                    taskMap_tb(:,k,bti) = mean(S(:, idxk), 2, 'omitnan');
                end
            end
        end

        % Per-cell corr across spatial bins per timebin vs template
        z_sc_with_tb = nan(Nc, B);
        for bti = 1:B
            for c = 1:Nc
                r = safe_corr(mapA(c,:).', taskMap_tb(c,:,bti).', 2, 1e-12); % minN=2 (not 3)
                if isfinite(r)
                    z_sc_with_tb(c,bti) = atanh(max(min(r,0.99999),-0.99999));
                end
            end
        end


        % collapse across timebins in Fisher-z per cell
        z_sc_with = mean(z_sc_with_tb, 2, 'omitnan');
        r_sc_with = tanh(z_sc_with);

        if ~isfield(C,'SC'), C.SC = struct(); end
        if ~isfield(C.SC,'byDay'), C.SC.byDay = struct(); end
        if ~isfield(C.SC.byDay,'withTask'), C.SC.byDay.withTask = struct(); end
        if ~isfield(C.SC.byDay,'withoutTask'), C.SC.byDay.withoutTask = struct(); end

        C.SC.byDay.withTask.z{d}     = z_sc_with;
        C.SC.byDay.withoutTask.z{d}  = z_sc_without;

        C.SC.byDay.withTask.zmean(d,1)    = mean(z_sc_with, 'omitnan');
        C.SC.byDay.withoutTask.zmean(d,1) = mean(z_sc_without, 'omitnan');

        % Optional: template-derange-bins null at SC level (pooled-task map)
        if ~isempty(o.NPerm_SC) && o.NPerm_SC >= 10
            validK = find(all(isfinite(mapA),1));
            if numel(validK) >= 2
                [taskMap_pooled, okTask] = build_taskMap_pooled_day(t, x, y, S, csd, edges_d, speed_ok, o);
                if okTask
                    delta_obs = mean(z_sc_with, 'omitnan') - mean(z_sc_without, 'omitnan');
                    perm_delta = nan(o.NPerm_SC,1);
                    for pidx = 1:o.NPerm_SC
                        kp = derange_indices(validK);
                        z_with_perm = nan(Nc,1);
                        for c = 1:Nc
                            r = safe_corr(mapA(c, kp).', taskMap_pooled(c, validK).', 3, 1e-12);
                            if isfinite(r), z_with_perm(c) = atanh(max(min(r,0.99999),-0.99999)); end
                        end
                        perm_delta(pidx) = mean(z_with_perm, 'omitnan') - mean(z_sc_without, 'omitnan');
                    end
                    pr = mean(perm_delta >= delta_obs);
                    pl = mean(perm_delta <= delta_obs);
                    pt = 2*min(pr, pl);

                    C.SC.delta_derange.byDay{d}   = perm_delta;
                    C.SC.delta_derange.obs_z(d,1) = delta_obs;
                    C.SC.delta_derange.p_right(d,1) = pr;
                    C.SC.delta_derange.p_left(d,1)  = pl;
                    C.SC.delta_derange.p_two(d,1)   = pt;
                end
            end
        end
    end
end

% After all days: build SC flip-cells null (paired within cell)
if isfield(C,'SC') && ~isempty(o.NPerm_SC) && o.NPerm_SC >= 10
    C.SC.delta_flipCells = build_flip_perm_singleCell_flipCells(C, o.NPerm_SC);
    fprintf('[%s][SC flip-cells] obs Δz=%.3f (r≈%.3f), p_two=%.3g (NPerm=%d)\n', ...
        L, C.SC.delta_flipCells.obs_z, tanh(C.SC.delta_flipCells.obs_z), ...
        C.SC.delta_flipCells.p_two, numel(C.SC.delta_flipCells.perm_z));
end

% Pick active SC delta
if isfield(C,'SC')
    C.SC.delta = [];
    switch lower(char(o.NullMode_SC))
        case 'flip-cells'
            if isfield(C.SC,'delta_flipCells') && ~isempty(C.SC.delta_flipCells) && ...
               isfield(C.SC.delta_flipCells,'perm_z') && ~isempty(C.SC.delta_flipCells.perm_z)
                C.SC.delta = C.SC.delta_flipCells;
            end
        case 'template-derange-bins'
            % Optional pooling of per-day derange; if absent, fall back
            if isfield(C.SC,'delta_derange') && ~isempty(C.SC.delta_derange)
                % If you later add pooling here, set C.SC.delta_derange.perm_z etc.
            end
    end
    if isempty(C.SC.delta) && isfield(C.SC,'delta_flipCells')
        C.SC.delta = C.SC.delta_flipCells;
    end
end

% Bin-level summaries
C.bin.byDay.r        = r_bin_byDay;
C.bin.byDay.z_mean   = z_day;
C.bin.byDay.n_valid  = n_valid_day;
C.bin.byDay.z_list   = perZ_list_byDay;
C.bin.group.z_mean   = mean(z_day(isfinite(z_day)),'omitnan');
C.bin.group.r_mean   = tanh(C.bin.group.z_mean);

% WITHOUT = template vs CTRL_B
with_byDay    = cell(1,D);
without_byDay = cell(1,D);
z_without_day = nan(D,1);

for d = 1:D
    r_task_bins = C.bin.byDay.r{d};
    if isempty(r_task_bins), continue; end

    if isfield(CTRL,'reliability_byDay') && numel(CTRL.reliability_byDay) >= d && ~isempty(CTRL.reliability_byDay{d})
        r_without = CTRL.reliability_byDay{d}(:);
    else
        r_without = nan(K,1);
    end

    mask = isfinite(r_task_bins) & isfinite(r_without);
    if ~any(mask), continue; end

    with_byDay{d}    = r_task_bins(mask);
    without_byDay{d} = r_without(mask);

    z_without_day(d) = mean(atanh(max(min(without_byDay{d},0.99999),-0.99999)), 'omitnan');
end

C.byDay.withTask.r    = with_byDay;
C.byDay.withoutTask.r = without_byDay;

delta_day = z_day - z_without_day; %#ok<NASGU>

% Flip null (within-day paired pairs)
C.delta_flip = build_flip_perm_withinDayPairs(C, o.NPerm);
if isfield(C,'delta_flip') && isfield(C.delta_flip,'p_two')
    fprintf('[%s] flip null: obs Δz=%.3f (r≈%.3f), p_two=%.3g (NPerm=%d)\n', ...
        L, C.delta_flip.obs_z, tanh(C.delta_flip.obs_z), C.delta_flip.p_two, numel(C.delta_flip.perm_z));
end

% Derange null (template bins)
C.delta_derange = build_derange_perm_templateBins(C, spikes, ts, pos, cs, CTRL, z_without_day, o);

% Pick active delta
which = lower(char(o.DeltaPermType));
if strcmp(which,'both')
    which = 'flip';
end
switch which
    case 'flip'
        C.delta = C.delta_flip;
    case 'derange'
        if ~isempty(C.delta_derange) && isfield(C.delta_derange,'perm_z') && ~isempty(C.delta_derange.perm_z)
            C.delta = C.delta_derange;
        else
            C.delta = C.delta_flip;
        end
    otherwise
        C.delta = C.delta_flip;
end
end

% ---------- Derange perm (template bins vs TASK bins) ----------
function dlt = build_derange_perm_templateBins(~, spikes, ts, pos, cs, CTRL, z_without_day, o)
dlt = struct('obs_z',NaN,'obs_r',NaN,'perm_z',[],'p_right',NaN,'p_left',NaN,'p_two',NaN);
P = o.NPerm;
if isempty(P) || P < 5
    return
end

D = numel(ts);
K = CTRL.grid.K;
perm_byDay = cell(D,1);
z_with_day = nan(D,1);

for d = 1:D
    if ~isfinite(z_without_day(d)), continue; end

    PVtempl = CTRL.PV_template{d};
    if isempty(PVtempl), continue; end
    edges_d = CTRL.grid.edges{d};

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

    csd = cs{d}(:);
    B = o.TimeBins;

    TASK_kb = cell(K,B);

    for tr = 1:numel(csd)
        t0 = csd(tr) + o.TraceWin(1);
        t1 = csd(tr) + o.TraceWin(2);
        idx_all = find(t >= t0 & t < t1);
        if isempty(idx_all), continue; end

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

        for bti = 1:(numel(edges_idx)-1)
            seg = idx_all(edges_idx(bti):edges_idx(bti+1)-1);
            if o.UseSpeedMask, seg = seg(speed_ok(seg)); end
            if numel(seg) < o.MinTraceFrames, continue; end

            [~, b_space] = pos2bin(x(seg), y(seg), edges_d);

            for k = 1:K
                hit = (b_space == k);
                if nnz(hit) >= o.MinTraceFrames && all(isfinite(PVtempl(:,k)))
                    vec = mean(S(:, seg(hit)), 2, 'omitnan');
                    if isempty(TASK_kb{k,bti})
                        TASK_kb{k,bti} = vec;
                    else
                        TASK_kb{k,bti} = nanmean([TASK_kb{k,bti}, vec], 2);
                    end
                end
            end
        end
    end

    validK = find(all(isfinite(PVtempl),1));
    if numel(validK) < 2, continue; end

    % Observed z_with_day from TASK_kb
    r_list_obs = [];
    for bti = 1:B
        for ii = 1:numel(validK)
            k = validK(ii);
            tv = TASK_kb{k,bti};
            if ~isempty(tv)
                r_list_obs(end+1,1) = safe_corr(PVtempl(:,k), tv); %#ok<AGROW>
            end
        end
    end
    if isempty(r_list_obs), continue; end
    z_with_day(d) = mean(atanh(max(min(r_list_obs,0.99999),-0.99999)), 'omitnan');

    perm_d = nan(P,1);
    for pidx = 1:P
        map = derange_indices(validK); % deranged mapping of valid bins
        r_list = [];
        for bti = 1:B
            for ii = 1:numel(validK)
                k  = validK(ii);
                kp = map(ii);
                tv = TASK_kb{k,bti};
                if ~isempty(tv)
                    r_list(end+1,1) = safe_corr(PVtempl(:,kp), tv); %#ok<AGROW>
                end
            end
        end
        if isempty(r_list), continue; end
        z_with_perm = mean(atanh(max(min(r_list,0.99999),-0.99999)), 'omitnan');
        perm_d(pidx) = z_with_perm - z_without_day(d);
    end
    perm_byDay{d} = perm_d;
end

mask = isfinite(z_with_day) & isfinite(z_without_day);
if nnz(mask) < 1
    return
end
delta_obs = mean(z_with_day(mask) - z_without_day(mask), 'omitnan');

valid = ~cellfun('isempty', perm_byDay);
if ~any(valid)
    return
end
Pcommon = min(cellfun(@(v) nnz(isfinite(v)), perm_byDay(valid)));
if Pcommon < 5
    return
end

delta_perm_anim = nan(Pcommon,1);
for pp = 1:Pcommon
    dvals = [];
    for d = 1:D
        if ~valid(d), continue; end
        if numel(perm_byDay{d}) >= pp && isfinite(perm_byDay{d}(pp))
            dvals(end+1,1) = perm_byDay{d}(pp); %#ok<AGROW>
        end
    end
    delta_perm_anim(pp) = mean(dvals,'omitnan');
end

pr = mean(delta_perm_anim >= delta_obs);
pl = mean(delta_perm_anim <= delta_obs);
pt = 2*min(pr, pl);

dlt.obs_z   = delta_obs;
dlt.obs_r   = tanh(delta_obs);
dlt.perm_z  = delta_perm_anim;
dlt.p_right = pr;
dlt.p_left  = pl;
dlt.p_two   = pt;
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
        [~, pval, ~, st] = ttest(zw, zu); dz = mean(delta)/std(delta);
        dayStats(d) = struct('nPairs',numel(delta),'t',st.tstat,'df',st.df,'p',pval,'dz',dz);
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

if ~isfield(C,'distStats'), C.distStats = struct(); end
C.distStats.perDay  = dayStats;
C.distStats.pooled  = struct('t',st_all.tstat,'df',st_all.df,'p',p_all,'dz',dz_all, 'nBins', nnz(mask));
C.distStats.cluster_perm = struct('NPerm',NPerm,'obs_meanDelta_z',obs_cluster_mean, ...
                                  'p_right',p_right,'p_left',p_left,'p_two',p_two);
end

% ========================== UTILITIES ==========================
function [taskMap, ok] = build_taskMap_pooled_day(t, x, y, S, csd, edges_d, speed_ok, o)
% Build pooled task map across ALL trace frames (ignores time bins).
% Output: taskMap [Nc×K] mean activity in each spatial bin within TraceWin.
K  = (numel(edges_d.x)-1) * (numel(edges_d.y)-1);
Nc = size(S,1);
taskMap = nan(Nc, K);
ok = false;
if isempty(csd), return; end

sumV = zeros(Nc, K);
cnt  = zeros(1, K);

for tr = 1:numel(csd)
    t0 = csd(tr) + o.TraceWin(1);
    t1 = csd(tr) + o.TraceWin(2);
    seg = find(t >= t0 & t < t1);
    if isempty(seg), continue; end
    if o.UseSpeedMask, seg = seg(speed_ok(seg)); end
    if numel(seg) < o.MinTraceFrames, continue; end

    [~, b_space] = pos2bin(x(seg), y(seg), edges_d);

    for k = 1:K
        hit = (b_space == k);
        if nnz(hit) >= o.MinTraceFrames
            v = mean(S(:, seg(hit)), 2, 'omitnan');
            if all(isfinite(v))
                sumV(:,k) = sumV(:,k) + v;
                cnt(k)    = cnt(k) + 1;
            end
        end
    end
end

for k = 1:K
    if cnt(k) > 0
        taskMap(:,k) = sumV(:,k) ./ cnt(k);
    end
end
ok = any(all(isfinite(taskMap),1));
end

function kp = derange_indices(validK)
% Return a derangement of validK (no fixed points), length = numel(validK).
validK = validK(:)';
n = numel(validK);
if n < 2
    kp = validK;
    return
end
for rep = 1:50
    perm = validK(randperm(n));
    if all(perm ~= validK)
        kp = perm;
        return
    end
end
% Fallback: cyclic shift
kp = validK([2:end,1]);
end

function OUT = build_flip_perm_withinDayPairs(C, NPerm)
OUT = struct('obs_z',NaN,'obs_r',NaN,'perm_z',[],'p_right',NaN,'p_left',NaN,'p_two',NaN, ...
             'nPairsByDay',[],'nDays',0);

if ~isfield(C,'byDay') || ~isfield(C.byDay,'withTask') || ~isfield(C.byDay,'withoutTask')
    return
end
if ~isfield(C.byDay.withTask,'r') || ~isfield(C.byDay.withoutTask,'r')
    return
end

rw = C.byDay.withTask.r;
ru = C.byDay.withoutTask.r;

D = min(numel(rw), numel(ru));
day_delta_mu  = nan(D,1);
day_delta_vec = cell(D,1);
nPairs        = zeros(D,1);

for d = 1:D
    if isempty(rw{d}) || isempty(ru{d}), continue; end
    a = rw{d}(:); b = ru{d}(:);
    m = isfinite(a) & isfinite(b);
    if nnz(m) < 2, continue; end

    zw = atanh(max(min(a(m),0.999999),-0.999999));
    zu = atanh(max(min(b(m),0.999999),-0.999999));
    dv = zw - zu;

    day_delta_vec{d} = dv;
    day_delta_mu(d)  = mean(dv,'omitnan');
    nPairs(d)        = numel(dv);
end

validDays = find(isfinite(day_delta_mu) & ~cellfun(@isempty, day_delta_vec));
OUT.nDays = numel(validDays);
OUT.nPairsByDay = nPairs;

if OUT.nDays < 1
    return
end

obs = mean(day_delta_mu(validDays), 'omitnan');

if nargin < 2 || isempty(NPerm) || NPerm < 5
    OUT.obs_z = obs;
    OUT.obs_r = tanh(obs);
    return
end

perm = nan(NPerm,1);
for pidx = 1:NPerm
    day_mu_perm = nan(OUT.nDays,1);
    for j = 1:OUT.nDays
        d = validDays(j);
        dv = day_delta_vec{d};
        s = (rand(size(dv)) > 0.5)*2 - 1;
        day_mu_perm(j) = mean(s .* dv, 'omitnan');
    end
    perm(pidx) = mean(day_mu_perm, 'omitnan');
end

pr = mean(perm >= obs);
pl = mean(perm <= obs);
pt = 2*min(pr, pl);

OUT.obs_z   = obs;
OUT.obs_r   = tanh(obs);
OUT.perm_z  = perm;
OUT.p_right = pr;
OUT.p_left  = pl;
OUT.p_two   = pt;
end

function OUT = build_flip_perm_singleCell_flipCells(C, NPerm)
OUT = struct('obs_z',NaN,'obs_r',NaN,'perm_z',[],'p_right',NaN,'p_left',NaN,'p_two',NaN, ...
             'nCellsByDay',[],'nDays',0);

if nargin < 2 || isempty(NPerm), NPerm = 500; end
if ~isfield(C,'SC') || ~isfield(C.SC,'byDay') || ...
   ~isfield(C.SC.byDay,'withTask') || ~isfield(C.SC.byDay,'withoutTask') || ...
   ~isfield(C.SC.byDay.withTask,'z') || ~isfield(C.SC.byDay.withoutTask,'z')
    return
end

zw = C.SC.byDay.withTask.z;
zu = C.SC.byDay.withoutTask.z;

D = min(numel(zw), numel(zu));
day_mu   = nan(D,1);
day_dv   = cell(D,1);
nCells   = zeros(D,1);

for d = 1:D
    if isempty(zw{d}) || isempty(zu{d}), continue; end
    a = zw{d}(:); b = zu{d}(:);
    m = isfinite(a) & isfinite(b);
    if nnz(m) < 2, continue; end

    dv = a(m) - b(m);
    day_dv{d} = dv;
    day_mu(d) = mean(dv,'omitnan');
    nCells(d) = numel(dv);
end

validDays = find(isfinite(day_mu) & ~cellfun(@isempty, day_dv));
OUT.nDays = numel(validDays);
OUT.nCellsByDay = nCells;

if OUT.nDays < 1
    return
end

obs = mean(day_mu(validDays), 'omitnan');

perm = nan(NPerm,1);
for pidx = 1:NPerm
    mu_perm = nan(OUT.nDays,1);
    for j = 1:OUT.nDays
        d  = validDays(j);
        dv = day_dv{d};
        s  = (rand(size(dv)) > 0.5)*2 - 1;
        mu_perm(j) = mean(s .* dv, 'omitnan');
    end
    perm(pidx) = mean(mu_perm, 'omitnan');
end

OUT.obs_z   = obs;
OUT.obs_r   = tanh(obs);
OUT.perm_z  = perm;
OUT.p_right = mean(perm >= obs);
OUT.p_left  = mean(perm <= obs);
OUT.p_two   = 2*min(OUT.p_right, OUT.p_left);
end

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

function cs_vec = coerce_cs_day(csd, ~)
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

% ======================= FALLBACK DAY HELPERS =======================
function dateList = autoDateList_fallback(rat)
% Prefer your existing helper if present on path.
if exist('autoDateList','file') == 2
    try
        dateList = autoDateList(rat);
        return
    catch
    end
end
% Fallback: infer from Ca_peaks fields or Ca_ts fields
dateList = {};
if isfield(rat,'Ca_peaks') && isstruct(rat.Ca_peaks)
    dateList = fieldnames(rat.Ca_peaks);
elseif isfield(rat,'Ca_ts') && isstruct(rat.Ca_ts)
    dateList = fieldnames(rat.Ca_ts);
end
dateList = dateList(:);
end

function S = filterFieldsByDay_fallback(Sin, daysToUse)
% Prefer your existing helper if present.
if exist('filterFieldsByDay','file') == 2
    try
        S = filterFieldsByDay(Sin, daysToUse);
        return
    catch
    end
end
% Fallback: pass through; standardize_day_format will pick matching fields.
S = Sin;
end

% ===================== PLOTS / SUMMARIES (same as your version) =====================

function summarize_taskSpacePV_clean(R)
titleStr = 'C: WITH vs WITHOUT task';
[Zw_all, Za_all, aid, did] = get_C_dayZ(R);

if isempty(Zw_all) || isempty(Za_all)
    warning('No WITH/WITHOUT day-level data found to summarize.');
    return
end

animals = {R.animal}';
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

[pA, tA, dfA, dzA] = paired_t_z(Zw_rats, Za_rats);

keepD = isfinite(Zw_all) & isfinite(Za_all);
Zw_days = Zw_all(keepD);
Za_days = Za_all(keepD);
[pD, tD, dfD, dzD] = paired_t_z(Zw_days, Za_days);

figure('Color','w','Position',[160 160 900 540]); hold on

cmap = lines(nA);
make_light = @(c,frac) (1-frac)*c + frac*[1 1 1];

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

rat_with_r    = tanh(rat_mean_with_z);
rat_without_r = tanh(rat_mean_without_z);

for i = 1:nA
    if isfinite(rat_with_r(i)) && isfinite(rat_without_r(i))
        plot([1 2], [rat_with_r(i) rat_without_r(i)], '-', 'Color', cmap(i,:), 'LineWidth', 2.5);
        plot(1, rat_with_r(i),    'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
        plot(2, rat_without_r(i), 'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
    end
end

bar(1, tanh(mean(Zw_rats,'omitnan')), 0.6, 'FaceColor',[0.35 0.70 1.00], 'EdgeColor','k');
bar(2, tanh(mean(Za_rats,'omitnan')), 0.6, 'FaceColor',[0.85 0.55 0.25], 'EdgeColor','k');

xlim([0.5 2.5]); xticks([1 2]);
xticklabels({'with task','without task'});
ylabel('Mean PV correlation (r)');
title(titleStr);
yline(0,'k:'); grid on; box on

yl = ylim; y = yl(2) + 0.06*range(yl);
line([1 2], [y y], 'Color','k', 'LineWidth', 1.2);
txt1 = sprintf('per-rat paired t (z): t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d', dfA, tA, pA, dzA, numel(Zw_rats));
text(1.5, y + 0.02*range(yl), txt1, 'HorizontalAlignment','center');
ylim([yl(1) y + 0.10*range(yl)]);

fprintf('\n== %s ==\n', titleStr);
fprintf('Primary (per-rat):  paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d\n', dfA, tA, pA, dzA, numel(Zw_rats));
fprintf('Secondary (per-day): paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nDays=%d\n', dfD, tD, pD, dzD, numel(Zw_days));
end

function [Zw, Za, aid, did] = get_C_dayZ(R)
Zw=[]; Za=[]; aid=[]; did=[];
for i = 1:numel(R)
    has_with = isfield(R(i),'C') && isfield(R(i).C,'bin') && isfield(R(i).C.bin,'byDay') ...
               && isfield(R(i).C.bin.byDay,'z_mean');
    has_without = isfield(R(i),'C') && isfield(R(i).C,'byDay') ...
                  && isfield(R(i).C.byDay,'withoutTask') && isfield(R(i).C.byDay.withoutTask,'r');
    if ~has_with || ~has_without, continue; end

    z_with = R(i).C.bin.byDay.z_mean(:);
    ur     = R(i).C.byDay.withoutTask.r;

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

    m = isfinite(z_with) & isfinite(z_without);
    Zw  = [Zw;  z_with(m)];
    Za  = [Za;  z_without(m)];
    aid = [aid; i*ones(nnz(m),1)];
    did_local = (1:numel(z_with)).';
    did = [did; did_local(m)];
end
end

function [p, t, df, dz] = paired_t_z(z1, z2)
if numel(z1) ~= numel(z2)
    error('paired_t_z: inputs must be same length.');
end
mask = isfinite(z1) & isfinite(z2);
z1 = z1(mask); z2 = z2(mask);

if numel(z1) < 2
    p = NaN; t = NaN; df = 0; dz = NaN; return
end

[~, p, ~, st] = ttest(z1, z2);
t  = st.tstat;
df = st.df;

d  = z1 - z2;
dz = mean(d,'omitnan') / std(d, 0, 'omitnan');
end

function plot_taskToSpace_deltaPermHist_perRat(R, varargin)
p = inputParser;
addParameter(p,'NBins',30, @(x) isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'ShowTwoSided',true, @(x) islogical(x)&&isscalar(x));
addParameter(p,'PermType','derange', @(s) any(strcmpi(s,{'derange','flip'})));
parse(p,varargin{:});
NBins = p.Results.NBins;
ShowTwoSided = p.Results.ShowTwoSided;
PermType = lower(char(p.Results.PermType));

nRats = numel(R);
figure('Color','w','Position',[200 200 1100 700]);

for ii = 1:nRats
    if isfield(R(ii),'animal') && ~isempty(R(ii).animal)
        ratName = R(ii).animal;
    else
        ratName = sprintf('rat%02d', ii);
    end
    if ~isfield(R(ii),'C') || isempty(R(ii).C), continue; end

    dlt = get_delta_struct(R(ii).C, PermType);
    if isempty(dlt) || ~isfield(dlt,'obs_z') || ~isfield(dlt,'perm_z'), continue; end

    obs  = dlt.obs_z;
    perm = dlt.perm_z;

    if ~isfinite(obs) || isempty(perm), continue; end
    perm = perm(isfinite(perm));
    if isempty(perm), continue; end

    subplot(3,2,ii); hold on;
    histogram(perm, NBins, 'Normalization','pdf', 'EdgeColor','none');
    yl = ylim;
    plot([obs obs], yl, 'k-', 'LineWidth', 2);

    xlabel('\Delta z = (WITH - WITHOUT) Fisher z');
    ylabel('Null density');
    title(sprintf('%s: %s null', ratName, PermType));
    grid on; box on

    pR = NaN; pL = NaN; pT = NaN;
    if isfield(dlt,'p_right'), pR = dlt.p_right; end
    if isfield(dlt,'p_left'),  pL = dlt.p_left;  end
    if isfield(dlt,'p_two'),   pT = dlt.p_two;   end

    txt = sprintf('obs \\Delta z = %.3f  (r \\approx %.3f)\n', obs, tanh(obs));
    if ShowTwoSided
        if isfinite(pT)
            txt = [txt, sprintf('p_{two} = %.3g', pT)];
        else
            pR_ = mean(perm >= obs);
            pL_ = mean(perm <= obs);
            txt = [txt, sprintf('p_{two} \\approx %.3g', 2*min(pR_,pL_))];
        end
    else
        if isfinite(pR)
            txt = [txt, sprintf('p_{right} = %.3g', pR)];
        else
            txt = [txt, sprintf('p_{right} \\approx %.3g', mean(perm >= obs))];
        end
    end

    text(obs, yl(2), ['  ', txt], 'VerticalAlignment','top', 'HorizontalAlignment','left');
    ylim(yl);
end
end

function dlt = get_delta_struct(C, permType)
dlt = [];
switch lower(permType)
  case 'flip'
      if isfield(C,'delta_flip') && ~isempty(C.delta_flip)
          dlt = C.delta_flip;
      end
  case 'derange'
      if isfield(C,'delta_derange') && ~isempty(C.delta_derange)
          dlt = C.delta_derange;
      elseif isfield(C,'delta') && ~isempty(C.delta)
          dlt = C.delta;
      end
end
end

function summarize_taskSpacePV_singleCell_clean(R)
titleStr = 'Single-cell: WITH vs WITHOUT task (corr across spatial bins per cell)';

[Zw_all, Za_all, aid, did] = get_SC_dayZ(R);

if isempty(Zw_all) || isempty(Za_all)
    warning('No single-cell WITH/WITHOUT day-level data found.');
    return
end

animals = {R.animal}';
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

[pA, tA, dfA, dzA] = paired_t_z(Zw_rats, Za_rats);

keepD = isfinite(Zw_all) & isfinite(Za_all);
Zw_days = Zw_all(keepD);
Za_days = Za_all(keepD);

[pD, tD, dfD, dzD] = paired_t_z(Zw_days, Za_days); %#ok<ASGLU>

perm = build_flip_perm_singleCell_byDay(Zw_all, Za_all, aid, did, 500);

figure('Color','w','Position',[160 160 920 560]); hold on

cmap = lines(nA);
make_light = @(c,frac) (1-frac)*c + frac*[1 1 1];

rw_day = tanh(Zw_all);
ru_day = tanh(Za_all);

for i = 1:nA
    c_base  = cmap(i,:);
    c_light = make_light(c_base, 0.60);

    m = (aid == i) & isfinite(rw_day) & isfinite(ru_day);
    idx = find(m);
    for j = 1:numel(idx)
        plot([1 2], [rw_day(idx(j)) ru_day(idx(j))], '-', 'Color', c_light, 'LineWidth', 1.0);
        plot(1, rw_day(idx(j)), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
        plot(2, ru_day(idx(j)), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
    end
end

rat_with_r    = tanh(rat_mean_with_z);
rat_without_r = tanh(rat_mean_without_z);

for i = 1:nA
    if isfinite(rat_with_r(i)) && isfinite(rat_without_r(i))
        plot([1 2], [rat_with_r(i) rat_without_r(i)], '-', 'Color', cmap(i,:), 'LineWidth', 2.5);
        plot(1, rat_with_r(i),    'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
        plot(2, rat_without_r(i), 'o', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
    end
end

bar(1, tanh(mean(Zw_rats,'omitnan')), 0.6, 'FaceColor',[0.35 0.70 1.00], 'EdgeColor','k');
bar(2, tanh(mean(Za_rats,'omitnan')), 0.6, 'FaceColor',[0.85 0.55 0.25], 'EdgeColor','k');

xlim([0.5 2.5]); xticks([1 2]);
xticklabels({'with task','without task'});
ylabel('Mean per-cell spatial-map correlation (r)');
title(titleStr);
yline(0,'k:'); grid on; box on

yl = ylim; y = yl(2) + 0.06*range(yl);
line([1 2], [y y], 'Color','k', 'LineWidth', 1.2);

txtA = sprintf('per-rat paired t (z): t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d', dfA, tA, pA, dzA, numel(Zw_rats));
txtP = sprintf('day-flip perm (z): obs \\Delta z=%.3f (r\\approx%.3f), p_{two}=%.3g', perm.obs_z, tanh(perm.obs_z), perm.p_two);

text(1.5, y + 0.02*range(yl), txtA, 'HorizontalAlignment','center');
text(1.5, y - 0.02*range(yl), txtP, 'HorizontalAlignment','center');

ylim([yl(1) y + 0.12*range(yl)]);

fprintf('\n== %s ==\n', titleStr);
fprintf('Primary (per-rat): paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d\n', dfA, tA, pA, dzA, numel(Zw_rats));
fprintf('Secondary (per-day): paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nDays=%d\n', dfD, tD, pD, dzD, numel(Zw_days));
fprintf('Perm (day-flip within rat): obs Δz=%.3f (r≈%.3f), p_two=%.3g (NPerm=%d)\n', perm.obs_z, tanh(perm.obs_z), perm.p_two, numel(perm.perm_z));
end

function [Zw, Za, aid, did] = get_SC_dayZ(R)
Zw=[]; Za=[]; aid=[]; did=[];
for i = 1:numel(R)
    ok = isfield(R(i),'C') && isfield(R(i).C,'SC') && isfield(R(i).C.SC,'byDay') && ...
         isfield(R(i).C.SC.byDay,'withTask') && isfield(R(i).C.SC.byDay.withTask,'zmean') && ...
         isfield(R(i).C.SC.byDay,'withoutTask') && isfield(R(i).C.SC.byDay.withoutTask,'zmean');
    if ~ok, continue; end

    z_with    = R(i).C.SC.byDay.withTask.zmean(:);
    z_without = R(i).C.SC.byDay.withoutTask.zmean(:);

    nD = min(numel(z_with), numel(z_without));
    z_with    = z_with(1:nD);
    z_without = z_without(1:nD);

    m = isfinite(z_with) & isfinite(z_without);
    Zw  = [Zw; z_with(m)]; %#ok<AGROW>
    Za  = [Za; z_without(m)]; %#ok<AGROW>
    aid = [aid; i*ones(nnz(m),1)]; %#ok<AGROW>
    dd = (1:nD).';
    did = [did; dd(m)]; %#ok<AGROW>
end
end

function OUT = build_flip_perm_singleCell_byDay(Zw_all, Za_all, aid, ~, NPerm)
OUT = struct('obs_z',NaN,'perm_z',[],'p_right',NaN,'p_left',NaN,'p_two',NaN);
if nargin < 5 || isempty(NPerm), NPerm = 500; end

nA = max(aid);
rat_mu = nan(nA,1);
rat_days = cell(nA,1);

for i = 1:nA
    m = (aid == i) & isfinite(Zw_all) & isfinite(Za_all);
    if ~any(m), continue; end
    dz = Zw_all(m) - Za_all(m);
    rat_mu(i) = mean(dz,'omitnan');
    rat_days{i} = dz(:);
end

keepA = isfinite(rat_mu);
if nnz(keepA) < 1
    return
end

obs = mean(rat_mu(keepA), 'omitnan');
OUT.obs_z = obs;

perm = nan(NPerm,1);
for p = 1:NPerm
    mu_perm = nan(nnz(keepA),1);
    jj = 1;
    for i = 1:nA
        if ~keepA(i), continue; end
        dz = rat_days{i};
        s  = (rand(size(dz)) > 0.5)*2 - 1;
        mu_perm(jj) = mean(s .* dz, 'omitnan');
        jj = jj + 1;
    end
    perm(p) = mean(mu_perm, 'omitnan');
end

OUT.perm_z  = perm;
OUT.p_right = mean(perm >= obs);
OUT.p_left  = mean(perm <= obs);
OUT.p_two   = 2*min(OUT.p_right, OUT.p_left);
end

function plot_taskToSpace_SC_deltaPermHist_perRat(R, varargin)
p = inputParser;
addParameter(p,'NBins',30, @(x) isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'ShowTwoSided',true, @(x) islogical(x)&&isscalar(x));
parse(p,varargin{:});
NBins = p.Results.NBins;
ShowTwoSided = p.Results.ShowTwoSided;

nRats = numel(R);
figure('Color','w','Position',[220 220 1100 700]);

for ii = 1:nRats
    if ~isfield(R(ii),'C') || isempty(R(ii).C) || ~isfield(R(ii).C,'SC') || isempty(R(ii).C.SC)
        continue;
    end
    if ~isfield(R(ii).C.SC,'delta') || isempty(R(ii).C.SC.delta)
        continue;
    end

    dlt = R(ii).C.SC.delta;
    if ~isfield(dlt,'obs_z') || ~isfield(dlt,'perm_z') || isempty(dlt.perm_z)
        continue;
    end

    obs  = dlt.obs_z;
    perm = dlt.perm_z;
    perm = perm(isfinite(perm));

    if ~isfinite(obs) || isempty(perm), continue; end

    subplot(3,2,ii); hold on;
    histogram(perm, NBins, 'Normalization','pdf', 'EdgeColor','none');
    yl = ylim;
    plot([obs obs], yl, 'k-', 'LineWidth', 2);

    ratName = '';
    if isfield(R(ii),'animal'), ratName = R(ii).animal; else, ratName = sprintf('rat%02d', ii); end

    xlabel('\Delta z (single-cell) = WITH - WITHOUT');
    ylabel('Null density');
    title(sprintf('%s: SC null', ratName));
    grid on; box on

    pR = NaN; pL = NaN; pT = NaN;
    if isfield(dlt,'p_right'), pR = dlt.p_right; end
    if isfield(dlt,'p_left'),  pL = dlt.p_left;  end
    if isfield(dlt,'p_two'),   pT = dlt.p_two;   end

    txt = sprintf('obs \\Delta z = %.3f (r\\approx%.3f)\n', obs, tanh(obs));
    if ShowTwoSided
        if isfinite(pT), txt = [txt, sprintf('p_{two} = %.3g', pT)];
        else, txt = [txt, sprintf('p_{two} \\approx %.3g', 2*min(mean(perm>=obs), mean(perm<=obs)))];
        end
    else
        if isfinite(pR), txt = [txt, sprintf('p_{right} = %.3g', pR)];
        else, txt = [txt, sprintf('p_{right} \\approx %.3g', mean(perm>=obs))];
        end
    end
    text(obs, yl(2), ['  ', txt], 'VerticalAlignment','top', 'HorizontalAlignment','left');
    ylim(yl);
end
end
