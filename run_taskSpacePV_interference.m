function R = run_taskSpacePV_interference(ratNames, varargin)
  % run_taskSpacePV_interference
  % Implements mentor's suggestion in two parts:
  %   (A) Task→Space interference: corr( PV_task(t), PV_control(same location) ) as a time course
  %   (B) Space→Task intrusion: within-trace within-block vs across-block PV similarity.
  %
  % Control spatial PVs are built from NON-task time, excluding [CS, CS+BufferPost] per trial.
  %
  % Usage:
  %   R = run_taskSpacePV_interference({'rat0222','rat0307',...}, 'DoPlots',false);
  %
  % Key options (name/value):
  %   'DaysMode'       : 'last3toAn' | 'all' (default 'last3toAn')
  %   'TraceWin'       : [0 2]                  % seconds from CS
  %   'PrePostWin'     : [-4 16]                % window for time-course (A)
  %   'BufferPost'     : 4                      % exclude Task+[0 BufferPost] from control
  %   'binSize'        : 1/7.5                  % Ca frame/bin (s)
  %   'GridRC'         : [2 3]                  % fallback rows×cols if NumBins=[]

  %   'NumBins'        : []                     % NEW: if set (e.g., 4), auto-grid trial ROI into this many bins
  %   'UseTrialROI'    : true                   % NEW: grid over trial-position ROI (else arena extents)
  %   'ROIPrc'         : [5 95]                 % NEW: percentile bounds for ROI from trace samples
  %   'ROIMarginFrac'  : 0.05                   % NEW: pad ROI by this fraction

  %   'MinOcc'         : 0.5                    % seconds of control occupancy required per bin
  %   'VelThresh'      : 4                      % cm/s minimum (optional kinematic filter)
  %   'DoPlots'        : true
  %   'NPerm'          : 500                    % permutations for nulls (B)
  %   'UseSpeedMask'   : true                   % if true, apply speed>=VelThresh
  %   'CtrlSplitMode'  : 'interleaved'          % even/odd split for reliability
  %
  % Output R: struct per animal...
  %   .A.t, .A.mean_r, .A.n, .A.ctrl_reliability
  %   .B.within, .B.across, .B.delta, .B.p_perm, .B.valid_blocks
  %   .meta.options, .meta.days, .meta.ctrl_reliability, .meta.grid

  % ---------------- Options ----------------
  p = inputParser;
  addParameter(p,'DaysMode','last3toAn');
  addParameter(p,'TraceWin',[0 2]);
  addParameter(p,'PrePostWin',[-4 6]);
  addParameter(p,'BufferPost',0);
  addParameter(p,'binSize',1/7.5);

  % grid / ROI controls
  addParameter(p,'GridRC',[4 4]);        % legacy fallback
  addParameter(p,'NumBins',[4]);          % NEW
  addParameter(p,'UseTrialROI',true);    % NEW
  addParameter(p,'ROIPrc',[10 90]);       % NEW
  addParameter(p,'ROIMarginFrac',0.03);  % NEW
  addParameter(p,'ROIByDay',true);   % NEW: build a day-specific trial ROI/grid


  addParameter(p,'MinOcc',1/7.5);
  addParameter(p,'VelThresh',4);
  addParameter(p,'DoPlots',true);
  addParameter(p,'NPerm',500);
  addParameter(p,'UseSpeedMask',true);
  addParameter(p,'CtrlSplitMode','interleaved');
  addParameter(p,'CellNorm','demean');   % 'zscore' | 'demean' | 'none'

  addParameter(p,'MinTraceSamples',3);        % total even+odd trace samples required per bin
  addParameter(p,'MinTraceEachHalf',1);       % require at least this many in each half (even & odd)
  addParameter(p,'BinWeight','traceCount');   % 'none' | 'traceCount' | 'ctrlReliability' | 'traceCount*ctrlReliability'


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

      % ----- figure out which days to use -----
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



      % --- pull and STANDARDIZE to cell-per-day format ---
      spikes_raw = filterFieldsByDay(rat.Ca_peaks, daysToUse);
      ts_raw     = filterFieldsByDay(rat.Ca_ts,    daysToUse);
      pos_raw    = filterFieldsByDay(rat.pos,      daysToUse);
      cs_raw     = filterFieldsByDay(rat.CS_times, daysToUse);
      [spikes, ts, pos, cs] = standardize_day_format(spikes_raw, ts_raw, pos_raw, cs_raw, daysToUse);

      % ----- build CONTROL spatial PV (exclude Task+[0 BufferPost]) -----
      ctrl = compute_control_spatialPV(spikes, ts, pos, cs, ...
          'TraceWin',        opt.TraceWin, ...
          'BufferPost',      opt.BufferPost, ...
          'binSize',         opt.binSize, ...
          'GridRC',          opt.GridRC, ...
          'NumBins',         opt.NumBins, ...
          'UseTrialROI',     opt.UseTrialROI, ...
          'ROIByDay',        opt.ROIByDay, ...        % << NEW
          'ROIPrc',          opt.ROIPrc, ...
          'ROIMarginFrac',   opt.ROIMarginFrac, ...
          'MinOcc',          opt.MinOcc, ...
          'VelThresh',       opt.VelThresh, ...
          'UseSpeedMask',    opt.UseSpeedMask, ...
          'SplitMode',       opt.CtrlSplitMode, ...
          'CellNorm', opt.CellNorm, ...
          'Label', ratVar);

      % ----- (A) Task→Space interference -----
      A = compute_task_to_space_timecourse(spikes, ts, pos, cs, ctrl, ...
          'PrePostWin',   opt.PrePostWin, ...
          'binSize',      opt.binSize, ...
          'VelThresh',    opt.VelThresh, ...
          'UseSpeedMask', opt.UseSpeedMask, ...
          'CellNorm',     opt.CellNorm, ...
          'TraceWin',     opt.TraceWin, ...      % <- add
          'PreTestWin',   [-2 0], ...            % <- add (adjust if you prefer)
          'Label',        ratVar);

          if isfield(A,'pretrace') && ~isempty(A.pretrace.ttest)
              tt = A.pretrace.ttest;
              fprintf('[%s] Pre (%.1f–%.1fs) vs Trace (%.1f–%.1fs): paired t(%d)=%.2f, p=%.3g, nTrials=%d\n', ...
                  ratVar, A.pretrace.win_pre(1), A.pretrace.win_pre(2), ...
                  A.pretrace.win_trace(1), A.pretrace.win_trace(2), ...
                  ifelse(isnan(tt.df), -1, tt.df), ifelse(isnan(tt.t), NaN, tt.t), ...
                  ifelse(isnan(tt.p), NaN, tt.p), A.pretrace.n_pairs);
          end


      % ----- (B) Space→Task intrusion -----
      B = compute_space_to_task_intrusion(spikes, ts, pos, cs, ctrl, ...
          'TraceWin',     opt.TraceWin, ...
          'binSize',      opt.binSize, ...
          'NPerm',        opt.NPerm, ...
          'VelThresh',    opt.VelThresh, ...
          'UseSpeedMask', opt.UseSpeedMask, ...
          'CellNorm',     opt.CellNorm, ...
          'MinTraceSamples',   opt.MinTraceSamples, ...
          'MinTraceEachHalf',  opt.MinTraceEachHalf, ...
          'BinWeight',         opt.BinWeight, ...
          'Label',        ratVar);


      R(ii).animal = ratVar;
      R(ii).A = A;
      R(ii).B = B;
      R(ii).meta.options = opt;
      R(ii).meta.days = daysToUse;
      R(ii).meta.ctrl_reliability = ctrl.reliability;
      R(ii).meta.grid = ctrl.grid;


  end
    plot_taskSpacePV_results(R);
end


function t = coerce_ts_day(td)
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

function C = to_cell(x)
if iscell(x), C=x; return; end
if isstruct(x) && numel(x)>1, C=arrayfun(@(k) x(k), 1:numel(x), 'uni',0); return; end
if isstruct(x) && numel(x)==1, C={x}; return; end
if istable(x) || isnumeric(x) || isstring(x) || ischar(x) || isempty(x), C={x}; return; end
error('Unsupported container type: %s', class(x));
end

function days = auto_or_all_days(rat)
if isfield(rat,'dates')
    days = rat.dates(:)';
elseif isfield(rat,'pos')
    if iscell(rat.pos), days = arrayfun(@(k) sprintf('day%02d',k), 1:numel(rat.pos), 'uni',0);
    elseif isstruct(rat.pos), days = arrayfun(@(k) sprintf('day%02d',k), 1:numel(rat.pos), 'uni',0);
    else, days = {'day01'};
    end
else
    days = {'day01'};
end
end

function CTRL = compute_control_spatialPV(spikes, ts, pos, cs, varargin)
p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'BufferPost',4);
addParameter(p,'binSize',1/7.5);          % not used for occupancy (we use true dt)
addParameter(p,'GridRC',[2 3]);           % fallback if NumBins=[]
addParameter(p,'NumBins',[]);             % if set, factor into rows×cols
addParameter(p,'UseTrialROI',true);
addParameter(p,'ROIByDay',true);          % per-day ROI/grids
addParameter(p,'ROIPrc',[5 95]);
addParameter(p,'ROIMarginFrac',0.05);
addParameter(p,'MinOcc',0.5);             % seconds of control occupancy required
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'SplitMode','interleaved');
addParameter(p,'Label','');
addParameter(p,'CellNorm','demean');   % 'zscore' | 'demean' | 'none'

parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);

% ---- Build grid edges ----
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

% ---- Control PVs per day ----
PV_by  = cell(1,D);
occ_by = cell(1,D);
rel_by = cell(1,D);
Nc = [];

for d = 1:D
    t = ts{d}(:);
    [x, y] = interp_pos(pos{d}, t);

    % exclude Task + BufferPost from control
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
    t_u = t(use);
    x_u = x(use);
    y_u = y(use);

    [~, b] = pos2bin(x_u, y_u, edgesByDay{d});

    % accumulate true seconds per bin (guard NaNs)
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

    % PV per bin (per-day cell count)
    Nc_d = size(S,1);
    PV   = nan(Nc_d, K);
    for k = 1:K
        sel = (b == k);
        if any(sel)
            PV(:,k) = mean(S(:,sel), 2, 'omitnan');
        end
    end
    % split-half reliability
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

    % per-day diagnostic
    kept = all(isfinite(PV),1);
    if o.UseTrialROI && o.ROIByDay && ~isempty(ROIinfo{d})
        R = ROIinfo{d};
        fprintf(['[%s][CONTROL d=%d] Grid=%dx%d (K=%d via trialROI x=[%.2f %.2f], y=[%.2f %.2f]) | ', ...
                 'MinOcc=%.2fs | kept bins=%d/%d | split-half median r=%.2f\n'], ...
            L, d, rc_eff(1), rc_eff(2), K, R.xmin, R.xmax, R.ymin, R.ymax, ...
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


function B = compute_space_to_task_intrusion(spikes, ts, pos, cs, CTRL, varargin)
% Compute within-bin vs across-bin task PV similarity and permutation tests.
% NullMode:
%   'derange-pairing'   : mild null; break even↔odd pairing across bins with a derangement
%   'frame-redistribute': strong null; reassign task frames across bins (preserving per-bin counts)
%
% RecomputeAcrossInPerm (only for 'derange-pairing'):
%   false (default) keeps ACROSS fixed to the observed value; true recomputes ACROSS per perm.

p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'binSize',1/7.5);
addParameter(p,'NPerm',500);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'CellNorm','demean');          % 'zscore'|'demean'|'none'
addParameter(p,'MinTraceSamples',6);
addParameter(p,'MinTraceEachHalf',1);
addParameter(p,'BinWeight','traceCount');     % 'none'|'traceCount'|'ctrlReliability'|'traceCount*ctrlReliability'
addParameter(p,'NullMode','derange-pairing'); % 'derange-pairing'|'frame-redistribute'
addParameter(p,'RecomputeAcrossInPerm',true);
addParameter(p,'AcrossMode','mean-of-halves') %'cross-halves' or 'mean-of-halves'
addParameter(p,'Label','');

parse(p,varargin{:});
o = p.Results; L = char(o.Label);

D = numel(ts);
K = CTRL.grid.K;

z_within_day = nan(D,1);
z_across_day = nan(D,1);
has_any      = false(D,1);

% caches for permutations
PV_even_cache = cell(1,D);
PV_odd_cache  = cell(1,D);
validK_cache  = cell(1,D);
w_cache       = cell(1,D);

% frame index caches (for 'frame-redistribute')
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

    % keep the exact frame indices per bin (init ONCE per day)
    even_idx_by_bin = cell(K,1);
    odd_idx_by_bin  = cell(K,1);

    csd = cs{d};
    for tr = 1:numel(csd)
        tidx = t >= csd(tr)+o.TraceWin(1) & t <= csd(tr)+o.TraceWin(2);
        tidx = tidx & mask_speed;
        idxs = find(tidx);
        if numel(idxs) < 2, continue; end

        idx_even = idxs(2:2:end);
        idx_odd  = idxs(1:2:end);

        % clip any stray indices (safety)
        idx_even = clamp_idx(idx_even, numel(t));
        idx_odd  = clamp_idx(idx_odd,  numel(t));

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

    % gating by trace samples + control PV availability
    n_tot = n_even + n_odd;
    gate_counts = (n_even >= o.MinTraceEachHalf) & (n_odd >= o.MinTraceEachHalf) & (n_tot >= o.MinTraceSamples);
    validK = find(has_k & gate_counts & all(isfinite(PVd),1)');
    if isempty(validK)
        fprintf('[%s][Intrusion d=%d] no valid bins.\n', L, d);
        continue;
    end
    has_any(d) = true;

    % weights per bin
    switch lower(o.BinWeight)
        case 'none'
            w = ones(K,1);
        case 'tracecount'
            w = n_tot;
        case 'ctrlreliability'
            r = CTRL.reliability_byDay{d};
            if isempty(r), r = zeros(K,1); end
            w = max(r(:),0);
        case {'tracecount*ctrlreliability','tracecount×ctrlreliability'}
            r = CTRL.reliability_byDay{d};
            if isempty(r), r = zeros(K,1); end
            w = n_tot .* max(r(:),0);
        otherwise
            error('Unknown BinWeight: %s', o.BinWeight);
    end
    w(~isfinite(w)) = 0;
    w_valid = w(validK);
    if ~any(w_valid>0), w_valid = ones(numel(validK),1); end

    % observed WITHIN: corr(E_k, O_k)
    zW = nan(numel(validK),1);
    for i = 1:numel(validK)
        k = validK(i);
        zW(i) = atanh(max(min(safe_corr(PV_even(:,k), PV_odd(:,k)),0.99999),-0.99999));
    end
    z_within_day(d) = nansum(w_valid .* zW) / max(1, nansum(w_valid));

    % observed ACROSS: corr( mean(E_k,O_k), mean(E_k',O_k') ), k≠k'
    if numel(validK) >= 2
        pairs = nchoosek(validK(:)',2);
        zA = nan(size(pairs,1),1);
        wA = nan(size(pairs,1),1);
        for i = 1:size(pairs,1)
            k1 = pairs(i,1); k2 = pairs(i,2);
            v1 = mean([PV_even(:,k1), PV_odd(:,k1)],2,'omitnan');
            v2 = mean([PV_even(:,k2), PV_odd(:,k2)],2,'omitnan');
            zA(i) = atanh(max(min(safe_corr(v1, v2),0.99999),-0.99999));
            wA(i) = sqrt(w(k1)*w(k2));
        end
        % after building PV_even, PV_odd, validK, w
        switch lower(o.AcrossMode)  % default
          case 'mean-of-halves'
            % current behavior
            pairs = nchoosek(validK(:)',2);
            zA = nan(size(pairs,1),1); wA = nan(size(pairs,1),1);
            for i = 1:size(pairs,1)
                k1 = pairs(i,1); k2 = pairs(i,2);
                v1 = mean([PV_even(:,k1), PV_odd(:,k1)],2,'omitnan');
                v2 = mean([PV_even(:,k2), PV_odd(:,k2)],2,'omitnan');
                zA(i) = atanh(bound_r(safe_corr(v1, v2)));
                wA(i) = sqrt(w(k1)*w(k2));
            end
          case 'cross-halves'  % symmetric-halves
            pairs = nchoosek(validK(:)',2);
            zA = nan(size(pairs,1),1); wA = nan(size(pairs,1),1);
            for i = 1:size(pairs,1)
                k1 = pairs(i,1); k2 = pairs(i,2);
                r1 = safe_corr(PV_even(:,k1), PV_odd(:,k2));
                r2 = safe_corr(PV_odd(:,k1),  PV_even(:,k2));
                zA(i) = atanh(bound_r(mean([r1 r2],'omitnan')));
                wA(i) = sqrt(w(k1)*w(k2));
            end
          otherwise
            error('Unknown AcrossMode');
        end
        wA(~isfinite(wA))=0; if ~any(wA>0), wA(:)=1; end
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

    pairs_count = max(numel(validK)*(numel(validK)-1)/2, 0);
    fprintf('[%s][Intrusion d=%d] valid bins=%s | (n_even+n_odd) median=%g | across pairs=%d\n', ...
        L, d, mat2str(validK(:)'), median(n_tot(validK)), pairs_count);
end

% average across days (Fisher z)
z_within_mean = mean(z_within_day(has_any), 'omitnan');
z_across_mean = mean(z_across_day(has_any), 'omitnan');
delta_z = z_within_mean - z_across_mean;

% --- permutations
Nperm = o.NPerm;
perm_delta = nan(Nperm,1);
zW_perm_mean_all = [];
zA_perm_mean_all = [];

if Nperm>0 && any(has_any)
    switch lower(o.NullMode)

        case 'derange-pairing'
            for pidx = 1:Nperm
                zW_perm_day = nan(D,1);
                zA_perm_day = nan(D,1);  % only used if RecomputeAcrossInPerm==true
                for d = 1:D
                    if ~has_any(d), continue; end
                    validK = validK_cache{d};
                    Kvalid = numel(validK);
                    if Kvalid < 1, continue; end

                    w = w_cache{d};
                    w_valid = w(validK); if ~any(w_valid>0), w_valid(:)=1; end

                    PV_even = PV_even_cache{d};
                    PV_odd  = PV_odd_cache{d};

                    % derangement over validK indices
                    pi_idx = rand_derangement_idx(Kvalid);
                    map_to = validK(pi_idx);  % O from different bins, no fixed points

                    % WITHIN under null: corr(E_k, O_{pi(k)})
                    zP = nan(Kvalid,1); wP = nan(Kvalid,1);
                    for j = 1:Kvalid
                        k1 = validK(j);
                        k2 = map_to(j);
                        zP(j) = atanh(max(min(safe_corr(PV_even(:,k1), PV_odd(:,k2)),0.99999),-0.99999));
                        wP(j) = sqrt(w(k1) * w(k2));
                    end
                    wP(~isfinite(wP))=0; if ~any(wP>0), wP(:)=1; end
                    zW_perm_day(d) = nansum(wP .* zP) / max(1,nansum(wP));

                    if o.RecomputeAcrossInPerm && Kvalid >= 2
                        pairs = nchoosek(validK(:)',2);
                        zA = nan(size(pairs,1),1);
                        wA = nan(size(pairs,1),1);
                        for i = 1:size(pairs,1)
                            k1 = pairs(i,1); k2 = pairs(i,2);
                            v1 = mean([PV_even(:,k1), PV_odd(:, map_to(validK==k1))],2,'omitnan');
                            v2 = mean([PV_even(:,k2), PV_odd(:, map_to(validK==k2))],2,'omitnan');
                            zA(i) = atanh(max(min(safe_corr(v1, v2),0.99999),-0.99999));
                            wA(i) = sqrt(w(k1)*w(k2));
                        end
                        wA(~isfinite(wA))=0; if ~any(wA>0), wA(:)=1; end
                        zA_perm_day(d) = nansum(wA .* zA) / max(1, nansum(wA));
                    end
                end
                zW_perm_mean = mean(zW_perm_day(has_any),'omitnan');
                zW_perm_mean_all(end+1) = zW_perm_mean;

                if o.RecomputeAcrossInPerm
                    zA_perm_mean = mean(zA_perm_day(has_any),'omitnan');
                    zA_perm_mean_all(end+1) = zA_perm_mean;

                    perm_delta(pidx) = zW_perm_mean - zA_perm_mean;
                else
                    perm_delta(pidx) = zW_perm_mean - z_across_mean;
                end
            end



        case 'frame-redistribute'
            for pidx = 1:Nperm
                zW_perm_day = nan(D,1);
                zA_perm_day = nan(D,1);

                for d = 1:D
                    if ~has_any(d), continue; end
                    validK = validK_cache{d};
                    Kvalid = numel(validK);
                    if Kvalid < 1, continue; end

                    w = w_cache{d};
                    w_valid = w(validK); if ~any(w_valid>0), w_valid(:)=1; end

                    % Build pools (deduplicate & sort for safety)
                    even_idx_by_bin = even_idx_cache{d};
                    odd_idx_by_bin  = odd_idx_cache{d};
                    ne = n_even_cache{d};
                    no = n_odd_cache{d};

                    E_pool = vertcat(even_idx_by_bin{validK});
                    O_pool = vertcat(odd_idx_by_bin{validK});
                    E_pool = unique(clamp_idx(E_pool, numel(ts{d}))); % clamp+dedupe
                    O_pool = unique(clamp_idx(O_pool, numel(ts{d})));
                    if isempty(E_pool) || isempty(O_pool)
                        continue;
                    end

                    % Randomize and split using pointer (prevents overrun)
                    E_perm = E_pool(randperm(numel(E_pool)));
                    O_perm = O_pool(randperm(numel(O_pool)));

                    ce = ne(validK); ce = max(0, round(ce(:)));        % ensure nonneg ints
                    co = no(validK); co = max(0, round(co(:)));

                    % If requested counts exceed pool size (rare), downscale last bin(s)
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

                    % rebuild permuted PV_even / PV_odd for the valid bins
                    % recompute S here to keep scope local & consistent
                    t_day = ts{d}(:);
                    S = spikes_to_matrix(spikes{d}, t_day);
                    S = normalize_cells(S, o.CellNorm);

                    PV_even_p = nan(size(S,1), K);
                    PV_odd_p  = nan(size(S,1), K);

                    for j = 1:Kvalid
                        k = validK(j);
                        segE = clamp_idx(E_chunks{j}, size(S,2));
                        segO = clamp_idx(O_chunks{j}, size(S,2));
                        if ~isempty(segE)
                            PV_even_p(:,k) = mean(S(:, segE), 2, 'omitnan');
                        end
                        if ~isempty(segO)
                            PV_odd_p(:,k)  = mean(S(:, segO),  2, 'omitnan');
                        end
                    end

                    % WITHIN (perm)
                    zW_vec = nan(Kvalid,1);
                    for j = 1:Kvalid
                        k = validK(j);
                        zW_vec(j) = atanh(max(min(safe_corr(PV_even_p(:,k), PV_odd_p(:,k)),0.99999),-0.99999));
                    end
                    zW_perm_day(d) = nansum(w_valid .* zW_vec) / max(1,nansum(w_valid));

                    % ACROSS (perm)
                    if Kvalid >= 2
                        pairs = nchoosek(validK(:)', 2);
                        zA_vec = nan(size(pairs,1),1);
                        wA     = nan(size(pairs,1),1);
                        for iPair = 1:size(pairs,1)
                            k1 = pairs(iPair,1); k2 = pairs(iPair,2);
                            v1 = mean([PV_even_p(:,k1), PV_odd_p(:,k1)],2,'omitnan');
                            v2 = mean([PV_even_p(:,k2), PV_odd_p(:,k2)],2,'omitnan');
                            zA_vec(iPair) = atanh(max(min(safe_corr(v1, v2),0.99999),-0.99999));
                            wA(iPair)     = sqrt(w(k1)*w(k2));
                        end
                        wA(~isfinite(wA))=0; if ~any(wA>0), wA(:)=1; end
                        zA_perm_day(d) = nansum(wA .* zA_vec) / max(1, nansum(wA));
                    else
                        zA_perm_day(d) = NaN;
                    end
                end

                zW_perm_mean = mean(zW_perm_day(has_any),'omitnan');
                zA_perm_mean = mean(zA_perm_day(has_any),'omitnan');
                perm_delta(pidx) = zW_perm_mean - zA_perm_mean;
                zW_perm_mean_all(end+1) = zW_perm_mean;
                zA_perm_mean_all(end+1) = zA_perm_mean;

            end

        otherwise
            error('Unknown NullMode: %s', o.NullMode);

    end



end

fprintf('Obs: within_z=%.3f, across_z=%.3f, Δ_z=%.3f\n', ...
        z_within_mean, z_across_mean, delta_z);

fprintf('Null (perm): mean(within_z)=%.3f (sd=%.3f), mean(across_z)=%.3f (sd=%.3f), mean(Δ_z)=%.3f (sd=%.3f)\n', ...
        mean(zW_perm_mean_all), std(zW_perm_mean_all), ...
        mean(zA_perm_mean_all), std(zA_perm_mean_all), ...
        mean(perm_delta),       std(perm_delta));

% Build p-values
perm_delta = perm_delta(isfinite(perm_delta));
delta_obs  = delta_z;
if isempty(perm_delta)
    p_right = NaN; p_left = NaN; p_two = NaN;
else
    p_right = mean(perm_delta >= delta_obs);
    p_left  = mean(perm_delta <= delta_obs)
    p_two   = 2*min(p_right, p_left);
end

% Pack outputs (convert z back to r)
B.within.r = tanh(z_within_mean);
B.within.z = z_within_mean;
B.across.r = tanh(z_across_mean);
B.across.z = z_across_mean;
B.delta.r  = tanh(delta_z);
B.delta.z  = delta_z;

B.p_perm_right = p_right;
B.p_perm_left  = p_left;
B.p_perm_two   = p_two;
B.perm_delta   = perm_delta;

B.null = struct('mode', o.NullMode, 'recomputeAcross', logical(o.RecomputeAcrossInPerm), ...
                'NPerm', o.NPerm);

B.notes = struct('MinTraceSamples',o.MinTraceSamples, ...
                 'MinTraceEachHalf',o.MinTraceEachHalf, ...
                 'BinWeight',o.BinWeight);
end




% ---------- local helpers (scoped to this file) ----------
function idx = clamp_idx(idx, N)
    if isempty(idx), return; end
    idx = idx(:);
    idx = idx(isfinite(idx));
    idx = idx(idx>=1 & idx<=N);
    idx = unique(round(idx)); % ensure integer, dedupe, sort asc
end

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



function pi = rand_derangement_idx(n)
% Return a random derangement of 1:n (no fixed points).
% Simple rejection sampler; efficient for small n (e.g., n<=100).
if n<=1
    pi = 1:n;  % degenerate; caller should avoid n<=1 cases
    return
end
ok = false;
while ~ok
    pi = randperm(n);
    ok = all(pi ~= 1:n);
end
end



function out = ifelse(cond, a, b), if cond, out=a; else, out=b; end, end


% -------------- small utilities --------------
function r = safe_corr(a,b)
if isempty(a) || isempty(b) || all(~isfinite(a)) || all(~isfinite(b)), r = NaN; return; end
m = isfinite(a) & isfinite(b);
if nnz(m) < 3, r = NaN; return; end
r = corr(a(m), b(m), 'type','Pearson');
end

function x = bound_r(x), x = max(min(x,0.99999), -0.99999); end

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


function [x_i, y_i] = interp_pos(posd, t)
tt = double(posd.t(:)); xx = double(posd.x(:)); yy = double(posd.y(:));
if ~isnumeric(t), t = coerce_ts_day(t); else, t = double(t(:)); end
[ttu, ia] = unique(tt, 'stable'); xxu = xx(ia); yyu = yy(ia);
x_i = interp1(ttu, xxu, t, 'linear','extrap');
y_i = interp1(ttu, yyu, t, 'linear','extrap');
end

function v = speed_cm_per_s(posd)
t  = double(posd.t(:)); x  = double(posd.x(:)); y  = double(posd.y(:));
n  = min([numel(t), numel(x), numel(y)]);
TXY = [t(1:n), x(1:n), y(1:n)];
try
    V = ca_velocity(TXY); v_times = double(V(2,:)).'; v_vals = double(V(1,:)).';
    [v_times_u, ia] = unique(v_times, 'stable'); v_vals_u = v_vals(ia);
    v = interp1(v_times_u, v_vals_u, t(1:n), 'linear', 'extrap');
    if any(~isfinite(v)), v = fillmissing(v,'nearest'); end
catch
    dt = diff(t(1:n)); dt(end+1,1) = median(dt(dt>0),'omitnan');
    dx = [diff(x(1:n)); 0]; dy = [diff(y(1:n)); 0]; v = hypot(dx,dy) ./ max(dt, eps);
end
end

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

function [rc_idx, k] = pos2bin(x, y, edges)
[~, cx] = histc(x, edges.x);  [~, cy] = histc(y, edges.y);
cx(cx<1 | cx>=numel(edges.x)) = NaN;
cy(cy<1 | cy>=numel(edges.y)) = NaN;
rc_idx = [cy, cx];
GridR = numel(edges.y)-1; GridC = numel(edges.x)-1;
k = nan(size(x));
m = isfinite(cx) & isfinite(cy);
k(m) = sub2ind([GridR, GridC], cy(m), cx(m));
end

% --------- plotting (optional) ---------
function plot_task_to_space(A, CTRL, animal)
figure('Color','w'); hold on;
plot(A.t, A.mean_r, 'k-', 'LineWidth', 2);
yline(nanmean(CTRL.reliability.r),'--','Color',[.5 .5 .5], 'LineWidth',1.5);
xlabel('Time from CS (s)'); ylabel('corr(PV_{task}(t), PV_{ctrl}(same place))');
title(sprintf('[%s] Task→Space interference', animal), 'Interpreter','none'); grid on;
end

function plot_space_to_task(B, CTRL, animal)
figure('Color','w'); hold on;
bar(1, B.within.r); bar(2, B.across.r);
xticks([1 2]); xticklabels({'within (trace)', 'across (trace)'});
ylabel('Mean PV correlation (r)');
title(sprintf('[%s] Space→Task intrusion (\\Delta=%.3f, p_{perm}=%.3g)', animal, B.delta.r, B.p_perm_right), 'Interpreter','none'); grid on;


end

function P = coerce_pos_day(pd, tday)
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
            C=candidates{c}; fC=fieldnames(C);
            for j=1:numel(fC)
                val=C.(fC{j});
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

function [t,x,y] = coerce_from_numeric(M, tday)
[nr,nc]=size(M);
if nc==3, t=M(:,1); x=M(:,2); y=M(:,3);
elseif nc==2, x=M(:,1); y=M(:,2); if ~isempty(tday)&&numel(tday)==nr, t=tday; else, t=(1:nr)'; end
else, error('numeric matrix must be [n×3] or [n×2], got [n×%d].', nc);
end
[t,x,y] = finalize_txy(t,x,y,tday);
end

function name = pick_name(names, options)
name=""; for k=1:numel(options), idx=find(names==options(k),1); if ~isempty(idx), name=names(idx); return; end, end
end

function val = get_field_if_exists(S, name)
if strlength(name)==0 || ~isstruct(S), val=[]; return; end
fn=fieldnames(S); idx=find(strcmpi(fn, char(name)),1);
if isempty(idx), val=[]; else, val=S.(fn{idx}); end
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

function plot_taskSpacePV_results(R)
if isempty(R), warning('Empty results.'); return; end
TraceWin = R(1).meta.options.TraceWin; tvec = R(1).A.t; nR = numel(R);

% per-rat
for i = 1:nR
    if ~isfield(R(i),'A') || isempty(R(i).A), continue; end
    A = R(i).A; B = R(i).B;
    animal = R(i).animal;
    ctrlRel = R(i).A.ctrl_reliability.r;
    ctrlMean = mean(ctrlRel,'omitnan');

    figure('Color','w','Position',[120 120 820 560]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile(1,[1 2]); hold on
    y = A.mean_r;
    yl = [-0.2 0.35]; if any(isfinite(y)), yl=[min(y,[],'omitnan')-0.05, max(y,[],'omitnan')+0.05]; end
    patch([TraceWin(1) TraceWin(2) TraceWin(2) TraceWin(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.95 0.95 0.95], 'EdgeColor','none','DisplayName','Trace window');
    plot(A.t, y, 'k-', 'LineWidth', 2, 'DisplayName','corr(PV_{task}(t), PV_{ctrl})');
    yline(ctrlMean,'--','Color',[.5 .5 .5],'LineWidth',1.5,'DisplayName','ctrl split-half');
    xline(0,':','Color',[.5 .5 .5],'LineWidth',1.0,'DisplayName','CS');
    xlabel('Time from CS (s)'); ylabel('Correlation r');
    title(sprintf('[%s] Task→Space interference', animal), 'Interpreter','none');
    grid on; ylim(yl); legend('Location','best');

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

% group summary
T = numel(tvec); M = nan(nR, T); ctrlAll = nan(nR,1);
for i=1:nR
    if isempty(R(i).A), continue; end
    if numel(R(i).A.t) == T && all(R(i).A.t(:)==tvec(:))
        M(i,:) = R(i).A.mean_r(:);
        ctrlAll(i) = mean(R(i).A.ctrl_reliability.r,'omitnan');
    end
end

% baseline vs trace, paired t across rats (Fisher-z)
preMask   = (tvec >= -2 & tvec < -0);
traceMask = (tvec >=  2 & tvec <= 4);

z_group = atanh(M);
pre_z   = mean(z_group(:,preMask),  2,'omitnan');
trace_z = mean(z_group(:,traceMask),2,'omitnan');
ok = isfinite(pre_z) & isfinite(trace_z);
[~, p_preTrace, ~, st] = ttest(trace_z(ok), pre_z(ok));     % paired
dz = mean(trace_z(ok)-pre_z(ok)) / std(trace_z(ok)-pre_z(ok)); % Cohen's dz
fprintf('Timecourse: trace vs pre (z): t(%d)=%.2f, p=%.3g, dz=%.2f\n', st.df, st.tstat, p_preTrace, dz);

okR = all(isfinite(M),2) | any(isfinite(M),2); okR = okR & isfinite(ctrlAll);
if any(okR)
    mu = nanmean(M(okR,:),1); se = nanstd(M(okR,:),0,1)./sqrt(nnz(okR));
    figure('Color','w','Position',[140 140 860 420]); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile(1); hold on
    yl = [min(mu-se,[],'omitnan')-0.05, max(mu+se,[],'omitnan')+0.05];
    patch([TraceWin(1) TraceWin(2) TraceWin(2) TraceWin(1)], [yl(1) yl(1) yl(2) yl(2)], [0.95 0.95 0.95], 'EdgeColor','none');
    xx = [tvec; flipud(tvec)]; yy = [mu(:)-se(:); flipud(mu(:)+se(:))];
    fill(xx, yy, [0 0 0], 'FaceAlpha',0.08, 'EdgeColor','none'); plot(tvec, mu, 'k-', 'LineWidth',2.5);
    yline(mean(ctrlAll(okR)),'--','Color',[.5 .5 .5],'LineWidth',1.5); xline(0,':','Color',[.5 .5 .5]);
    xlabel('Time from CS (s)'); ylabel('Mean correlation r'); title(sprintf('Task→Space (group mean ± SEM; n=%d)', nnz(okR)));
    grid on; ylim(yl);

    nexttile(2); hold on
    w = nan(1,nR); a = nan(1,nR);
    for i=1:nR, if ~isempty(R(i).B), w(i)=R(i).B.within.r; a(i)=R(i).B.across.r; end, end
    mask = isfinite(w) & isfinite(a);
    plot([ones(1,nnz(mask)); 2*ones(1,nnz(mask))], [w(mask); a(mask)], '-o', 'LineWidth',1.4);
    bar(1, mean(w(mask),'omitnan'), 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k');
    bar(2, mean(a(mask),'omitnan'), 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k');
    xticks([1 2]); xticklabels({'within','across'}); ylabel('Correlation r'); box on
    title(sprintf('Space→Task (per-rat & mean; n=%d)', nnz(mask)));

    % Paired t-test across rats on Fisher-z (within vs across)
    wz = arrayfun(@(s) s.B.within.z, R);
    az = arrayfun(@(s) s.B.across.z, R);
    maskZ = isfinite(wz) & isfinite(az);
    if any(maskZ)
        [~, p, ~, stats] = ttest(wz(maskZ), az(maskZ));      % paired t on z
        ax = gca;
        yl = ylim(ax); y = yl(2) + 0.05*range(yl);
        line(ax, [1 2], [y y], 'Color','k','LineWidth',1.2);
        text(1.5, y + 0.02*range(yl), ...
            sprintf('paired t: t(%d)=%.2f, p=%.3g', stats.df, stats.tstat, p), ...
            'HorizontalAlignment','center');
        ylim(ax, [yl(1) y + 0.08*range(yl)]);

        dz = mean(wz(maskZ)-az(maskZ)) / std(wz(maskZ)-az(maskZ));
        fprintf('Paired t-test (Fisher-z): t(%d)=%.2f, p=%.3g, dz=%.2f\n', ...
                stats.df, stats.tstat, p, dz);
    end
end
end


function A = compute_task_to_space_timecourse(spikes, ts, pos, cs, CTRL, varargin)
  p = inputParser;
  addParameter(p,'PrePostWin',[-4 16]);
  addParameter(p,'TraceWin',[0 2]);      % << new (used for the per-trial test)
  addParameter(p,'PreTestWin',[-2 0]);   % << new (window to compare vs TraceWin)
  addParameter(p,'binSize',1/7.5);
  addParameter(p,'VelThresh',4);
  addParameter(p,'UseSpeedMask',false);
  addParameter(p,'CellNorm','demean');
  addParameter(p,'Label','');
  parse(p,varargin{:});
  o = p.Results;

  % ---- Build time axis for pooled time course (unchanged) ----
  t_vec = (o.PrePostWin(1):o.binSize:o.PrePostWin(2))';
  T = numel(t_vec);
  sum_z = zeros(T,1);  cnt = zeros(T,1);

  % masks for the per-trial windows (relative to CS)
  isPreMask   = (t_vec >= o.PreTestWin(1) & t_vec < o.PreTestWin(2));
  isTraceMask = (t_vec >= o.TraceWin(1)   & t_vec < o.TraceWin(2));

  % per-trial accumulators (Fisher-z means)
  z_pre_trials   = [];   % one number per trial (mean of z’s inside PreTestWin)
  z_trace_trials = [];   % one number per trial (mean of z’s inside TraceWin)

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

          % per-trial temporary containers
          z_pre_tp   = [];   % z’s that land in PreTestWin for this trial
          z_trace_tp = [];   % z’s that land in TraceWin for this trial

          % loop over relative time bins
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

              % pooled time course (unchanged)
              sum_z(i) = sum_z(i) + z;
              cnt(i)   = cnt(i) + 1;
              dbg.n_ok = dbg.n_ok + 1;

              % per-trial window assignment
              if isPreMask(i)
                  z_pre_tp(end+1) = z; %#ok<AGROW>
              elseif isTraceMask(i)
                  z_trace_tp(end+1) = z; %#ok<AGROW>
              end
          end

          % finalize per-trial means if both windows have data
          if ~isempty(z_pre_tp) && ~isempty(z_trace_tp)
              z_pre_trials(end+1)   = mean(z_pre_tp,   'omitnan'); %#ok<AGROW>
              z_trace_trials(end+1) = mean(z_trace_tp, 'omitnan'); %#ok<AGROW>
          end
      end
  end

  % ---- pack pooled time course (same as before) ----
  mean_z = sum_z ./ max(1,cnt);
  A.t = t_vec;
  A.mean_z = mean_z;
  A.mean_r = tanh(mean_z);
  A.n = cnt;

  % reliability passthrough (unchanged)
  A.ctrl_reliability.r = CTRL.reliability.r;
  A.ctrl_reliability.z = CTRL.reliability.z;

  % ---- per-rat paired t-test across trials (pre vs trace) on Fisher-z ----
  A.pretrace.win_pre   = o.PreTestWin;
  A.pretrace.win_trace = o.TraceWin;
  A.pretrace.z_pre     = z_pre_trials(:);
  A.pretrace.z_trace   = z_trace_trials(:);
  A.pretrace.r_pre     = tanh(z_pre_trials(:));
  A.pretrace.r_trace   = tanh(z_trace_trials(:));
  A.pretrace.n_pairs   = numel(z_pre_trials);

  if A.pretrace.n_pairs >= 2
      [h,p,~,stats] = ttest(A.pretrace.z_pre, A.pretrace.z_trace);  % paired on z
      A.pretrace.ttest = struct('h',h,'p',p,'t',stats.tstat,'df',stats.df);
  else
      A.pretrace.ttest = struct('h',NaN,'p',NaN,'t',NaN,'df',NaN);
  end

  % debug
  A.debug = dbg;
  A.debug.ctrl_valid_bins = arrayfun(@(d) find(all(isfinite(CTRL.PV{d}),1)), 1:D, 'uni',0);
  A.debug.ctrl_occ_sec    = CTRL.occ;
  A.debug.ctrl_MinOcc     = CTRL.params.MinOcc;
  A.debug.time_binSize    = o.binSize;
end



function [edges, K, GridRC, ROI] = build_grid_edges_from_trialROI(pos, ts, cs, TraceWin, NumBins, prc, marginFrac)
% Build axis-aligned edges covering where the animal is during the trace window
% across all days, then tile that ROI into NumBins bins (rows×cols).

% 1) collect trial positions
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

% Fallback to arena extents if no trial samples found
if isempty(X) || isempty(Y)
    allx = []; ally = [];
    for d = 1:numel(pos), allx = [allx; pos{d}.x(:)]; ally = [ally; pos{d}.y(:)]; end
    allx = allx(isfinite(allx)); ally = ally(isfinite(ally));
    if isempty(allx) || isempty(ally)
        % degenerate fallback
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

% 2) pad ROI a bit
dx = xmax - xmin; dy = ymax - ymin;
xmin = xmin - marginFrac*dx; xmax = xmax + marginFrac*dx;
ymin = ymin - marginFrac*dy; ymax = ymax + marginFrac*dy;

% 3) choose rows×cols to exactly factor NumBins (near-square)
GridRC = best_factors(NumBins);
rows = GridRC(1); cols = GridRC(2);
edges.x = linspace(xmin, xmax, cols+1);
edges.y = linspace(ymin, ymax, rows+1);
K = rows*cols;

ROI = struct('xmin',xmin,'xmax',xmax,'ymin',ymin,'ymax',ymax,'prc',prc,'marginFrac',marginFrac);
end

function rc = best_factors(N)
% Return integer [rows cols] such that rows*cols == N and rows≈sqrt(N)
% Always succeeds (worst case: [1 N]).
r = floor(sqrt(N));
while r > 1 && mod(N,r) ~= 0
    r = r - 1;
end
c = N / r;
rc = [r, c];
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
