function ATTS = run_space_to_task_interference(ratNames, varargin)
% compute_task_temporal_sameVsDiffSpace
% Self-contained across-trials, time-locked (bin b ↔ bin b) PV similarity.
% Splits TRACE into B temporal bins; for each b:
%   SAME-space: (trial1,b,k) vs (trial2,b,k)
%   DIFF-space: (trial1,b,k1) vs (trial2,b,k2), k1≠k2
%
% Call:
%   ATTS = compute_task_temporal_sameVsDiffSpace({'rat0222','rat0307',...});
%
% Params (name/value):
%   'TraceWin'     [0 2]      % seconds relative to CS
%   'TimeBins'     10         % number of temporal bins in TraceWin
%   'GridN'        [7 7]      % spatial grid (nx,ny) for x/y arena binning
%   'UseSpeedMask' true       % only keep frames with speed >= VelThresh
%   'VelThresh'    4          % cm/s
%   'CellNorm'     'demean'   % 'none'|'demean'|'zscore'
%   'MinFrames'    2          % min frames in (temporal-bin ∩ spatial-bin)
%   'MinNCorr'     10         % min overlapping cells for corr
%   'MinStdBin'    1e-3       % min std for each PV vector
%   'MinTrials'    2          % need at least this many trials/day
%
% Output: ATTS (struct)
%   .rats(i).days(d).byBin(b).sameZ / diffZ  (vectors of Fisher-z)
%   .rats(i).days(d).same.z(b) / diff.z(b)   (day means per bin, Fisher-z)
%   .rats(i).days(d).same.r(b) / diff.r(b)   (Fisher backtransform)
%   .rats(i).delta.obs_z / .delta.perm_z / p-values  <-- NEW Δ-based test
%   .group.same.z(b) / diff.z(b) / delta.z(b)
%   .group.same.r(b) / diff.r(b) / delta.r(b)
%   .distStats.perm.* (Δ-based nulls)
%   .params

% ------------ options -------------
p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'TimeBins',15);
addParameter(p,'GridN',[2 2]);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'VelThresh',4);
addParameter(p,'CellNorm','demean');
addParameter(p,'MinFrames',2);
addParameter(p,'MinNCorr',0);
addParameter(p,'MinStdBin',0);
addParameter(p,'MinTrials',0);
addParameter(p,'NPerm',500);           % <-- NEW: expose permutation count
parse(p,varargin{:});
o = p.Results;

if ischar(ratNames) || isstring(ratNames), ratNames = cellstr(ratNames); end
nR = numel(ratNames);

% containers for group summaries
B = o.TimeBins;
all_same_z = cell(B,1);
all_diff_z = cell(B,1);

ATTS = struct();
ATTS.params = o;
ATTS.rats   = repmat(struct('name','','days',[],'delta',struct()), nR, 1);

for ri = 1:nR
    rname = ratNames{ri};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', rname))
        warning('Variable %s not found in base workspace. Skipping.', rname);
        continue;
    end
    RAT   = evalin('base', rname);              % pull rat struct from base
    ATTS.rats(ri).name = rname;

    % pick 3 days around rat.An
    dates = autoDateList(RAT);
    idx   = [];
    if isfield(RAT,'An')
        idx = find(strcmp(dates, RAT.An), 1);
    end
    if isempty(idx)
        dayIdx = max(1, numel(dates)-2):numel(dates);
    elseif idx < 3
        dayIdx = max(1, idx-2):idx;
    else
        dayIdx = (idx-2):idx;
    end
    days = dates(dayIdx);

    ATTS.rats(ri).days = repmat(struct('label','', 'same',[], 'diff',[], 'byBin',[]), numel(days), 1);

    for di = 1:numel(days)
        dlabel = days{di};
        ATTS.rats(ri).days(di).label = dlabel;

        % ---- extract day data (matches your field layout) ----
        spk  = RAT.Ca_peaks.(sprintf('CA_peaks_%s', dlabel));
        posM = RAT.pos.(sprintf('pos_%s', dlabel));
        csT  = RAT.CS_times.(sprintf('CS_%s', dlabel));

        % ratemask is optional; if missing, keep all
        keepCells = [];
        if isfield(RAT,'ratemask') && isfield(RAT.ratemask, sprintf('ratemask_%s', dlabel))
            rm = RAT.ratemask.(sprintf('ratemask_%s', dlabel));
            keepCells = (rm == 1);
        end
        if ~isempty(keepCells)
            spk = spk(keepCells, :);
        end

        % coerce position and timebase
        [ts, x, y] = coerce_pos(posM);  % ts in seconds
        if numel(ts) < 5 || size(spk,1) < 2, continue; end

        % coerce CS vector
        cs_vec = coerce_cs(csT);
        if numel(cs_vec) < o.MinTrials, continue; end

        % spikes → per-cell spike-time cell array
        spikeCell = to_cell_spikes(spk);

        % ---- speed mask (optional) ----
        if o.UseSpeedMask
            v = local_speed(ts, x, y); % cm/s
            speed_ok = v >= o.VelThresh;
        else
            speed_ok = true(size(ts));
        end

        % ---- build spatial grid on the fly (nx × ny) ----
        nx = o.GridN(1); ny = o.GridN(2);
        xedges = linspace(min(x,[],'omitnan'), max(x,[],'omitnan'), nx+1);
        yedges = linspace(min(y,[],'omitnan'), max(y,[],'omitnan'), ny+1);

        % discretize returns 1..nx/ny or NaN for out-of-range
        bX = discretize(x, xedges);
        bY = discretize(y, yedges);

        % make a linear spatial-bin index only where both are valid
        bSpace = nan(size(x));
        m = isfinite(bX) & isfinite(bY);
        bSpace(m) = sub2ind([nx, ny], bX(m), bY(m));
        K = nx * ny;

        % ---- split TRACE into temporal bins aligned to frames ----
        tw0 = o.TraceWin(1); tw1 = o.TraceWin(2);
        nT  = numel(cs_vec);
        framesByTrial = cell(nT, B, K);

        for tr = 1:nT
            t0 = cs_vec(tr) + tw0;
            t1 = cs_vec(tr) + tw1;

            % all frames in the TRACE window (no speed filter yet)
            idx_all = find(ts >= t0 & ts < t1);
            if isempty(idx_all), continue; end

            % split the frame indices into B consecutive blocks, exactly on frame edges
            e = round(linspace(1, numel(idx_all)+1, B+1));  % 1..(len+1)
            e(1) = 1; e(end) = numel(idx_all)+1;

            for bti = 1:B
                if e(bti) >= e(bti+1), continue; end
                seg = idx_all(e(bti):e(bti+1)-1);
                if o.UseSpeedMask, seg = seg(speed_ok(seg)); end
                if isempty(seg), continue; end

                % spatial bin per frame (already computed as bSpace on full ts)
                bs_b = bSpace(seg);

                % stash per-(trial, temporal-bin, spatial-bin) frame lists
                for k = 1:K
                    sel = (bs_b == k);
                    if any(sel)
                        framesByTrial{tr, bti, k} = seg(sel);
                    end
                end
            end
        end

        % ---- per-(trial,b,k) PVs (mean over frames in the exact intersection) ----
        PV = cell(nT, B, K);  % each entry: [nCells x 1]
        S  = spikes_to_matrix(spikeCell, ts);  % frames on ts
        S  = normalize_cells(S, o.CellNorm);
        for tr = 1:nT
            for bti = 1:B
                for k = 1:K
                    f = framesByTrial{tr,bti,k};
                    if numel(f) >= o.MinFrames
                        PV{tr,bti,k} = mean(S(:, f), 2, 'omitnan');
                    end
                end
            end
        end

        % ---- SAME vs DIFF per temporal bin (strict b↔b across trials) ----
        same_z_mu = nan(B,1); diff_z_mu = nan(B,1);
        day_byBin = repmat(struct('sameZ',[],'diffZ',[]), B, 1);

        for bti = 1:B
            z_same = []; z_diff = [];

            % SAME: (t1,b,k) vs (t2,b,k)
            for k = 1:K
                trs = find(~cellfun(@isempty, squeeze(PV(:,bti,k))));
                if numel(trs) < 2, continue; end
                for i1 = 1:numel(trs)-1
                    t1 = trs(i1); v1 = PV{t1,bti,k};
                    for i2 = i1+1:numel(trs)
                        t2 = trs(i2); v2 = PV{t2,bti,k};
                        [r,~] = safe_corr(v1, v2, o.MinNCorr, o.MinStdBin);
                        if isfinite(r)
                            z_same(end+1,1) = atanh(max(min(r,0.999999),-0.999999)); %#ok<AGROW>
                        end
                    end
                end
            end

            % DIFF: (t1,b,k1) vs (t2,b,k2), k1≠k2
            for k1 = 1:K
                tr1 = find(~cellfun(@isempty, squeeze(PV(:,bti,k1))));
                if isempty(tr1), continue; end
                for k2 = 1:K
                    if k2 == k1, continue; end
                    tr2 = find(~cellfun(@isempty, squeeze(PV(:,bti,k2))));
                    if isempty(tr2), continue; end
                    for i1 = 1:numel(tr1)
                        t1 = tr1(i1); v1 = PV{t1,bti,k1};
                        for i2 = 1:numel(tr2)
                            t2 = tr2(i2); if t2 == t1, continue; end
                            v2 = PV{t2,bti,k2};
                            [r,~] = safe_corr(v1, v2, o.MinNCorr, o.MinStdBin);
                            if isfinite(r)
                                z_diff(end+1,1) = atanh(max(min(r,0.999999),-0.999999)); %#ok<AGROW>
                            end
                        end
                    end
                end
            end

            day_byBin(bti).sameZ = z_same;
            day_byBin(bti).diffZ = z_diff;

            if ~isempty(z_same), same_z_mu(bti) = mean(z_same,'omitnan'); end
            if ~isempty(z_diff), diff_z_mu(bti) = mean(z_diff,'omitnan'); end

            all_same_z{bti} = [all_same_z{bti}; z_same];
            all_diff_z{bti} = [all_diff_z{bti}; z_diff];
        end

        % stash day summary (z + r)
        ATTS.rats(ri).days(di).byBin = day_byBin;
        ATTS.rats(ri).days(di).same = struct('z',same_z_mu, 'r',tanh(same_z_mu));
        ATTS.rats(ri).days(di).diff = struct('z',diff_z_mu, 'r',tanh(diff_z_mu));
    end

    % ===== Δ-based permutation (SAME − DIFF) at the ANIMAL level =====
    % We build a *day-level* Δz by averaging SAME/DIFF in z across bins,
    % then the ANIMAL statistic is the mean Δ across this animal’s days.
    % Null: within each day, pool SAME and DIFF z’s per bin, randomly
    %       reassign n_same / n_diff counts, compute Δ^perm_b, then average
    %       across bins; aggregate to animal mean across its days.

    P = o.NPerm;
    day_delta_obs = [];              % store each day’s Δz (observed)
    day_delta_perm = {};             % each cell: [P x 1] Δz^perm for that day

    Ddays = ATTS.rats(ri).days;
    for di = 1:numel(Ddays)
        Db = Ddays(di).byBin;
        if isempty(Db), continue; end

        % observed Δz for the day (mean over bins)
        z_same_mu = nan(B,1); z_diff_mu = nan(B,1);
        for bti = 1:B
            zs = Db(bti).sameZ; zs = zs(isfinite(zs));
            zd = Db(bti).diffZ; zd = zd(isfinite(zd));
            if ~isempty(zs), z_same_mu(bti) = mean(zs,'omitnan'); end
            if ~isempty(zd), z_diff_mu(bti) = mean(zd,'omitnan'); end
        end
        z_same_mu = z_same_mu(isfinite(z_same_mu));
        z_diff_mu = z_diff_mu(isfinite(z_diff_mu));
        if isempty(z_same_mu) || isempty(z_diff_mu), continue; end

        day_delta_obs(end+1,1) = mean(z_same_mu,'omitnan') - mean(z_diff_mu,'omitnan'); %#ok<AGROW>

        % per-day permutation via label-shuffling within each bin
        perm_d = nan(P,1);
        if P > 0
            for pidx = 1:P
                deltas_b = nan(B,1);
                for bti = 1:B
                    zs = Db(bti).sameZ; zs = zs(isfinite(zs));
                    zd = Db(bti).diffZ; zd = zd(isfinite(zd));
                    if isempty(zs) || isempty(zd), continue; end
                    pool = [zs(:); zd(:)];
                    nS = numel(zs); nD = numel(zd);
                    if numel(pool) < 2 || nS<1 || nD<1, continue; end
                    perm_idx = randperm(numel(pool));
                    same_p = pool(perm_idx(1:nS));
                    diff_p = pool(perm_idx(nS+1:nS+nD));
                    deltas_b(bti) = mean(same_p,'omitnan') - mean(diff_p,'omitnan');
                end
                perm_d(pidx) = mean(deltas_b(isfinite(deltas_b)),'omitnan');
            end
        end
        day_delta_perm{end+1,1} = perm_d; %#ok<AGROW>
    end

    % animal-level statistic + null
    if isempty(day_delta_obs)
        ATTS.rats(ri).delta = struct('obs_z',NaN,'obs_r',NaN,'perm_z',[],'p_left',NaN,'p_right',NaN,'p_two',NaN);
    else
        delta_obs_anim = mean(day_delta_obs,'omitnan');
        % align permutations across days by common P
        Pcommon = min(cellfun(@(v) sum(isfinite(v)), day_delta_perm(~cellfun('isempty',day_delta_perm))));
        if isempty(Pcommon) || ~isfinite(Pcommon) || Pcommon < 5
            perm_anim = [];
            pL = NaN; pR = NaN; pT = NaN;
        else
            perm_anim = nan(Pcommon,1);
            for p = 1:Pcommon
                vals = [];
                for d = 1:numel(day_delta_perm)
                    v = day_delta_perm{d};
                    if ~isempty(v) && numel(v) >= p && isfinite(v(p))
                        vals(end+1,1) = v(p); %#ok<AGROW>
                    end
                end
                perm_anim(p) = mean(vals,'omitnan');
            end
            pR = mean(perm_anim >= delta_obs_anim);
            pL = mean(perm_anim <= delta_obs_anim);
            pT = 2*min(pR,pL);
        end
        pT
        ATTS.rats(ri).delta = struct( ...
            'obs_z',  delta_obs_anim, ...
            'obs_r',  tanh(delta_obs_anim), ...
            'perm_z', perm_anim, ...
            'p_left', pL, 'p_right', pR, 'p_two', pT ...
        );
    end
end

% ---- group summaries (mean across all days from all rats) ----
% gather day means per bin
Zs = []; Zd = [];
for ri = 1:nR
    for di = 1:numel(ATTS.rats(ri).days)
        s = ATTS.rats(ri).days(di).same;
        d = ATTS.rats(ri).days(di).diff;
        if isempty(s) || isempty(d), continue; end
        Zs = [Zs; s.z(:)']; %#ok<AGROW>
        Zd = [Zd; d.z(:)']; %#ok<AGROW>
    end
end
group_same_z = mean(Zs,1,'omitnan'); group_diff_z = mean(Zd,1,'omitnan');

ATTS.group = struct( ...
    'same',  struct('z', group_same_z(:),          'r', tanh(group_same_z(:))), ...
    'diff',  struct('z', group_diff_z(:),          'r', tanh(group_diff_z(:))), ...
    'delta', struct('z', group_same_z(:)-group_diff_z(:), ...
                    'r', tanh(group_same_z(:)) - tanh(group_diff_z(:))) ...
);
ATTS.allPairs = struct('same_z',{all_same_z}, 'diff_z',{all_diff_z});

% stash Δ-based permutation collections for convenience (per-animal)
ATTS.distStats = struct();
ATTS.distStats.delta_byAnimal = arrayfun(@(r) r.delta, ATTS.rats, 'uni', 0);

fprintf('[Temporal b-matched] B=%d | mean Δz over bins = %.3f (r≈%.3f)\n', ...
    B, mean(ATTS.group.delta.z,'omitnan'), mean(ATTS.group.delta.r,'omitnan'));

plot_taskTemporal_sameVsDiffSpace_group(ATTS);

end

% -------------------- helpers (self-contained) --------------------

function [ts,x,y] = coerce_pos(pos)
% Accepts table or numeric; returns column vectors ts,x,y (double)
if istable(pos)
    vn = lower(string(pos.Properties.VariableNames));
    tname = pick(vn, ["t","time","ts"]);
    xname = pick(vn, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"]);
    yname = pick(vn, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"]);
    ts = getcol(pos, tname); x = getcol(pos, xname); y = getcol(pos, yname);
else
    pos = double(pos);
    if size(pos,2) >= 3
        ts = pos(:,1); x = pos(:,2); y = pos(:,3);
    elseif size(pos,2) == 2
        x = pos(:,1); y = pos(:,2); ts = (1:size(pos,1))';
    else
        error('pos must be [n×3] or [n×2]');
    end
end
ts = double(ts(:)); x = double(x(:)); y = double(y(:));
end

function v = local_speed(ts,x,y)
% cm/s from (x,y) samples at times ts
dt = diff(ts); dt(~isfinite(dt) | dt<=0) = NaN;
dx = diff(x);  dy = diff(y);
spd = sqrt(dx.^2 + dy.^2) ./ dt;         % units of pos per second (assume cm)
v = [spd; spd(end)];
v(~isfinite(v)) = 0;
end

function cs = coerce_cs(csIn)
% table or numeric vector -> seconds (col)
if istable(csIn)
    vn = lower(string(csIn.Properties.VariableNames));
    cname = pick(vn, ["cs_ms","cs_time","cstime","cs","onset","onsets","cs_onset","cs_time_ms","cue_onset","time","ts"]);
    cs = getcol(csIn, cname);
else
    cs = double(csIn(:));
end
% convert ms to s if it looks like ms
if max(cs,[],'omitnan') > 1e4, cs = cs/1000; end
end

function cellSpk = to_cell_spikes(spk)
% Accepts [nCells x ?] matrix of spike times (NaNs), or cell-of-vectors
if iscell(spk), cellSpk = spk; return; end
[nc, ~] = size(spk);
cellSpk = cell(nc,1);
for c = 1:nc
    st = double(spk(c,:));
    st = st(isfinite(st) & st > 0);
    cellSpk{c} = st(:);
end
end

function S = spikes_to_matrix(spikeCell, ts)
% Histogram spikes into frames on ts (uniform edges using median dt)
ts = double(ts(:));
if numel(ts) < 2
    S = zeros(numel(spikeCell), numel(ts), 'single'); return;
end
dt = median(diff(ts),'omitnan'); if ~isfinite(dt) || dt<=0, dt = max(eps, mean(diff(ts),'omitnan')); end
edges = [ts - dt/2; ts(end) + dt/2];
nc = numel(spikeCell);
S  = zeros(nc, numel(ts), 'single');
for c = 1:nc
    st = spikeCell{c};
    if isempty(st), continue; end
    S(c,:) = histcounts(st, edges);
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
        error('CellNorm must be ''none'',''demean'',''zscore''.');
end
end

function [r, nC] = safe_corr(a,b,minN,minStd)
% Robust to tiny mismatches or degenerate vectors; prevents mask-size errors.
a = a(:); b = b(:);
if numel(a) ~= numel(b)
    n = min(numel(a), numel(b)); a = a(1:n); b = b(1:n);
end
if nargin < 3 || isempty(minN),   minN = 3;      end
if nargin < 4 || isempty(minStd), minStd = 1e-12; end
m = isfinite(a) & isfinite(b);
nC = nnz(m);
if nC < minN, r = NaN; return, end
a = a(m); b = b(m);
sa = std(a); sb = std(b);
if ~isfinite(sa) || ~isfinite(sb) || sa < minStd || sb < minStd
    r = 0; return
end
r = corr(a,b,'Type','Pearson');
end

function v = getcol(T, name)
if strlength(name)==0, v = []; return; end
v = T{:, find(strcmpi(T.Properties.VariableNames, char(name)),1)};
end

function name = pick(names, options)
name = "";
for k = 1:numel(options)
    idx = find(names == options(k), 1);
    if ~isempty(idx), name = names(idx); return; end
end
end

function plot_taskTemporal_sameVsDiffSpace_group(ATTS, varargin)
% Plot SAME-space vs DIFF-space (task, b-matched) with subject-level stats.
%
% Usage:
%   plot_taskTemporal_sameVsDiffSpace_group(ATTS)             % mean over bins (primary)
%   plot_taskTemporal_sameVsDiffSpace_group(ATTS,'Bin',3)     % a specific temporal bin
%
% Options:
%   'Bin'  : integer bin index, or 'all-mean' (default) to average across bins in z.

p = inputParser;
addParameter(p,'Bin','all-mean');  % 'all-mean' or a scalar bin index
parse(p,varargin{:});
Bin = p.Results.Bin;

% ---------- harvest per-day Fisher-z for SAME and DIFF ----------
Zw_day = [];  % SAME (z) per day
Zd_day = [];  % DIFF (z) per day
aid    = [];  % animal index per day

nR = numel(ATTS.rats);
for ri = 1:nR
    D = ATTS.rats(ri).days;
    if isempty(D), continue; end
    for di = 1:numel(D)
        if isempty(D(di)) || ~isfield(D(di),'same') || ~isfield(D(di),'diff'), continue; end

        if isnumeric(Bin)
            z_same = safe_pick(D(di).same.z, Bin);
            z_diff = safe_pick(D(di).diff.z, Bin);
        else
            % mean across bins in z
            z_same = mean(D(di).same.z(isfinite(D(di).same.z)), 'omitnan');
            z_diff = mean(D(di).diff.z(isfinite(D(di).diff.z)), 'omitnan');
        end

        if isfinite(z_same) && isfinite(z_diff)
            Zw_day(end+1,1) = z_same; %#ok<AGROW>
            Zd_day(end+1,1) = z_diff; %#ok<AGROW>
            aid(end+1,1)    = ri;     %#ok<AGROW>
        end
    end
end

if isempty(Zw_day)
    warning('No SAME/DIFF day-level data found to plot.'); return
end

% ---------- PRIMARY STATS: subject-level (per-rat) ----------
nA = nR;
rat_mean_same_z = nan(nA,1);
rat_mean_diff_z = nan(nA,1);
for ri = 1:nA
    m = (aid==ri);
    if any(m)
        rat_mean_same_z(ri) = mean(Zw_day(m), 'omitnan');
        rat_mean_diff_z(ri) = mean(Zd_day(m), 'omitnan');
    end
end
keepA = isfinite(rat_mean_same_z) & isfinite(rat_mean_diff_z);
Zs_rats = rat_mean_same_z(keepA);
Zd_rats = rat_mean_diff_z(keepA);

% paired t across rats in z
[pA, tA, dfA, dzA] = paired_t_z(Zs_rats, Zd_rats);

% ---------- SECONDARY: pooled per-day paired t in z ----------
keepD = isfinite(Zw_day) & isfinite(Zd_day);
[pD, tD, dfD, dzD] = paired_t_z(Zw_day(keepD), Zd_day(keepD));

% ---------- PLOT ----------
figure('Color','w','Position',[160 160 900 540]); hold on
cmap = lines(nA);
make_light = @(c,frac) (1-frac)*c + frac*[1 1 1];

% 1) per-day light lines (in r)
rw_day = tanh(Zw_day);
rd_day = tanh(Zd_day);
for ri = 1:nA
    c_base  = cmap(ri,:);
    c_light = make_light(c_base, 0.60);
    idx = find(aid==ri & isfinite(rw_day) & isfinite(rd_day));
    for j = 1:numel(idx)
        plot([1 2], [rw_day(idx(j)) rd_day(idx(j))], '-', 'Color', c_light, 'LineWidth', 1.0);
        plot(1, rw_day(idx(j)), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
        plot(2, rd_day(idx(j)), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
    end
end

% 2) per-rat dark means (in r)
rat_same_r = tanh(rat_mean_same_z);
rat_diff_r = tanh(rat_mean_diff_z);
for ri = 1:nA
    if isfinite(rat_same_r(ri)) && isfinite(rat_diff_r(ri))
        plot([1 2], [rat_same_r(ri) rat_diff_r(ri)], '-', 'Color', cmap(ri,:), 'LineWidth', 2.5);
        plot(1, rat_same_r(ri), 'o', 'MarkerFaceColor', cmap(ri,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
        plot(2, rat_diff_r(ri), 'o', 'MarkerFaceColor', cmap(ri,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
    end
end

% 3) group bars (mean across rats in z, back-transform)
bar(1, tanh(mean(Zs_rats,'omitnan')), 0.6, 'FaceColor',[0.30 0.60 1.00], 'EdgeColor','k');
bar(2, tanh(mean(Zd_rats,'omitnan')), 0.6, 'FaceColor',[0.85 0.40 0.20], 'EdgeColor','k');

xlim([0.5 2.5]); xticks([1 2]);
xticklabels({'SAME space','DIFF space'});
ylabel('PV correlation (r)');
title('Task (temporal b-matched): SAME vs DIFF space');
yline(0,'k:'); grid on; box on

% annotate stats (primary + secondary)
yl = ylim; y = yl(2) + 0.06*range(yl);
line([1 2], [y y], 'Color','k','LineWidth',1.2);
txt1 = sprintf('per-rat paired t (z): t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d', dfA, tA, pA, dzA, numel(Zs_rats));
text(1.5, y + 0.02*range(yl), txt1, 'HorizontalAlignment','center');
ylim([yl(1) y + 0.10*range(yl)]);

% also print secondary per-day result to console
fprintf('\n== SAME vs DIFF (b-matched) ==\n');
fprintf('Primary (per-rat):  paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d\n', dfA, tA, pA, dzA, numel(Zs_rats));
fprintf('Secondary (per-day): paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nDays=%d\n', dfD, tD, pD, dzD, nnz(keepD));

% ------------- helpers -------------
function v = safe_pick(vec, k)
    v = NaN;
    if ~isvector(vec) || k<1 || k>numel(vec), return; end
    v = vec(k);
end
end

function [p, t, df, dz] = paired_t_z(z1, z2)
% Paired t-test in z space + Cohen's dz.
mask = isfinite(z1) & isfinite(z2);
z1 = z1(mask); z2 = z2(mask);
if numel(z1) < 2
    p = NaN; t = NaN; df = 0; dz = NaN; return
end
[~, p, ~, st] = ttest(z1, z2);   % two-sided
t  = st.tstat;
df = st.df;
d  = z1 - z2;
dz = mean(d,'omitnan') / std(d, 0, 'omitnan');
end
