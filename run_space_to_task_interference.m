function ATTS = run_space_to_task_interference(ratNames, varargin)
% run_space_to_task_interference
% Self-contained across-trials, time-locked (bin b ↔ bin b) PV similarity.
% Splits TRACE into B temporal bins; for each b:
%   SAME-space: (trial1,b,k) vs (trial2,b,k)
%   DIFF-space: (trial1,b,k1) vs (trial2,b,k2), k1≠k2
%
% Call:
%   ATTS = run_space_to_task_interference({'rat0222','rat0307',...});
%
% Params (name/value):
%   'TraceWin'       [0 2]      % seconds relative to CS
%   'TimeBins'       15         % number of temporal bins in TraceWin
%   'GridN'          [2 2]      % spatial grid [nx ny] for x/y arena binning
%   'UseSpeedMask'   false      % only keep frames with speed >= VelThresh
%   'VelThresh'      4          % cm/s
%   'CellNorm'       'demean'   % 'none'|'demean'|'zscore'
%   'MinFrames'      2          % min frames in (temporal-bin ∩ spatial-bin)
%   'MinNCorr'       0          % min overlapping cells for corr
%   'MinStdBin'      0          % min std for each PV vector
%   'MinTrials'      0          % need at least this many trials/day
%   'NPerm'          500        % permutation count
%   'DeltaPermType'  'flip'     % 'flip' | 'derange' | 'both'
%
% Outputs:
%   ATTS.rats(ri).delta_flip    (if requested)
%   ATTS.rats(ri).delta_derange (if requested)
%   ATTS.rats(ri).delta         (the "active" one, based on DeltaPermType)
%   plus group summaries + plots

% ------------ options -------------
p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'TimeBins',15);
addParameter(p,'GridN',[2 2]);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'VelThresh',4);
addParameter(p,'CellNorm','none');
addParameter(p,'MinFrames',2);
addParameter(p,'MinNCorr',0);
addParameter(p,'MinStdBin',0);
addParameter(p,'MinTrials',0);
addParameter(p,'NPerm',500);
addParameter(p,'DeltaPermType','flip', @(s) any(strcmpi(s,{'flip','derange','both'})));
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
ATTS.rats   = repmat(struct('name','','days',[],'delta',struct(), ...
                            'delta_flip',struct(), 'delta_derange',struct()), nR, 1);

for ri = 1:nR
    rname = ratNames{ri};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', rname))
        warning('Variable %s not found in base workspace. Skipping.', rname);
        continue;
    end
    RAT   = evalin('base', rname);
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

        % ---- extract day data ----
        spk  = RAT.Ca_peaks.(sprintf('CA_peaks_%s', dlabel));
        posM = RAT.pos.(sprintf('pos_%s', dlabel));
        csT  = RAT.CS_times.(sprintf('CS_%s', dlabel));

        % ratemask optional
        keepCells = [];
        if isfield(RAT,'ratemask') && isfield(RAT.ratemask, sprintf('ratemask_%s', dlabel))
            rm = RAT.ratemask.(sprintf('ratemask_%s', dlabel));
            keepCells = (rm == 1);
        end
        if ~isempty(keepCells)
            spk = spk(keepCells, :);
        end

        % coerce position and timebase
        [ts, x, y] = coerce_pos(posM);
        if numel(ts) < 5 || size(spk,1) < 2, continue; end

        % coerce CS vector
        cs_vec = coerce_cs(csT);
        if numel(cs_vec) < o.MinTrials, continue; end

        % spikes -> cell array
        spikeCell = to_cell_spikes(spk);

        % speed mask
        if o.UseSpeedMask
            v = local_speed(ts, x, y);
            speed_ok = v >= o.VelThresh;
        else
            speed_ok = true(size(ts));
        end

        % ---- spatial grid ----
        nx = o.GridN(1); ny = o.GridN(2);
        xedges = linspace(min(x,[],'omitnan'), max(x,[],'omitnan'), nx+1);
        yedges = linspace(min(y,[],'omitnan'), max(y,[],'omitnan'), ny+1);

        bX = discretize(x, xedges);
        bY = discretize(y, yedges);

        bSpace = nan(size(x));
        msp = isfinite(bX) & isfinite(bY);
        % NOTE: sub2ind expects [rows cols]. Here we want K = nx*ny with bX in 1..nx and bY in 1..ny.
        % Use [nx ny] and index (bX,bY) consistently.
        bSpace(msp) = sub2ind([nx, ny], bX(msp), bY(msp));
        K = nx * ny;

        % ---- split TRACE into temporal bins (frame-aligned) ----
        tw0 = o.TraceWin(1); tw1 = o.TraceWin(2);
        nT  = numel(cs_vec);
        framesByTrial = cell(nT, B, K);

        for tr = 1:nT
            t0 = cs_vec(tr) + tw0;
            t1 = cs_vec(tr) + tw1;

            idx_all = find(ts >= t0 & ts < t1);
            if isempty(idx_all), continue; end

            e = round(linspace(1, numel(idx_all)+1, B+1));
            e(1) = 1; e(end) = numel(idx_all)+1;

            for bti = 1:B
                if e(bti) >= e(bti+1), continue; end
                seg = idx_all(e(bti):e(bti+1)-1);
                if o.UseSpeedMask, seg = seg(speed_ok(seg)); end
                if isempty(seg), continue; end

                bs_b = bSpace(seg);

                for k = 1:K
                    sel = (bs_b == k);
                    if any(sel)
                        framesByTrial{tr, bti, k} = seg(sel);
                    end
                end
            end
        end

        % ---- per-(trial,b,k) PVs ----
        PV = cell(nT, B, K);  % each: [nCells x 1]
        S  = spikes_to_matrix(spikeCell, ts);
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

        % ---- SAME vs DIFF per temporal bin ----
        same_z_mu = nan(B,1); diff_z_mu = nan(B,1);
        day_byBin = repmat(struct('sameZ',[],'diffZ',[]), B, 1);

        for bti = 1:B
            z_same = []; z_diff = [];

            % SAME: (t1,b,k) vs (t2,b,k)
            for k = 1:K
                trs = find(~cellfun(@isempty, squeeze(PV(:,bti,k))));
                if numel(trs) < 2, continue; end
                for i1 = 1:numel(trs)-1
                    t1i = trs(i1); v1 = PV{t1i,bti,k};
                    for i2 = i1+1:numel(trs)
                        t2i = trs(i2); v2 = PV{t2i,bti,k};
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
                        t1i = tr1(i1); v1 = PV{t1i,bti,k1};
                        for i2 = 1:numel(tr2)
                            t2i = tr2(i2); if t2i == t1i, continue; end
                            v2  = PV{t2i,bti,k2};
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

        ATTS.rats(ri).days(di).byBin = day_byBin;
        ATTS.rats(ri).days(di).same  = struct('z',same_z_mu, 'r',tanh(same_z_mu));
        ATTS.rats(ri).days(di).diff  = struct('z',diff_z_mu, 'r',tanh(diff_z_mu));
    end

    % ===================== Δ-based stats per rat =====================
    % Observed per-day Δz = mean_b( mean(z_same_b) - mean(z_diff_b) )
    % Then animal stat = mean_d(Δz_day)

    [day_delta_obs, day_delta_perm_derange, day_deltas_byDay] = build_dayDelta_and_derangeNull(ATTS.rats(ri).days, B, o.NPerm);

    % FLIP null: within-day label swaps across temporal bins (cluster-aware by day)
    ATTS.rats(ri).delta_flip = build_flip_perm_withinDayBins(day_deltas_byDay, o.NPerm);

    % DERANGE null: your existing within-bin label shuffle aggregated to animal
    ATTS.rats(ri).delta_derange = build_derange_perm_from_dayPerms(day_delta_obs, day_delta_perm_derange);

    % pick active delta based on flag
    which = lower(char(o.DeltaPermType));
    switch which
        case 'flip'
            ATTS.rats(ri).delta = ATTS.rats(ri).delta_flip;
        case 'derange'
            ATTS.rats(ri).delta = ATTS.rats(ri).delta_derange;
        otherwise % 'both'
            % keep .delta as flip by default (change if you want)
            ATTS.rats(ri).delta = ATTS.rats(ri).delta_flip;
    end

    fprintf('[%s] Δz_obs=%.3f (r≈%.3f) | flip p_two=%.3g | derange p_two=%.3g\n', ...
        ATTS.rats(ri).name, ATTS.rats(ri).delta.obs_z, tanh(ATTS.rats(ri).delta.obs_z), ...
        ATTS.rats(ri).delta_flip.p_two, ATTS.rats(ri).delta_derange.p_two);
end

% ---- group summaries (mean across all days from all rats) ----
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
group_same_z = mean(Zs,1,'omitnan');
group_diff_z = mean(Zd,1,'omitnan');

ATTS.group = struct( ...
    'same',  struct('z', group_same_z(:),          'r', tanh(group_same_z(:))), ...
    'diff',  struct('z', group_diff_z(:),          'r', tanh(group_diff_z(:))), ...
    'delta', struct('z', group_same_z(:)-group_diff_z(:), ...
                    'r', tanh(group_same_z(:)) - tanh(group_diff_z(:))) ...
);
ATTS.allPairs = struct('same_z',{all_same_z}, 'diff_z',{all_diff_z});

ATTS.distStats = struct();
ATTS.distStats.delta_byAnimal = arrayfun(@(r) r.delta, ATTS.rats, 'uni', 0);

fprintf('[Temporal b-matched] B=%d | mean Δz over bins = %.3f (r≈%.3f)\n', ...
    B, mean(ATTS.group.delta.z,'omitnan'), mean(ATTS.group.delta.r,'omitnan'));

plot_taskTemporal_sameVsDiffSpace_group(ATTS);
plot_deltaPermHist_perRat(ATTS, 'PermType', o.DeltaPermType);

end

% ======================================================================
% ============================ helpers ==================================
% ======================================================================

function [day_delta_obs, day_delta_perm_derange, day_deltas_byDay] = build_dayDelta_and_derangeNull(Ddays, B, P)
% day_delta_obs: [nDays x 1] observed day-level Δz (mean over bins)
% day_delta_perm_derange: {nDays x 1}, each is [P x 1] day-level perm Δz (within-bin label shuffle)
% day_deltas_byDay: {nDays x 1}, each is [B x 1] observed per-bin Δz for that day

day_delta_obs = [];
day_delta_perm_derange = {};
day_deltas_byDay = {};

for di = 1:numel(Ddays)
    Db = Ddays(di).byBin;
    if isempty(Db), continue; end

    % observed per-bin deltas for the day
    deltas_b_obs = nan(B,1);
    for bti = 1:B
        zs = Db(bti).sameZ; zs = zs(isfinite(zs));
        zd = Db(bti).diffZ; zd = zd(isfinite(zd));
        if isempty(zs) || isempty(zd), continue; end
        deltas_b_obs(bti) = mean(zs,'omitnan') - mean(zd,'omitnan');
    end

    if ~any(isfinite(deltas_b_obs)), continue; end

    day_deltas_byDay{end+1,1} = deltas_b_obs; %#ok<AGROW>
    day_delta_obs(end+1,1)    = mean(deltas_b_obs(isfinite(deltas_b_obs)),'omitnan'); %#ok<AGROW>

    % derange-perm: within each bin, pool z and reassign counts
    perm_d = nan(P,1);
    if P >= 1
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
    day_delta_perm_derange{end+1,1} = perm_d; %#ok<AGROW>
end
end

function OUT = build_flip_perm_withinDayBins(day_deltas_byDay, NPerm)
% Within-day label-swap null using per-bin Δz within each day.
% For each permutation:
%   - within each day, flip sign per temporal bin (equiv. swapping SAME/DIFF labels)
%   - average across bins -> day stat
%   - average across days -> animal stat

OUT = struct('obs_z',NaN,'obs_r',NaN,'perm_z',[],'p_left',NaN,'p_right',NaN,'p_two',NaN, ...
             'nBinsByDay',[],'nDays',0);

if isempty(day_deltas_byDay), return; end

% keep valid days with at least 2 finite bins
valid = false(numel(day_deltas_byDay),1);
nBins = zeros(numel(day_deltas_byDay),1);
for d = 1:numel(day_deltas_byDay)
    dv = day_deltas_byDay{d};
    if isempty(dv), continue; end
    dv = dv(:);
    nBins(d) = nnz(isfinite(dv));
    valid(d) = (nBins(d) >= 2);
end
day_deltas_byDay = day_deltas_byDay(valid);
nBins = nBins(valid);

OUT.nDays = numel(day_deltas_byDay);
OUT.nBinsByDay = nBins;

if OUT.nDays < 1, return; end

% observed animal stat = mean over days of (mean over bins of Δz)
day_mu = nan(OUT.nDays,1);
for j = 1:OUT.nDays
    dv = day_deltas_byDay{j};
    day_mu(j) = mean(dv(isfinite(dv)), 'omitnan');
end
obs = mean(day_mu(isfinite(day_mu)), 'omitnan');

OUT.obs_z = obs;
OUT.obs_r = tanh(obs);

if nargin < 2 || isempty(NPerm) || NPerm < 5, return; end

perm = nan(NPerm,1);
for p = 1:NPerm
    day_mu_perm = nan(OUT.nDays,1);
    for j = 1:OUT.nDays
        dv = day_deltas_byDay{j};
        dv = dv(:);
        m  = isfinite(dv);
        if nnz(m) < 2, continue; end
        s = (rand(nnz(m),1) > 0.5)*2 - 1;   % +/-1 per bin (within day)
        day_mu_perm(j) = mean(s .* dv(m), 'omitnan');
    end
    perm(p) = mean(day_mu_perm(isfinite(day_mu_perm)), 'omitnan');
end

pR = mean(perm >= obs);
pL = mean(perm <= obs);
pT = 2*min(pR,pL);

OUT.perm_z  = perm;
OUT.p_right = pR;
OUT.p_left  = pL;
OUT.p_two   = pT;
end

function OUT = build_derange_perm_from_dayPerms(day_delta_obs, day_delta_perm_derange)
% Aggregate the within-bin label-shuffle day nulls to an animal-level null.
OUT = struct('obs_z',NaN,'obs_r',NaN,'perm_z',[],'p_left',NaN,'p_right',NaN,'p_two',NaN);

dd = day_delta_obs(:);
dd = dd(isfinite(dd));

if isempty(dd)
    return
end

obs = mean(dd,'omitnan');
OUT.obs_z = obs;
OUT.obs_r = tanh(obs);

valid = ~cellfun('isempty', day_delta_perm_derange);
if ~any(valid)
    return
end

Pcommon = min(cellfun(@(v) sum(isfinite(v)), day_delta_perm_derange(valid)));
if isempty(Pcommon) || ~isfinite(Pcommon) || Pcommon < 5
    return
end

perm_anim = nan(Pcommon,1);
for p = 1:Pcommon
    vals = nan(sum(valid),1);
    idx = 1;
    for d = 1:numel(day_delta_perm_derange)
        if ~valid(d), continue; end
        v = day_delta_perm_derange{d};
        if numel(v) >= p && isfinite(v(p))
            vals(idx) = v(p); idx = idx + 1;
        end
    end
    perm_anim(p) = mean(vals(isfinite(vals)),'omitnan');
end

pR = mean(perm_anim >= obs);
pL = mean(perm_anim <= obs);
pT = 2*min(pR,pL);

OUT.perm_z  = perm_anim;
OUT.p_right = pR;
OUT.p_left  = pL;
OUT.p_two   = pT;
end

% -------------------- helpers: data coercion --------------------

function [ts,x,y] = coerce_pos(pos)
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
dt = diff(ts); dt(~isfinite(dt) | dt<=0) = NaN;
dx = diff(x);  dy = diff(y);
spd = sqrt(dx.^2 + dy.^2) ./ dt;
v = [spd; spd(end)];
v(~isfinite(v)) = 0;
end

function cs = coerce_cs(csIn)
if istable(csIn)
    vn = lower(string(csIn.Properties.VariableNames));
    cname = pick(vn, ["cs_ms","cs_time","cstime","cs","onset","onsets","cs_onset","cs_time_ms","cue_onset","time","ts"]);
    cs = getcol(csIn, cname);
else
    cs = double(csIn(:));
end
if max(cs,[],'omitnan') > 1e4, cs = cs/1000; end
end

function cellSpk = to_cell_spikes(spk)
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

% -------------------- plotting --------------------

function plot_taskTemporal_sameVsDiffSpace_group(ATTS, varargin)
p = inputParser;
addParameter(p,'Bin','all-mean');
parse(p,varargin{:});
Bin = p.Results.Bin;

Zw_day = [];
Zd_day = [];
aid    = [];

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

rat_mean_same_z = nan(nR,1);
rat_mean_diff_z = nan(nR,1);
for ri = 1:nR
    m = (aid==ri);
    if any(m)
        rat_mean_same_z(ri) = mean(Zw_day(m), 'omitnan');
        rat_mean_diff_z(ri) = mean(Zd_day(m), 'omitnan');
    end
end
keepA = isfinite(rat_mean_same_z) & isfinite(rat_mean_diff_z);
Zs_rats = rat_mean_same_z(keepA);
Zd_rats = rat_mean_diff_z(keepA);

[pA, tA, dfA, dzA] = paired_t_z(Zs_rats, Zd_rats);

keepD = isfinite(Zw_day) & isfinite(Zd_day);
[pD, tD, dfD, dzD] = paired_t_z(Zw_day(keepD), Zd_day(keepD));

figure('Color','w','Position',[160 160 900 540]); hold on
cmap = lines(nR);
make_light = @(c,frac) (1-frac)*c + frac*[1 1 1];

rw_day = tanh(Zw_day);
rd_day = tanh(Zd_day);
for ri = 1:nR
    c_base  = cmap(ri,:);
    c_light = make_light(c_base, 0.60);
    idx = find(aid==ri & isfinite(rw_day) & isfinite(rd_day));
    for j = 1:numel(idx)
        plot([1 2], [rw_day(idx(j)) rd_day(idx(j))], '-', 'Color', c_light, 'LineWidth', 1.0);
        plot(1, rw_day(idx(j)), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
        plot(2, rd_day(idx(j)), 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
    end
end

rat_same_r = tanh(rat_mean_same_z);
rat_diff_r = tanh(rat_mean_diff_z);
for ri = 1:nR
    if isfinite(rat_same_r(ri)) && isfinite(rat_diff_r(ri))
        plot([1 2], [rat_same_r(ri) rat_diff_r(ri)], '-', 'Color', cmap(ri,:), 'LineWidth', 2.5);
        plot(1, rat_same_r(ri), 'o', 'MarkerFaceColor', cmap(ri,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
        plot(2, rat_diff_r(ri), 'o', 'MarkerFaceColor', cmap(ri,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
    end
end

bar(1, tanh(mean(Zs_rats,'omitnan')), 0.6, 'FaceColor',[0.30 0.60 1.00], 'EdgeColor','k');
bar(2, tanh(mean(Zd_rats,'omitnan')), 0.6, 'FaceColor',[0.85 0.40 0.20], 'EdgeColor','k');

xlim([0.5 2.5]); xticks([1 2]);
xticklabels({'SAME space','DIFF space'});
ylabel('PV correlation (r)');
title('Task (temporal b-matched): SAME vs DIFF space');
yline(0,'k:'); grid on; box on

yl = ylim; y = yl(2) + 0.06*range(yl);
line([1 2], [y y], 'Color','k','LineWidth',1.2);
txt1 = sprintf('per-rat paired t (z): t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d', dfA, tA, pA, dzA, numel(Zs_rats));
text(1.5, y + 0.02*range(yl), txt1, 'HorizontalAlignment','center');
ylim([yl(1) y + 0.10*range(yl)]);

fprintf('\n== SAME vs DIFF (b-matched) ==\n');
fprintf('Primary (per-rat):  paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nRats=%d\n', dfA, tA, pA, dzA, numel(Zs_rats));
fprintf('Secondary (per-day): paired t in z: t(%d)=%.2f, p=%.3g, dz=%.2f, nDays=%d\n', dfD, tD, pD, dzD, nnz(keepD));

function v = safe_pick(vec, k)
v = NaN;
if ~isvector(vec) || k<1 || k>numel(vec), return; end
v = vec(k);
end
end

function [p, t, df, dz] = paired_t_z(z1, z2)
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

function plot_deltaPermHist_perRat(ATTS, varargin)
% Single figure (3x2) with one subplot per rat.
% PermType: 'flip' | 'derange' | 'both'

p = inputParser;
addParameter(p,'NBins',30, @(x) isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'ShowTwoSided',true, @(x) islogical(x)&&isscalar(x));
addParameter(p,'PermType','flip', @(s) any(strcmpi(s,{'flip','derange','both'})));
parse(p,varargin{:});
NBins = p.Results.NBins;
ShowTwoSided = p.Results.ShowTwoSided;
permType = lower(char(p.Results.PermType));

nR = numel(ATTS.rats);

toPlot = {};
if strcmpi(permType,'both')
    toPlot = {'flip','derange'};
else
    toPlot = {permType};
end

for jj = 1:numel(toPlot)
    pt = toPlot{jj};

    figure('Color','w','Position',[200 200 1100 700]);
    for ri = 1:nR
        r = ATTS.rats(ri);
        if ~isfield(r,'name') || isempty(r.name)
            rname = sprintf('rat%02d',ri);
        else
            rname = r.name;
        end

        if strcmpi(pt,'flip')
            dlt = r.delta_flip;
        else
            dlt = r.delta_derange;
        end

        if isempty(dlt) || ~isfield(dlt,'obs_z') || ~isfield(dlt,'perm_z')
            continue
        end

        obs = dlt.obs_z;
        perm = dlt.perm_z;

        if ~isfinite(obs) || isempty(perm), continue; end
        perm = perm(isfinite(perm));
        if isempty(perm), continue; end

        subplot(3,2,ri); hold on
        histogram(perm, NBins, 'Normalization','pdf', 'EdgeColor','none');
        yl = ylim;
        plot([obs obs], yl, 'k-', 'LineWidth', 2);

        xlabel('\Delta z = (SAME - DIFF) Fisher z');
        ylabel('Null density');
        title(sprintf('%s: %s', rname, upper(pt)));
        grid on; box on

        % p-values
        pR = NaN; pL = NaN; pT = NaN;
        if isfield(dlt,'p_right'), pR = dlt.p_right; end
        if isfield(dlt,'p_left'),  pL = dlt.p_left;  end
        if isfield(dlt,'p_two'),   pT = dlt.p_two;   end

        txt = sprintf('obs \\Delta z = %.3f (r\\approx%.3f)\n', obs, tanh(obs));
        if ShowTwoSided
            if isfinite(pT), txt = [txt, sprintf('p_{two}=%.3g', pT)];
            else
                pR_ = mean(perm >= obs); pL_ = mean(perm <= obs);
                txt = [txt, sprintf('p_{two}\\approx%.3g', 2*min(pR_,pL_))];
            end
        else
            if isfinite(pR), txt = [txt, sprintf('p_{right}=%.3g', pR)];
            else, txt = [txt, sprintf('p_{right}\\approx%.3g', mean(perm>=obs))];
            end
        end

        text(obs, yl(2), ['  ' txt], 'VerticalAlignment','top', 'HorizontalAlignment','left');
        ylim(yl);
    end

    sgtitle(sprintf('Space→Task Δz permutation null (%s)', upper(pt)));
end
end
