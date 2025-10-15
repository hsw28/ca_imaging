function ATTS = compute_task_temporal_sameVsDiffSpace(ratNames, varargin)
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
%   .group.same.z(b) / diff.z(b) / delta.z(b)
%   .group.same.r(b) / diff.r(b) / delta.r(b)
%   .params

% ------------ options -------------
p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'TimeBins',15);
addParameter(p,'GridN',[3 2]);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'VelThresh',4);
addParameter(p,'CellNorm','demean');
addParameter(p,'MinFrames',2);
addParameter(p,'MinNCorr',0);
addParameter(p,'MinStdBin',0);
addParameter(p,'MinTrials',0);
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
ATTS.rats   = repmat(struct('name','','days',[]), nR, 1);

for ri = 1:nR
    rname = ratNames{ri};
    RAT   = evalin('base', rname);              % pull rat struct from base
    ATTS.rats(ri).name = rname;

    % pick 3 days around rat.An
    dates = autoDateList(RAT);
    idx   = find(strcmp(dates, RAT.An), 1);
    if isempty(idx) || idx < 3
        dayIdx = max(1, idx-2) : idx;
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
        rm   = RAT.ratemask.(sprintf('ratemask_%s', dlabel));

        keepCells = (rm == 1);
        spk = spk(keepCells, :);   % keep only accepted cells

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


        % ---- split TRACE into temporal bins ----
        % ---- split TRACE into temporal bins aligned to frames (7.5 Hz → 15 bins over 2s) ----
        tw0 = o.TraceWin(1); tw1 = o.TraceWin(2);
        nT  = numel(cs_vec);
        framesByTrial = cell(nT, B, K);

        % (use dt measured on ts; we will bin by frame index anyway)
        % NOTE: we bin using absolute frame indices between [t0,t1),
        %       then (optionally) intersect with the speed mask inside each bin.

        for tr = 1:nT
            t0 = cs_vec(tr) + tw0;
            t1 = cs_vec(tr) + tw1;

            % all frames in the TRACE window (no speed filter yet)
            idx_all = find(ts >= t0 & ts < t1);
            if isempty(idx_all), continue; end

            % split the frame indices into B consecutive blocks, exactly on frame edges
            % edges are in index-space so they align perfectly with the sampled frames
            e = round(linspace(1, numel(idx_all)+1, B+1));  % 1..(len+1)
            e(1) = 1; e(end) = numel(idx_all)+1;            % guard exact endpoints

            for bti = 1:B
                if e(bti) >= e(bti+1), continue; end
                seg = idx_all(e(bti):e(bti+1)-1);           % frames for temporal bin bti

                % now apply the speed mask *within* the temporal bin (if enabled)
                if o.UseSpeedMask
                    seg = seg(speed_ok(seg));
                end
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
if nargin < 3 || isempty(minN), minN = 3; end
if nargin < 4 || isempty(minStd), minStd = 1e-12; end
m = isfinite(a) & isfinite(b);
nC = nnz(m);
if nC < minN
    r = NaN; return
end
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
% Plot SAME-space vs DIFF-space (task, b-matched) in your "per-day lines,
% per-rat mean" style. Works on ATTS (not ATTS.group).
%
% Usage:
%   plot_taskTemporal_sameVsDiffSpace_group(ATTS)             % mean over all bins
%   plot_taskTemporal_sameVsDiffSpace_group(ATTS,'Bin',3)     % temporal bin #3
%
% Options:
%   'Bin'  : integer bin index, or 'all-mean' (default) to average across bins (Fisher-z then tanh).

p = inputParser;
addParameter(p,'Bin','all-mean');  % 'all-mean' or a scalar bin index
parse(p,varargin{:});
Bin = p.Results.Bin;

% ---------- harvest per-day r-values ----------
nR = numel(ATTS.rats);
cmap = lines(nR);
make_light = @(c,frac) (1-frac)*c + frac*[1 1 1];

perday_same = [];     % r per day (row)
perday_diff = [];
perday_rat  = [];     % rat index per day

figure('Color','w','Position',[160 160 860 520]); hold on
% legend bars (for colors)
bar(1, 0, 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k');
bar(2, 0, 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k');
xticks([1 2]); xticklabels({'SAME space','DIFF space'});
ylabel('Correlation r'); box on
title('Task (temporal b-matched): SAME vs DIFF space');

rat_means_same = [];
rat_means_diff = [];

for ri = 1:nR
    if isempty(ATTS.rats(ri).days), continue; end
    c_base  = cmap(ri,:);
    c_light = make_light(c_base, 0.60);

    w_days = []; a_days = [];

    for di = 1:numel(ATTS.rats(ri).days)
        D = ATTS.rats(ri).days(di);
        if isempty(D) || ~isfield(D,'same') || ~isfield(D,'diff'), continue; end

        % pull r for requested bin (or mean over bins via Fisher z)
        if isnumeric(Bin)
            r_same = getfieldsafe(D,'same','r',Bin);
            r_diff = getfieldsafe(D,'diff','r',Bin);
        else
            z_s = D.same.z; z_d = D.diff.z;
            r_same = tanh(mean(z_s(isfinite(z_s)),'omitnan'));
            r_diff = tanh(mean(z_d(isfinite(z_d)),'omitnan'));
        end

        if isfinite(r_same) && isfinite(r_diff)
            % per-day line
            plot([1 2], [r_same r_diff], '-', 'Color', c_light, 'LineWidth', 1);
            plot(1, r_same, 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
            plot(2, r_diff, 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);

            w_days(end+1) = r_same; %#ok<AGROW>
            a_days(end+1) = r_diff; %#ok<AGROW>

            perday_same(end+1,1) = r_same; %#ok<AGROW>
            perday_diff(end+1,1) = r_diff; %#ok<AGROW>
            perday_rat(end+1,1)  = ri;     %#ok<AGROW>
        end
    end

    % per-rat mean (dark)
    if ~isempty(w_days) && ~isempty(a_days)
        w_mu = mean(w_days,'omitnan');
        a_mu = mean(a_days,'omitnan');
        plot([1 2], [w_mu a_mu], '-', 'Color', c_base, 'LineWidth', 2.5);
        plot(1, w_mu, 'o', 'MarkerFaceColor', c_base, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 6);
        plot(2, a_mu, 'o', 'MarkerFaceColor', c_base, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 6);
        rat_means_same(end+1) = w_mu; %#ok<AGROW>
        rat_means_diff(end+1) = a_mu; %#ok<AGROW>
    end
end

if isempty(perday_same)
    warning('No per-day same/diff data found to plot.'); return
end

% overlay group means (mean of per-rat means, like your other plots)
if ~isempty(rat_means_same), bar(1, mean(rat_means_same,'omitnan'), 0.6, 'FaceColor',[0.3 0.6 1], 'EdgeColor','k'); end
if ~isempty(rat_means_diff), bar(2, mean(rat_means_diff,'omitnan'), 0.6, 'FaceColor',[0.85 0.4 0.2], 'EdgeColor','k'); end
yline(0,'k:'); ylim auto

% ---------- paired t-test across all days ----------
m = isfinite(perday_same) & isfinite(perday_diff);
s_all = perday_same(m); d_all = perday_diff(m);
[~, p, ~, st] = ttest(s_all, d_all);
dz = mean(atanh(s_all) - atanh(d_all)) / std(atanh(s_all) - atanh(d_all)); % Cohen’s dz on z

% annotate
yl = ylim; y = yl(2) + 0.06*range(yl);
line([1 2], [y y], 'Color','k','LineWidth',1.2);
txt = sprintf('paired t across days: t(%d)=%.2f, p=%.3g, dz=%.2f, n=%d', st.df, st.tstat, p, dz, numel(s_all));
text(1.5, y + 0.02*range(yl), txt, 'HorizontalAlignment','center');
ylim([yl(1) y + 0.10*range(yl)]);

% helper to safely pull field/element
function v = getfieldsafe(S, f1, f2, idx)
    v = NaN;
    if ~isfield(S,f1) || ~isfield(S.(f1),f2), return; end
    V = S.(f1).(f2);
    if isempty(V) || ~isvector(V) || idx<1 || idx>numel(V), return; end
    v = V(idx);
end
end
