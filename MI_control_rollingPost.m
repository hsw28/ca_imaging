function rolling_struct = MI_control_rollingPost(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, CS_times_struct, nIter, winEdges, ControlMatch, NSpeedBins, SpaceBinSize)
% MI_control_rollingPost  (CS-aligned rolling removal/control; negatives allowed)
%
% Test  : remove spikes ONLY in the CURRENT CS-aligned window from the movement MI pool.
% Control: remove the SAME number of spikes from a WINDOW-LOCAL pool that EXCLUDES:
%           (1) trial time: [CS, CS+2]  (hard-coded here as TrialWinSecs = [0 2])
%           (2) and the current window ± guard.
%          with optional matching of SPEED and/or SPACE.
%
% ControlMatch  : 'none' | 'speed' | 'space' | 'speedspace'|'speed_space'|'both'
% NSpeedBins    : integer (# bins for speed matching; default 5)
% SpaceBinSize  : spatial bin size for space matching (units of pos; default uses 'dim')
%
% Output:
%   rolling_struct.(sprintf('MI_%s', date)) with fields:
%     .winEdges, .MI_win (W×Nc), .MI_rand (W×Nc), .MI_base (1×Nc), .dMI
%     .counts: struct of spike-count diagnostics (see below)

% ---- defaults ----
if nargin < 7 || isempty(nIter),       nIter   = 100; end
if nargin < 8 || isempty(winEdges),    winEdges = [0 2; 1 3; 2 4; 3 5; 4 6; 5 7; 8 10; 9 11; 10 12; 11 13]; end
if nargin < 9 || isempty(ControlMatch), ControlMatch = 'none'; end
if nargin < 10 || isempty(NSpeedBins),  NSpeedBins = 5; end
if nargin < 11 || isempty(SpaceBinSize),SpaceBinSize = dim; end

% ---- NEW: local control pool definition ----
TrialWinSecs = [0 2];     % exclude [CS, CS+2] from control sampling (your request)
GuardSecs    = 0.5;      % exclude current window ± guard from control sampling

% normalize ControlMatch string
cm = lower(strrep(ControlMatch,'_',''));
if strcmp(cm,'both'), cm = 'speedspace'; end

set(0,'DefaultFigureVisible','off');

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_cats   = fieldnames(CA_timestamps);
fields_CS     = fieldnames(CS_times_struct);

% last 3 days
fields_spikes = fields_spikes(max(1,end-2):end);
fields_pos    = fields_pos(   max(1,end-2):end);
fields_cats   = fields_cats(  max(1,end-2):end);
fields_CS     = fields_CS(    max(1,end-2):end);

rolling_struct = struct();

for iDay = 1:numel(fields_spikes)
    iDay
    fieldName_sp  = fields_spikes{iDay};
    fieldName_pos = fields_pos{iDay};
    fieldName_ts  = fields_cats{iDay}; %#ok<NASGU>
    fieldName_CS  = fields_CS{iDay};

    peaks_time = spike_structure.(fieldName_sp);   % Ncells x NspikeTimes
    pos        = pos_structure.(fieldName_pos);    % [t x y]
    cs_times   = CS_times_struct.(fieldName_CS);   % vector of CS onset times (sec)

    spikes_date = extract_date_suffix(fieldName_sp);

    % ---- velocity & interpolation ----
    [~, uniqPos] = unique(pos(:,1), 'stable'); pos = pos(uniqPos,:);
    vel      = ca_velocity(pos);         % [speed; time]
    vel_time = vel(2,:)'; vel_mag = vel(1,:)';

    % interpolate pos onto vel_time
    ix = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    iy = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    if isempty(cs_times)
        warning('No CS times for %s; skipping day.', fieldName_CS);
        out = struct('winEdges',winEdges,'MI_win',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_rand',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_base',nan(1,size(peaks_time,1)),'dMI',nan(size(winEdges,1),size(peaks_time,1)));
        out.counts = struct('N_move',nan(1,size(peaks_time,1)), 'N_inWin',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_win',nan(size(winEdges,1),size(peaks_time,1)), 'N_ctrl_idx',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_ctrl',nan(size(winEdges,1),size(peaks_time,1)), 'Iter_OK',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_removed_ctrl',nan(size(winEdges,1),size(peaks_time,1)));
        rolling_struct.(sprintf('MI_%s', spikes_date)) = out; continue;
    end

    % movement occupancy samples
    validIdx = (vel_mag >= velthreshold) & ~isnan(ix) & ~isnan(iy);
    posDat   = [vel_time(validIdx), ix(validIdx), iy(validIdx)];

    nCells = size(peaks_time,1);
    nWins  = size(winEdges,1);

    MI_base = nan(1, nCells);
    MI_win  = nan(nWins, nCells);
    MI_rand = nan(nWins, nCells);

    % ---- spike-count diagnostics ----
    N_inWin        = nan(nWins, nCells);   % spikes removed by TEST (nRemove)
    N_keep_win     = nan(nWins, nCells);   % spikes kept for MI_win
    N_ctrl_idx     = nan(nWins, nCells);   % control candidate pool size (pre-matching)
    N_keep_ctrl    = nan(nWins, nCells);   % spikes kept for MI_rand (should equal move-nRemove when success)
    N_removed_ctrl = nan(nWins, nCells);   % spikes removed by CONTROL (should equal nRemove when success)
    N_move         = nan(1,    nCells);    % movement spikes available before any removal
    Iter_OK        = nan(nWins, nCells);   % # successful control iterations used (for MI_rand)

    % ---------- MI_base over movement ----------
    for k = 1:nCells
        spk = peaks_time(k,:); spk = spk(~isnan(spk) & spk > 0);
        if isempty(spk), continue; end
        spike_vel = interp1(vel_time, vel_mag, spk, 'linear', NaN);
        move_spikes = spk(spike_vel >= velthreshold);
        if isempty(move_spikes) || size(posDat,1) < 10, continue; end

        [~,~,~,~,spikeprob,occprob] = CA_normalizePosData(move_spikes, posDat, dim, 1.0);
        if size(spikeprob,1) < size(spikeprob,2), spikeprob = spikeprob'; end
        if size(occprob,1)  < size(occprob,2),  occprob  = occprob';  end
        MI_base(k) = mutualinfo([spikeprob, occprob]);
    end

    % ---------- Rolling windows ----------
    parfor k = 1:nCells
        spk = peaks_time(k,:); spk = spk(~isnan(spk) & spk > 0);

        MI_win_k  = nan(nWins,1);
        MI_rand_k = nan(nWins,1);

        % local copies for diagnostics (parfor-safe)
        N_inWin_k        = nan(nWins,1);
        N_keep_win_k     = nan(nWins,1);
        N_ctrl_idx_k     = nan(nWins,1);
        N_keep_ctrl_k    = nan(nWins,1);
        N_removed_ctrl_k = nan(nWins,1);
        Iter_OK_k        = nan(nWins,1);
        N_move_k         = NaN;

        if isempty(spk) || isnan(MI_base(k)) || isempty(posDat) || size(posDat,1) < 10
            MI_win(:,k)  = MI_win_k;
            MI_rand(:,k) = MI_rand_k;
            N_inWin(:,k)        = N_inWin_k;
            N_keep_win(:,k)     = N_keep_win_k;
            N_ctrl_idx(:,k)     = N_ctrl_idx_k;
            N_keep_ctrl(:,k)    = N_keep_ctrl_k;
            N_removed_ctrl(:,k) = N_removed_ctrl_k;
            Iter_OK(:,k)        = Iter_OK_k;
            N_move(1,k)         = N_move_k;
            continue;
        end

        spike_vel   = interp1(vel_time, vel_mag, spk, 'linear', NaN);
        move_spikes = spk(spike_vel >= velthreshold);
        N_move_k    = numel(move_spikes);
        if numel(move_spikes) < 2
            MI_win(:,k)  = MI_win_k;
            MI_rand(:,k) = MI_rand_k;
            N_inWin(:,k)        = N_inWin_k;
            N_keep_win(:,k)     = N_keep_win_k;
            N_ctrl_idx(:,k)     = N_ctrl_idx_k;
            N_keep_ctrl(:,k)    = N_keep_ctrl_k;
            N_removed_ctrl(:,k) = N_removed_ctrl_k;
            Iter_OK(:,k)        = Iter_OK_k;
            N_move(1,k)         = N_move_k;
            continue;
        end

        % relative times
        nearest_cs_for_spike   = interp1(cs_times, cs_times, spk, 'nearest', NaN);
        rel_t_all              = spk - nearest_cs_for_spike;

        nearest_cs_for_move    = interp1(cs_times, cs_times, move_spikes, 'nearest', NaN);
        rel_t_move             = move_spikes - nearest_cs_for_move;

        % attributes for matching
        move_speed = interp1(vel_time, vel_mag, move_spikes, 'linear', NaN);
        move_x = interp1(pos(:,1), pos(:,2), move_spikes, 'linear', NaN);
        move_y = interp1(pos(:,1), pos(:,3), move_spikes, 'linear', NaN);

        % spatial bin edges for matching (based on movement cloud)
        if any(strcmp(cm, {'space','speedspace'}))
            minX = nanmin(move_x); maxX = nanmax(move_x);
            minY = nanmin(move_y); maxY = nanmax(move_y);
            if ~isfinite(minX) || ~isfinite(minY) || ~isfinite(maxX) || ~isfinite(maxY)
                edgesX = []; edgesY = [];
            else
                edgesX = make_lin_edges(minX, maxX, SpaceBinSize);
                edgesY = make_lin_edges(minY, maxY, SpaceBinSize);
            end
        else
            edgesX = []; edgesY = [];
        end

        for w = 1:nWins
            wstart = winEdges(w,1);
            wend   = winEdges(w,2);

            % ---- define the test window ON THE MOVEMENT SPIKES ----
            w_mask_on_move = (rel_t_move > wstart) & (rel_t_move <= wend);   % logical idx over move_spikes
            w_idx          = find(w_mask_on_move);                            % indices into move_spikes for this window
            nRemove        = numel(w_idx);

            % record counts early
            N_inWin_k(w)    = nRemove;

            if nRemove < 2
                N_keep_win_k(w)     = NaN;
                % control pool size for visibility (local ctrl definition)
                in_trial = (rel_t_move >= TrialWinSecs(1)) & (rel_t_move <= TrialWinSecs(2));
                in_win_guard = (rel_t_move > (wstart-GuardSecs)) & (rel_t_move <= (wend+GuardSecs));
                ctrl_mask_tmp = ~in_trial & ~in_win_guard;
                N_ctrl_idx_k(w) = sum(ctrl_mask_tmp);
                N_keep_ctrl_k(w)    = NaN;
                N_removed_ctrl_k(w) = NaN;
                Iter_OK_k(w)        = 0;
                MI_win_k(w)         = NaN;
                MI_rand_k(w)        = NaN;
                continue;
            end

            % ---------------- TEST (deterministic): remove those indices from move_spikes ----------------
            keep_mask_win = true(size(move_spikes));
            keep_mask_win(w_idx) = false;                      % drop all window spikes
            keep_win = move_spikes(keep_mask_win);

            % counts for test
            N_keep_win_k(w) = numel(keep_win);

            % compute MI_win if enough remain
            if numel(keep_win) >= 2
                [~,~,~,~,spikeprobW,occprobW] = CA_normalizePosData(keep_win, posDat, dim, 1.0);
                if size(spikeprobW,1) < size(spikeprobW,2), spikeprobW = spikeprobW'; end
                if size(occprobW,1)  < size(occprobW,2),  occprobW  = occprobW';  end
                MI_win_k(w) = mutualinfo([spikeprobW, occprobW]);
            else
                MI_win_k(w) = NaN;
            end

            % ---------------- CONTROL (UPDATED): window-local, trial-excluded ----------------
            % Exclude:
            %   - trial time [0,2]
            %   - current window (± GuardSecs) so control isn’t drawn from the same time band
            in_trial = (rel_t_move >= TrialWinSecs(1)) & (rel_t_move <= TrialWinSecs(2));
            in_win_guard = (rel_t_move > (wstart-GuardSecs)) & (rel_t_move <= (wend+GuardSecs));

            ctrl_mask = ~in_trial & ~in_win_guard;

            ctrl_idx  = find(ctrl_mask);                        % candidate indices (into move_spikes)
            N_ctrl_idx_k(w) = numel(ctrl_idx);

            % prepare window attributes for matching options
            w_speed = move_speed(w_idx);
            w_x     = move_x(w_idx);
            w_y     = move_y(w_idx);

            mi_accum = 0;
            n_ok = 0;

            for it = 1:nIter
                switch cm
                    case 'none'
                        if numel(ctrl_idx) < nRemove
                            sample_idx = [];
                        else
                            sample_idx = randsample(ctrl_idx, nRemove, false);
                        end

                    case 'speed'
                        ctrl_speed = move_speed(ctrl_idx);
                        edges = make_quant_edges([w_speed(:); ctrl_speed(:)], NSpeedBins);
                        wbins    = discretize(w_speed,    edges);
                        ctrlbins = discretize(ctrl_speed, edges);
                        sample_idx = stratified_sample(ctrl_idx, ctrlbins, wbins, nRemove);

                    case 'space'
                        if isempty(edgesX) || isempty(edgesY)
                            sample_idx = (numel(ctrl_idx) >= nRemove) * randsample(ctrl_idx, min(nRemove,numel(ctrl_idx)), false);
                            if isscalar(sample_idx) && sample_idx==0, sample_idx = []; end
                        else
                            [wBX,wBY]  = deal(discretize(w_x, edgesX), discretize(w_y, edgesY));
                            [cBX,cBY]  = deal(discretize(move_x(ctrl_idx), edgesX), discretize(move_y(ctrl_idx), edgesY));
                            w_joint    = wBX.*1e6 + wBY;
                            ctrl_joint = cBX.*1e6 + cBY;
                            sample_idx = stratified_sample(ctrl_idx, ctrl_joint, w_joint, nRemove);
                        end

                    case 'speedspace'
                        ctrl_speed = move_speed(ctrl_idx);
                        edges = make_quant_edges([w_speed(:); ctrl_speed(:)], NSpeedBins);
                        wbins    = discretize(w_speed,    edges);
                        ctrlbins = discretize(ctrl_speed, edges);
                        if isempty(edgesX) || isempty(edgesY)
                            sample_idx = stratified_sample(ctrl_idx, ctrlbins, wbins, nRemove);   % speed-only fallback
                        else
                            [wBX,wBY]  = deal(discretize(w_x, edgesX), discretize(w_y, edgesY));
                            [cBX,cBY]  = deal(discretize(move_x(ctrl_idx), edgesX), discretize(move_y(ctrl_idx), edgesY));
                            wjoint     = wbins.*1e8 + wBX.*1e4 + wBY;           % joint (speed, xBin, yBin)
                            cjoint     = ctrlbins.*1e8 + cBX.*1e4 + cBY;
                            sample_idx = stratified_sample(ctrl_idx, cjoint, wjoint, nRemove);
                        end

                    otherwise
                        if numel(ctrl_idx) < nRemove
                            sample_idx = [];
                        else
                            sample_idx = randsample(ctrl_idx, nRemove, false);
                        end
                end

                % must remove EXACTLY nRemove; otherwise skip this iteration
                if isempty(sample_idx) || numel(sample_idx) ~= nRemove, continue; end

                keep_rd = move_spikes;
                keep_rd(sample_idx) = [];         % remove by index (no set/value issues)

                if numel(keep_rd) < 2, continue; end

                [~,~,~,~,spikeprobR,occprobR] = CA_normalizePosData(keep_rd, posDat, dim, 1.0);
                if size(spikeprobR,1) < size(spikeprobR,2), spikeprobR = spikeprobR'; end
                if size(occprobR,1)  < size(occprobR,2),  occprobR  = occprobR';  end
                mi_accum = mi_accum + mutualinfo([spikeprobR, occprobR]);
                n_ok = n_ok + 1;
            end

            Iter_OK_k(w) = n_ok;
            if n_ok > 0
                MI_rand_k(w)        = mi_accum / n_ok;
                N_removed_ctrl_k(w) = nRemove;
                N_keep_ctrl_k(w)    = numel(move_spikes) - nRemove;   % identical to test by construction
            else
                MI_rand_k(w)        = NaN;
                N_removed_ctrl_k(w) = NaN;
                N_keep_ctrl_k(w)    = NaN;
            end
        end

        MI_win(:,k)  = MI_win_k;
        MI_rand(:,k) = MI_rand_k;

        N_inWin(:,k)        = N_inWin_k;
        N_keep_win(:,k)     = N_keep_win_k;
        N_ctrl_idx(:,k)     = N_ctrl_idx_k;
        N_keep_ctrl(:,k)    = N_keep_ctrl_k;
        N_removed_ctrl(:,k) = N_removed_ctrl_k;
        Iter_OK(:,k)        = Iter_OK_k;
        N_move(1,k)         = N_move_k;
    end

    out.winEdges = winEdges;
    out.MI_win   = MI_win;
    out.MI_rand  = MI_rand;
    out.MI_base  = MI_base;
    out.dMI      = MI_win - MI_rand;

    % spike-count diagnostics attached
    out.counts = struct( ...
        'N_move',         N_move, ...
        'N_inWin',        N_inWin, ...
        'N_keep_win',     N_keep_win, ...
        'N_ctrl_idx',     N_ctrl_idx, ...
        'N_keep_ctrl',    N_keep_ctrl, ...
        'N_removed_ctrl', N_removed_ctrl, ...
        'Iter_OK',        Iter_OK);

    rolling_struct.(sprintf('MI_%s', spikes_date)) = out;
end
end


% ------------------------------- helpers -------------------------------

function edges = make_lin_edges(a, b, step)
if ~isfinite(a) || ~isfinite(b) || ~(b>a) || step<=0
    edges = [];
    return;
end
start = floor(a/step)*step;  stop = ceil(b/step)*step;
edges = start:step:stop;
if numel(edges)<2, edges = [a b]; end
end

function edges = make_quant_edges(x, NB)
% quantile-ish edges with fallback to linspace if duplicates
x = x(isfinite(x));
if isempty(x) || NB<1
    edges = [];
    return;
end
qs = linspace(0,1,NB+1);
edges = quantile(x, qs);
edges = unique(edges);
if numel(edges) < 2
    mn = min(x); mx = max(x);
    edges = linspace(mn, mx, NB+1);
elseif numel(edges) <= NB
    mn = min(x); mx = max(x);
    edges = linspace(mn, mx, NB+1);
end
end

function sample_idx = stratified_sample(ctrl_idx, ctrl_bins, w_bins, nRemove)
% Nearest-bin stratified sampling (for ordered integer bins like speed 1..K).
% Returns EXACTLY nRemove indices from ctrl_idx whenever possible.
% If impossible (not enough total candidates), returns [].

ctrl_idx  = ctrl_idx(:);
ctrl_bins = ctrl_bins(:);
w_bins    = w_bins(:);

sample_idx = zeros(0,1);
if isempty(ctrl_idx) || nRemove <= 0
    sample_idx = [];
    return;
end

all_bins = [ctrl_bins; w_bins];
all_bins = all_bins(isfinite(all_bins));
if isempty(all_bins)
    if numel(ctrl_idx) >= nRemove
        sample_idx = randsample(ctrl_idx, nRemove, false); sample_idx = sample_idx(:);
    else
        sample_idx = [];
    end
    return;
end

ub = unique(all_bins);
map = containers.Map(num2cell(ub.'), num2cell(1:numel(ub)));
remap = @(v) cell2mat(values(map, num2cell(v(isfinite(v)).')));
ctrl_isfinite = isfinite(ctrl_bins);
w_isfinite    = isfinite(w_bins);

ctrl_b = nan(size(ctrl_bins));
w_b    = nan(size(w_bins));
if any(ctrl_isfinite)
    ctrl_b(ctrl_isfinite) = remap(ctrl_bins(ctrl_isfinite));
end
if any(w_isfinite)
    w_b(w_isfinite) = remap(w_bins(w_isfinite));
end
K = numel(ub);

hist_w = accumarray(w_b(~isnan(w_b)), 1, [K,1], @sum, 0);
if sum(hist_w) == 0
    if numel(ctrl_idx) < nRemove, sample_idx = []; return; end
    sample_idx = randsample(ctrl_idx, nRemove, false); sample_idx = sample_idx(:);
    return;
end
raw   = nRemove * (hist_w / sum(hist_w));
tgt   = floor(raw);
left  = nRemove - sum(tgt);
if left > 0
    frac = raw - tgt;
    [~, ord] = sort(frac, 'descend');
    tgt(ord(1:left)) = tgt(ord(1:left)) + 1;
end

cands = cell(K,1);
avail = zeros(K,1);
for b = 1:K
    cands{b} = find(ctrl_b == b);
    avail(b) = numel(cands{b});
end
if sum(avail) < nRemove
    sample_idx = []; return;
end

taken  = zeros(K,1);
chosen = cell(K,1);
for b = 1:K
    want = tgt(b);
    if want <= 0 || avail(b) == 0, continue; end
    take = min(want, avail(b));
    if take > 0
        pick = randsample(cands{b}, take, false);
        chosen{b} = pick(:);
        cands{b} = setdiff(cands{b}, pick, 'stable');
        avail(b) = numel(cands{b});
        taken(b) = take;
    end
end

deficit = tgt - taken;

maxRadius = K-1;
for radius = 1:maxRadius
    todo = find(deficit > 0);
    if isempty(todo), break; end
    for i = 1:numel(todo)
        b = todo(i);
        need = deficit(b);
        if need <= 0, continue; end

        neigh = [];
        L = b - radius; if L >= 1, neigh = [neigh, L]; end %#ok<AGROW>
        R = b + radius; if R <= K, neigh = [neigh, R]; end %#ok<AGROW>
        if isempty(neigh), continue; end

        neigh_av = sum(avail(neigh));
        if neigh_av <= 0, continue; end

        take_tot = min(need, neigh_av);
        take_each = zeros(size(neigh));
        if numel(neigh) == 1
            take_each(1) = take_tot;
        else
            prop = avail(neigh) / neigh_av;
            take_each = floor(take_tot * prop);
            rem = take_tot - sum(take_each);
            if rem > 0
                [~, nnord] = sort(prop, 'descend');
                take_each(nnord(1:rem)) = take_each(nnord(1:rem)) + 1;
            end
        end

        for j = 1:numel(neigh)
            nb = neigh(j);
            if take_each(j) <= 0 || avail(nb) <= 0, continue; end
            t = min(take_each(j), avail(nb));
            if t <= 0, continue; end
            pick = randsample(cands{nb}, t, false);
            chosen{b} = [chosen{b}; pick(:)]; %#ok<AGROW>
            cands{nb} = setdiff(cands{nb}, pick, 'stable');
            avail(nb) = numel(cands{nb});
            deficit(b) = deficit(b) - t;
            if deficit(b) <= 0, break; end
        end
    end
end

if any(deficit > 0)
    sample_idx = [];
    return;
end

pick_pos_within_ctrl = vertcat(chosen{:});
sample_idx = ctrl_idx(pick_pos_within_ctrl);
sample_idx = sample_idx(:);

if numel(sample_idx) ~= nRemove
    sample_idx = [];
end
end

function suf = extract_date_suffix(fieldName)
us = strfind(fieldName, '_');
if numel(us) >= 2
    suf = fieldName(us(end-2)+1:end);
else
    suf = fieldName(max(1,end-10):end);
end
end
