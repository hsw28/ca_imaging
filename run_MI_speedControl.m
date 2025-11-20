function run_MI_speedControl(ratVar, velthreshold, dim, NSpeedBins)
% Updates base workspace variable <ratVar>:
%   <ratVar>.MI_speedControl.MI_<YYYY_MM_DD> = Nc×1 MI vector
%
% Args:
%   ratVar      : name of the rat variable in base workspace (e.g., 'rat0222')
%   velthreshold: speed threshold (e.g., 4)
%   dim         : bin size for CA_normalizePosData (e.g., 2.5)
%   NSpeedBins  : # bins for speed matching (e.g., 200)

if nargin < 2 || isempty(velthreshold), velthreshold = 4; end
if nargin < 3 || isempty(dim),         dim = 2.5;        end
if nargin < 4 || isempty(NSpeedBins),  NSpeedBins = 200; end

% pull the rat
assert(evalin('base', sprintf('exist(''%s'',''var'')', ratVar))==1, ...
    'Variable %s not found in base workspace.', ratVar);
rat = evalin('base', ratVar);

need = {'Ca_peaks','pos','Ca_ts','CS_times'};
missing = need(~cellfun(@(f)isfield(rat,f), need));
assert(isempty(missing), 'Missing fields in %s: %s', ratVar, strjoin(missing,', '));

% ---- pick 3 days around rat.An using autoDateList(rat) ----
dates = autoDateList(rat);  % cellstr of date strings like 'YYYY_MM_DD'
idx = [];
if isfield(rat,'An') && ~isempty(rat.An)
    idx = find(strcmp(dates, rat.An), 1);
end
if isempty(idx)
    dayIdx = max(1, numel(dates)-2):numel(dates);
elseif idx < 3
    dayIdx = max(1, idx-2):idx;
else
    dayIdx = (idx-2):idx;
end
days = dates(dayIdx);  % 1×3 (or fewer if not enough days)

% ---- filter each per-day struct to those dates ----
spikes3 = filterFieldsByDay(rat.Ca_peaks, days);
pos3    = filterFieldsByDay(rat.pos,     days);
ts3     = filterFieldsByDay(rat.Ca_ts,   days);
cs3     = filterFieldsByDay(rat.CS_times,days);

% compute per-day MI vectors (after control removal) on these 3 days
out = MI_taskWindow02_speedControl_days(spikes3, pos3, velthreshold, dim, ts3, cs3, NSpeedBins);

% attach to rat and write back
if ~isfield(rat, 'MI_speedControl') || ~isstruct(rat.MI_speedControl)
    rat.MI_speedControl = struct();
end
dayFns = fieldnames(out);
for i = 1:numel(dayFns)
    rat.MI_speedControl.(dayFns{i}) = out.(dayFns{i});   % Nc×1 vector
end
assignin('base', ratVar, rat);

fprintf('[%s] wrote %d day(s) into %s.MI_speedControl\n', mfilename, numel(dayFns), ratVar);
end

function rolling_struct = MI_taskWindow02_speedControl_days(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, CS_times_struct, NSpeedBins)
% For each selected day (fields in the provided structs):
%   1) Find spikes in [0,2] s relative to nearest CS (task window).
%   2) Remove SAME COUNT of spikes from times outside [0,2], matched by speed.
%   3) Compute MI on remaining movement spikes (same code path you use).
%   4) Store an Nc×1 vector of MI values in rolling_struct.MI_<date>.

if nargin < 7 || isempty(NSpeedBins), NSpeedBins = 200; end
set(0,'DefaultFigureVisible','off');

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_cats   = fieldnames(CA_timestamps);
fields_CS     = fieldnames(CS_times_struct);

% Align lists to common length by date suffix matching
% (we assume filterFieldsByDay already aligned, but double-check ordering)
date_sp = cellfun(@extract_date_suffix, fields_spikes, 'UniformOutput', false);
date_po = cellfun(@extract_date_suffix, fields_pos,    'UniformOutput', false);
date_ts = cellfun(@extract_date_suffix, fields_cats,   'UniformOutput', false);
date_cs = cellfun(@extract_date_suffix, fields_CS,     'UniformOutput', false);
common = intersect(intersect(date_sp, date_po), intersect(date_ts, date_cs), 'stable');
rolling_struct = struct();

for i = 1:numel(common)
    d = common{i};
    f_sp = fields_spikes(strcmp(date_sp, d));
    f_po = fields_pos(   strcmp(date_po, d));
    f_ts = fields_cats(  strcmp(date_ts, d));
    f_cs = fields_CS(    strcmp(date_cs, d));
    if isempty(f_sp) || isempty(f_po) || isempty(f_ts) || isempty(f_cs)
        continue;
    end
    f_sp = f_sp{1}; f_po = f_po{1}; f_ts = f_ts{1}; f_cs = f_cs{1};

    peaks_time = spike_structure.(f_sp);   % Nc x NspikeTimes
    pos        = pos_structure.(f_po);     % [t x y]
    cs_times   = CS_times_struct.(f_cs);   % vector of CS onset times (sec)
    date_suf   = d;
    Nc         = size(peaks_time,1);

    % default output for this day
    per_cell_MI = nan(Nc,1);

    if isempty(cs_times) || all(~isfinite(cs_times))
        rolling_struct.(sprintf('MI_%s',date_suf)) = per_cell_MI;
        continue;
    end

    % ---- velocity/position resampling grid (movement only) ----
    [~, uniqPos] = unique(pos(:,1), 'stable'); pos = pos(uniqPos,:);
    vel      = ca_velocity(pos);        % [speed; time]
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

    ix = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    iy = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    validIdx = (vel_mag >= velthreshold) & isfinite(ix) & isfinite(iy);
    posDat   = [vel_time(validIdx), ix(validIdx), iy(validIdx)];
    if size(posDat,1) < 10
        rolling_struct.(sprintf('MI_%s',date_suf)) = per_cell_MI;
        continue;
    end

    % ---- per cell ----
    parfor k = 1:Nc
        spk = peaks_time(k,:); spk = spk(isfinite(spk) & spk > 0);
        if isempty(spk)
            per_cell_MI(k) = NaN; %#ok<PFBNS>
            continue;
        end

        % movement spikes
        spike_vel   = interp1(vel_time, vel_mag, spk, 'linear', NaN);
        move_spikes = spk(spike_vel >= velthreshold);
        if numel(move_spikes) < 2
            per_cell_MI(k) = NaN; %#ok<PFBNS>
            continue;
        end

        % relative time to nearest CS (for movement spikes)
        nearest_cs_for_move = interp1(cs_times, cs_times, move_spikes, 'nearest', NaN);
        rel_t_move          = move_spikes - nearest_cs_for_move;

        % task window [0,2] s
        in_task  = (rel_t_move >= 0) & (rel_t_move <= 2);
        nTask    = nnz(in_task);

        % if no task spikes, MI = baseline movement MI
        if nTask == 0
            try
                [~,~,~,~,spikeprob,occprob] = CA_normalizePosData(move_spikes, posDat, dim, 1.0);
                if size(spikeprob,1) < size(spikeprob,2), spikeprob = spikeprob'; end
                if size(occprob,1)  < size(occprob,2),  occprob  = occprob';  end
                per_cell_MI(k) = mutualinfo([spikeprob, occprob]); %#ok<PFBNS>
            catch
                per_cell_MI(k) = NaN; %#ok<PFBNS>
            end
            continue;
        end

        % control pool: all movement spikes *outside* [0,2]
        ctrl_mask = ~in_task;
        ctrl_idx  = find(ctrl_mask);
        ctrl_idx  = ctrl_idx(:);  % shape safety
        if numel(ctrl_idx) < nTask
            per_cell_MI(k) = NaN; %#ok<PFBNS>
            continue;
        end

        % speed matching: stratified by quantile edges from union(task, control speeds)
        move_speed = interp1(vel_time, vel_mag, move_spikes, 'linear', NaN);
        w_speed    = move_speed(in_task);
        c_speed    = move_speed(ctrl_idx);

        edges = make_quant_edges([w_speed(:); c_speed(:)], NSpeedBins);
        wbins = discretize(w_speed,    edges);
        cbins = discretize(c_speed,    edges);
        wbins = wbins(:); cbins = cbins(:);  % shape safety

        % target counts per bin = histogram of window (task) speeds
        K = max([wbins(:); cbins(:)], [], 'omitnan');
        if isempty(K) || ~isfinite(K), K = 1; end

        subs = wbins(~isnan(wbins));  % <-- FIX: accumarray subs must be col
        subs = subs(:);
        hist_w = accumarray(subs, 1, [K, 1], @sum, 0);

        % candidates per bin from control pool
        cand_idx_within_ctrl = cell(K,1);
        for b = 1:K
            cand_idx_within_ctrl{b} = find(cbins == b);  % positions within ctrl_idx
        end

        % assemble a matched sample of EXACT nTask spikes
        chosen_ctrl_positions = [];
        left = nTask;
        % first pass: take up to the needed count per bin (cap by availability)
        for b = 1:K
            need_b = hist_w(b);
            if need_b <= 0, continue; end
            avail_b = numel(cand_idx_within_ctrl{b});
            take_b  = min(need_b, avail_b);
            if take_b > 0
                pick = randsample(cand_idx_within_ctrl{b}, take_b, false);
                chosen_ctrl_positions = [chosen_ctrl_positions; pick(:)]; %#ok<AGROW>
                % remove picked from availability
                cand_idx_within_ctrl{b} = setdiff(cand_idx_within_ctrl{b}, pick, 'stable');
                left = left - take_b;
            end
        end
        % top-up uniformly from whatever remains if we still need more
        if left > 0
            remaining = vertcat(cand_idx_within_ctrl{:});
            if numel(remaining) >= left
                extra = randsample(remaining, left, false);
                chosen_ctrl_positions = [chosen_ctrl_positions; extra(:)]; %#ok<AGROW>
                left = 0;
            end
        end

        if left > 0  % still short → fail for this cell
            per_cell_MI(k) = NaN; %#ok<PFBNS>
            continue;
        end

        % convert positions-within-ctrl to indices into move_spikes
        sample_idx = ctrl_idx(chosen_ctrl_positions);

        % remove those control spikes
        keep_idx = true(size(move_spikes));
        keep_idx(sample_idx) = false;
        keep_spikes = move_spikes(keep_idx);

        if numel(keep_spikes) < 2
            per_cell_MI(k) = NaN; %#ok<PFBNS>
            continue;
        end

        % --- MI on remaining movement spikes (same structure as your code) ---
        try
            [~,~,~,~,spikeprob,occprob] = CA_normalizePosData(keep_spikes, posDat, dim, 1.0);
            if size(spikeprob,1) < size(spikeprob,2), spikeprob = spikeprob'; end
            if size(occprob,1)  < size(occprob,2),  occprob  = occprob';  end
            per_cell_MI(k) = mutualinfo([spikeprob, occprob]); %#ok<PFBNS>
        catch
            per_cell_MI(k) = NaN; %#ok<PFBNS>
        end
    end

    rolling_struct.(sprintf('MI_%s', date_suf)) = per_cell_MI;
end
end

function Sout = filterFieldsByDay(Sin, days)
% Keep only fields whose date suffix matches any in 'days' (cellstr)
% and preserve input order of 'days'.
Sout = struct();
if isempty(Sin) || ~isstruct(Sin) || isempty(days), return; end
fns = fieldnames(Sin);
fn_dates = cellfun(@extract_date_suffix, fns, 'UniformOutput', false);
for i = 1:numel(days)
    idx = find(strcmp(fn_dates, days{i}), 1);
    if ~isempty(idx)
        Sout.(fns{idx}) = Sin.(fns{idx});
    end
end
end

function edges = make_quant_edges(x, NB)
x = x(isfinite(x));
if isempty(x) || NB < 1
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

function suf = extract_date_suffix(fieldName)
% expects something like 'peaks_2023_05_03' or 'pos_2023_05_03'
us = strfind(fieldName, '_');
if numel(us) >= 2
    suf = fieldName(us(end-2)+1:end);
else
    suf = fieldName(max(1,end-10):end);
end
end
