function rolling_struct = MI_control_rollingPost( ...
    spike_structure, pos_structure, velthreshold, dim, ...
    CA_timestamps, CS_times_struct, nIter, winEdges, ...
    ControlMatch, NSpeedBins, SpaceBinSize, RM_structure, SpeedBinSize)
% MI_control_rollingPost
%
% CS-aligned rolling removal/control using count-based joint MI:
%
%   I(X;K) = sum_x sum_k p(x,k) * log2( p(x,k) / (p(x)*p(k)) )
%
% where:
%   X = spatial bin
%   K = event count per included sample
%
% Test  : remove spikes ONLY in the CURRENT CS-aligned window from the movement MI pool.
% Control: remove the SAME number of spikes from a WINDOW-LOCAL pool that EXCLUDES:
%           (1) trial time: [CS, CS+2]
%           (2) and the current window ± guard
%          with optional matching of SPEED and/or SPACE.
%
% Output:
%   rolling_struct.(sprintf('MI_%s', date)) with fields:
%     .winEdges, .MI_win (W×Nc), .MI_rand (W×Nc), .MI_base (1×Nc), .dMI
%     .counts: struct of spike-count diagnostics

if nargin < 7 || isempty(nIter),        nIter = 5; end
if nargin < 8 || isempty(winEdges),     winEdges = [0 2; 1 3; 2 4; 3 5; 4 6; 5 7; 8 10; 9 11; 10 12; 11 13]; end
if nargin < 9 || isempty(ControlMatch), ControlMatch = 'none'; end
if nargin < 10,                         NSpeedBins = 5; end
if nargin < 11 || isempty(SpaceBinSize),SpaceBinSize = dim; end
if nargin < 13,                         SpeedBinSize = []; end

if isempty(NSpeedBins) && isempty(SpeedBinSize)
    NSpeedBins = 5;
elseif ~isempty(NSpeedBins) && ~isempty(SpeedBinSize)
    error('Choose either NSpeedBins or SpeedBinSize; leave the other empty.');
end
if ~isempty(NSpeedBins) && (~isscalar(NSpeedBins) || NSpeedBins < 1 || NSpeedBins ~= floor(NSpeedBins))
    error('NSpeedBins must be a positive integer.');
end
if ~isempty(SpeedBinSize) && (~isscalar(SpeedBinSize) || ~isfinite(SpeedBinSize) || SpeedBinSize <= 0)
    error('SpeedBinSize must be a positive scalar in cm/s.');
end

TrialWinSecs = [0 2];
GuardSecs    = 0.5;
maxK         = 1;   % 0,1,2,3+

cm = lower(strrep(ControlMatch,'_',''));
if strcmp(cm,'both'), cm = 'speedspace'; end

set(0,'DefaultFigureVisible','off');

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_cats   = fieldnames(CA_timestamps);
fields_CS     = fieldnames(CS_times_struct);
fields_RM     = fieldnames(RM_structure);

fields_spikes = fields_spikes(max(1,end-2):end);
fields_pos    = fields_pos(   max(1,end-2):end);
fields_cats   = fields_cats(  max(1,end-2):end);
fields_CS     = fields_CS(    max(1,end-2):end);
fields_RM     = fields_RM(    max(1,end-2):end);

rolling_struct = struct();

for iDay = 1:numel(fields_spikes)
    fieldName_sp  = fields_spikes{iDay}
    fieldName_pos = fields_pos{iDay};
    fieldName_ts  = fields_cats{iDay};
    fieldName_CS  = fields_CS{iDay};
    fieldName_RM  = fields_RM{iDay};

    peaks_time = spike_structure.(fieldName_sp);
    pos        = pos_structure.(fieldName_pos);
    curr_CA_timestamps = CA_timestamps.(fieldName_ts);
    cs_times   = CS_times_struct.(fieldName_CS);
    RM         = RM_structure.(fieldName_RM);

    spikes_date = extract_date_suffix(fieldName_sp);

    if isempty(pos)
        warning('Empty pos for %s; skipping day.', fieldName_pos);
        out = struct('winEdges',winEdges,'MI_win',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_rand',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_base',nan(1,size(peaks_time,1)),'dMI',nan(size(winEdges,1),size(peaks_time,1)));
        out.counts = struct('N_move',nan(1,size(peaks_time,1)), 'N_inWin',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_win',nan(size(winEdges,1),size(peaks_time,1)), 'N_ctrl_idx',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_ctrl',nan(size(winEdges,1),size(peaks_time,1)), 'Iter_OK',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_removed_ctrl',nan(size(winEdges,1),size(peaks_time,1)));
        rolling_struct.(sprintf('MI_%s', spikes_date)) = out;
        continue;
    end

    if size(pos,2) > 3
        error('Position data should have columns [time x y]');
    end
    if size(pos,2) > size(pos,1)
        pos = pos';
    end

    if ~isempty(pos) && numel(curr_CA_timestamps) > 1
        dt_pos = median(diff(pos(:,1)), 'omitnan');
        dt_ca  = median(diff(curr_CA_timestamps), 'omitnan');
        if ~isnan(dt_pos) & ~isnan(dt_ca) & dt_pos < 0.1 * dt_ca
            pos = convertpostoframe(pos, curr_CA_timestamps);
        end
    end

    [~, uniqPos] = unique(pos(:,1), 'stable');
    pos = pos(uniqPos,:);

    vel      = ca_velocity(pos);
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

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
        rolling_struct.(sprintf('MI_%s', spikes_date)) = out;
        continue;
    end

    validIdx = (vel_mag >= velthreshold) & ~isnan(ix) & ~isnan(iy);
    sample_time = vel_time(validIdx);
    sample_x    = ix(validIdx);
    sample_y    = iy(validIdx);

    if isempty(sample_time)
        warning('No valid movement samples for %s; skipping day.', spikes_date);
        out = struct('winEdges',winEdges,'MI_win',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_rand',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_base',nan(1,size(peaks_time,1)),'dMI',nan(size(winEdges,1),size(peaks_time,1)));
        out.counts = struct('N_move',nan(1,size(peaks_time,1)), 'N_inWin',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_win',nan(size(winEdges,1),size(peaks_time,1)), 'N_ctrl_idx',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_ctrl',nan(size(winEdges,1),size(peaks_time,1)), 'Iter_OK',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_removed_ctrl',nan(size(winEdges,1),size(peaks_time,1)));
        rolling_struct.(sprintf('MI_%s', spikes_date)) = out;
        continue;
    end

    xmin = min(sample_x);
    xmax = max(sample_x);
    ymin = min(sample_y);
    ymax = max(sample_y);

    if xmax == xmin || ymax == ymin
        warning('Degenerate position range for %s; skipping day.', spikes_date);
        out = struct('winEdges',winEdges,'MI_win',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_rand',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_base',nan(1,size(peaks_time,1)),'dMI',nan(size(winEdges,1),size(peaks_time,1)));
        out.counts = struct('N_move',nan(1,size(peaks_time,1)), 'N_inWin',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_win',nan(size(winEdges,1),size(peaks_time,1)), 'N_ctrl_idx',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_ctrl',nan(size(winEdges,1),size(peaks_time,1)), 'Iter_OK',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_removed_ctrl',nan(size(winEdges,1),size(peaks_time,1)));
        rolling_struct.(sprintf('MI_%s', spikes_date)) = out;
        continue;
    end

    xEdges = xmin:dim:(xmax + dim);
    yEdges = ymin:dim:(ymax + dim);

    xBin_true = discretize(sample_x, xEdges);
    yBin_true = discretize(sample_y, yEdges);

    validBinSamples = ~isnan(xBin_true) & ~isnan(yBin_true);
    sample_time_valid = sample_time(validBinSamples);
    xBin_true         = xBin_true(validBinSamples);
    yBin_true         = yBin_true(validBinSamples);

    if isempty(sample_time_valid)
        warning('No valid binned movement samples for %s; skipping day.', spikes_date);
        out = struct('winEdges',winEdges,'MI_win',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_rand',nan(size(winEdges,1),size(peaks_time,1)), ...
            'MI_base',nan(1,size(peaks_time,1)),'dMI',nan(size(winEdges,1),size(peaks_time,1)));
        out.counts = struct('N_move',nan(1,size(peaks_time,1)), 'N_inWin',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_win',nan(size(winEdges,1),size(peaks_time,1)), 'N_ctrl_idx',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_keep_ctrl',nan(size(winEdges,1),size(peaks_time,1)), 'Iter_OK',nan(size(winEdges,1),size(peaks_time,1)), ...
            'N_removed_ctrl',nan(size(winEdges,1),size(peaks_time,1)));
        rolling_struct.(sprintf('MI_%s', spikes_date)) = out;
        continue;
    end

    nXBins = numel(xEdges) - 1;
    nYBins = numel(yEdges) - 1;
    linBin_true = sub2ind([nXBins, nYBins], xBin_true, yBin_true);

    if numel(sample_time_valid) > 1
        dt_samp = median(diff(sample_time_valid), 'omitnan');
    else
        dt_samp = NaN;
    end
    if isnan(dt_samp) || dt_samp <= 0
        dt_samp = 1/15;
    end
    maxAssignDist = dt_samp / 2;

    nCells = size(peaks_time,1);
    nWins  = size(winEdges,1);

    MI_base = nan(1, nCells);
    MI_win  = nan(nWins, nCells);
    MI_rand = nan(nWins, nCells);

    N_inWin        = nan(nWins, nCells);
    N_keep_win     = nan(nWins, nCells);
    N_ctrl_idx     = nan(nWins, nCells);
    N_keep_ctrl    = nan(nWins, nCells);
    N_removed_ctrl = nan(nWins, nCells);
    N_move         = nan(1,    nCells);
    Iter_OK        = nan(nWins, nCells);

    for k = 1:nCells
        spk = peaks_time(k,:);
        spk = spk(~isnan(spk) & spk > 0);

        if RM(k) == 0 || isempty(spk)
            continue;
        end

        spike_vel = interp1(vel_time, vel_mag, spk, 'linear', NaN);
        move_spikes = spk(~isnan(spike_vel) & (spike_vel >= velthreshold));

        if isempty(move_spikes) || numel(sample_time_valid) < 10
            continue;
        end

        MI_base(k) = count_mi_from_spikes(move_spikes, sample_time_valid, linBin_true, maxAssignDist, nXBins, nYBins, maxK);
    end

    parfor k = 1:nCells
        spk = peaks_time(k,:);
        spk = spk(~isnan(spk) & spk > 0);

        MI_win_k  = nan(nWins,1);
        MI_rand_k = nan(nWins,1);

        N_inWin_k        = nan(nWins,1);
        N_keep_win_k     = nan(nWins,1);
        N_ctrl_idx_k     = nan(nWins,1);
        N_keep_ctrl_k    = nan(nWins,1);
        N_removed_ctrl_k = nan(nWins,1);
        Iter_OK_k        = nan(nWins,1);
        N_move_k         = NaN;

        if isempty(spk) || isnan(MI_base(k)) || numel(sample_time_valid) < 10
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
        move_spikes = spk(~isnan(spike_vel) & (spike_vel >= velthreshold));
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

        % Assign each movement spike to a position sample once.  All rolling
        % windows and control iterations below reuse this assignment instead
        % of repeating a full nearest-sample search for every MI calculation.
        [base_count_per_sample, assigned_sample_idx] = assign_spikes_to_samples( ...
            move_spikes, sample_time_valid, maxAssignDist);

        nearest_cs_for_move = interp1(cs_times, cs_times, move_spikes, 'nearest', NaN);
        rel_t_move          = move_spikes - nearest_cs_for_move;

        move_speed = interp1(vel_time, vel_mag, move_spikes, 'linear', NaN);
        move_x = interp1(pos(:,1), pos(:,2), move_spikes, 'linear', NaN);
        move_y = interp1(pos(:,1), pos(:,3), move_spikes, 'linear', NaN);

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

            w_mask_on_move = (rel_t_move > wstart) & (rel_t_move <= wend);
            w_idx          = find(w_mask_on_move);
            nRemove        = numel(w_idx);

            N_inWin_k(w) = nRemove;

            if nRemove < 2
                in_trial = (rel_t_move >= TrialWinSecs(1)) & (rel_t_move <= TrialWinSecs(2));
                in_win_guard = (rel_t_move > (wstart-GuardSecs)) & (rel_t_move <= (wend+GuardSecs));
                ctrl_mask_tmp = ~in_trial & ~in_win_guard;
                N_ctrl_idx_k(w) = sum(ctrl_mask_tmp);
                Iter_OK_k(w) = 0;
                continue;
            end

            nKeepWin = numel(move_spikes) - nRemove;
            N_keep_win_k(w) = nKeepWin;

            if nKeepWin >= 2
                win_counts = remove_assigned_spikes( ...
                    base_count_per_sample, assigned_sample_idx, w_idx);
                MI_win_k(w) = count_mi_from_counts( ...
                    win_counts, linBin_true, nXBins, nYBins, maxK);
            end

            in_trial = (rel_t_move >= TrialWinSecs(1)) & (rel_t_move <= TrialWinSecs(2));
            in_win_guard = (rel_t_move > (wstart-GuardSecs)) & (rel_t_move <= (wend+GuardSecs));

            ctrl_mask = ~in_trial & ~in_win_guard;
            ctrl_idx  = find(ctrl_mask);
            N_ctrl_idx_k(w) = numel(ctrl_idx);

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
                        edges = make_speed_edges( ...
                            [w_speed(:); ctrl_speed(:)], NSpeedBins, SpeedBinSize);
                        wbins    = discretize(w_speed,    edges);
                        ctrlbins = discretize(ctrl_speed, edges);
                        sample_idx = stratified_sample(ctrl_idx, ctrlbins, wbins, nRemove);

                    case 'space'
                        if isempty(edgesX) || isempty(edgesY)
                            if numel(ctrl_idx) < nRemove
                                sample_idx = [];
                            else
                                sample_idx = randsample(ctrl_idx, nRemove, false);
                            end
                        else
                            [wBX,wBY]  = deal(discretize(w_x, edgesX), discretize(w_y, edgesY));
                            [cBX,cBY]  = deal(discretize(move_x(ctrl_idx), edgesX), discretize(move_y(ctrl_idx), edgesY));
                            w_joint    = wBX.*1e6 + wBY;
                            ctrl_joint = cBX.*1e6 + cBY;
                            sample_idx = stratified_sample(ctrl_idx, ctrl_joint, w_joint, nRemove);
                        end

                    case 'speedspace'
                        ctrl_speed = move_speed(ctrl_idx);
                        edges = make_speed_edges( ...
                            [w_speed(:); ctrl_speed(:)], NSpeedBins, SpeedBinSize);
                        wbins    = discretize(w_speed,    edges);
                        ctrlbins = discretize(ctrl_speed, edges);
                        if isempty(edgesX) || isempty(edgesY)
                            sample_idx = stratified_sample(ctrl_idx, ctrlbins, wbins, nRemove);
                        else
                            [wBX,wBY]  = deal(discretize(w_x, edgesX), discretize(w_y, edgesY));
                            [cBX,cBY]  = deal(discretize(move_x(ctrl_idx), edgesX), discretize(move_y(ctrl_idx), edgesY));
                            wjoint     = wbins.*1e8 + wBX.*1e4 + wBY;
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

                if isempty(sample_idx) || numel(sample_idx) ~= nRemove
                    continue;
                end

                if (numel(move_spikes) - numel(sample_idx)) < 2
                    continue;
                end

                ctrl_counts = remove_assigned_spikes( ...
                    base_count_per_sample, assigned_sample_idx, sample_idx);
                thisMI = count_mi_from_counts( ...
                    ctrl_counts, linBin_true, nXBins, nYBins, maxK);
                if isnan(thisMI)
                    continue;
                end

                mi_accum = mi_accum + thisMI;
                n_ok = n_ok + 1;
            end

            Iter_OK_k(w) = n_ok;
            if n_ok > 0
                MI_rand_k(w)        = mi_accum / n_ok;
                N_removed_ctrl_k(w) = nRemove;
                N_keep_ctrl_k(w)    = numel(move_spikes) - nRemove;
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


function mi = count_mi_from_spikes(spike_times, sample_time_valid, linBin_true, maxAssignDist, nXBins, nYBins, maxK)

if isempty(spike_times) || isempty(sample_time_valid) || isempty(linBin_true)
    mi = NaN;
    return;
end

[countPerSample, ~] = assign_spikes_to_samples( ...
    spike_times, sample_time_valid, maxAssignDist);

mi = count_mi_from_counts(countPerSample, linBin_true, nXBins, nYBins, maxK);
end


function [countPerSample, assignedSampleIdx] = assign_spikes_to_samples(spike_times, sample_time_valid, maxAssignDist)

countPerSample = zeros(numel(sample_time_valid),1);
assignedSampleIdx = nan(numel(spike_times),1);

for s = 1:numel(spike_times)
    [d, idxNearest] = min(abs(sample_time_valid - spike_times(s)));
    if d <= maxAssignDist
        assignedSampleIdx(s) = idxNearest;
        countPerSample(idxNearest) = countPerSample(idxNearest) + 1;
    end
end
end


function countPerSample = remove_assigned_spikes(baseCountPerSample, assignedSampleIdx, spikeIdx)

countPerSample = baseCountPerSample;
sampleIdx = assignedSampleIdx(spikeIdx);
sampleIdx = sampleIdx(isfinite(sampleIdx));

if isempty(sampleIdx)
    return;
end

removedPerSample = accumarray(sampleIdx, 1, size(countPerSample), @sum, 0);
countPerSample = countPerSample - removedPerSample;
end


function mi = count_mi_from_counts(countPerSample, linBin_true, nXBins, nYBins, maxK)

if isempty(countPerSample) || isempty(linBin_true)
    mi = NaN;
    return;
end

countPerSample = countPerSample(:);
linBin_true = linBin_true(:);

if ~isempty(maxK)
    countPerSample(countPerSample > maxK) = maxK;
end

[countVals, ~, kIdx] = unique(countPerSample); %#ok<ASGLU>
nK = numel(countVals);

jointCounts = accumarray([linBin_true, kIdx], 1, [nXBins*nYBins, nK], @sum, 0);

occCounts = sum(jointCounts, 2);
validSpace = occCounts > 0;
jointCounts = jointCounts(validSpace, :);

N = sum(jointCounts(:));
if N <= 0
    mi = NaN;
    return;
end

P_xk = jointCounts / N;
P_x  = sum(P_xk, 2);
P_k  = sum(P_xk, 1);

mi = 0;
for ix = 1:size(P_xk,1)
    for ik = 1:size(P_xk,2)
        p = P_xk(ix,ik);
        if p > 0 && P_x(ix) > 0 && P_k(ik) > 0
            mi = mi + p * log2(p / (P_x(ix) * P_k(ik)));
        end
    end
end
end


function edges = make_lin_edges(a, b, step)
if ~isfinite(a) || ~isfinite(b) || ~(b>a) || step<=0
    edges = [];
    return;
end
start = floor(a/step)*step;
stop  = ceil(b/step)*step;
edges = start:step:stop;
if numel(edges)<2
    edges = [a b];
end
end

function edges = make_quant_edges(x, NB)
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


function edges = make_speed_edges(x, NSpeedBins, SpeedBinSize)
x = x(isfinite(x));
if isempty(x)
    edges = [];
elseif ~isempty(SpeedBinSize)
    edges = make_lin_edges(min(x), max(x), SpeedBinSize);
else
    edges = make_quant_edges(x, NSpeedBins);
end
end


function sample_idx = stratified_sample(ctrl_idx, ctrl_bins, w_bins, nRemove)
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
        sample_idx = randsample(ctrl_idx, nRemove, false);
        sample_idx = sample_idx(:);
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
    if numel(ctrl_idx) < nRemove
        sample_idx = [];
        return;
    end
    sample_idx = randsample(ctrl_idx, nRemove, false);
    sample_idx = sample_idx(:);
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
    sample_idx = [];
    return;
end

taken  = zeros(K,1);
chosen = cell(K,1);
for b = 1:K
    want = tgt(b);
    if want <= 0 || avail(b) == 0
        continue;
    end
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
    if isempty(todo)
        break;
    end
    for i = 1:numel(todo)
        b = todo(i);
        need = deficit(b);
        if need <= 0
            continue;
        end

        neigh = [];
        L = b - radius; if L >= 1, neigh = [neigh, L]; end
        R = b + radius; if R <= K, neigh = [neigh, R]; end
        if isempty(neigh)
            continue;
        end

        neigh_av = sum(avail(neigh));
        if neigh_av <= 0
            continue;
        end

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
            if take_each(j) <= 0 || avail(nb) <= 0
                continue;
            end
            t = min(take_each(j), avail(nb));
            if t <= 0
                continue;
            end
            pick = randsample(cands{nb}, t, false);
            chosen{b} = [chosen{b}; pick(:)];
            cands{nb} = setdiff(cands{nb}, pick, 'stable');
            avail(nb) = numel(cands{nb});
            deficit(b) = deficit(b) - t;
            if deficit(b) <= 0
                break;
            end
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
