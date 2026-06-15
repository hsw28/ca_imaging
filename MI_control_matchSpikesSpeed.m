function mutualinfo_struct = MI_control_matchSpikesSpeed( ...
    spike_structure, pos_structure, velthreshold, dim, ...
    CA_timestamps, CSUS_id_struct, ca_MI, nIter, velEdges)
% MI_control_matchSpikesSpeed
% Control MI by removing spikes OUTSIDE CSUS, with removal chosen to match
% the velocity distribution of CSUS-period spikes (per cell, per iteration).
%
% Uses the same count-based joint-distribution MI formulation as
% mutualinfo_openfield_wCSUS and mutualinfo_openfield_shuff_noCSUS.
%
% Output:
%   mutualinfo_struct.(sprintf('MIspeedmatch_%s', spikes_date)) = [nCells x 3]
%
% Columns:
%   1 = 95th percentile of control distribution
%   2 = mean of control distribution
%   3 = percentile rank of observed MI relative to control

if nargin < 8 || isempty(nIter), nIter = 100; end
if nargin < 9, velEdges = []; end

maxK = 1;   % cap counts per sample: 0,1,2,3+

set(0,'DefaultFigureVisible','off');

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_MI     = fieldnames(ca_MI);
fields_cats   = fieldnames(CA_timestamps);
fields_CSUS   = fieldnames(CSUS_id_struct);

if numel(fields_spikes) ~= numel(fields_pos)
    error('spike_structure and pos_structure have different # days (pad if needed).');
end
if numel(fields_spikes) ~= numel(fields_MI)
    error('spike_structure and ca_MI have different # days.');
end
if numel(fields_spikes) ~= numel(fields_cats)
    error('spike_structure and CA_timestamps have different # days.');
end
if numel(fields_spikes) ~= numel(fields_CSUS)
    error('spike_structure and CSUS_id_struct have different # days.');
end

mutualinfo_struct = struct();

% only use last 3 entries (An, An-1, An-2)
fields_spikes = fields_spikes(end-2:end);
fields_pos    = fields_pos(end-2:end);
fields_cats   = fields_cats(end-2:end);
fields_CSUS   = fields_CSUS(end-2:end);
fields_MI     = fields_MI(end-2:end);

for iDay = 1:numel(fields_spikes)

    fieldName_MI = fields_MI{iDay}
    MI = ca_MI.(fieldName_MI);

    fieldName_CSUS = fields_CSUS{iDay};
    CSUS_id = CSUS_id_struct.(fieldName_CSUS);

    fieldName_spikes = fields_spikes{iDay};
    peaks_time = spike_structure.(fieldName_spikes);

    idxu = strfind(fieldName_spikes,'_');
    spikes_date = fieldName_spikes(idxu(2)+1:end);

    fieldName_pos = fields_pos{iDay};
    pos = pos_structure.(fieldName_pos);

    fieldName_cats = fields_cats{iDay};
    curr_CA_timestamps = CA_timestamps.(fieldName_cats);

    fprintf('\nProcessing %s\n', spikes_date);

    if isempty(pos)
        warning('Empty pos for %s', spikes_date);
        mutualinfo_struct.(sprintf('MIspeedmatch_%s', spikes_date)) = NaN;
        continue;
    end

    if size(pos,2) > 3
        error('Position data should have columns [time x y]');
    end
    if size(pos,2) > size(pos,1)
        pos = pos';
    end

    % Convert to timestamp units if needed
    if ~isempty(pos) && numel(curr_CA_timestamps) > 1
        dt_pos = median(diff(pos(:,1)), 'omitnan');
        dt_ca  = median(diff(curr_CA_timestamps), 'omitnan');
        if ~isnan(dt_pos) & ~isnan(dt_ca) & dt_pos < 0.1 * dt_ca
            pos = convertpostoframe(pos, curr_CA_timestamps);
        end
    end

    % Remove duplicate pos timestamps
    [~, uniqueIdxPos] = unique(pos(:,1), 'stable');
    pos = pos(uniqueIdxPos, :);

    % Remove duplicate CSUS timestamps
    [csus_t_unique, idxCSUS] = unique(CSUS_id(2,:), 'stable');
    csus_val_unique = CSUS_id(1, idxCSUS);

    % velocity
    vel = ca_velocity(pos);
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

    % Interpolate CSUS labels to velocity timestamps (0 = non-task)
    interp_CSUS = interp1(csus_t_unique, csus_val_unique, vel_time, 'nearest', 0);

    % Interpolate x/y onto vel_time
    interp_x = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    interp_y = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    % Non-task, running, valid position samples
    validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ...
               ~isnan(interp_x) & ~isnan(interp_y);

    sample_time = vel_time(validIdx);
    sample_x    = interp_x(validIdx);
    sample_y    = interp_y(validIdx);

    if isempty(sample_time)
        warning('No valid position samples for %s', spikes_date);
        mutualinfo_struct.(sprintf('MIspeedmatch_%s', spikes_date)) = NaN;
        continue;
    end

    xmin = min(sample_x);
    xmax = max(sample_x);
    ymin = min(sample_y);
    ymax = max(sample_y);

    if xmax == xmin || ymax == ymin
        warning('Degenerate position range for %s', spikes_date);
        mutualinfo_struct.(sprintf('MIspeedmatch_%s', spikes_date)) = NaN;
        continue;
    end

    xEdges = xmin:dim:(xmax + dim);
    yEdges = ymin:dim:(ymax + dim);

    % Assign included samples to bins
    xBin_true = discretize(sample_x, xEdges);
    yBin_true = discretize(sample_y, yEdges);

    validBinSamples = ~isnan(xBin_true) & ~isnan(yBin_true);
    sample_time_valid = sample_time(validBinSamples);
    xBin_true         = xBin_true(validBinSamples);
    yBin_true         = yBin_true(validBinSamples);

    if isempty(sample_time_valid)
        warning('No valid binned samples for %s', spikes_date);
        mutualinfo_struct.(sprintf('MIspeedmatch_%s', spikes_date)) = NaN;
        continue;
    end

    nXBins = numel(xEdges) - 1;
    nYBins = numel(yEdges) - 1;
    linBin_true = sub2ind([nXBins, nYBins], xBin_true, yBin_true);

    % Auto velocity edges if not provided
    if isempty(velEdges)
        vUse = vel_mag(validIdx);
        vUse = vUse(isfinite(vUse));
        if isempty(vUse)
            velEdges_day = linspace(velthreshold, velthreshold+50, 10);
        else
            hi = prctile(vUse, 99);
            hi = max(hi, velthreshold + 1);
            velEdges_day = linspace(velthreshold, hi, 11);
        end
    else
        velEdges_day = velEdges;
    end

    % Sample timing tolerance
    if numel(sample_time_valid) > 1
        dt_samp = median(diff(sample_time_valid), 'omitnan');
    else
        dt_samp = NaN;
    end
    if isnan(dt_samp) || dt_samp <= 0
        dt_samp = 1/7.5;
          warning('funky business')
    end
    maxAssignDist = dt_samp / 2;

    nCells = size(peaks_time,1);
    controlMI = nan(nCells, nIter);

    parfor iIter = 1:nIter
        s = RandStream('Threefry','Seed', iIter + 1000*iDay);
        RandStream.setGlobalStream(s);

        iterVals = nan(nCells,1);

        for k = 1:nCells

            if numel(MI) < k || isnan(MI(k))
                continue;
            end

            spk = peaks_time(k,:);
            spk = spk(~isnan(spk) & spk > 0);
            if isempty(spk)
                continue;
            end

            spk_vel  = interp1(vel_time, vel_mag, spk, 'linear', NaN);
            csus_spk = interp1(csus_t_unique, csus_val_unique, spk, 'nearest', 0);

            runMask = isfinite(spk_vel) & (spk_vel >= velthreshold);

            % spikes during CSUS (task) while running
            csusMask = runMask & (csus_spk > 0);
            v_csus = spk_vel(csusMask);
            nRemove = numel(v_csus);

            % candidate pool: NON-task spikes while running
            candMask = runMask & (csus_spk == 0);
            cand_spk = spk(candMask);
            cand_vel = spk_vel(candMask);

            if isempty(cand_spk)
                continue;
            end

            % remove same # from non-task spikes, matched to velocity profile of task spikes
            if nRemove > 0
                keep_spk = remove_velProfileMatched(cand_spk, cand_vel, v_csus, velEdges_day, s);
            else
                keep_spk = cand_spk;
            end

            if numel(keep_spk) < 1
                continue;
            end

            % Count per included non-task sample
            countPerSample = zeros(numel(sample_time_valid),1);

            for ss = 1:numel(keep_spk)
                [d, idxNearest] = min(abs(sample_time_valid - keep_spk(ss)));
                if d <= maxAssignDist
                    countPerSample(idxNearest) = countPerSample(idxNearest) + 1;
                end
            end

            iterVals(k) = count_mi_from_counts(countPerSample, linBin_true, nXBins, nYBins, maxK);
        end

        controlMI(:, iIter) = iterVals;
    end

    % Summarize shuffles per cell
    mutinfo = nan(3, nCells);
    for k = 1:nCells
        shuffVals = controlMI(k,:);
        shuffVals = shuffVals(~isnan(shuffVals));

        if isempty(shuffVals) || numel(MI) < k || isnan(MI(k))
            continue;
        end

        p95cut = prctile(shuffVals, 95);
        muShuff = mean(shuffVals);
        perc = mean(shuffVals <= MI(k));

        mutinfo(:,k) = [p95cut; muShuff; perc];
    end

    mutualinfo_struct.(sprintf('MIspeedmatch_%s', spikes_date)) = mutinfo';
end
end

% -------------------------------------------------------------------------
function keep_spk = remove_velProfileMatched(cand_spk, cand_vel, v_csus, velEdges, s)
% Remove N spikes from cand_spk where N=numel(v_csus),
% choosing removals so their velocity histogram matches v_csus across velEdges.

cand_spk = cand_spk(:);
cand_vel = cand_vel(:);
v_csus   = v_csus(:);

nRemove = numel(v_csus);
nCand   = numel(cand_spk);

if nRemove <= 0
    keep_spk = cand_spk;
    return;
end

if nRemove >= nCand
    keep_spk = [];
    return;
end

csusCounts = histcounts(v_csus, velEdges);

candBin = discretize(cand_vel, velEdges);
candBin = candBin(:);

removeMask = false(nCand, 1);
removed = 0;
nBins = numel(velEdges) - 1;

for b = 1:nBins
    want = csusCounts(b);
    if want <= 0, continue; end

    idx = find((candBin == b) & (~removeMask));
    if isempty(idx), continue; end

    take = min(want, numel(idx));
    if take > 0
        pick = randsample(s, idx, take, false);
        removeMask(pick) = true;
        removed = removed + take;
    end
end

remain = nRemove - removed;

if remain > 0
    targetBins = find(csusCounts > 0);

    if isempty(targetBins)
        idxLeft = find(~removeMask & ~isnan(candBin));
        if numel(idxLeft) < remain
            idxLeft = find(~removeMask);
        end
        take = min(remain, numel(idxLeft));
        if take > 0
            pick = randsample(s, idxLeft, take, false);
            removeMask(pick) = true;
        end
    else
        binDist = inf(nBins,1);
        for bb = 1:nBins
            binDist(bb) = min(abs(bb - targetBins));
        end

        dmax = max(binDist(isfinite(binDist)));
        if isempty(dmax) || ~isfinite(dmax), dmax = nBins; end

        for d = 0:dmax
            if remain <= 0, break; end

            binsAtD = find(binDist == d);
            if isempty(binsAtD), continue; end

            idx = find(~removeMask & ismember(candBin, binsAtD));
            if isempty(idx), continue; end

            take = min(remain, numel(idx));
            pick = randsample(s, idx, take, false);
            removeMask(pick) = true;
            remain = remain - take;
        end

        if remain > 0
            idxLeft = find(~removeMask);
            take = min(remain, numel(idxLeft));
            if take > 0
                pick = randsample(s, idxLeft, take, false);
                removeMask(pick) = true;
            end
        end
    end
end

keep_spk = cand_spk(~removeMask);
end

function mi = count_mi_from_counts(countPerSample, linBin_true, nXBins, nYBins, maxK)
% Full joint-distribution MI:
% I(X;K) = sum_x sum_k p(x,k) log2( p(x,k) / (p(x)p(k)) )

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
