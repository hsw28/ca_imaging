function mutualinfo_struct = MI_control_matchSpikesSpeed(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, CSUS_id_struct, ca_MI, nIter, velEdges)
% MI_control_matchSpikes_velMatched
% Control MI by removing spikes OUTSIDE CSUS, with removal chosen to match
% the velocity distribution of CSUS-period spikes (per cell, per iteration).
%
% Inputs same as MI_control_matchSpikes, plus optional:
%   velEdges – velocity bin edges for matching (e.g., linspace(4,40,13)).
%              If empty/not provided, it will auto-define edges from data.

if nargin < 8 || isempty(nIter), nIter = 100; end
if nargin < 9, velEdges = []; end

set(0,'DefaultFigureVisible','off');

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_MI     = fieldnames(ca_MI);
fields_cats   = fieldnames(CA_timestamps);
fields_CSUS   = fieldnames(CSUS_id_struct);

if numel(fields_spikes) ~= numel(fields_pos)
    error('spike_structure and pos_structure have different # days (pad if needed).');
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
    curr_CA_timestamps = CA_timestamps.(fieldName_cats); %#ok<NASGU>  % kept for parity

    % velocity (your helper)
    vel = ca_velocity(pos);
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

    % Interpolate CSUS labels to velocity timestamps (0 = non-task)
    interp_CSUS = interp1(CSUS_id(2,:), CSUS_id(1,:), vel_time, 'nearest', 0);

    % de-dup pos timestamps, then interpolate x/y onto vel_time
    [~, uniqueIdx] = unique(pos(:,1), 'stable');
    pos = pos(uniqueIdx, :);

    interp_x = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    interp_y = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    % Non-task, running, valid position frames
    validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ~isnan(interp_x) & ~isnan(interp_y);
    posDat = [vel_time(validIdx), interp_x(validIdx), interp_y(validIdx)];


    % Auto velocity edges if not provided
    if isempty(velEdges)
        vUse = vel_mag(validIdx);
        vUse = vUse(isfinite(vUse));
        if isempty(vUse)
            % fallback
            velEdges = linspace(velthreshold, velthreshold+50, 10);
        else
            hi = prctile(vUse, 99);
            hi = max(hi, velthreshold + 1);
            velEdges = linspace(velthreshold, hi, 11); % 10 bins by default
        end
    end

    nCells = size(peaks_time,1);
    controlMI = nan(nCells, nIter);

    parfor iIter = 1:nIter
        % stable per-iteration randomness inside parfor
        s = RandStream('Threefry','Seed', iIter + 1000*iDay);
        RandStream.setGlobalStream(s);

        for k = 1:nCells
            spk = peaks_time(k,:);
            spk = spk(~isnan(spk) & spk > 0);
            if isempty(spk), continue; end

            spk_vel = interp1(vel_time, vel_mag, spk, 'linear', NaN);
            csus_spk = interp1(CSUS_id(2,:), CSUS_id(1,:), spk, 'nearest', 0);

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

            % If there are task spikes, remove same # from cand_spk,
            % matching the velocity profile of v_csus.
            if nRemove > 0
                keep_spk = remove_velProfileMatched(cand_spk, cand_vel, v_csus, velEdges, s);
            else
                keep_spk = cand_spk;
            end

            if numel(keep_spk) < 1
                continue;
            end

            [~, ~, ~, ~, spikeprob, occprob] = CA_normalizePosData(keep_spk, posDat, dim, 1.0);

            if size(spikeprob,1) < size(spikeprob,2), spikeprob = spikeprob'; end
            if size(occprob,1)   < size(occprob,2),   occprob   = occprob';   end

            controlMI(k, iIter) = mutualinfo([spikeprob, occprob]);
        end
    end

    % Summarize shuffles per cell
    mutinfo = nan(3, nCells);
    for k = 1:nCells
        shuffVals = controlMI(k,:);
        shuffVals = shuffVals(~isnan(shuffVals));
        if isempty(shuffVals), continue; end

        p95cut = prctile(shuffVals, 95);
        muShuff = mean(shuffVals);
        perc = sum(MI(k) > shuffVals) / numel(shuffVals);

        mutinfo(:,k) = [p95cut; muShuff; perc];
    end

    mutualinfo_struct.(sprintf('MIspeedmatch_%s', spikes_date)) = mutinfo';

end
end

% -------------------------------------------------------------------------
function keep_spk = remove_velProfileMatched(cand_spk, cand_vel, v_csus, velEdges, s)
% Remove N spikes from cand_spk where N=numel(v_csus),
% choosing removals so their velocity histogram matches v_csus across velEdges.

% ---- force column vectors (prevents implicit expansion N×N) ----
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

% Desired removals per velocity bin (from CSUS spikes)
csusCounts = histcounts(v_csus, velEdges);

% Candidate bin membership (column)
candBin = discretize(cand_vel, velEdges);
candBin = candBin(:);

removeMask = false(nCand, 1);
removed = 0;
nBins = numel(velEdges) - 1;

for b = 1:nBins
    want = csusCounts(b);
    if want <= 0, continue; end

    idx = find((candBin == b) & (~removeMask));   % now 1×N & 1×N (columns) -> OK
    if isempty(idx), continue; end

    take = min(want, numel(idx));
    if take > 0
        pick = randsample(s, idx, take, false);
        removeMask(pick) = true;
        removed = removed + take;
    end
end

% If bins underfilled, remove remaining from "close" bins (nearest in velocity-bin space)
remain = nRemove - removed;

if remain > 0

    targetBins = find(csusCounts > 0);

    % If somehow no target bins, just remove uniformly from remaining candidates
    if isempty(targetBins)
        idxLeft = find(~removeMask & ~isnan(candBin));
        if numel(idxLeft) < remain
            idxLeft = find(~removeMask); % last resort
        end
        take = min(remain, numel(idxLeft));
        if take > 0
            pick = randsample(s, idxLeft, take, false);
            removeMask(pick) = true;
        end
    else
        % distance of each bin to nearest target bin
        binDist = inf(nBins,1);
        for bb = 1:nBins
            binDist(bb) = min(abs(bb - targetBins));
        end

        % iterate outward by distance: 0,1,2,... taking uniformly within bins at that distance
        dmax = max(binDist(isfinite(binDist)));
        if isempty(dmax) || ~isfinite(dmax), dmax = nBins; end

        for d = 0:dmax
            if remain <= 0, break; end

            binsAtD = find(binDist == d);
            if isempty(binsAtD), continue; end

            idx = find(~removeMask & ismember(candBin, binsAtD));
            if isempty(idx), continue; end

            take = min(remain, numel(idx));
            pick = randsample(s, idx, take, false);  % uniform within closest bins
            removeMask(pick) = true;
            remain = remain - take;
        end

        % If still need removals (e.g., lots of NaN bins), remove uniformly from anything left
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
