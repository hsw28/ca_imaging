function f = mutualinfo_openfield_noCSUS_DEP( ...
    spike_structure, pos_structure, velthreshold, dim, ...
    CA_timestamps, CSUS_id_struct)
% mutualinfo_openfield_noCSUS
%
% Computes spatial mutual information outside CSUS periods and above a
% velocity threshold using a binary event / no-event formulation:
%
%   I_pos(x_i) = sum_{k=0}^1 P(k|x_i) * log2( P(k|x_i) / P(k) )
%   SI         = sum_i P(x_i) * I_pos(x_i)
%
% where:
%   P(x_i)   = probability of being in spatial bin i
%   P(k)     = probability of observing k events in a sample (k = 0 or 1)
%   P(k|x_i) = conditional probability of observing k events in bin i
%
% Output per day:
%   f.MI_YYYY_MM_DD = [numUnits x 1]
%
% Notes:
%   - Position is sampled on the included velocity timestamps.
%   - Calcium events are binarized per included sample:
%         1 = at least one event assigned to that sample
%         0 = no event
%   - Multiple spikes mapping to the same sample still count as 1.

fprintf('running mutualinfo_openfield_noCSUS\n');

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_cats   = fieldnames(CA_timestamps);
fields_CSUS   = fieldnames(CSUS_id_struct);

if numel(fields_spikes) ~= numel(fields_pos)
    error('spike and pos structures do not have the same number of fields');
end
if numel(fields_spikes) ~= numel(fields_cats)
    error('pos and timestamp structures do not have the same number of fields');
end
if numel(fields_spikes) ~= numel(fields_CSUS)
    error('spike and CSUS structures do not have the same number of fields');
end

mutualinfo_struct = struct();

for i = 1:numel(fields_spikes)

    fieldName_spikes = fields_spikes{i};
    peaks_time = spike_structure.(fieldName_spikes);

    fieldName_pos = fields_pos{i};
    pos = pos_structure.(fieldName_pos);

    fieldName_cats = fields_cats{i};
    curr_CA_timestamps = CA_timestamps.(fieldName_cats);

    fieldName_CSUS = fields_CSUS{i};
    CSUS_id = CSUS_id_struct.(fieldName_CSUS);

    idx = strfind(fieldName_spikes, '_');
    spikes_date = fieldName_spikes(idx(2)+1:end);

    fprintf('\nProcessing %s\n', spikes_date);

    if isempty(pos)
        warning('Empty pos for %s', spikes_date);
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
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
    pos = pos(uniqueIdxPos,:);

    % Remove duplicate CSUS timestamps
    [csus_t_unique, idxCSUS] = unique(CSUS_id(2,:), 'stable');
    csus_val_unique = CSUS_id(1, idxCSUS);

    % Velocity
    vel = ca_velocity(pos);
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

    % Interpolate CSUS labels and XY onto velocity timestamps
    interp_CSUS = interp1(csus_t_unique, csus_val_unique, vel_time, 'nearest', 0);
    interp_x = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    interp_y = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    % Keep only valid moving, non-CSUS samples
    validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ...
               ~isnan(interp_x) & ~isnan(interp_y);

    sample_time = vel_time(validIdx);
    sample_x    = interp_x(validIdx);
    sample_y    = interp_y(validIdx);

    goodpos = [sample_time, sample_x, sample_y];

    if isempty(goodpos)
        warning('No valid position samples for %s', spikes_date);
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue;
    end

    mintime = sample_time(1);
    maxtime = sample_time(end);

    if isempty(peaks_time) || size(peaks_time,1) < 1
        warning('No units found for %s', spikes_date);
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue;
    end

    % Spatial bins
    xmin = min(sample_x);
    xmax = max(sample_x);
    ymin = min(sample_y);
    ymax = max(sample_y);

    if xmax == xmin || ymax == ymin
        warning('Degenerate position range for %s', spikes_date);
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue;
    end

    xEdges = xmin:dim:(xmax + dim);
    yEdges = ymin:dim:(ymax + dim);

    % Assign every included sample to a spatial bin
    xBin = discretize(sample_x, xEdges);
    yBin = discretize(sample_y, yEdges);

    validBinSamples = ~isnan(xBin) & ~isnan(yBin);
    xBin = xBin(validBinSamples);
    yBin = yBin(validBinSamples);
    sample_time_valid = sample_time(validBinSamples);

    if isempty(sample_time_valid)
        warning('No valid binned samples for %s', spikes_date);
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue;
    end

    nXBins = numel(xEdges) - 1;
    nYBins = numel(yEdges) - 1;
    linBin = sub2ind([nXBins, nYBins], xBin, yBin);

    numunits = size(peaks_time,1);
    mutinfo = NaN(numunits,1);

    % Frame tolerance for assigning spikes to included samples
    if numel(sample_time_valid) > 1
        dt_samp = median(diff(sample_time_valid), 'omitnan');
    else
        dt_samp = NaN;
    end
    if isnan(dt_samp) || dt_samp <= 0
        dt_samp = 1/15; % fallback
    end
    maxAssignDist = dt_samp / 2;

    for k = 1:numunits

        currspikes = peaks_time(k,:);
        currspikes = currspikes(~isnan(currspikes));
        currspikes = currspikes(currspikes >= mintime & currspikes <= maxtime);

        if isempty(currspikes)
            mutinfo(k) = NaN;
            continue;
        end

        % Filter spikes by movement and CSUS exclusion
        spike_vel = interp1(vel_time, vel_mag, currspikes, 'linear', NaN);
        csus_currspikes = interp1(csus_t_unique, csus_val_unique, currspikes, 'nearest', 0);

        keepSpikes = ~isnan(spike_vel) & (spike_vel >= velthreshold) & (csus_currspikes == 0);
        highspeedspikes = currspikes(keepSpikes);

        if isempty(highspeedspikes)
            mutinfo(k) = NaN;
            continue;
        end

        % Assign spikes to nearest included sample, then binarize event/no-event per sample
        eventPerSample = false(numel(sample_time_valid),1);

        for s = 1:numel(highspeedspikes)
            [d, idxNearest] = min(abs(sample_time_valid - highspeedspikes(s)));
            if d <= maxAssignDist
                eventPerSample(idxNearest) = true;
            end
        end

        % Occupancy counts and event counts per spatial bin
        occCounts   = accumarray(linBin, 1, [nXBins*nYBins, 1], @sum, 0);
        eventCounts = accumarray(linBin, double(eventPerSample), [nXBins*nYBins, 1], @sum, 0);

        validOccBins = occCounts > 0;
        if ~any(validOccBins)
            mutinfo(k) = NaN;
            continue;
        end

        occCounts   = occCounts(validOccBins);
        eventCounts = eventCounts(validOccBins);

        totalSamples = sum(occCounts);
        totalEventSamples = sum(eventCounts);

        if totalSamples <= 0
            mutinfo(k) = NaN;
            continue;
        end

        Pxi = occCounts / totalSamples;              % P(x_i)
        P1_given_x = eventCounts ./ occCounts;       % P(event | x_i)
        P0_given_x = 1 - P1_given_x;                 % P(no event | x_i)

        P1 = totalEventSamples / totalSamples;       % P(event)
        P0 = 1 - P1;                                 % P(no event)

        % Compute binary mutual information
        Ipos = zeros(size(Pxi));

        use1 = (P1_given_x > 0) & (P1 > 0);
        use0 = (P0_given_x > 0) & (P0 > 0);

        Ipos(use1) = Ipos(use1) + P1_given_x(use1) .* log2(P1_given_x(use1) ./ P1);
        Ipos(use0) = Ipos(use0) + P0_given_x(use0) .* log2(P0_given_x(use0) ./ P0);

        mutinfo(k) = sum(Pxi .* Ipos);

    end

    mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
end

f = mutualinfo_struct;
end
