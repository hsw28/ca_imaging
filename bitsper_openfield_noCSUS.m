function f = bitsper_openfield_noCSUS( ...
    spike_structure, pos_structure, velthreshold, dim, ...
    CA_timestamps, CSUS_id_struct)

fprintf('running bitsper_openfield_noCSUS\n');

f = struct();

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_cats   = fieldnames(CA_timestamps);
fields_CSUS   = fieldnames(CSUS_id_struct);

if numel(fields_spikes) ~= numel(fields_pos)
    error('spike and pos structures do not have the same number of fields');
end
if numel(fields_spikes) ~= numel(fields_cats)
    error('spike and timestamp structures do not have the same number of fields');
end
if numel(fields_spikes) ~= numel(fields_CSUS)
    error('spike and CSUS structures do not have the same number of fields');
end

set(0,'DefaultFigureVisible','off');

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
        f.(sprintf('MI_%s', spikes_date)) = NaN;
        continue;
    end

    if size(pos,2) > 3
        error('you are not using a fixed position');
    end

    if size(pos,2) > size(pos,1)
        pos = pos';
    end

    % same heuristic as shuffle version
    if ~isempty(pos) && numel(curr_CA_timestamps) > 1
        dt_pos = median(diff(pos(:,1)), 'omitnan');
        dt_ca  = median(diff(curr_CA_timestamps), 'omitnan');
        if ~isnan(dt_pos) & ~isnan(dt_ca) & dt_pos < 0.1 * dt_ca
            pos = convertpostoframe(pos, curr_CA_timestamps);
        end
    end

    % remove duplicate pos timestamps
    [~, uniqueIdxPos] = unique(pos(:,1), 'stable');
    pos = pos(uniqueIdxPos,:);

    % remove duplicate CSUS timestamps
    [csus_t_unique, idxCSUS] = unique(CSUS_id(2,:), 'stable');
    csus_val_unique = CSUS_id(1, idxCSUS);

    vel = ca_velocity(pos);
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

    interp_CSUS = interp1(csus_t_unique, csus_val_unique, vel_time, 'nearest', 0);
    interp_x = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    interp_y = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ...
               ~isnan(interp_x) & ~isnan(interp_y);

    goodpos = [vel_time(validIdx), interp_x(validIdx), interp_y(validIdx)];

    mintime = vel_time(1);
    maxtime = vel_time(end);

    if isempty(peaks_time) || size(peaks_time,1) < 1
        warning('No units found for %s', spikes_date);
        f.(sprintf('MI_%s', spikes_date)) = NaN;
        continue;
    end

    numunits = size(peaks_time,1);
    bitsper_info = NaN(2, numunits);

    for k = 1:numunits
        currspikes = peaks_time(k,:);
        currspikes = currspikes(~isnan(currspikes));
        currspikes = currspikes(currspikes >= mintime & currspikes <= maxtime);

        if isempty(currspikes) || isempty(goodpos)
            continue;
        end

        spike_vel = interp1(vel_time, vel_mag, currspikes, 'linear', NaN);
        csus_currspikes = interp1(csus_t_unique, csus_val_unique, currspikes, 'nearest', 0);

        keepSpikes = ~isnan(spike_vel) & (spike_vel >= velthreshold) & (csus_currspikes == 0);
        highspeedspikes = currspikes(keepSpikes);

        if numel(highspeedspikes) <= 1
            continue;
        end

        [rate, ~, ~, ~, ~, occprob] = CA_normalizePosData(highspeedspikes, goodpos, dim, 1.000);

        if size(occprob,1) < size(occprob,2)
            occprob = occprob';
        end
        if size(rate,1) < size(rate,2)
            rate = rate';
        end

        [bitsPerSpike, bitsPerSecond] = bits_per(rate, occprob);

        bitsper_info(1,k) = bitsPerSpike;
        bitsper_info(2,k) = bitsPerSecond;
    end

    f.(sprintf('MI_%s', spikes_date)) = bitsper_info;
end
end
