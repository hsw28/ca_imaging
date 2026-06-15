function [per_spike, per_sec] = bitsper_openfield_shuff_noCSUS( ...
    spike_structure, pos_structure, velthreshold, dim, ...
    CA_timestamps, CSUS_id_struct, ca_bitsper, num_times_to_run)
% bitsper_openfield_shuff_noCSUS
% Computes shuffled bits/spike and bits/sec outside CSUS periods and while
% above velocity threshold.
%
% Output per day:
%   per_spike.bitsPerSpike_YYYY_MM_DD = [numUnits x 4]
%   per_sec.bitsPerSec_YYYY_MM_DD     = [numUnits x 4]
%
% Columns:
%   1 = 95th percentile of shuffle
%   2 = mean of shuffle
%   3 = rank of observed value relative to shuffle
%   4 = firing rate (Hz) in included periods

fprintf('running bitsper_openfield_shuff_noCSUS\n');

per_spike = struct();
per_sec   = struct();

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_MI     = fieldnames(ca_bitsper);
fields_cats   = fieldnames(CA_timestamps);
fields_CSUS   = fieldnames(CSUS_id_struct);

if numel(fields_spikes) ~= numel(fields_pos)
    error('spike and pos structures do not have the same number of fields');
end

if numel(fields_spikes) ~= numel(fields_MI)
    error('spike and ca_bitsper structures do not have the same number of fields');
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

    fieldName_MI = fields_MI{i};
    MI = ca_bitsper.(fieldName_MI);

    fieldName_cats = fields_cats{i};
    curr_CA_timestamps = CA_timestamps.(fieldName_cats);

    fieldName_CSUS = fields_CSUS{i};
    CSUS_id = CSUS_id_struct.(fieldName_CSUS);

    % Extract date from spike field name
    idx = strfind(fieldName_spikes, '_');
    spikes_date = fieldName_spikes(idx(2)+1:end);

    fprintf('\nProcessing %s\n', spikes_date);

    if isempty(pos)
        warning('Empty pos for %s', spikes_date);
        per_spike.(sprintf('bitsPerSpike_%s', spikes_date)) = NaN;
        per_sec.(sprintf('bitsPerSec_%s', spikes_date)) = NaN;
        continue;
    end

    if size(pos,2) > 3
        error('you are not using a fixed position');
    end

    if size(pos,2) > size(pos,1)
        pos = pos';
    end

    % Optional conversion if pos is not in timestamp units
    % This condition is still heuristic; replace if you have a better test
    if ~isempty(pos) && numel(curr_CA_timestamps) > 1
        dt_pos = median(diff(pos(:,1)), 'omitnan');
        dt_ca  = median(diff(curr_CA_timestamps), 'omitnan');
        if ~isnan(dt_pos) & ~isnan(dt_ca) & dt_pos < 0.1 * dt_ca
            pos = convertpostoframe(pos, curr_CA_timestamps);
        end
    end

    % Remove duplicate pos timestamps before interpolation
    [~, uniqueIdxPos] = unique(pos(:,1), 'stable');
    pos = pos(uniqueIdxPos,:);

    % Remove duplicate CSUS timestamps before interpolation
    [csus_t_unique, idxCSUS] = unique(CSUS_id(2,:), 'stable');
    csus_val_unique = CSUS_id(1, idxCSUS);

    % Velocity from position
    vel = ca_velocity(pos);
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

    % Interpolate CSUS labels and position to velocity timestamps
    interp_CSUS = interp1(csus_t_unique, csus_val_unique, vel_time, 'nearest', 0);
    interp_x = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    interp_y = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    % Keep only high-velocity, non-CSUS, valid-position samples
    validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ...
               ~isnan(interp_x) & ~isnan(interp_y);

    goodpos = [vel_time(validIdx), interp_x(validIdx), interp_y(validIdx)];

    mintime = vel_time(1);
    maxtime = vel_time(end);

    if isempty(peaks_time) || size(peaks_time,1) < 1
        warning('No units found for %s', spikes_date);
        per_spike.(sprintf('bitsPerSpike_%s', spikes_date)) = NaN;
        per_sec.(sprintf('bitsPerSec_%s', spikes_date)) = NaN;
        continue;
    end

    numunits = size(peaks_time,1);

    % Keep same output format as your original code:
    % rows = metric type, cols = units, then transpose at assignment
    bitsPspike = NaN(4, numunits);
    bitsPsec   = NaN(4, numunits);

    fprintf('about to go through units\n');

    for k = 1:numunits

        if size(MI,1) < 2 || size(MI,2) < k
            warning('MI dimensions do not match expected 2 x numUnits for %s', spikes_date);
            continue;
        end

        if isnan(MI(1,k)) || isnan(MI(2,k))
            bitsPspike(:,k) = NaN;
            bitsPsec(:,k)   = NaN;
            continue;
        end

        % Clean spike vector
        currspikes = peaks_time(k,:);
        currspikes = currspikes(~isnan(currspikes));
        currspikes = currspikes(currspikes >= mintime & currspikes <= maxtime);

        if isempty(currspikes) || isempty(goodpos)
            bitsPspike(:,k) = NaN;
            bitsPsec(:,k)   = NaN;
            continue;
        end

        % Filter spikes by velocity and CSUS exclusion
        spike_vel = interp1(vel_time, vel_mag, currspikes, 'linear', NaN);
        csus_currspikes = interp1(csus_t_unique, csus_val_unique, currspikes, 'nearest', 0);

        keepSpikes = ~isnan(spike_vel) & (spike_vel >= velthreshold) & (csus_currspikes == 0);
        highspeedspikes = currspikes(keepSpikes);

        if numel(highspeedspikes) <= 1
            bitsPspike(:,k) = NaN;
            bitsPsec(:,k)   = NaN;
            continue;
        end

        % Firing rate in included period
        occDur = size(goodpos,1) / 15;   % keep your original 15 Hz assumption
        if occDur <= 0
            hertz = NaN;
        else
            hertz = numel(highspeedspikes) / occDur;
        end

        shuf = NaN(num_times_to_run, 2);

        pos_only = goodpos(:,2:3);
        timevec  = goodpos(:,1);
        nPos = size(pos_only,1);

        parfor l = 1:num_times_to_run

            if nPos < 9
                shuf(l,:) = [NaN, NaN];
                continue;
            end

            % Circular shift by at least 8 samples and at most nPos-1
            shift = randi([8, nPos-1], 1);
            if rand < 0.5
                shift = -shift;
            end

            shiftedData = circshift(pos_only, shift, 1);
            shuff_pos   = [timevec, shiftedData];

            [rate, ~, ~, ~, ~, occprob] = CA_normalizePosData(highspeedspikes, shuff_pos, dim, 1.000);

            if size(occprob,1) < size(occprob,2)
                occprob = occprob';
            end
            if size(rate,1) < size(rate,2)
                rate = rate';
            end

            [bitsPerSpike, bitsPerSecond] = bits_per(rate, occprob);
            shuf(l,:) = [bitsPerSpike, bitsPerSecond];
        end

        shOK1 = shuf(~isnan(shuf(:,1)), 1);
        shOK2 = shuf(~isnan(shuf(:,2)), 2);

        if isempty(shOK1) || isempty(shOK2)
            bitsPspike(:,k) = NaN;
            bitsPsec(:,k)   = NaN;
            continue;
        end

        obscurrspikes = MI(1,k);
        obsSec        = MI(2,k);

        bitsPspike(1,k) = prctile(shOK1, 95);
        bitsPsec(1,k)   = prctile(shOK2, 95);

        bitsPspike(2,k) = mean(shOK1);
        bitsPsec(2,k)   = mean(shOK2);

        bitsPspike(3,k) = mean(shOK1 >= obscurrspikes);
        bitsPsec(3,k)   = mean(shOK2 >= obsSec);

        bitsPspike(4,k) = hertz;
        bitsPsec(4,k)   = hertz;
    end

    per_spike.(sprintf('bitsPerSpike_%s', spikes_date)) = bitsPspike';
    per_sec.(sprintf('bitsPerSec_%s', spikes_date))     = bitsPsec';
end
end
