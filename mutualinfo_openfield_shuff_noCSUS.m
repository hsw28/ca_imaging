function f = mutualinfo_openfield_shuff_noCSUS( ...
    spike_structure, pos_structure, velthreshold, dim, ...
    CA_timestamps, CSUS_id_struct, ca_MI, num_times_to_run)
% mutualinfo_openfield_shuff_noCSUS
%
% Computes shuffled spatial mutual information outside CSUS periods and
% above a velocity threshold using the full joint-distribution count MI:
%
%   I(X;K) = sum_x sum_k p(x,k) * log2( p(x,k) / (p(x)*p(k)) )
%
% where:
%   X = spatial bin
%   K = event count per included sample
%
% Output per day:
%   f.MI_YYYY_MM_DD = [numUnits x 4]
%
% Columns:
%   1 = 95th percentile of shuffle
%   2 = mean of shuffle
%   3 = percentile rank of observed MI relative to shuffle
%   4 = firing rate (Hz) in included periods

fprintf('running mutualinfo_openfield_shuff_noCSUS\n');

maxK = 1;   % cap counts per sample: 0,1,2,3+

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_MI     = fieldnames(ca_MI);
fields_cats   = fieldnames(CA_timestamps);
fields_CSUS   = fieldnames(CSUS_id_struct);

if numel(fields_spikes) ~= numel(fields_pos)
    error('spike and pos structures do not have the same number of fields');
end
if numel(fields_spikes) ~= numel(fields_MI)
    error('spike and ca_MI structures do not have the same number of fields');
end
if numel(fields_spikes) ~= numel(fields_cats)
    error('spike and timestamp structures do not have the same number of fields');
end
if numel(fields_spikes) ~= numel(fields_CSUS)
    error('spike and CSUS structures do not have the same number of fields');
end

mutualinfo_struct = struct();

for i = 1:numel(fields_spikes)

    fieldName_spikes = fields_spikes{i}
    peaks_time = spike_structure.(fieldName_spikes);

    fieldName_pos = fields_pos{i};
    pos = pos_structure.(fieldName_pos);

    fieldName_MI = fields_MI{i};
    MI = ca_MI.(fieldName_MI);

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

    % Interpolate CSUS and XY onto velocity timestamps
    interp_CSUS = interp1(csus_t_unique, csus_val_unique, vel_time, 'nearest', 0);
    interp_x = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    interp_y = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    % Keep only moving, non-CSUS, valid XY samples
    validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ...
               ~isnan(interp_x) & ~isnan(interp_y);

    sample_time = vel_time(validIdx);
    sample_x    = interp_x(validIdx);
    sample_y    = interp_y(validIdx);

    if isempty(sample_time)
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

    % Bin edges from included true positions
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

    % Assign every included sample to a true spatial bin
    xBin_true = discretize(sample_x, xEdges);
    yBin_true = discretize(sample_y, yEdges);

    validBinSamples = ~isnan(xBin_true) & ~isnan(yBin_true);
    sample_time_valid = sample_time(validBinSamples);
    sample_x_valid    = sample_x(validBinSamples);
    sample_y_valid    = sample_y(validBinSamples);
    xBin_true         = xBin_true(validBinSamples);
    yBin_true         = yBin_true(validBinSamples);

    if isempty(sample_time_valid)
        warning('No valid binned samples for %s', spikes_date);
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue;
    end

    nXBins = numel(xEdges) - 1;
    nYBins = numel(yEdges) - 1;
    linBin_true = sub2ind([nXBins, nYBins], xBin_true, yBin_true);

    numunits = size(peaks_time,1);
    mutinfo = NaN(4, numunits);

    % Tolerance for assigning spikes to included samples
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

    for k = 1:numunits

        if numel(MI) < k || isnan(MI(k))
            mutinfo(:,k) = NaN;
            continue;
        end

        currspikes = peaks_time(k,:);
        currspikes = currspikes(~isnan(currspikes));
        currspikes = currspikes(currspikes >= mintime & currspikes <= maxtime);

        if isempty(currspikes)
            mutinfo(:,k) = NaN;
            continue;
        end

        % Filter spikes by movement and CSUS exclusion
        spike_vel = interp1(vel_time, vel_mag, currspikes, 'linear', NaN);
        csus_currspikes = interp1(csus_t_unique, csus_val_unique, currspikes, 'nearest', 0);

        keepSpikes = ~isnan(spike_vel) & (spike_vel >= velthreshold) & (csus_currspikes == 0);
        highspeedspikes = currspikes(keepSpikes);

        if isempty(highspeedspikes)
            mutinfo(:,k) = NaN;
            continue;
        end

        % Count per included sample using true sample times
        countPerSample = zeros(numel(sample_time_valid),1);

        for s = 1:numel(highspeedspikes)
            [d, idxNearest] = min(abs(sample_time_valid - highspeedspikes(s)));
            if d <= maxAssignDist
                countPerSample(idxNearest) = countPerSample(idxNearest) + 1;
            end
        end

        % firing rate in included periods
        if numel(sample_time_valid) > 1
            occDur = median(diff(sample_time_valid), 'omitnan') * numel(sample_time_valid);
        else
            occDur = NaN;
        end
        if isnan(occDur) || occDur <= 0
            hertz = NaN;
        else
            hertz = numel(highspeedspikes) / occDur;
        end

        % Observed MI
        observedMI = count_mi_from_counts(countPerSample, linBin_true, nXBins, nYBins, maxK);
        if isnan(observedMI)
            mutinfo(:,k) = NaN;
            continue;
        end

        % Shuffle by circularly shifting XY across included samples
        shuf = NaN(num_times_to_run,1);

        pos_only = [sample_x_valid, sample_y_valid];
        nPos = size(pos_only,1);

        parfor l = 1:num_times_to_run

            if nPos < 9
                shuf(l) = NaN;
                continue;
            end

            shift = randi([8, nPos-1], 1);
            if rand < 0.5
                shift = -shift;
            end

            shiftedXY = circshift(pos_only, shift, 1);

            xBin_shuf = discretize(shiftedXY(:,1), xEdges);
            yBin_shuf = discretize(shiftedXY(:,2), yEdges);

            validShuf = ~isnan(xBin_shuf) & ~isnan(yBin_shuf);
            if ~any(validShuf)
                shuf(l) = NaN;
                continue;
            end

            linBin_shuf = sub2ind([nXBins, nYBins], xBin_shuf(validShuf), yBin_shuf(validShuf));
            countShuf = countPerSample(validShuf);

            shuf(l) = count_mi_from_counts(countShuf, linBin_shuf, nXBins, nYBins, maxK);
        end

        shOK = shuf(~isnan(shuf));
        if isempty(shOK)
            mutinfo(:,k) = NaN;
            continue;
        end

        mutinfo(1,k) = prctile(shOK, 95);
        mutinfo(2,k) = mean(shOK);
        mutinfo(3,k) = mean(shOK <= observedMI);
        mutinfo(4,k) = hertz;
    end

    mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
end

f = mutualinfo_struct;
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
