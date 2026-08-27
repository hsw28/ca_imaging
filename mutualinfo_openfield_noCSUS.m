function f = mutualinfo_openfield_noCSUS(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, CSUS_id_struct)
% Computes spatial MI outside CSUS periods using the full joint distribution:
%
%   I(X;K) = sum_x sum_k p(x,k) * log2( p(x,k) / (p(x)*p(k)) )
%
% where:
%   X = spatial bin
%   K = spike/event count per included sample
%
% Count states are capped at 3 by default: 0,1,2,3+

tic

maxK = 1;   % cap counts per sample to reduce sparsity

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_cats   = fieldnames(CA_timestamps);
fields_CSUS   = fieldnames(CSUS_id_struct);

if numel(fields_spikes) ~= numel(fields_pos)
    error('your spike and pos structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

if numel(fields_pos) ~= numel(fields_cats)
    error('your pos and timestamp structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

mutualinfo_struct = struct();

for i = 1:numel(fields_spikes)

    fieldName_spikes = fields_spikes{i}
    peaks_time = spike_structure.(fieldName_spikes);

    fieldName_CSUS = fields_CSUS{i};
    CSUS_id = CSUS_id_struct.(fieldName_CSUS);

    index = strfind(fieldName_spikes, '_');
    spikes_date = fieldName_spikes(index(2)+1:end);

    fieldName_pos = fields_pos{i};
    pos = pos_structure.(fieldName_pos);

    if size(pos,2) > 3
        error('you are not using a fixed position')
    end

    if size(pos,2) > size(pos,1)
        pos = pos';
    end

    fieldName_cats = fields_cats{i};
    curr_CA_timestamps = CA_timestamps.(fieldName_cats);

    if (pos(1,1)-pos(end,1))./length(pos) < 1
        pos = convertpostoframe(pos, curr_CA_timestamps);
    end

    mutinfo = NaN(size(peaks_time,1),1);

    % velocity
    vel = ca_velocity(pos);
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

    % Interpolate CSUS labels to velocity timestamps
    interp_CSUS = interp1(CSUS_id(2,:), CSUS_id(1,:), vel_time, 'nearest', 0);

    [~, uniqueIdx] = unique(pos(:,1), 'stable');
    pos = pos(uniqueIdx, :);

    % Interpolate X and Y position to velocity timestamps
    interp_x = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    interp_y = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    % Keep only moving, non-CSUS, valid-position samples
    validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ...
               ~isnan(interp_x) & ~isnan(interp_y);

    sample_time = vel_time(validIdx);
    sample_x    = interp_x(validIdx);
    sample_y    = interp_y(validIdx);

    if isempty(sample_time)
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
        continue
    end

    % spatial bins
    xmin = min(sample_x);
    xmax = max(sample_x);
    ymin = min(sample_y);
    ymax = max(sample_y);

    if xmax == xmin || ymax == ymin
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
        continue
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
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
        continue
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
        dt_samp = 1/7.5;
    end
    maxAssignDist = dt_samp / 2;

    mintime = sample_time_valid(1);
    maxtime = sample_time_valid(end);

    numunits = size(peaks_time,1);

    if numunits <= 1
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        warning('you have no spikes')
    else
        for k = 1:numunits

            currspikes = peaks_time(k,:);
            currspikes = currspikes(~isnan(currspikes));
            currspikes = currspikes(currspikes >= mintime & currspikes <= maxtime);

            spike_vel = interp1(vel_time, vel_mag, currspikes, 'linear', NaN);
            csus_currspikes = interp1(CSUS_id(2,:), CSUS_id(1,:), currspikes, 'nearest', 0);

            highspeedspikes = currspikes((spike_vel >= velthreshold) & ...
                                         (csus_currspikes == 0) & ...
                                         ~isnan(spike_vel));

            if ~isempty(highspeedspikes)
                mutinfo(k) = count_mi_from_spikes( ...
                    highspeedspikes, sample_time_valid, linBin_true, ...
                    maxAssignDist, nXBins, nYBins, maxK);
            else
                mutinfo(k) = NaN;
            end
        end

        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
    end
end

f = mutualinfo_struct;
end


function mi = count_mi_from_spikes(spike_times, sample_time_valid, linBin_true, maxAssignDist, nXBins, nYBins, maxK)
% Full joint-distribution MI with count per sample

if isempty(spike_times) || isempty(sample_time_valid) || isempty(linBin_true)
    mi = NaN;
    return;
end

nSamples = numel(sample_time_valid);
countPerSample = zeros(nSamples,1);

for s = 1:numel(spike_times)
    [d, idxNearest] = min(abs(sample_time_valid - spike_times(s)));
    if d <= maxAssignDist
        countPerSample(idxNearest) = countPerSample(idxNearest) + 1;
    end
end

% cap counts: 0,1,2,3+
if ~isempty(maxK)
    countPerSample(countPerSample > maxK) = maxK;
end

[countVals, ~, kIdx] = unique(countPerSample); %#ok<ASGLU>
nK = numel(countVals);

jointCounts = accumarray([linBin_true(:), kIdx(:)], 1, [nXBins*nYBins, nK], @sum, 0);

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
