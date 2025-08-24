function mutualinfo_struct = MI_control_matchSpikes_post(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, CSUS_id_struct, ca_MI, nIter)
% MI_control_matchSpikes
% Computes shuffled mutual information (MI) as a control  excluding task periods and
% removing matched number of cells to task.
% Returns a struct with per-cell control stats.
%
% Inputs:
%   spike_structure – struct of spike timestamps per day
%   pos_structure   – struct of position data per day
%   velthreshold    – speed threshold (e.g., 4 cm/s)
%   dim             – bin resolution for CA_normalizePosData (e.g., 2.5)
%   CA_timestamps   – imaging frame timestamps per day
%   CSUS_id_struct  – struct of trial indicator (2×N; 1,: > 0 = task, 2,: = time)
%   ca_MI           – struct of actual MI values per day (1×nCells)
%   nIter           – number of shuffles (default = 500)
%
% Output:
%   mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = [p95; mean; percentile] × nCells

if nargin < 8, nIter = 500; end

  set(0,'DefaultFigureVisible', 'off');


fields_spikes = fieldnames(spike_structure);
fields_pos = fieldnames(pos_structure);
try
    fprintf('loading MI\n')
    fields_MI = fieldnames(ca_MI);
catch
    fprintf('issue loading, what it is nobody knows\n')
    error('problem loading')
end
fields_cats = fieldnames(CA_timestamps);
fields_CSUS = fieldnames(CSUS_id_struct);

if numel(fields_spikes) ~= numel(fields_pos)
    fprintf('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days\n')
    error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

fprintf('starting loop\n')
mutualinfo_struct = struct();

% only use the last 3 entries (An, An-1, An-2)
fields_spikes = fields_spikes(end-2:end);
fields_pos    = fields_pos(end-2:end);
fields_cats   = fields_cats(end-2:end);
fields_CSUS   = fields_CSUS(end-2:end);
fields_MI     = fields_MI(end-2:end);

for i = 1:numel(fields_spikes)
    fieldName_MI = fields_MI{i}
    MI = ca_MI.(fieldName_MI);

    fieldName_CSUS = fields_CSUS{i};
    CSUS_id = CSUS_id_struct.(fieldName_CSUS);

    fieldName_spikes = fields_spikes{i};
    peaks_time = spike_structure.(fieldName_spikes);

    index = strfind(fieldName_spikes, '_');
    spikes_date = fieldName_spikes(index(2)+1:end);

    fieldName_pos = fields_pos{i};
    pos = pos_structure.(fieldName_pos);

    fieldName_cats = fields_cats{i};
    curr_CA_timestamps = CA_timestamps.(fieldName_cats);

    %get vel
    vel = ca_velocity(pos);
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

    % Interpolate CSUS labels to velocity timestamps
    interp_CSUS = interp1(CSUS_id(2,:), CSUS_id(1,:), vel_time, 'nearest', 0);


    [~, uniqueIdx] = unique(pos(:,1), 'stable');
    pos = pos(uniqueIdx, :);

    % Interpolate x and y to velocity timestamps
    interp_x = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    interp_y = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);


    % Build mask based on velocity, CSUS period, and valid position values
    validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ...
               ~isnan(interp_x) & ~isnan(interp_y);

    % Now build the posDat used downstream
    posDat = [vel_time(validIdx), interp_x(validIdx), interp_y(validIdx)];




    % Initialize shuffle matrix
    nCells = size(peaks_time,1);
    controlMI = nan(nCells, nIter);

    %for iIter = 1:nIter
    parfor iIter = 1:nIter
        for k = 1:nCells
            spk = peaks_time(k,:);
            spk = spk(~isnan(spk) & spk > 0);

            spike_vel = interp1(vel_time, vel_mag, spk, 'linear');
            csus_spk = interp1(CSUS_id(2,:), CSUS_id(1,:), spk, 'nearest', 0);

            valid_spikes = spk((spike_vel >= velthreshold) & (csus_spk == 0));

            all_spikes = spk(spike_vel >= velthreshold);

            CSUS_spikes = spk((spike_vel >= velthreshold) & (csus_spk <= 6) & (csus_spk > 0)); %always elim these

            nonCSUS_spikes = setdiff(all_spikes, CSUS_spikes);

            non_post = spk((spike_vel >= velthreshold) & (csus_spk <= 6));

            post_spikes = spk((spike_vel >= velthreshold) & (csus_spk > 6));



            %want to elim CSUS spikes for all
            %select spikes from any time
            %number of spikes should be total-post
            if numel(valid_spikes) >= 1

                resampled = randsample(nonCSUS_spikes, numel(all_spikes)-numel(post_spikes)-numel(CSUS_spikes), false);


                [~, ~, ~, ~, spikeprob, occprob] = CA_normalizePosData(resampled, posDat, dim, 1.0);

                if size(spikeprob,1) < size(spikeprob,2)
                    spikeprob = spikeprob';
                end
                if size(occprob,1) < size(occprob,2)
                    occprob = occprob';
                end

                controlMI(k, iIter) = mutualinfo([spikeprob, occprob]);
            end
        end
    end

    % Summarize across iterations
    mutinfo = nan(3, nCells);
    for k = 1:nCells
        shuffVals = controlMI(k, :);
        shuffVals = shuffVals(~isnan(shuffVals));
        if isempty(shuffVals), continue; end

        p95cut = prctile(shuffVals, 95);
        muShuff = mean(shuffVals);
        perc = sum(MI(k) > shuffVals) / numel(shuffVals);

        mutinfo(:,k) = [p95cut; muShuff; perc];
    end


    mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';

    % Wait briefly to ensure parfor finishes cleanly
  %  pause(2)
  %  p = gcp('nocreate');
  %  if ~isempty(p)
  %      delete(p)
  %  end


end
end
