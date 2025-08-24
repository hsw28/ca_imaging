function f = mutualinfo_openfield_shuff_noCSUS(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, CSUS_id_struct, ca_MI, num_times_to_run)
%finds mutual info for a bunch of cells
%little did I know i already had code for this: ca_mutualinfo_openfield.m
%returns 95% cutoff, average MI, and rank of actual MI



fprintf('running mutualinfo_CSUS_trace_shuff')
fprintf('starting task')
fprintf('loading spikes')
fields_spikes = fieldnames(spike_structure);
fprintf('loading pos')
fields_pos = fieldnames(pos_structure);
try
fprintf('loading MI')
fields_MI = fieldnames(ca_MI);
catch
  fprintf('issue loading, what it is nobody knows')
  error('problem loading')
end
fprintf('loading TS')
fields_cats = fieldnames(CA_timestamps);
fprintf('all loaded')

fields_CSUS = fieldnames(CSUS_id_struct);


if numel(fields_spikes) ~= numel(fields_pos)
  fprintf('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

fprintf('starting loop')
for i = 1:numel(fields_spikes)

      fieldName_MI = fields_MI{i};
      fieldValue_MI = ca_MI.(fieldName_MI);
      MI = fieldValue_MI;

      fieldName_CSUS = fields_CSUS{i};
      CSUS_id = CSUS_id_struct.(fieldName_CSUS);

      fieldName_spikes = fields_spikes{i};
      fieldValue_spikes = spike_structure.(fieldName_spikes);
      peaks_time = fieldValue_spikes;

      index = strfind(fieldName_spikes, '_');
      spikes_date = fieldName_spikes(index(2)+1:end)

      fieldName_pos = fields_pos{i};
      fieldValue_pos = pos_structure.(fieldName_pos);
      pos = fieldValue_pos;


      fieldName_cats = fields_cats{i};
      curr_CA_timestamps = CA_timestamps.(fieldName_cats);



      fprintf('trimming date')
      tm = pos(:, 1);
      biggest = max(peaks_time(:));
      [minValue,closestIndex] = min(abs(biggest-tm));


      index = strfind(fieldName_spikes, '_');
      pos_date = fieldName_spikes(index(2)+1:end)

      if length(peaks_time) <5
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue
      end


            if (pos(1,1)-pos(end,1))./length(pos) < 1
              pos = convertpostoframe(pos, curr_CA_timestamps);
            end

      fprintf('trimming velocity')

      if abs(length(CSUS_id)-length(pos))>3
        error('your ID and pos have different lengths')
      end


      % Find all time points where CSUS_id > 0
      taskIdx = find(CSUS_id(1,:) > 0);
      CSUS_time = pos(taskIdx, 1);

      % Find all time points where CSUS_id <=0
      goodCSUS = find(CSUS_id(1,:) <= 0);

      % Now find the full range of indices to keep
      %get vel
      vel = ca_velocity(pos);
      vel_time = vel(2,:)';
      vel_mag  = vel(1,:)';

      % Interpolate CSUS labels to velocity timestamps
      interp_CSUS = interp1(CSUS_id(2,:), CSUS_id(1,:), vel_time, 'nearest', 0);


      [~, uniqueIdx] = unique(pos(:,1), 'stable');
      pos = pos(uniqueIdx, :);

      % Interpolate X and Y position to velocity timestamps too
      interp_x = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
      interp_y = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

      % Build mask based on velocity, CSUS period, and valid position values
      validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ...
                 ~isnan(interp_x) & ~isnan(interp_y);

      % Now build the posDat used downstream
        goodpos = [vel_time(validIdx), interp_x(validIdx), interp_y(validIdx)];


      mintime = vel(2,1);
      maxtime = vel(2,end);

      numunits = size(peaks_time,1);
      mutinfo = NaN(3,numunits);


      fprintf('done loading')

      if numunits<=1
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        warning('you have no spikes')
      else
          fprintf('about to go through units')
          for k=1:numunits

                currspikes = peaks_time(k,:);
                [c indexmin] = (min(abs(peaks_time(k,:)-mintime))); %
                [c indexmax] = (min(abs(peaks_time(k,:)-maxtime))); %
                currspikes = peaks_time(k,indexmin:indexmax);


                if isnan(MI(k))==1
                  mutinfo(1, k) = NaN;
                  mutinfo(2, k) = NaN;
                  continue
                end

                spike_vel = interp1(vel_time, vel_mag, currspikes, 'linear');
                csus_currspikes = interp1(CSUS_id(2,:), CSUS_id(1,:), currspikes, 'nearest', 0);
                highspeedspikes = currspikes((spike_vel >= velthreshold) & (csus_currspikes == 0));

                set(0,'DefaultFigureVisible', 'off');


                hertz = length(highspeedspikes)./(length(goodpos)/15); %% is this right or should it be 7.5??

                shuf = NaN(num_times_to_run,1);
                %for l = 1:num_times_to_run
                parfor l = 1:num_times_to_run
                      %fprintf('survived the great parfor loop trauma of jan 10')
                      if isnan(MI(k))==0 && length(highspeedspikes)>1

                        %code for random pos shuffle
                        %shuff_pos = goodpos;
                        %shuffled_indices = randperm(size(shuff_pos, 1));
                        %shuff_pos(:, 2:3) = shuff_pos(shuffled_indices, 2:3);
                        %end random post shuffle


                        % code for circular shift
                        pos_only = goodpos(:, 2:3);
                        time = goodpos(:, 1);
                        shift = randi([8 length(pos_only)], 1);
                        if rand < 0.5
                          shift = -shift;
                        end
                        shiftedData = circshift(pos_only, shift);
                        shuff_pos = [time, shiftedData];
                        % end circular shift

                        [rate totspikes totstime colorbar spikeprob occprob] = CA_normalizePosData(highspeedspikes,shuff_pos,dim, 1.000);

                        if (size(spikeprob,1)) < (size(spikeprob,2))
                          spikeprob = spikeprob';
                        end
                        if (size(occprob,1)) < (size(occprob,2))
                          occprob = occprob';
                        end

                        shuf(l) = mutualinfo([spikeprob, occprob]);
                      else
                        shuf(l) = NaN;
                      end
                  end


                  % remove NaNs from the shuffle distribution
                shOK = shuf(~isnan(shuf));

                % 95-th percentile and mean of shuffle
                p95cut  = prctile(shOK,95);
                muShuff = mean(shOK);

                % p-value (upper-tail) and percentile rank
                pVal = mean(shOK >= MI(k));
                perc = mean(shOK <= MI(k));     % 0 = worst, 1 = best

                % store
                mutinfo(1,k) = p95cut;     % 95-th-percentile threshold
                mutinfo(2,k) = muShuff;    % shuffle mean
                mutinfo(3,k) = perc;       % percentile of actual MI
                mutinfo(4,k) = hertz;      % firing rate

          end

    fprintf('assigning MI')
    mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
  end
end


%{
  results_MI_shuff = mutualinfo_struct;
  fprintf('saving\n');
  MI_trace_shuff = mutualinfo_struct;
  fprintf('Get the current date and time as a string\n');
  currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');
  fprintf('Create a filename with the timestamp\n');
  filename = ['results_MI_shuff_', currentDateTime, '.mat'];
  fprintf('Save the output to the .mat file with the timestamped filename\n');
  save(filename, 'results_MI_shuff');
  fprintf('Save is a success\n');



results_MI_shift = mutualinfo_struct;
fprintf('saving\n');
MI_trace_shift = mutualinfo_struct;
fprintf('Get the current date and time as a string\n');
currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');
fprintf('Create a filename with the timestamp\n');
filename = ['results_MI_shift_', currentDateTime, '.mat'];
fprintf('Save the output to the .mat file with the timestamped filename\n');
save(filename, 'results_MI_shift');
fprintf('Save is a success\n');
%}

  f = mutualinfo_struct;

end
