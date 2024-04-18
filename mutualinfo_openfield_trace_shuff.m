function f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, velthreshold, dim, CA_timestamps, num_times_to_run, ca_MI)
%finds mutual info for a bunch of cells


fprintf('running mutualinfo_openfield_trace_shuff')

fprintf('starting function');
cct = (class(calcium_traces) == 'struct');
if cct(1) == 0
  allvariables = load('allvariables2.mat');
  calcium_traces = allvariables.all_traces;
end


spike_structure = calcium_traces;

fields_spikes = fieldnames(spike_structure);
fields_pos = fieldnames(pos_structure);
fields_MI = fieldnames(ca_MI);
fields_cats = fieldnames(CA_timestamps);


if numel(fields_spikes) ~= numel(fields_pos)
  warning('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end


for i = 1:numel(fields_spikes)
      fieldName_spikes = fields_spikes{i};
      fieldValue_spikes = spike_structure.(fieldName_spikes);
      peaks_time = fieldValue_spikes;

      fieldName_MI = fields_MI{i};
      fieldValue_MI = ca_MI.(fieldName_MI);
      MI = fieldValue_MI;


      fieldName_cats = fields_cats{i};
      curr_CA_timestamps = CA_timestamps.(fieldName_cats);

      index = strfind(fieldName_spikes, '_');
      spikes_date = fieldName_spikes(index(2)+1:end)

      fieldName_pos = fields_pos{i};
      fieldValue_pos = pos_structure.(fieldName_pos);
      pos = fieldValue_pos;

      index = strfind(fieldName_spikes, '_');

      if all(isnan(MI(1,:))) == 1
          fprintf('All MI are NaN');
          mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
          continue
      end


      if length(peaks_time) <3
        fprintf('no spikes, continuing')
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue
      end

      fprintf('fixing positions') %%%
      if (pos(1,1)-pos(end,1))./length(pos) < 1
        pos = convertpostoframe(pos, curr_CA_timestamps);
      end

      if length(peaks_time)>length(pos)
        peaks_time = peaks_time(1:length(pos));
      elseif length(peaks_time)<length(pos)
        pos = pos(1:length(peaks_time),:);
      end

      pos = smoothpos(pos);

      vel = ca_velocity(pos);
      times = vel(2,:);
      velocities = vel(1,:);


      fprintf('thresholding spikes for vel') %%%
      %want highspeedspikes
      % Thresholds
      velThreshold = velthreshold; % cm/s
      timeThreshold = 1/15; % second
      % Find indices where velocity is greater than the threshold
      highVelIndices = find(velocities >= velThreshold);
      % Find indices where velocity is less than or equal to the threshold
      lowVelIndices = find(velocities < velThreshold);
      % Filter out high velocity indices that are too close to low velocities
      validHighVelIndices = [];
      for ii = 1:length(highVelIndices)
          highVelTime = times(highVelIndices(ii));
          % Find the closest low velocity time
          [~, closestLowVelIndex] = min(abs(highVelTime - times(lowVelIndices)));
          closestLowVelTime = times(lowVelIndices(closestLowVelIndex));

          % Check if the high velocity time is more than 1 second away from the closest low velocity time
          if abs(highVelTime - closestLowVelTime) > timeThreshold
              validHighVelIndices = [validHighVelIndices, highVelIndices(ii)];
          end
      end

      goodtime = pos(validHighVelIndices, 1);
      goodpos = pos(validHighVelIndices,:);
      all_highspeedspikes = peaks_time(:,validHighVelIndices);

      numunits = size(peaks_time,1);
      mutinfo = NaN(3,numunits);



      if (numunits)<=3 | isnan(numunits)== 1
        fprintf('you have no cells and no spikes')
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
      else
        fprintf('going through units') %%%
          for k=1:numunits
                    if isnan(MI(k))==1
                      fprintf('check 1') %%%
                      mutinfo(1, k) = NaN;
                      mutinfo(2, k) = NaN;
                      continue
                    else
                      fprintf('check 2') %%%
                      highspeedspikes = all_highspeedspikes(k,:);
                    end


                    currspikes = peaks_time(k,:);

                    shuf = NaN(num_times_to_run,1);
                    %for l = 1:num_times_to_run
                    fprintf('check 3') %%%
                    parfor l = 1:num_times_to_run

                          if isnan(MI(k))==0 && length(highspeedspikes)>1
                            fprintf('check 4') %%%
                            %shuff_pos = goodpos;
                            %shuffled_indices = randperm(size(shuff_pos, 1));
                            % Apply the shuffled indices to the first two columns
                            %shuff_pos(:, 2:3) = shuff_pos(shuffled_indices, 1:2);

                            % code for circular shift
                            pos_only = goodpos(:, 2:3);
                            time = goodpos(:, 1);
                            shift = randi([8 length(pos_only)], 1);
                            if rand < 0.5
                              shift = -shift;
                            end
                            shiftedData = circshift(pos_only, shift);
                            shuff_pos = [time, shiftedData];
                            %% end circular shift

                            fprintf('check 5') %%%
                            [trace_mean occprob] = CA_normalizePosData_trace(highspeedspikes, shuff_pos, dim, 1.000);
                            if (size(trace_mean,1)) < (size(trace_mean,2))
                              trace_mean = trace_mean';
                            end
                            if (size(occprob,1)) < (size(occprob,2))
                              occprob = occprob';
                            end
                            fprintf('check 6') %%%
                            shuf(l) = mutualinfo([trace_mean, occprob]);
                            fprintf('check 7') %%%

                          else
                            shuf(l) = NaN;
                            fprintf('check 8') %%%
                          end
                          fprintf('check 9') %%%

                      end

                      fprintf('check 10') %%%
                      topMI5 = floor(num_times_to_run*.95);
                      topMI1 = floor(num_times_to_run*.99);
                      fprintf('check 11') %%%
                      shuf = sort(shuf);
                      if isnan(topMI5)==0
                        mutinfo(1, k) = shuf(topMI5);
                      else
                        mutinfo(1, k) = NaN;
                      end
                      if isnan(topMI1)==0
                        mutinfo(2, k) = nanmean(shuf);
                      else
                        mutinfo(2, k) = NaN;
                      end
                      fprintf('check 12') %%%
                      [c index] = (min(abs(MI(k)-shuf)));
                      if isnan(index)==0
                        rank = index./length(shuf);
                        mutinfo(3, k) = rank;
                      else
                        mutinfo(3,k) = NaN;
                      end
                      fprintf('check 13') %%%
                  end %ending the for loop for units

                  %fprintf('assigning MI')
    mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
  end % ending if numunits<=1



end % ending for i = 1:numel(fields_spikes)


fprintf('saving')
results_MI_trace_shuff = mutualinfo_struct;
fprintf('Get the current date and time as a string')
currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');
fprintf('Create a filename with the timestamp')
filename = ['results_MI_trace_shuff_', currentDateTime, '.mat'];
fprintf('Save the output to the .mat file with the timestamped filename')
save(filename, 'results_MI_trace_shuff');
print('save is a success')



f = mutualinfo_struct
