function f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, velthreshold, dim, CA_timestamps, num_times_to_run, ca_MI)
%finds mutual info for a bunch of cells


tic
cct = (class(calcium_traces) == 'struct')
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

      if length(peaks_time)>1

      fieldName_cats = fields_cats{i};
      curr_CA_timestamps = CA_timestamps.(fieldName_cats);

      index = strfind(fieldName_spikes, '_');
      spikes_date = fieldName_spikes(index(2)+1:end)

      fieldName_pos = fields_pos{i};
      fieldValue_pos = pos_structure.(fieldName_pos);
      pos = fieldValue_pos;

      index = strfind(fieldName_spikes, '_');
      pos_date = fieldName_spikes(index(2)+1:end)





      if (pos(1,1)-pos(end,1))./length(pos) < 1
        pos = convertpostoframe(pos, curr_CA_timestamps);
      end


      if length(peaks_time)>length(pos)
        peaks_time = peaks_time(1:length(pos));
      elseif length(peaks_time)<length(pos)
        pos = pos(1:length(peaks_time),:);
      end



      velthreshold = 2;
      vel = ca_velocity(pos);
      times = vel(2,:);
      velocities = vel(1,:);

      %want highspeedspikes
      % Thresholds
      velThreshold = 2; % cm/s
      timeThreshold = 1; % second
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
              validHighVelIndices = [validHighVelIndices, highVelIndices(i)];
          end
      end

      goodtime = pos(validHighVelIndices, 1);
      goodpos = pos(validHighVelIndices,:);
      all_highspeedspikes = peaks_time(:,validHighVelIndices);

      numunits = size(peaks_time,1);
      mutinfo = NaN(3,numunits);


      if numunits<=1
        warning('you have no cells and no spikes')
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
      else
          for k=1:numunits
                    currspikes = peaks_time(k,:);
                    if isnan(MI(k))==1
                      mutinfo(1, k) = NaN;
                      mutinfo(2, k) = NaN;
                      continue
                    else
                      highspeedspikes = all_highspeedspikes(k,:);
                    end

                    shuf = NaN(num_times_to_run,1);
                    %for l = 1:num_times_to_run
                    parfor l = 1:num_times_to_run

                          if isnan(MI(k))==0 && length(highspeedspikes)>1

                            shufff = highspeedspikes(randperm(length(highspeedspikes)))
                            [trace_mean occprob] = CA_normalizePosData_trace(shufff, goodpos, dim, 1.000);
                            shuf(l) = mutualinfo([trace_mean', occprob']);

                          else
                            shuf(l) = NaN;
                          end


                      end


                      topMI5 = floor(num_times_to_run*.95);
                      topMI1 = floor(num_times_to_run*.99);
                      shuf = sort(shuf);
                      if isnan(topMI5)==0
                        mutinfo(1, k) = shuf(topMI5);
                      else
                        mutinfo(1, k) = NaN;
                      end
                      if isnan(topMI1)==0
                        mutinfo(2, k) = shuf(topMI1);
                      else
                        mutinfo(2, k) = NaN;
                      end

                      [c index] = (min(abs(MI(k)-shuf)));
                      if isnan(index)==0
                        rank = index./length(shuf);
                        mutinfo(3, k) = rank;
                      else
                        mutinfo(3,k) = NaN;
                      end

                  end %ending the for loop for units


    mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
  end % ending if numunits<=1
end % ending   if length(peaks_time)>1

end % ending for i = 1:numel(fields_spikes)

 MI_trace_shuff = mutualinfo_struct;
% Get the current date and time as a string
currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');
% Create a filename with the timestamp
filename = ['results_MI_trace_shuff', currentDateTime, '.mat'];
% Save the output to the .mat file with the timestamped filename
save(filename, 'results_MI_trace_shuff');


f = mutualinfo_struct;
