function f = mutualinfo_openfield_trace(spike_structure, pos_structure, velthreshold, dim, CA_timestamps)
%finds mutual info for a bunch of cells
%idea modified from  HippoBellum: acute cerebellar modulation alters hippocampal dynamics and function
%The spatial information content of a cell's activity (in bits) was defined as: I= Σ୧(λ୧⁄λ) x logଶ(λ୧⁄λ) x p୧, where λ is the mean calcium activity ( Σ୧ p୧ x λ୧), λ୧  is the mean
%calcium activity in the i-th pixel, and p୧ is the probability of being in bin I (Skaggs et al., 1997;
%Rochefort et al., 2011; Shuman et al., 2020).

tic


          set(0,'DefaultFigureVisible', 'off');

fields_spikes = fieldnames(spike_structure);
fields_pos = fieldnames(pos_structure);
fields_cats = fieldnames(CA_timestamps);

if numel(fields_spikes) ~= numel(fields_pos)
%  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end


for i = 1:numel(fields_spikes)

      fieldName_spikes = fields_spikes{i}
      fieldValue_spikes = spike_structure.(fieldName_spikes);
      peaks_time = fieldValue_spikes;
      if length(peaks_time)>1


      fieldName_cats = fields_cats{i};
      curr_CA_timestamps = CA_timestamps.(fieldName_cats);

      index = strfind(fieldName_spikes, '_');
      spikes_date = fieldName_spikes(index(2)+1:end);

      fieldName_pos = fields_pos{i};
      fieldValue_pos = pos_structure.(fieldName_pos);
      pos = fieldValue_pos;

      index = strfind(fieldName_spikes, '_');
      if length(peaks_time) <5
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        fprintf('fewer than five spikes')
        continue
      end

      mutinfo = NaN(size(peaks_time,1),1);


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

      fprintf('got velocities')
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
              validHighVelIndices = [validHighVelIndices, highVelIndices(ii)];
          end
      end


      goodpos = pos(validHighVelIndices,:);

      all_highspeedspikes = peaks_time(:,validHighVelIndices);

  fprintf('starting units')
      numunits = size(peaks_time,1);
      if numunits<=1
        mutinfo = NaN;
        warning('you have no cells and no spikes')
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue
      else
          for k=1:numunits

          highspeedspikes = all_highspeedspikes(k,:);



            if length(highspeedspikes)>0


              [trace_mean occprob] = CA_normalizePosData_trace(highspeedspikes, goodpos, dim, 1.000);

                if (size(trace_mean,1)) < (size(trace_mean,2))
                  trace_mean = trace_mean';
                end
                if (size(occprob,1)) < (size(occprob,2))
                  occprob = occprob';
                end

              mutinfo(k) = mutualinfo([trace_mean, occprob]);
            else
              mutinfo(k) = NaN;
              fprintf('not enough high speed spikes')
            end
          end
      end


mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
else
  mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
end

end

f = mutualinfo_struct;
