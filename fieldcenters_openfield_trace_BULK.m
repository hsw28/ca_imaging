function field_centers  = fieldcenters_openfield_trace_BULK(trace_time_struct, pos_struct, dim, velthreshold, CA_timestamps)
  %good cells IS AN OPTIONAL INPUT and are indices of the cells you know have fields
  % field ceenters are the highest spiking point, not the geometric center
  %rates returns max rate, av rate, min rate



  fields_peaks_time = fieldnames(trace_time_struct);
  fields_pos = fieldnames(pos_struct);
  fields_cats = fieldnames(CA_timestamps);

  if numel(fields_peaks_time) ~= numel(fields_pos)
    error('your spike and pos structures do not have the same number of values. you may need to pad your US structure for exploration days')
  end


  for i = 1:numel(fields_pos)

        fieldName_spikes = fields_peaks_time{i};
        fieldValue_spikes = trace_time_struct.(fieldName_spikes);
        peaks_time = fieldValue_spikes;


        index = strfind(fieldName_spikes, '_');
        spikes_date = fieldName_spikes(index(2)+1:end)

        fieldName_pos = fields_pos{i};
        fieldValue_pos = pos_struct.(fieldName_pos);
        pos = fieldValue_pos;

        fieldName_cats = fields_cats{i};
        curr_CA_timestamps = CA_timestamps.(fieldName_cats);

        if length(peaks_time) < 5
          field_centers.(sprintf('MI_%s', spikes_date)) = NaN;
          continue
        end

        if (pos(1,1)-pos(end,1))./length(pos) < 1
          pos = convertpostoframe(pos, curr_CA_timestamps);
        end
        if length(peaks_time)>length(pos)
          peaks_time = peaks_time(1:length(pos));
        elseif length(peaks_time)<length(pos)
          pos = pos(1:length(peaks_time),:);
        end

        vel = ca_velocity(pos);
        times = vel(2,:);
        highVelIndices = find(vel(1,:)>=velthreshold);

        % Find indices where velocity is less than or equal to the threshold
        lowVelIndices = find(vel(1,:) < velthreshold);
        % Filter out high velocity indices that are too close to low velocities
        validHighVelIndices = [];

        timeThreshold = 1/7.5;
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



        numunits = size(peaks_time,1);
        for k=1:numunits

          currspikes = all_highspeedspikes(k,:);
          set(0,'DefaultFigureVisible', 'off');
          fr = ca_firingrate(currspikes, pos);

          if fr > .0000000001 && length(currspikes)>0

            [trace_mean occprob] = CA_normalizePosData_trace(currspikes,goodpos,dim, 1.000);
            rate=trace_mean;
            [maxval, maxindex] = max(rate(:));
            [x,y] = ind2sub(size(rate), maxindex);
            maxrate(1, k) = x*dim;
            maxrate(2, k) = y*dim;
          else
            maxrate(1, k) = NaN;
            maxrate(2, k) = NaN;
          end

        end

        field_centers.(sprintf('centers_%s', spikes_date)) = maxrate';
  end
