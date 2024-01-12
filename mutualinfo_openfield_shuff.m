function f = mutualinfo_openfield_shuff(spike_structure, pos_structure, velthreshold, dim, num_times_to_run, ca_MI)
%finds mutual info for a bunch of cells
%little did I know i already had code for this: ca_mutualinfo_openfield.m

tic




fields_spikes = fieldnames(spike_structure);
fields_pos = fieldnames(pos_structure);
fields_MI = fieldnames(pos_structure);

if numel(fields_spikes) ~= numel(fields_pos)
  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end


for i = 1:numel(fields_spikes)
      fieldName_spikes = fields_spikes{i};
      fieldValue_spikes = spike_structure.(fieldName_spikes);
      peaks_time = fieldValue_spikes;

      index = strfind(fieldName_spikes, '_');
      spikes_date = fieldName_spikes(index(2)+1:end)

      fieldName_pos = fields_pos{i};
      fieldValue_pos = pos_structure.(fieldName_pos);
      pos = fieldValue_pos;

      index = strfind(fieldName_spikes, '_');
      pos_date = fieldName_spikes(index(2)+1:end)


      size(pos)
      vel = ca_velocity(pos);
      %vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30.0005); %originally had this at 30, trying with 15 now
      goodvel = find(vel(1,:)>=velthreshold);
      goodtime = pos(goodvel, 1);
      goodpos = pos(goodvel,:);

      mintime = vel(2,1);
      maxtime = vel(2,end);

      numunits = size(peaks_time,1);
      mutinfo = NaN(3,numunits);

      fieldName_MI = fields_MI{i};
      fieldValue_MI = pos_structure.(fieldName_MI);
      MI = fieldValue_MI;

      if numunits<=1
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        warning('you have no spikes')
      else
          for k=1:numunits

                currspikes = peaks_time(k,:);
                if isnan(MI)==1
                  mutinfo(1, k) = NaN;
                  mutinfo(2, k) = NaN;
                  continue
                else
                  highspeedspikes = [];
                end

                for i=1:length(currspikes) %finding if in good vel
                    [minValue,closestIndex] = min(abs(currspikes(i)-goodtime));
                    if minValue <= 1 %if spike is within 1 second of moving. no idea if good time
                      highspeedspikes(end+1) = currspikes(i);
                    end
                end


                shuf = NaN(num_times_to_run,1);
                %for l = 1:num_times_to_run
                parfor l = 1:num_times_to_run
                      %fprintf('survived the great parfor loop trauma of jan 10')
                      if isnan(MI(k))==0 && length(highspeedspikes)>1
                        shufff = randsample(goodtime, length(highspeedspikes));
                        shufff = sort(shufff);

                        [rate totspikes totstime colorbar spikeprob occprob] = CA_normalizePosData(shufff,goodpos,dim, 1.000);

                        shuf(l) = mutualinfo([spikeprob', occprob']);
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
                    rank = index./length(shuf)
                    mutinfo(3, k) = rank;
                  else
                    mutinfo(3,k) = NaN;
                  end



              end

    mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
    end
  end


f = mutualinfo_struct;
end
