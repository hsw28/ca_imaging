function f = mutualinfo_openfield_shuff(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, num_times_to_run, ca_MI)
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
size(ca_MI)
class(ca_MI)
fields_MI = fieldnames(ca_MI);
catch
  fprintf('issue loading, what it is nobody knows')
  error('problem loading')
end
fprintf('loading TS')
fields_cats = fieldnames(CA_timestamps);
fprintf('all loaded')

if numel(fields_spikes) ~= numel(fields_pos)
  fprintf('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

fprintf('starting loop')
for i = 1:numel(fields_spikes)

      fieldName_MI = fields_MI{i};
      fieldValue_MI = ca_MI.(fieldName_MI);
      MI = fieldValue_MI;

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
      pos = pos(1:closestIndex, :);

      index = strfind(fieldName_spikes, '_');
      pos_date = fieldName_spikes(index(2)+1:end)

      if length(peaks_time) <5
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue
      end


            if (pos(1,1)-pos(end,1))./length(pos) < 1
              pos = convertpostoframe(pos, curr_CA_timestamps);
            end

      pos = smoothpos(pos);

      fprintf('trimming velocity')
      vel = ca_velocity(pos);
      goodvel = find(vel(1,:)>=velthreshold);
      goodtime = pos(goodvel, 1);
      goodpos = pos(goodvel,:);

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
                else
                  highspeedspikes = [];
                end

                for ii=1:length(currspikes) %finding if in good vel
                    [minValue,closestIndex] = min(abs(currspikes(ii)-goodtime));
                    if minValue <= 1/15 %if spike is within 1 second of moving. no idea if good time
                      highspeedspikes(end+1) = currspikes(ii);
                    end
                end


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


                topMI5 = floor(num_times_to_run*.95);
                topMI1 = floor(num_times_to_run*.99);
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

                  [c index] = (min(abs(MI(k)-shuf)));
                  if isnan(index)==0
                    rank = index./length(shuf);
                    mutinfo(3, k) = rank;
                  else
                    mutinfo(3,k) = NaN;
                  end



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
%}


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


  f = mutualinfo_struct
