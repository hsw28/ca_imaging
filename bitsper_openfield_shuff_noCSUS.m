function f = bitsper_openfield_shuff_noCSUS(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, CSUS_id_struct, ca_bitsper, num_times_to_run)
%finds mutual info for a bunch of cells
%little did I know i already had code for this: ca_mutualinfo_openfield.m
%returns 95% cutoff, average MI, and rank of actual MI

%IN PREP
%need to input ca_bitsper instead of MI


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

      % Find all time points where CSUS_id > 0
      taskIdx = find(CSUS_id(1,:) > 0);

      % Initialize logical mask the same size as CSUS_id
      keepIdx = false(size(CSUS_id));

      % For each task time point, mark the current and next 15 points
      for i = 1:length(taskIdx)
          idx = taskIdx(i);
          maxIdx = min(idx + 5, length(CSUS_id));  % avoid going out of bounds
          keepIdx(idx:maxIdx) = true;
      end

      % Now find the full range of indices to keep
      goodCSUS = find(keepIdx);
      good_CSUStime = pos(goodCSUS,1);
      good_CSUSpos = pos(goodCSUS,:);


      vel = ca_velocity(pos);
      goodvel = find(vel(1,:)>=velthreshold);
      goodtime = pos(goodvel, 1);
      goodpos = pos(goodvel,:);
      goodvel = setdiff(goodvel, goodCSUS);

      mintime = vel(2,1);
      maxtime = vel(2,end);

      numunits = size(peaks_time,1);
      bitsPspike = NaN(4,numunits);
      bitsPsec = NaN(4,numunits);


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
                  bitsPspike(1, k) = NaN;
                  bitsPspike(2, k) = NaN;
                  bitsPsec(1, k) = NaN;
                  bitsPsec(2, k) = NaN;
                  continue
                else
                  highspeedspikes = [];
                end

                for ii=1:length(currspikes) %finding if in good vel
                  [minValue_CSUS,closestIndex] = min(abs(currspikes(ii)-good_CSUStime));
                  [minValue_vel,closestIndex] = min(abs(currspikes(ii)-goodtime));
                  if minValue_CSUS <= 1/15 & isnan(currspikes(ii))==0 %being CSUS takes precedence
                    continue;
                  elseif minValue_vel <= 1/15 & isnan(currspikes(ii))==0
                    highspeedspikes(end+1) = currspikes(ii);
                  end
                end

                hertz = length(highspeedspikes)./(length(goodpos)/15);

                shuf = NaN(num_times_to_run,2);
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

                        [bitsPerSpike, bitsPerSecond] = bits_per([spikeprob, occprob]);

                        shuf(l,1) = bitsPerSpike;
                        shuf(l,2) = bitsPerSecond;

                      else
                        shuf(l,1) = NaN;
                        shuf(l,2) = NaN;
                      end
                  end



                topMI5 = floor(num_times_to_run*.95);
                topMI1 = floor(num_times_to_run*.99);
                shuf = sort(shuf);
                    if isnan(topMI5)==0
                      bitsPspike(1, k) = shuf(topMI5,1);
                      bitsPsec(1, k) = shuf(topMI5,2);
                    else
                      bitsPspike(1, k) = NaN;
                      bitsPsec(1, k) = NaN;

                    end
                    if isnan(topMI1)==0
                      bitsPspike(2, k) = nanmean(shuf(:,1));
                      bitsPsec(2, k) = nanmean(shuf(:,2));

                    else
                      bitsPspike(2, k) = NaN;
                      bitsPsec(2, k) = NaN;

                    end

                  [c index1] = (min(abs(bitsPspike_OG(k)-shuf(:,1))));
                  [c index2] = (min(abs(bitsPsec_OG(k)-shuf(:,2))));
                  if isnan(index)==0
                    rank1 = index1./length(shuf);
                    bitsPspike(3, k) = rank1;
                    rank1 = index2./length(shuf);
                    bitsPsec(3, k) = rank;

                  else
                    bitsPspike(3,k) = NaN;
                    bitsPsec(3,k) = NaN;

                  end

                  bitsPspike(4,k) = hertz;
                  bitsPsec(4,k) = hertz;



              end
    fprintf('assigning MI')
    mutualinfo_struct.(sprintf('bitsPerSpike_%s', spikes_date)) = bitsPspike';
        mutualinfo_struct.(sprintf('bitsPerSec_%s', spikes_date)) = bitsPsec';
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


  f = mutualinfo_struct;
