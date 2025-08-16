
function [per_spike per_sec] = bitsper_CSUS15_shuff(spike_structure, CSUS_structure, do_you_want_pretrial, num_times_to_run, ca_bitsper)
%finds 'mutual info' for CS/US/ non CS/US
%CSUS_structure should come from BULKconverttoframe.m
%do_you_want_pretrial: 0 for only cs us, 1 for cs us pretrial
%how many divisions you wanted-- for ex,
    % do_you_want_pretrial = 1
    % how_many_divisions = 2 will just split between cs and us
                        %= 10 will split CS and US each into 5
%right now because im lazy how_many_divisions must be a factor of 10


fprintf('running mutualinfo_CSUS_shuff')
divisions = 15;
fields_spikes = fieldnames(spike_structure);
fields_CSUS = fieldnames(CSUS_structure);
fields_MI = fieldnames(ca_bitsper)

if numel(fields_spikes) ~= numel(fields_CSUS)
  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

shuf_all = [];
for i = 1:numel(fields_spikes)
      fieldName_spikes = fields_spikes{i};
      fieldValue_spikes = spike_structure.(fieldName_spikes);
      peaks_time = fieldValue_spikes;

      index = strfind(fieldName_spikes, '_');
      spikes_date = fieldName_spikes(index(2)+1:end);

      fieldName_CSUS = fields_CSUS{i};
      fieldValue_CSUS = CSUS_structure.(fieldName_CSUS);

      if startsWith(spikes_date, '2022')
          fieldValue_CSUS(1,1:15) = zeros(1,15);
      end

      numtrials = sum(fieldValue_CSUS(1,:)==1);
      CSUS = fieldValue_CSUS;

      fieldName_MI = fields_MI{i};
      fieldValue_MI = ca_bitsper.(fieldName_MI);
      MI = fieldValue_MI;

      index = strfind(fieldName_spikes, '_');
      CSUS_date = fieldName_spikes(index(2)+1:end)

      if length(peaks_time) <1
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue
      end

      numunits = size(peaks_time,1);

      time = CSUS(2,:);
      CSUS = CSUS(1,:);

      biggest = max([peaks_time(:)]);
      [minValue,closestIndex] = min(abs(biggest-time));
      CSUS = CSUS(1:closestIndex);
      time = time(1:closestIndex);

      biggest = max(time);
      [I,J] = find(peaks_time>biggest);
      peaks_time(I,J) = NaN;

      occ_in_CS_US = zeros(1,15);
      occ_intertrial = zeros(1,1);
      occ_pretrial = zeros(1,1);

        numbering = 15./divisions;
        previousz = 0;
          for z=0:numbering:15
            if z==0
              occ_pretrial = length(find(CSUS ==-1));
              occ_intertrial = length(find(CSUS ==0));
            else
              wanted1 = find(CSUS > previousz);
              wanted2 = find(CSUS <= z);
              wanted = intersect(wanted1,wanted2);
              occ_in_CS_US(z) = length(wanted);
              CSUS(wanted) = z;
              previousz = z;
            end
          end



      mutinfo = NaN(3,numunits);
      if numunits<=1 || length(unique(CSUS))<3
          mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
          warning('you have no cells or no spikes or no CS/US')
      else
        for k=1:size(peaks_time,1) %selecting cell
                currspikes = peaks_time(k,:);
                if isnan(MI(k))==1
                    mutinfo(1, k) = NaN;
                    mutinfo(2, k) = NaN;
                    mutinfo(3, k) = NaN;
                    continue
                end


              set(0,'DefaultFigureVisible', 'off');
              shuf = NaN(num_times_to_run,2);
            if length(currspikes)>0  %finding how many spikes in each time bin
                parfor l = 1:num_times_to_run
                  spikes_in_CS_US = zeros(1,15);
                  spikes_intertrial = zeros(1,1);
                  spikes_pretial = zeros(1,1);

                  if do_you_want_pretrial == 0
                    wanted = find(CSUS > 0);
                    shuffCSUS = CSUS;
                    shufff = wanted(randperm(length(wanted)));
                    shuffCSUS(wanted) = CSUS(shufff);
                  else
                    wanted = find(CSUS > 0 | CSUS == -1);
                    shuffCSUS = CSUS;
                    shufff = wanted(randperm(length(wanted)));
                    shuffCSUS(wanted) = CSUS(shufff);
                  end

                                for q =1:length(currspikes)
                                  if isnan(currspikes(q))==1
                                    continue
                                  end
                                [c index] = (min(abs(currspikes(q)-time))); %
                                spikebin = shuffCSUS(index);
                                        if spikebin == 0
                                          spikes_intertrial = spikes_intertrial+1;
                                        elseif spikebin == -1
                                          spikes_pretial = spikes_pretial+1;
                                        else
                                          spikes_in_CS_US(spikebin) = spikes_in_CS_US(spikebin)+1;
                                        end
                                  end
                              if do_you_want_pretrial == 1
                                    pretrial_occprob = occ_pretrial*(1/7.5);
                                    spikes_occprob = occ_in_CS_US.*(1/7.5);
                                    occprob = [pretrial_occprob, spikes_occprob];
                                    occprob = occprob./nansum(occprob);
                                    spikeprob =  [spikes_pretial, spikes_in_CS_US];
                                    if (size(spikeprob,1)) < (size(spikeprob,2))
                                      spikeprob = spikeprob';
                                    end
                                    if (size(occprob,1)) < (size(occprob,2))
                                      occprob = occprob';
                                    end
                                    [bitsPerSpike, bitsPerSecond] = bits_perCSUS(spikeprob./(numtrials/7.5), occprob);
                                    shuf(l,:) = [bitsPerSpike, bitsPerSecond];
                              else
                                    occprob = occ_in_CS_US.*(1/7.5);
                                    occprob = occprob./nansum(occprob);
                                    spikeprob =  [spikes_in_CS_US];
                                    if (size(spikeprob,1)) < (size(spikeprob,2))
                                      spikeprob = spikeprob';
                                    end
                                    if (size(occprob,1)) < (size(occprob,2))
                                      occprob = occprob';
                                    end
                                    [bitsPerSpike, bitsPerSecond] = bits_perCSUS(spikeprob./(numtrials/7.5), occprob);
                                    shuf(l,:) = [bitsPerSpike, bitsPerSecond];
                              end
                      end



                      shSort1 = sort(shuf(:,1),'ascend');
                      shSort2 = sort(shuf(:,2),'ascend');

                      % 95-th percentile and shuffle mean
                      p95_idx = round(0.95 * num_times_to_run);
                      bitsPspike(1,k) = shSort1(p95_idx);
                      bitsPsec  (1,k) = shSort2(p95_idx);

                      bitsPspike(2,k) = nanmean(shSort1);
                      bitsPsec  (2,k) = nanmean(shSort2);

                      % Rank (percent of shuffle ≥ observed)
                      % Rank (percent of shuffle ≥ observed)
                      obsSpk = MI(1,k);
                      obsSec = MI(2,k);
                      rankSpk = mean(shSort1 >= obsSpk);
                      rankSec = mean(shSort2 >= obsSec);
                      bitsPspike(3,k) = rankSpk;
                      bitsPsec  (3,k) = rankSec;

            end %for k=1:size(peaks_time,1)

        end %if numunits<=1

        fprintf('assigning MI')
        per_spike.(sprintf('bitsPerSpike_%s', spikes_date)) = bitsPspike';
        per_sec.(sprintf('bitsPerSec_%s', spikes_date)) = bitsPsec';
    end

%{
      f = mutualinfo_struct;
      fprintf('saving\n');

      MI_CSUS5_shuff = f;

      % Determine the suffix based on do_you_want_pretrial
      if do_you_want_pretrial == 0
          suffix = '';
      elseif do_you_want_pretrial == 1
          suffix = '_pt';
      end


      % Create the dynamic variable name with the suffix
      variableName = sprintf('MI_CSUS5_shuff%s', suffix);

      % Assign the structure to the new variable name
      eval([variableName ' = MI_CSUS5_shuff;']);

      % Get the current date and time as a string
      currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');

      % Create a filename with the timestamp
      filename = sprintf('results_%s_%s.mat', variableName, currentDateTime);

      % Save the output to the .mat file with the timestamped filename
      save(filename, variableName);
      fprintf('File saved successfully as %s\n', filename);

%}

    end
