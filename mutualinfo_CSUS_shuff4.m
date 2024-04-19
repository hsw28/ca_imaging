function [mutualinfo_struct, shuf_all] = mutualinfo_CSUS_shuff4(spike_structure, CSUS_structure, do_you_want_pretrial, num_times_to_run, MI_CSUS)
%finds 'mutual info' for CS/US/ non CS/US
%CSUS_structure should come from BULKconverttoframe.m
%do_you_want_pretrial: 0 for only cs us, 1 for cs us pretrial
%how many divisions you wanted-- for ex,
    % do_you_want_pretrial = 1
    % how_many_divisions = 2 will just split between cs and us
                        %= 10 will split CS and US each into 5
%right now because im lazy how_many_divisions must be a factor of 10


fprintf('running mutualinfo_CSUS_shuff')
divisions = 4;
fields_spikes = fieldnames(spike_structure);
fields_CSUS = fieldnames(CSUS_structure);
fields_MI = fieldnames(MI_CSUS);

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
      CSUS = fieldValue_CSUS;

      fieldName_MI = fields_MI{i};
      fieldValue_MI = MI_CSUS.(fieldName_MI);
      MI = fieldValue_MI;

      index = strfind(fieldName_spikes, '_');
      CSUS_date = fieldName_spikes(index(2)+1:end)

      if length(peaks_time) <5
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

      biggest = max(time)
      [I,J] = find(peaks_time>biggest);
      peaks_time(I,J) = NaN;

      occ_in_CS_US = zeros(1,8);
      occ_intertrial = zeros(1,1);
      occ_pretrial = zeros(1,1);

        numbering = 8./divisions;
        previousz = 0;
          for z=0:numbering:8
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
              shuf = NaN(num_times_to_run,1);
            if length(currspikes)>0  %finding how many spikes in each time bin
                parfor l = 1:num_times_to_run
                  spikes_in_CS_US = zeros(1,8);
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
                                    shuf(l) = mutualinfo([spikeprob, occprob]); %is this oriented the right way
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
                                    shuf(l) = mutualinfo([spikeprob, occprob]); %is this oriented the right way
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

                  else
                          mutinfo(1,k) = NaN;
                          mutinfo(2,k) = NaN;
                          mutinfo(3,k) = NaN;
                  end


            shuf_all = [shuf_all, shuf];

            end %for k=1:size(peaks_time,1)

        end %if numunits<=1


            mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
    end


      f = mutualinfo_struct;
      fprintf('saving\n');

      MI_CSUS_shuff = f;

      % Determine the suffix based on do_you_want_pretrial
      if do_you_want_pretrial == 0
          suffix = '';
      elseif do_you_want_pretrial == 1
          suffix = 'pretrial';
      end
      % Create the dynamic variable name
      variableName = sprintf('MI_CSUS%d_%s_shuff', how_many_divisions, suffix);

      % Assign the structure to the new variable name
      eval([variableName ' = MI_CSUS_shuff;']);

      % Get the current date and time as a string
      currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');

      % Create a filename with the timestamp
      filename = sprintf('results_%s_%s.mat', variableName, currentDateTime);

      % Save the output to the .mat file with the timestamped filename
      save(filename, variableName);
      fprintf('File saved successfully as %s\n', filename);


    end
