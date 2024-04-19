function f = mutualinfo_CSUS_trace_shuff(spike_structure, CSUS_structure, do_you_want_pretrial, num_times_to_run, MI_CSUS)
%SHUFFLES CSUS
%finds 'mutual info' for CS/US/ non CS/US
%CSUS_structure should come from BULKconverttoframe.m
%do_you_want_pretrial: 0 for only cs us, 1 for cs us pretrial
%how many divisions you wanted-- for ex,
    % do_you_want_pretrial = 1
    % how_many_divisions = 2 will just split between cs and us
                        %= 10 will split CS and US each into 5
%right now because im lazy how_many_divisions must be a factor of 10


fprintf('running mutualinfo_CSUS_trace_shuff')

divisions = 2;
fields_spikes = fieldnames(spike_structure);
fields_CSUS = fieldnames(CSUS_structure);
fields_MI = fieldnames(MI_CSUS);
          set(0,'DefaultFigureVisible', 'off');

fprintf('error message')
if numel(fields_spikes) ~= numel(fields_CSUS)
  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

fprintf('going through days')
for i = 1:numel(fields_spikes)
      fieldName_spikes = fields_spikes{i};
      fieldValue_spikes = spike_structure.(fieldName_spikes);
      peaks_time = fieldValue_spikes;

      index = strfind(fieldName_spikes, '_');
      spikes_date = fieldName_spikes(index(2)+1:end);

      fieldName_CSUS = fields_CSUS{i};
      fieldValue_CSUS = CSUS_structure.(fieldName_CSUS);
      CSUS = fieldValue_CSUS;

      index = strfind(fieldName_spikes, '_');
      CSUS_date = fieldName_spikes(index(2)+1:end)

      fieldName_MI = fields_MI{i};
      fieldValue_MI = MI_CSUS.(fieldName_MI);
      MI = fieldValue_MI;

      if length(peaks_time) <5
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue
      end


      time = CSUS(2,:);
      CSUS = CSUS(1,:);

      if size(peaks_time,2)>length(CSUS)
        peaks_time = peaks_time(:,1:length(CSUS));
      else
        CSUS = CSUS(1:size(peaks_time,2));
      end

      occ_in_CS_US = zeros(1,2);
      occ_pretrial = zeros(1,1);
      occ_intertrial = zeros(1,1);

        previousz = 0;

        occ_pretrial = length(find(CSUS ==-1));
        occ_intertrial = length(find(CSUS ==0));

        CS_time1 = find(CSUS <= 6);
        CS_time2 = find(CSUS > 0);
        CS_time = intersect(CS_time1,CS_time2);
        US_time = find(CSUS > 6);
        CSUS(CS_time) = 1;
        CSUS(US_time) = 2;

        occ_in_CS_US(1) = length(CS_time);
        occ_in_CS_US(2) = length(US_time);



        fprintf('going through spikes')
      if length(peaks_time)<3 || length(unique(CSUS))<3
          mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
          warning('you have no spikes')
      else
        mutinfo = NaN(3,size(peaks_time,1));
        for k=1:size(peaks_time,1)
          currspikes = peaks_time(k,:);

              if length(currspikes)>0 &&  isnan(MI(k))==0 && length(unique(CSUS))>=3 %finding how many spikes in each time bin


                parfor d = 1:num_times_to_run
                  if do_you_want_pretrial == 0
                    wantedindex = find(CSUS>0);
                    shuffCSUS = CSUS;
                    shufff = wantedindex(randperm(length(wantedindex)));
                    shuffCSUS(wantedindex) = CSUS(shufff);
                  else
                    wantedindex = find(CSUS>0 | CSUS == -1);
                    shuffCSUS = CSUS;
                    shufff = wantedindex(randperm(length(wantedindex)));
                    shuffCSUS(wantedindex) = CSUS(shufff);
                  end

                  wanted = currspikes(wantedindex);

                  uni = unique(CSUS);
                  occ_in_CS_US = NaN(length(uni)-1,1);
                  spikes_in_CS_US = NaN(length(uni)-1,1);

                  occ_pretrial = zeros(1,1);
                  spikes_pretrial = NaN(1,1);

                  index=1;
                  for q=0:2
                        if q == 0
                          wanted= (find(shuffCSUS == -1));
                          occ_pretrial = length(wanted);
                          spikes_pretrial = nanmean(currspikes(wanted));
                          index=1;
                        else
                          wantedtimes = find(shuffCSUS == q);
                          wantedtimes = wantedtimes(find(wantedtimes<=length(currspikes)));
                          trace_mean = nanmean(currspikes(wantedtimes));
                          occ_in_CS_US(index) = length(wantedtimes);
                          spikes_in_CS_US(index) = trace_mean;
                          index = index+1;
                        end
                  end
                  if do_you_want_pretrial == 1
                        pretrial_occprob = occ_pretrial*(1/7.5);
                        spikes_occprob = occ_in_CS_US.*(1/7.5);
                        occprob = [pretrial_occprob, spikes_occprob'];
                        occprob = occprob./nansum(occprob);
                        spikeprob =  [spikes_pretrial, spikes_in_CS_US'];
                        if (size(spikeprob,1)) < (size(spikeprob,2))
                          spikeprob = spikeprob';
                        end
                        if (size(occprob,1)) < (size(occprob,2))
                          occprob = occprob';
                        end
                        shuf(d) = mutualinfo([spikeprob, occprob]); %is this oriented the right way
                  else
                    occprob = occ_in_CS_US.*(1/7.5);
                    occprob = occprob./nansum(occprob);
                    spikeprob =  spikes_in_CS_US;
                    if (size(spikeprob,1)) < (size(spikeprob,2))
                      spikeprob = spikeprob';
                    end
                    if (size(occprob,1)) < (size(occprob,2))
                      occprob = occprob';
                    end
                    compMI = mutualinfo([spikeprob, occprob]);
                    shuf(d) = compMI; %is this oriented the right way
                  end


                end
                topMI5 = floor(num_times_to_run*.95);;
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

                else
                fprintf('throw')
                mutinfo(1,k) = NaN;
                mutinfo(2,k) = NaN;
                mutinfo(3,k) = NaN;
                end


                end
                mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';

                end

                end



                f = mutualinfo_struct;


                MI_CSUS2_trace_shuff = f;

                % Determine the suffix based on do_you_want_pretrial
                if do_you_want_pretrial == 0
                suffix = '';
                elseif do_you_want_pretrial == 1
                suffix = 'pretrial';
                end


                % Create the dynamic variable name
                variableName = sprintf('MI_CSUS2%d_%s_trace_shuff', suffix);

                % Assign the structure to the new variable name
                eval([variableName ' = MI_CSUS2_trace_shuff;']);

                % Get the current date and time as a string
                currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');

                % Create a filename with the timestamp
                filename = sprintf('results_%s_%s.mat', variableName, currentDateTime);

                % Save the output to the .mat file with the timestamped filename
                save(filename, variableName);
                fprintf('File saved successfully as %s\n', filename);


                end
