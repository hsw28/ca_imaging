function f = mutualinfo_CSUS_trace_shuff(spike_structure, CSUS_structure, do_you_want_CSUS_or_CSUSnone, how_many_divisions, num_times_to_run, MI_CSUS)
%finds 'mutual info' for CS/US/ non CS/US
%CSUS_structure should come from BULKconverttoframe.m
%do_you_want_CSUS_or_CSUSnone: 1 for only cs us, 0 for cs us pretrial
%how many divisions you wanted-- for ex,
    % do_you_want_CSUS_or_CSUSnone = 1
    % how_many_divisions = 2 will just split between cs and us
                        %= 10 will split CS and US each into 5
%right now because im lazy how_many_divisions must be a factor of 10


fprintf('running mutualinfo_CSUS_trace_shuff')

divisions = how_many_divisions;
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

      occ_in_CS_US = zeros(1,10);
      occ_pretrial = zeros(1,1);
      occ_intertrial = zeros(1,1);

        numbering = 10./divisions;
        previousz = 0;

        fprintf('assigning IDs')
        for z=0:numbering:10
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


        fprintf('going through spikes')
      if length(peaks_time)<3 || length(unique(CSUS))<3
          mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
          warning('you have no spikes')
      else
        mutinfo = NaN(3,size(peaks_time,1));
        for k=1:size(peaks_time,1)
          currspikes = peaks_time(k,:);

              if length(currspikes)>0 &&  isnan(MI(k))==0 && length(unique(CSUS)>=3) %finding how many spikes in each time bin
                shuf = NaN(num_times_to_run,1);
                parfor d = 1:num_times_to_run
                          if do_you_want_CSUS_or_CSUSnone == 1
                            wantedindex = find(CSUS>0);
                          else
                            wantedindex = find(CSUS>0 | CSUS == -1);
                          end
                            wanted = currspikes(wantedindex);
                            shuff = currspikes;
                            shuff(wantedindex) = wanted(randperm(length(wanted)));

                          uni = unique(CSUS);
                          occ_in_CS_US = NaN(length(uni)-1,1);
                          spikes_in_CS_US = NaN(length(uni)-1,1);

                          occ_pretrial = zeros(1,1);
                          spikes_pretrial = NaN(1,1);

                          index=1;
                                for q=0:numbering:10
                                      if q == 0
                                        wanted= (find(CSUS == -1));
                                        occ_pretrial = length(wanted);
                                        spikes_pretrial = nanmean(currspikes(wanted));
                                        index=1;

                                      else
                                        wantedtimes = find(CSUS == q);
                                        wantedtimes = wantedtimes(find(wantedtimes<=length(currspikes)));
                                        trace_mean = nanmean(shuff(wantedtimes));
                                        occ_in_CS_US(index) = length(wantedtimes);
                                        spikes_in_CS_US(index) = trace_mean;
                                        index = index+1;
                                      end
                                end


                                if do_you_want_CSUS_or_CSUSnone == 0
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
                  fprintf('throw')
                  mutinfo(1,k) = NaN;
                  mutinfo(2,k) = NaN;
                  mutinfo(3,k) = NaN;
                end


            end
            mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';

          end

      end



      f = MI_CSUS_trace_shuff;

    end
