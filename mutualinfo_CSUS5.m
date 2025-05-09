function f = mutualinfo_CSUS5(spike_structure, CSUS_structure, do_you_want_pretrial)
%finds 'mutual info' for CS/US/ non CS/US
%CSUS_structure should come from BULKconverttoframe.m
%do_you_want_pretrial: 0 for only cs us, 1 for cs us pretrial
%how many divisions you wanted-- for ex,
    % do_you_want_pretrial = 1
    % how_many_divisions = 2 will just split between cs and us
                        %= 10 will split CS and US each into 5
%right now because im lazy how_many_divisions must be a factor of 10


divisions = 5;
fields_spikes = fieldnames(spike_structure);
fields_CSUS = fieldnames(CSUS_structure);

if numel(fields_spikes) ~= numel(fields_CSUS)
  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end


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

      if length(peaks_time) <5
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue
      end

      mutinfo = NaN(size(peaks_time,1),1);

      numunits = size(peaks_time,1);

      time = CSUS(2,:);
      CSUS = CSUS(1,:);

      biggest = max([peaks_time(:)]);
      [minValue,closestIndex] = min(abs(biggest-time));
      CSUS = CSUS(1:closestIndex);
      time = time(1:closestIndex);


      [I,J] = find(peaks_time>biggest);
      peaks_time(I,J) = NaN;



      occ_in_CS_US = zeros(1,10);
      occ_intertrial = zeros(1,1);
      occ_pretrial = zeros(1,1);

        numbering = 10./divisions;
        previousz = 0;
        for z=0:numbering:10 %%%%%%fix
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

      if numunits<=1
          mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
          warning('you have no spikes')
      else

          wanted = find(CSUS > 0);
          for k=1:size(peaks_time,1)
            k
          currspikes = peaks_time(k,:);


          set(0,'DefaultFigureVisible', 'off');
          spikes_in_CS_US = zeros(1,10);
          spikes_intertrial = zeros(1,1);
          spikes_pretial = zeros(1,1);
            if length(currspikes)>0  && length(unique(CSUS))>=3 %finding how many spikes in each time bin
                for q =1:length(currspikes)
                  if isnan(currspikes(q))==1
                    continue
                  end
                  [c index] = (min(abs(currspikes(q)-time))); %
                  spikebin = CSUS(index);
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
                mutinfo(k) = mutualinfo([spikeprob, occprob]) %is this oriented the right way
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
                mutinfo(k) = mutualinfo([spikeprob, occprob]); %is this oriented the right way
              end

            else
              mutinfo(k) = NaN;
            end
          end
      end

      mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
      end


      f = mutualinfo_struct;
