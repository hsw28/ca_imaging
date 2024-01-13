function f = mutualinfo_CSUS_trace(spike_structure, CSUS_structure, do_you_want_CSUS_or_CSUSnone, how_many_divisions)
%finds 'mutual info' for CS/US/ non CS/US
%CSUS_structure should come from BULKconverttoframe.m
%do_you_want_CSUS_or_CSUSnone: 1 for only cs us, 0 for cs us none
%how many divisions you wanted-- for ex,
    % do_you_want_CSUS_or_CSUSnone = 1
    % how_many_divisions = 2 will just split between cs and us
                        %= 10 will split CS and US each into 5
%right now because im lazy how_many_divisions must be a factor of 10


divisions = how_many_divisions;
fields_spikes = fieldnames(spike_structure);
fields_CSUS = fieldnames(CSUS_structure);
          set(0,'DefaultFigureVisible', 'off');

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

      mutinfo = NaN(size(peaks_time,1),1);

      numunits = size(peaks_time,1);

      time = CSUS(2,:);
      CSUS = CSUS(1,:);

      occ_in_CS_US = zeros(1,10);
      occ_intertrial = zeros(1,1);
      if do_you_want_CSUS_or_CSUSnone==1
        numbering = 10./divisions;
        previousz = 0;
        for z=0:numbering:10 %%%%%%fix

            if z==0
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

      elseif do_you_want_CSUS_or_CSUSnone==0
          error(' have not written this code yet') %%%%%%%%%%%%%%%%%%%%%%%%
          %should be easy to impliment, already tagging all those times with 0s and counting them
          %just need to decide if i want all the times or just a little before cs/us
          %probably want the latter so need to tag them, prob want to do that outside this code
      end

      if numunits<=1
          mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
          warning('you have no spikes')
      else
          wanted = find(CSUS > 0);
          currspikes = peaks_time(k,:);

          if length(currspikes)>0  %finding how many spikes in each time bin
            occ_in_CS_US = Nan(length(unique(CSUS))-1,1);
            spikes_in_CS_US = NaN(length(unique(CSUS))-1,1);
            for q=0:numbering:10
                  if q == 0 %%%%%%%%%%ADD THIS LATER
                    continue
                  else
                    wantedtimes = find(CSUS == q);
                    trace_mean = nanmean(currspikes(wantedtimes));
                    occ_in_CS_US(q) = length(wantedtimes);
                    spikes_in_CS_US(q) = trace_mean;
                  end
              occprob = occ_in_CS_US.*(1/7.5);
              spikeprob =  spikes_in_CS_US;
              mutinfo(k) = mutualinfo([spikeprob', occprob']); %is this oriented the right way
              end
            else
              mutinfo(k) = NaN;
            end
          end
      end

      mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
      end


      f = mutualinfo_struct;
