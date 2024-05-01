function ratios = eyeblink_ratios(US_time_structure, spike_structure, cutoff)
% input times_US structure and Ca_peaks structure
%outputs number of spikes during pretrial, CS, US periods and the ratios of CS:pretrial and US:pretrial
%cuttoff is spikes cell must have (>=) in CS, US, or intertrial period in order to count

  ratios = struct();
  fields_US = fieldnames(US_time_structure);
  fields_spikes = fieldnames(spike_structure);

if numel(fields_US) ~= numel(fields_spikes)
  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

for i = 1:numel(fields_US) %going through each day
    output = [];

    fieldName_US = fields_US{i};
    fieldValue_US = US_time_structure.(fieldName_US);
    US = fieldValue_US;

    index = strfind(fieldName_US, '_');
    US_date = fieldName_US(index(1)+1:end);

    fieldName_spikes = fields_spikes{i};
    fieldValue_spikes = spike_structure.(fieldName_spikes);
    spikes = fieldValue_spikes;

    index = strfind(fieldName_spikes, '_');
    spikes_date = fieldName_spikes(index(2)+1:end);

    if length(US)<5
      fieldName_US
      warning('there are no USs for this day')
      ratios.(sprintf('ratios_%s', spikes_date)) = NaN;
      continue
    end


    if strcmp(US_date, spikes_date)==1
      US_end = US+.5;
      US_start = US;
      CS_start = US_start-.750; %this is start of CS, -.5 is end of CS
      pretrial_start=CS_start-.750;
    else
      US_date
      spikes_date
      error('your spike name does not match US name')
    end

  for k = 1:size(spikes,1) %going through each spike
      currentspike = spikes(k,:);
      inPretrial = [];
      inCS = [];
      inUS = [];


      for z = 1:length(CS_start) %go through different intervals
        inPretrial(end+1) = length(find(currentspike>=pretrial_start(z) & currentspike<CS_start(z)));
        inCS(end+1) = length(find(currentspike>=CS_start(z) & currentspike<US_start(z)));
        inUS(end+1) = length(find(currentspike>=US_start(z) & currentspike<US_end(z)));
      end

      pretrial_sum = sum(inPretrial);


      CS_sum = sum(inCS);
      US_sum = sum(inUS);

      if CS_sum >= cutoff | US_sum >= cutoff | CS_sum >= cutoff
        if pretrial_sum+CS_sum>0
          CS_change = CS_sum-pretrial_sum;
        else
          CS_change = NaN;
        end

        if pretrial_sum+US_sum>0
          US_change = US_sum-pretrial_sum;
        else
          US_change = NaN;
        end
      else
        CS_change = NaN;
        US_change = NaN;
        pretrial_sum = NaN;
        CS_sum = NaN;
        US_sum = NaN;
      end

      newdata = [pretrial_sum, CS_sum, US_sum, CS_change, US_change];
      output = vertcat(output, newdata);
    end

    ratios.(sprintf('ratios_%s', spikes_date)) = output;
  end
