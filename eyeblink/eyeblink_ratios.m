function f = eyeblink_ratios(CS_time_structure, spike_structure)

  fields_CS = fieldnames(CS_time_structure);
  fields_spikes = fieldnames(spike_structure);

if numel(fields_CS) ~= numel(fields_spikes)
  error('your spike and CS structures do not have the same number of values')
end

for i = 1:numel(fields_CS)
    fieldName_TS = fields_CS{i}
    fieldValue_CS = CS_struct.(fieldName_CS);
    CS = fieldValue_CS;

    index = strfind(fieldName_CS, '_');
    CS_date = fieldName_CS(index(2)+1:end);


    fieldName_spikes = fields_spikes{i};
    fieldValue_spikes = spikes_struct.(fieldName_spikes);
    spikes = fieldValue_spikes;

    index = strfind(fieldName_spikes, '_');
    spikes_date = fieldName_spikes(index(2)+1:end);

    if strcmp(CS_date, spikes_date)==1
      CS = CS;
      US = CS-.750; %this is start of CS, -.5 is end of CS
      pretrial=US-.750;
    else
      CS_date
      spikes_date
      error('your spike name does not match CS name')
    end

    for k = 1:size(spikes,1) %check this
      currentspike = spike(1,:); %check this

      %%%%%%%%%%%%GO FROM HERE
      spikesIntertrial = intersect(currentspike, intertrial);
      spikesCue = intersect(currentspike, cueOnly);
      spikesReward = intersect(currentspike, reward);

      %finding rates
      rateIntertrial = length(spikesIntertrial)/timeintertrial;
      rateCue = length(spikesCue)/timecueOnly;
      rateReward = length(spikesReward)/timereward;

      %finding difference from baseline
      changeCue = rateCue/rateIntertrial;
      changeReward = rateReward/rateIntertrial;
