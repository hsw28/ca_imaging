function f = choose_eyeblink_PSTH(spike_structure, US_time_structure, alignment_matrix, aligned_ratio_matrix, number_to_be_greater_than, days_for_condition_in_vector_form);

%makes a psth for all days only using a subset of cells
%ex choose_eyeblink_PSTH(Ca_peaks, times_US, alignment, align_US_ratio, 5, [2,3,4]);
%this plots psth's for cells that have an entry greater than 5  in align_US_ratio for days 2, 3, and 4

fields_US = fieldnames(US_time_structure);
fields_spikes = fieldnames(spike_structure);
ratios = aligned_ratio_matrix;

if numel(fields_US) ~= numel(fields_spikes)
  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

wanted_index = [];
for i = min(days_for_condition_in_vector_form):max(days_for_condition_in_vector_form)

      fieldName_spikes = fields_spikes{i};
      fieldValue_spikes = spike_structure.(fieldName_spikes);
      spikes = fieldValue_spikes;

      index = strfind(fieldName_spikes, '_');
      spikes_date = fieldName_spikes(index(2)+1:end);


      %if length(US)<5
      %  fieldName_US
      %  warning('there are no USs for this day')
      %  continue
      %if strcmp(US_date, spikes_date)==1 | strcmp(ratios_date, spikes_date)==1
        wanted_cells = find(ratios(:,i)<number_to_be_greater_than);
        wanted_index = [wanted_index; wanted_cells];
      %else
      %  US_date
      %  spikes_date
      %  error('your spike name does not match US name')
      %end

end

for i=1:numel(fields_US)
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
  day_index = alignment_matrix(wanted_index, i);
  day_index = day_index(find(day_index>0));
  spike_structure.(fieldName_spikes) = spikes(day_index,:);
end

f = eyeblink_PSTH(spike_structure, US_time_structure);
