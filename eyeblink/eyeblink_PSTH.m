function bars = eyeblink_PSTH(spike_structure, US_time_structure)


fields_US = fieldnames(US_time_structure);
fields_spikes = fieldnames(spike_structure);

if numel(fields_US) ~= numel(fields_spikes)
  error('your spike and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end

figure
width = 5;
height = ceil(length(fields_US)./5);
%height = 3;

bars = [];

for i = 1:length(fields_US)
%  for i = 1:15

  fieldName_US = fields_US{i};
  fieldValue_US = US_time_structure.(fieldName_US);
  US = fieldValue_US;


  if length(US)<5
    fieldName_US
    warning('there are no USs for this day')
    continue
  end

  index = strfind(fieldName_US, '_');
  US_date = fieldName_US(index(1)+1:end);

  fieldName_spikes = fields_spikes{i};
  fieldValue_spikes = spike_structure.(fieldName_spikes);
  spikes = fieldValue_spikes;

  index = strfind(fieldName_spikes, '_');
  spikes_date = fieldName_spikes(index(2)+1:end);

  if strcmp(US_date, spikes_date)==1

    spikes = sort(spikes(:));
    center = US;
    center = center';

    center = sort(center);
    spikes = sort(spikes);
    spikes = spikes(~isnan(spikes));
    %psth_bars = psth(center, spikes);
    psth_bars = psth(center, spikes, 'lags', [-1.5:.1:1.5]);
    bars = [bars; psth_bars];

    subplot(height, width, i)
    bar(psth_bars)
    hold on

    set(gca,'XTick',0:5:30)
    vline(15, 'r')
    vline(15, 'r', 'US')
    vline(7.5, 'k')
    vline(7.5, 'k', 'CS on')
    set(gca,'XTickLabel',-1.5:.5:1.5)




  else
    US_date
    spikes_date
    error('your spike name does not match US name')
  end

end
