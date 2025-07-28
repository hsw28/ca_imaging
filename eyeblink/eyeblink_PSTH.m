function bars = eyeblink_PSTH(spike_structure, CS_time_structure)


fields_CS = fieldnames(CS_time_structure);
fields_spikes = fieldnames(spike_structure);

if numel(fields_CS) ~= numel(fields_spikes)
  error('your spike and CS structures do not have the same number of values. you may need to pad your CS structure for exploration days')
end

figure
width = 5;
height = ceil(length(fields_CS)./5);
%height = 3;

bars = [];

for i = 1:length(fields_CS)
%  for i = 1:15

  fieldName_CS = fields_CS{i};
  fieldValue_CS = CS_time_structure.(fieldName_CS);
  CS = fieldValue_CS;


  if length(CS)<5
    fieldName_CS
    warning('there are no CSs for this day')
    continue
  end

  index = strfind(fieldName_CS, '_');
  CS_date = fieldName_CS(index(1)+1:end);

  fieldName_spikes = fields_spikes{i};
  fieldValue_spikes = spike_structure.(fieldName_spikes);
  spikes = fieldValue_spikes;

  index = strfind(fieldName_spikes, '_');
  spikes_date = fieldName_spikes(index(2)+1:end);

  if strcmp(CS_date, spikes_date)==1

    spikes = sort(spikes(:));
    center = CS;
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
    vline(15, 'r', 'CS')
    vline(7.5, 'k')
    vline(7.5, 'k', 'CS on')
    set(gca,'XTickLabel',-1.5:.5:1.5)




  else
    CS_date
    spikes_date
    error('your spike name does not match CS name')
  end

end
