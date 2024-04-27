function bars = eyeblink_PSTH_single(spike_structure, US_time_structure, sp)


if max(US_time_structure)>100000
  US_time_structure = US_time_structure./1000000;
  max(US_time_structure)
end






bars = [];


    spikes = sort(spike_structure(:));
    center = US_time_structure;
    center = center';

    center = sort(center);
    spikes = sort(spikes);
    spikes = spikes(~isnan(spikes));
    %psth_bars = psth(center, spikes);
    b = 1.5;
    spikes = spikes';
    psth_bars = psth(center, spikes, 'lags', [-b:.1:b]);
    bars = [bars; psth_bars];

    subplot(4,2,sp)
    bar(psth_bars)
    hold on

    %set(gca,'XTick',0:5:30)
    vline((b*4)./b+8, 'r')
    %vline((b*4)./b+8, 'r', 'US')
    vline(16/2, 'k')
    %vline(16/2, 'k', 'CS on')
    %set(gca,'XTickLabel',0:3/15:3)
