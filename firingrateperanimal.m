function firingrateperanimal(ratNames)
% ratNames: e.g., {'rat0314', 'rat0222'}

colors = [0.1 0.6 0.9]; % match previous plot palette
dayLabels = {'An-2', 'An-1', 'An'};

allMeans = [];

figure('Color','w', 'Position', [300 300 800 400]);



for r = 1:length(ratNames)
    rat = evalin('base', ratNames{r});
    dateList = autoDateList(rat);
    An_day = rat.An;
    idx = find(strcmp(dateList, An_day));
    theseDays = dateList(idx-2:idx);

    for d = 1:3
      field = ['rates_' theseDays{d}];

      spikenum = rat.Ca_peaks.(['CA_peaks_' theseDays{d}]);
      times = rat.CS_times.(['CS_' theseDays{d}]);
      time = times(end)-times(1);
      spikenum = sum(~isnan(spikenum'));

      rate = spikenum./time;
      rat.rates.(field) = rate;
    end
assignin('base', ratNames{r}, rat);        

end
