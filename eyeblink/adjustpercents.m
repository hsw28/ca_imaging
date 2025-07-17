function f = adjustpercents(ratNames)

tot = [];
for r = 1:length(ratNames)
    rat = evalin('base', ratNames{r});
    dateList = autoDateList(rat);
    An_day = rat.An;
    idx = find(strcmp(dateList, An_day));
    theseDays = dateList(idx-2:idx);

    for d = 1:3
      miPlaceShuff = rat.MI_noCSUS_shuff.(['MI_' theseDays{d}]);
      miCSUS = rat.MI_CSUS15.(['MI_' theseDays{d}]);
      nanny = isnan(miCSUS);
      sum(~isnan(miCSUS))
      miPlaceShuff(nanny,3) = NaN;
      tot(end+1) = sum((miPlaceShuff(:,3)>=.95))./sum(~isnan(miPlaceShuff(:,3)));
    end
end

mean(tot);
std(tot);
