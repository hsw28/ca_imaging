function f = adjustpercents(ratNames)

  tot = [];

  for r = 1:numel(ratNames)
      % Pull the rat struct from the base workspace
      rat      = evalin('base', ratNames{r});

      % --- figure out the three “day-2, day-1, day-0” strings -------------
      dateList = autoDateList(rat);
      An_day   = rat.An;                                 % “analysis day”
      idx      = find(strcmp(dateList, An_day));
      theseDays = dateList(idx-2:idx);                   % cell array (1×3)

      % --- loop over the three days ---------------------------------------
      for d = 1:3
          field = ['MI_' theseDays{d}];
          field2 = ['rates_' theseDays{d}];

          % pull the matrices
          miPlace      = rat.MI_noCSUS.(field);
          miPlaceShuff = rat.MI_noCSUS_shuff.(field);
          rates       = rat.rates.(field2);

          % mask bins where CS-US MI is NaN
          mask = rates<0.05;
          miPlaceShuff(mask, 3) = NaN;                   % overwrite col-3
          miPlace(mask) = NaN;

          % *** SAVE THE UPDATED MATRIX BACK INTO THE STRUCT ***
          rat.MI_noCSUS_shuff.(field) = miPlaceShuff;
          rat.MI_noCSUS.(field) = miPlace;

          % compute your summary metric
          tot(end+1) = sum(miPlaceShuff(:,3) >= .95) / sum(~isnan(miPlaceShuff(:,3)));
      end

      % --- push the updated struct back to the caller workspace -----------
      assignin('base', ratNames{r}, rat);                % overwrites original


  end


mean(tot);
std(tot);
