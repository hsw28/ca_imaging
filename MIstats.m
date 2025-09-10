function MIstats

ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};


meantot = [];
meantot5 = [];
stdtot = [];
stdtot5 = [];
allval = [];
for i = 1:numel(ratNames)
  tot = [];
  tot5 = [];
    ratVar = ratNames{i};
    rat = evalin('base', ratVar);

    % Get target days
    dateList = autoDateList(rat);
    idx = find(strcmp(dateList, rat.An));
    if idx < 3
        warning('%s does not have enough days before An. Skipping...', ratVar);
        continue;
    end
    daysToUse = dateList(idx-2:idx);

    % Filter MI struct
    mi_struct = filterFieldsByDay(rat.MI_noCSUS15_shuff, daysToUse);
    %mi_struct = filterFieldsByDay(rat.MI_wCSUS_shuff, daysToUse);
  %  mi_struct = filterFieldsByDay(rat.MI_noCSUS15_8cm_shuff, daysToUse); %%%

    ratemask_struct = filterFieldsByDay(rat.ratemask, daysToUse);

    fprintf('\n---- %s ----\n', ratVar);

    for d = 1:numel(daysToUse)
        dateStr = daysToUse{d};
        M = mi_struct.(['MI_' dateStr]);  % Assuming 4-column format: [x x percentile firingRate]

        ratemask = ratemask_struct.(['ratemask_' dateStr]);

        percent = M(:,3);
        percent(ratemask==0) = NaN;

        allOver95 = sum(percent > 0.95) ./ sum(ratemask);

    %    frOver02 = sum(percent(fr > 0.02) > 0.95) / sum(fr > 0.02);
    %    frOver03 = sum(percent(fr > 0.03) > 0.95) / sum(fr > 0.03);
    %    frOver05 = sum(percent(fr > 0.05) > 0.95) / sum(fr > 0.05);
        tot(end+1) = allOver95;


      %  fprintf('%s: total=%.2f%%, fr>.02=%.2f%%, fr>.03=%.2f%%, fr>.05=%.2f%%\n', ...
      %      dateStr, 100*allOver95, 100*frOver02, 100*frOver03, 100*frOver05);

      fprintf('%s: total=%.2f%%\n', dateStr, 100*allOver95);

      allval(end+1) = allOver95;
    end

    fprintf('all')

    meantot(end+1) = mean(tot)
    stdtot = std(tot)

  %  fprintf('0.05hz')
  %  meantot5(end+1) = mean(tot5)
  %  stdtot5 = std(tot5)


end

fprintf('ALL RATS -- all cells')
mean(allval)
allval = reshape(allval, 3,5);

std(mean(allval))

%fprintf('ALL RATS -- 0.05 cells')
%mean(meantot5)
%std(meantot5)

%allval = reshape(allval, 3,5);
%figure
%bar(allval')
end
