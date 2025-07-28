function MIstatsCSUS

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
    mi_struct_shuff = filterFieldsByDay(rat.MI_CSUS15_shuff, daysToUse);
    mi_struct = filterFieldsByDay(rat.MI_CSUS15, daysToUse);

    mi_struct2 = filterFieldsByDay(rat.MI_noCSUS_shuff, daysToUse)
    ratemask_struct = filterFieldsByDay(rat.ratemask, daysToUse);


    fprintf('\n---- %s ----\n', ratVar);

    for d = 1:numel(daysToUse)
        dateStr = daysToUse{d};

        ratemask = ratemask_struct.(['ratemask_' dateStr]);

        dateList = daysToUse{d};
        M = mi_struct.(['MI_' dateStr]);  % Assuming 4-column format: [x x percentile firingRate]
        Mper = mi_struct_shuff.(['MI_' dateStr]);  % Assuming 4-column format: [x x percentile firingRate]


        percent = Mper(:,3);
        percent(ratemask==0) = NaN;

        allOver95 = sum(percent > 0.95) ./ sum(ratemask);

        tot(end+1) = allOver95;

        fprintf('%s: total=%.2f%%, fr>.05=%.2f%%\n', ...
            dateList, 100*allOver95);

        allval(end+1) = allOver95;
    end

        fprintf('all')
        meantot(end+1) = mean(tot)
        stdtot = std(tot)

end
fprintf('ALL RATS -- all cells')
mean(meantot)
std(meantot)

end
