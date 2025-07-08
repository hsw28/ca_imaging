function MIstats

ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};

for i = 1:numel(ratNames)
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
    mi_struct = filterFieldsByDay(rat.MI_noCSUS_shuff, daysToUse);

    fprintf('\n---- %s ----\n', ratVar);

    for d = 1:numel(daysToUse)
        dateStr = daysToUse{d};
        M = mi_struct.(['MI_' dateStr]);  % Assuming 4-column format: [x x percentile firingRate]

        percent = M(:,3);
        fr = M(:,4);

        allOver95 = sum(percent > 0.95) / length(percent);
        frOver02 = sum(percent(fr > 0.02) > 0.95) / sum(fr > 0.02);
        frOver03 = sum(percent(fr > 0.03) > 0.95) / sum(fr > 0.03);
        frOver05 = sum(percent(fr > 0.05) > 0.95) / sum(fr > 0.05);


        fprintf('%s: total=%.2f%%, fr>.02=%.2f%%, fr>.03=%.2f%%, fr>.05=%.2f%%\n', ...
            dateStr, 100*allOver95, 100*frOver02, 100*frOver03, 100*frOver05);
    end
end
end
