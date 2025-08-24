function bitsperstats

ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};


meantotSPIKES = [];
meantotSECS = [];
stdtot = [];
stdtot5 = [];
allvalspike = [];
allvalsec = [];

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

    % Filter MI struct -- FOR PLACE
    bitsperspike = filterFieldsByDay(rat.per_spike_shuff, daysToUse);
    bitspersec = filterFieldsByDay(rat.per_sec_shuff, daysToUse);

    % Filter MI struct -- FOR CSUS
%    bitsperspike = filterFieldsByDay(rat.per_spike_shuffCSUS, daysToUse);
%    bitspersec = filterFieldsByDay(rat.per_sec_shuffCSUS, daysToUse);


    ratemask_struct = filterFieldsByDay(rat.ratemask, daysToUse);

    fprintf('\n---- %s ----\n', ratVar);

    for d = 1:numel(daysToUse)
        dateStr = daysToUse{d};

        Mspike = bitsperspike.(['bitsPerSpike_' dateStr]);  % Assuming 4-column format: [x x percentile firingRate]
        Msec = bitspersec.(['bitsPerSec_' dateStr]);  % Assuming 4-column format: [x x percentile firingRate]

        ratemask = ratemask_struct.(['ratemask_' dateStr]);

        percent = Mspike(:,3);
        percent(ratemask==0) = NaN;
        allOver95 = sum(percent > 0.95) ./ sum(ratemask);
        allvalspike(end+1) = allOver95*100;

        percent = Msec(:,3);
        percent(ratemask==0) = NaN;
        allOver95 = sum(percent > 0.95) ./ sum(ratemask);
        allvalsec(end+1) = allOver95*100;

    %    frOver02 = sum(percent(fr > 0.02) > 0.95) / sum(fr > 0.02);
    %    frOver03 = sum(percent(fr > 0.03) > 0.95) / sum(fr > 0.03);
    %    frOver05 = sum(percent(fr > 0.05) > 0.95) / sum(fr > 0.05);


      %  fprintf('%s: total=%.2f%%, fr>.02=%.2f%%, fr>.03=%.2f%%, fr>.05=%.2f%%\n', ...
      %      dateStr, 100*allOver95, 100*frOver02, 100*frOver03, 100*frOver05);

      fprintf('%s: spike total=%.2f%%\n', dateStr, allvalspike(end));
      fprintf('%s: sec total=%.2f%%\n', dateStr, allvalsec(end));

    end

    %fprintf('all spike')

    meantotSPIKES(end+1) = mean(allvalspike);
    stdtot = std(allvalspike);

    %fprintf('all sec')

    meantotSECS(end+1) = mean(allvalsec);
    stdtot = std(allvalsec);

  %  fprintf('0.05hz')
  %  meantot5(end+1) = mean(tot5)
  %  stdtot5 = std(tot5)


end

fprintf('ALL RATS -- spikes')
mean(allvalspike)
allvalspike = reshape(allvalspike, 3,5);
std(mean(allvalspike))

fprintf('ALL RATS -- secs')
mean(allvalsec)
allvalsec = reshape(allvalsec, 3,5);
std(mean(allvalsec))


%allval = reshape(allval, 3,5);
%figure
%bar(allval')
end
