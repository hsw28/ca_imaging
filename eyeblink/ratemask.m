function f = ratemask
% gives a 1 for any cell woth EITHER task or nontask rate >=0.05, gives a 0 if neither

ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};

nRats = numel(ratNames);

% Preallocate stats
rate_task_mask = nan(nRats,1);
rate_nontask_mask = nan(nRats,1);
rate_task_mask_non = nan(nRats,1);
rate_nontask_mask_non = nan(nRats,1);
sem_task_mask = nan(nRats,1);
sem_nontask_mask = nan(nRats,1);
sem_task_mask_non = nan(nRats,1);
sem_nontask_mask_non = nan(nRats,1);
pvals_taskmask = nan(nRats,1);
pvals_nonmask = nan(nRats,1);


for r = 1:nRats
    rat = evalin('base', ratNames{r});
    dateList = autoDateList(rat);
    idx = find(strcmp(dateList, rat.An));
    theseDays = dateList(idx-2:idx);

    rt1_all = []; rt2_all = []; rp1_all = []; rp2_all = [];

    for d = 1:3
        dateStr = theseDays{d};
        spikeMat = rat.Ca_peaks.(['CA_peaks_' dateStr]);
        nNeurons = size(spikeMat,1);
        pos = rat.pos.(['pos_' dateStr]);
        cs_times = rat.CS_times.(['CS_' dateStr]);
        ts = pos(:,1);

        % Task mask
        taskMask = false(size(ts));
        for i = 1:numel(cs_times)
            taskMask = taskMask | (ts >= cs_times(i) & ts <= cs_times(i) + 2);
        end
        taskTimes = ts(taskMask);
        nonTaskTimes = ts(~taskMask);

        rate = NaN(nNeurons,1);
        for ni = 1:nNeurons
            spikes = spikeMat(ni,:);
            spikes = spikes(~isnan(spikes));

            % Split spikes
            taskSpikes = spikes(ismembertol(spikes, taskTimes, 1e-3));
            nonTaskSpikes = setdiff(spikes, taskSpikes, 'stable');

            trial_rate = length(taskSpikes)./(length(cs_times)*2);
            nontrial_rate = length(nonTaskSpikes)./((cs_times(end)+2)-(length(cs_times)*2));

            if trial_rate<0.05 && nontrial_rate<0.05
              rate(ni) = 0;
            else
              rate(ni) = 1;
            end
        end

        field = ['ratemask_' theseDays{d}];
        rat.ratemask.(field) = rate;

      end

      assignin('base', ratNames{r}, rat);
end
