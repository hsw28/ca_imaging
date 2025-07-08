function perDayStruct = alignEMGandVelocity_perDay(ratName)

    ratName = char(ratName);
    animal = evalin('base', ratName);
    emg = evalin('base', [ratName 'emg']);

    dateList = autoDateList(animal);
    finalDay = animal.An;
    finalIdx = find(strcmp(dateList, finalDay));

    if finalIdx < 3
        warning('%s has fewer than 3 learned sessions', ratName);
        perDayStruct = [];
        return;
    end

    learnedDays = dateList(finalIdx-2:finalIdx);
    win = [0 2];
    Fs_interp = 7.5;
    t_interp = win(1):1/Fs_interp:win(2);
    t_interp = round(t_interp, 6);
    trace_mask = t_interp >= 0 & t_interp < 2;

    dcount = 0;
    for d = 1:length(learnedDays)
        day = learnedDays{d};
        dateStr = day;

        try
            EMG   = emg.EMG.(['EMG_' dateStr]);
            EMGts = emg.EMG_ts.(['EMGts_' dateStr]);
            pos   = animal.pos.(['pos_' dateStr]);
            vel   = ca_velocity(pos);  % [2 × n]
            cs_times = animal.CS_times.(['CS_' dateStr]);
        catch
            warning('Missing data for %s on %s', ratName, dateStr);
            continue;
        end

        vel_trials = [];
        used_windows = [];

        for i = 1:length(cs_times)
          if strcmp(ratName, 'rat0816') && i == 1
              continue;
          end

            cs = cs_times(i);
            t_win = cs + win;

            vel_idx = vel(2,:) >= t_win(1) & vel(2,:) <= t_win(2);
            if sum(vel_idx) < 10, continue; end
            vel_interp = interp1(vel(2, vel_idx), vel(1, vel_idx), cs + t_interp, 'linear', 'extrap');
            vel_trials(end+1, :) = vel_interp;
            used_windows = [used_windows; t_win];
        end

        if isempty(vel_trials), continue; end

        % Get all trace values
        trace_vals = vel_trials(:, trace_mask);
        trace_vals = trace_vals(:);
        trace_vals = trace_vals(~isnan(trace_vals));

        % Get all non-trial values from full velocity trace
        all_times = vel(2,:);
        all_vals  = vel(1,:);
        trial_mask = false(size(all_times));
        for w = 1:size(used_windows,1)
            trial_mask = trial_mask | (all_times >= used_windows(w,1) & all_times <= used_windows(w,2));
        end
        nontrace_vals = all_vals(~trial_mask);
        nontrace_vals = nontrace_vals(~isnan(nontrace_vals));

        % Run unpaired t-test
        [~, p] = ttest2(trace_vals, nontrace_vals);

        % Store in struct
        dcount = dcount + 1;
        perDayStruct(dcount).rat = ratName;
        perDayStruct(dcount).date = dateStr;
        perDayStruct(dcount).trace_vel_vals = trace_vals;
        perDayStruct(dcount).nontrace_vel_vals = nontrace_vals;
        perDayStruct(dcount).mean_trace = mean(trace_vals);
        perDayStruct(dcount).mean_nontrace = mean(nontrace_vals);
        perDayStruct(dcount).std_trace = std(trace_vals);
        perDayStruct(dcount).std_nontrace = std(nontrace_vals);
        perDayStruct(dcount).p = p;
    end
end
