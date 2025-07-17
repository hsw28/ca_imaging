function [t_interp, mean_emg, spread_emg, mean_vel, spread_vel] = alignEMGandVelocity(ratName, spreadType)
    % spreadType: 'sem' (default) or 'std'

    if nargin < 2
        spreadType = 'sem';
    end

    ratName = char(ratName);
    animal = evalin('base', ratName);
    emg = evalin('base', [ratName 'emg']);

    dateList = autoDateList(animal);
    finalDay = animal.An;
    finalIdx = find(strcmp(dateList, finalDay));

    if finalIdx < 3
        warning('%s has fewer than 3 learned sessions', ratName);
        t_interp = []; mean_emg = []; spread_emg = [];
        mean_vel = []; spread_vel = []; vel_trials = [];
        return;
    end

    learnedDays = dateList(finalIdx-2:finalIdx);

    win = [-1 2];
  %  Fs_interp = 7.5;
    Fs_interp = 2000;
    t_interp = win(1):1/Fs_interp:win(2);
    t_interp = round(t_interp, 6);
    nTimepoints = length(t_interp);

    all_emg = [];
    all_vel = [];
    emg_trials = [];

    vel_trials = [];  % NEW: trial-by-trial velocity traces

    for d = 1:length(learnedDays)
        day = learnedDays{d};
        dateStr = day;

        EMG   = emg.EMG.(['EMG_' dateStr]);
        EMGts = emg.EMG_ts.(['EMGts_' dateStr]);
        pos   = animal.pos.(['pos_' dateStr]);
        vel   = ca_velocity(pos);  % [2 × n]

        cs_times = animal.CS_times.(['CS_' dateStr]);

        for i = 1:length(cs_times)
          % Skip first CS for rat0816
            if strcmp(ratName, 'rat0816') && i == 1
                continue;
            end

            cs = cs_times(i);
            t_win = cs + win;

            emg_idx = EMGts >= t_win(1) & EMGts <= t_win(2);
            if sum(emg_idx) < 10, continue; end
            emg_interp = interp1(EMGts(emg_idx), EMG(emg_idx), cs + t_interp, 'linear', 'extrap');

            vel_idx = vel(2,:) >= t_win(1) & vel(2,:) <= t_win(2);
            if sum(vel_idx) < 10, continue; end
            vel_interp = interp1(vel(2, vel_idx), vel(1, vel_idx), cs + t_interp, 'linear', 'extrap');

            all_emg(end+1,:) = emg_interp;
            all_vel(end+1,:) = vel_interp;
            emg_trials(end+1,:) = emg_interp;
            vel_trials(end+1,:) = vel_interp;  % Store trial-wise
        end
    end

    if isempty(all_emg)
        warning('%s: no valid trials found', ratName);
        t_interp = []; mean_emg = []; spread_emg = [];
        mean_vel = []; spread_vel = []; vel_trials = [];
        return;
    end

    mean_emg = mean(all_emg, 1);
    mean_vel = mean(all_vel, 1);

    switch lower(spreadType)
        case 'sem'
            spread_emg = std(all_emg, [], 1) / sqrt(size(all_emg,1));
            spread_vel = std(all_vel, [], 1) / sqrt(size(all_vel,1));
        case 'std'
            spread_emg = std(all_emg, [], 1);
            spread_vel = std(all_vel, [], 1);
        otherwise
            error('spreadType must be "sem" or "std"');
    end

size(t_interp)
end
