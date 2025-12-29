function [t_interp, mean_emg, spread_emg, mean_vel, spread_vel, perDayOut] = alignEMGandVelocity(ratName, spreadType, perDay)
% alignEMGandVelocity
%   Align eyelid EMG + running velocity to CS onset across the last 3 learned sessions.
%
% Inputs
%   ratName    : e.g., 'rat0314'
%   spreadType : 'sem' (default) or 'std'
%   perDay     : false (default) pool across days; true also return per-day outputs
%
% Outputs (pooled across the selected days)
%   t_interp, mean_emg, spread_emg, mean_vel, spread_vel
%
% Output (optional per-day)
%   perDayOut(d).day
%   perDayOut(d).mean_emg, perDayOut(d).spread_emg, perDayOut(d).mean_vel, perDayOut(d).spread_vel
%   perDayOut(d).emg_trials, perDayOut(d).vel_trials

    if nargin < 2 || isempty(spreadType), spreadType = 'sem'; end
    if nargin < 3 || isempty(perDay),     perDay = false;     end

    ratName = char(ratName);
    animal = evalin('base', ratName);
    emg    = evalin('base', [ratName 'emg']);

    dateList = autoDateList(animal);
    finalDay = animal.An;
    finalIdx = find(strcmp(dateList, finalDay));

    perDayOut = struct('day',{},'mean_emg',{},'spread_emg',{},'mean_vel',{},'spread_vel',{}, ...
                       'emg_trials',{},'vel_trials',{});

    if isempty(finalIdx) || finalIdx < 3
        warning('%s has fewer than 3 learned sessions', ratName);
        t_interp = []; mean_emg = []; spread_emg = [];
        mean_vel = []; spread_vel = [];
        return;
    end

    learnedDays = dateList(finalIdx-2:finalIdx);

    win = [-0.5 4];
    Fs_interp = 2000;
    t_interp = win(1):1/Fs_interp:win(2);
    t_interp = round(t_interp, 6);

    % --- pooled containers (across days) ---
    all_emg = [];
    all_vel = [];

    % --- loop days ---
    for d = 1:numel(learnedDays)
        dateStr = learnedDays{d};

        EMG   = emg.EMG.(['EMG_' dateStr]);
        EMGts = emg.EMG_ts.(['EMGts_' dateStr]);
        pos   = animal.pos.(['pos_' dateStr]);
        vel   = ca_velocity(pos);  % [2 x n], vel(1,:)=speed, vel(2,:)=time

        cs_times = animal.CS_times.(['CS_' dateStr]);

        day_emg_trials = [];
        day_vel_trials = [];

        for i = 1:numel(cs_times)
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

            % pooled
            all_emg(end+1,:) = emg_interp; %#ok<AGROW>
            all_vel(end+1,:) = vel_interp; %#ok<AGROW>

            % per-day
            day_emg_trials(end+1,:) = emg_interp; %#ok<AGROW>
            day_vel_trials(end+1,:) = vel_interp; %#ok<AGROW>
        end

        % save per-day if requested
        if perDay
            if isempty(day_emg_trials)
                perDayOut(d).day = dateStr;
                perDayOut(d).mean_emg = [];
                perDayOut(d).spread_emg = [];
                perDayOut(d).mean_vel = [];
                perDayOut(d).spread_vel = [];
                perDayOut(d).emg_trials = [];
                perDayOut(d).vel_trials = [];
            else
                [mE, sE] = local_mean_spread(day_emg_trials, spreadType);
                [mV, sV] = local_mean_spread(day_vel_trials, spreadType);

                perDayOut(d).day = dateStr;
                perDayOut(d).mean_emg = mE;
                perDayOut(d).spread_emg = sE;
                perDayOut(d).mean_vel = mV;
                perDayOut(d).spread_vel = sV;
                perDayOut(d).emg_trials = day_emg_trials;
                perDayOut(d).vel_trials = day_vel_trials;
            end
        end
    end

    if isempty(all_emg)
        warning('%s: no valid trials found', ratName);
        t_interp = []; mean_emg = []; spread_emg = [];
        mean_vel = []; spread_vel = [];
        return;
    end

    % --- pooled outputs (original behavior) ---
    [mean_emg, spread_emg] = local_mean_spread(all_emg, spreadType);
    [mean_vel, spread_vel] = local_mean_spread(all_vel, spreadType);
end

function [mu, sp] = local_mean_spread(X, spreadType)
    mu = mean(X, 1, 'omitnan');
    switch lower(spreadType)
        case 'sem'
            sp = std(X, [], 1, 'omitnan') ./ sqrt(size(X,1));
        case 'std'
            sp = std(X, [], 1, 'omitnan');
        otherwise
            error('spreadType must be "sem" or "std"');
    end
end
