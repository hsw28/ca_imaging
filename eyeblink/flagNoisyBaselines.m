function noisy_trials_struct = flagNoisyBaselines(EMG_struct, EMG_ts_struct, cs_times_struct, baseline_duration, noise_threshold)
% FLAGNOISYBASELINES Identifies trials with excessive baseline noise in EMG recordings.
%
% Inputs:
%   EMG_struct        - Struct containing EMG signals for each session.
%   EMG_ts_struct     - Struct containing timestamps corresponding to EMG samples.
%   cs_times_struct   - Struct containing CS onset times for each trial.
%   baseline_duration - Duration of the baseline period before CS onset (in seconds).
%   noise_threshold   - Threshold for flagging noise (e.g., 1.5 for 1.5x mean amplitude).
%
% Output:
%   noisy_trials_struct - Struct indicating noisy trials (true for noisy, false otherwise).

    fields = fieldnames(EMG_struct);
    noisy_trials_struct = struct();

    for i = 1:numel(fields)
        field = fields{i};
        emg_data = EMG_struct.(field);
        timestamps = EMG_ts_struct.(field);
        cs_times = cs_times_struct.(field);

        num_trials = numel(cs_times);
        noisy_trials = false(num_trials, 1);

        for t = 1:num_trials
            cs_time = cs_times(t);
            baseline_start = cs_time - baseline_duration;
            baseline_end = cs_time;

            % Identify indices corresponding to the baseline period
            idx_baseline = timestamps >= baseline_start & timestamps < baseline_end;

            if any(idx_baseline)
                baseline_signal = emg_data(idx_baseline);
                baseline_mean = mean(abs(baseline_signal));
                baseline_std = std(baseline_signal);

                % Flag trial if baseline noise exceeds threshold
                if baseline_std > noise_threshold * baseline_mean
                    noisy_trials(t) = true;
                end
            else
                % If baseline period is not found, consider the trial noisy
                noisy_trials(t) = true;
            end
        end

        noisy_trials_struct.(field) = noisy_trials;
    end
end
