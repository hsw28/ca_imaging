function [US_times_noise, noisy_trials] = flagNoisyBaselines(emg_data, timestamps, US_times)
% FLAGNOISYBASELINES Flags trials with outlier baseline noise compared to group.
%
% Inputs:
%   emg_data       - full EMG trace (vector)
%   timestamps     - timestamps for EMG data
%   US_times       - array of US onset times (vector)
%   noise_threshold - threshold in z-score units
%
% Outputs:
%   US_times_noise - US_times with noisy trials NaN'ed out
%   noisy_trials   - logical array, true = noisy

noise_threshold = .5;
baseline_duration = 0.5;
cs_times = US_times - 0.750;

num_trials = numel(cs_times);
baseline_stds = NaN(num_trials, 1);

% --- Collect baseline stds
for t = 1:num_trials
    cs_time = cs_times(t);
    idx = timestamps >= (cs_time - baseline_duration) & timestamps < cs_time;
    if any(idx)
        baseline_stds(t) = std(emg_data(idx));
    end
end

% --- Z-score baseline stds across trials
z_stds = (baseline_stds - nanmean(baseline_stds)) ./ nanstd(baseline_stds);
noisy_trials = z_stds > noise_threshold;

% --- Output
US_times_noise = US_times;
US_times_noise(noisy_trials) = NaN;

end
