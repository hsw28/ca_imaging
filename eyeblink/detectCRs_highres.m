function [EMG_cr_detected, EMG_cr_onset_times, EMG_cr_types] = detectCRs_fromUS_struct_highres(EMG_struct, EMG_ts_struct, us_times_struct)

%%only use if you are using the original 32k emg signal

% Parameters
fs_emg = 32000;
STD_THRESHOLD = 4;
MIN_CR_DURATION = 0.010;   % 10 ms
MIN_LATENCY = 0.035;       % 35 ms
MAX_LATENCY = 0.740;       % 740 ms after CS
ADAPTIVE_WINDOW = 0.020;   % 20 ms before US
BASELINE_DURATION = 0.5;   % 500 ms
US_DELAY = 0.75;           % CS to US delay (fixed)
dt = 1 / fs_emg;

% Struct outputs
EMG_cr_detected = struct();
EMG_cr_onset_times = struct();
EMG_cr_types = struct();

fields_EMG = fieldnames(EMG_struct);
fields_TS = fieldnames(EMG_ts_struct);
fields_US = fieldnames(us_times_struct);

for i = 1:numel(fields_EMG)
    field_EMG = fields_EMG{i};
    field_TS = fields_TS{i};
    field_US = fields_US{i};

    emg_data = EMG_struct.(field_EMG);
    ts_low = EMG_ts_struct.(field_TS);
    us_times = us_times_struct.(field_US);
    cs_times = us_times - US_DELAY;

    % Interpolate low-res timestamps to EMG resolution
    valid_idx = ~isnan(ts_low);
    ts_low_clean = ts_low(valid_idx);
    sample_idx_low = linspace(1, numel(emg_data), numel(ts_low_clean));
    ts_interp = interp1(sample_idx_low, ts_low_clean, 1:numel(emg_data), 'linear', 'extrap');
    ts_high = ts_interp(:);

    num_trials = numel(cs_times);
    cr_detected = false(num_trials, 1);
    cr_onsets = nan(num_trials, 1);
    cr_types = strings(num_trials, 1);

    for t = 1:num_trials
        cs = cs_times(t);
        us = us_times(t);

        idx_baseline = ts_high >= (cs - BASELINE_DURATION) & ts_high < cs;
        idx_response = ts_high >= (cs + MIN_LATENCY) & ts_high < (cs + MAX_LATENCY);
        idx_adaptive = ts_high >= (us - ADAPTIVE_WINDOW) & ts_high < us;

        baseline = emg_data(idx_baseline);
        threshold = nanmean(baseline) + STD_THRESHOLD * nanstd(baseline);

        response = emg_data(idx_response);
        t_response = ts_high(idx_response);
        above_thresh = response > threshold;

        [starts, lengths, times] = findRuns(above_thresh, t_response);
        durations = lengths * dt;
        valid = durations >= MIN_CR_DURATION;

        found = false;

        for j = find(valid)'
            onset_time = times(j);
            rel_onset = onset_time - cs;

            if rel_onset < MIN_LATENCY
                continue;  % Alpha response
            end

            if any(emg_data(idx_adaptive) > threshold)
                cr_detected(t) = true;
                cr_onsets(t) = onset_time;
                cr_types(t) = "adaptive";
                found = true;
                break;
            end
        end
    end

    EMG_cr_detected.(field_TS) = cr_detected;
    EMG_cr_onset_times.(field_EMG) = cr_onsets;
    EMG_cr_types.(field_TS) = cr_types;
end
end

function [starts, lengths, times] = findRuns(binary, t)
d = diff([0; binary(:); 0]);
starts = find(d == 1);
ends = find(d == -1) - 1;
lengths = ends - starts + 1;
times = t(starts);
end
