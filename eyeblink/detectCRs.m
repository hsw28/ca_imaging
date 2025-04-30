function [EMG_cr_detected, EMG_cr_onset_times] = detectCRs(EMG_struct, EMG_ts_struct, cs_times_struct)
% DETECTCRS Detects conditioned responses in EMG data
%   EMG_struct: Struct containing raw EMG signals
%   EMG_ts_struct: Struct containing timestamps corresponding to each EMG sample (in seconds)
%   cs_times_struct: Struct containing CS onset times (in seconds)
%
% Returns:
%   EMG_cr_detected: Struct indicating presence of CR for each trial
%   EMG_cr_onset_times: Struct of CR onset times (in seconds) for each trial

% Define timing parameters (in seconds)
BASELINE_DURATION = .5;    % 500ms baseline
CS_DURATION = 0.25;         % 250ms CS
TRACE_DURATION = 0.5;       % 500ms trace period
US_DURATION = 0.1;          % 100ms US
MIN_CR_DURATION = 0.015;    % 15ms minimum duration
MIN_LATENCY = 0.05;         % 50ms minimum latency after CS
STD_THRESHOLD = 4;          % 4 SD above baseline mean
adaptive = 0.2;

fields_EMG = fieldnames(EMG_struct);
fields_CS = fieldnames(cs_times_struct);
fields_TS = fieldnames(EMG_ts_struct);

if numel(fields_EMG) ~= numel(fields_CS)
    error('EMG and CS structures do not have the same number of fields.');
end

if numel(fields_TS) ~= numel(fields_CS)
    error('TS and CS structures do not have the same number of fields.');
end

EMG_cr_detected = struct();
EMG_cr_onset_times = struct();

for i = 1:numel(fields_TS)
    fieldName_TS = fields_TS{i};
    timestamps = EMG_ts_struct.(fieldName_TS);

    index = strfind(fieldName_TS, '_');
    TS_date = fieldName_TS(index(2)+1:end);

    fieldName_EMG = fields_EMG{i};
    emg_data = EMG_struct.(fieldName_EMG);

    index = strfind(fieldName_EMG, '_');
    EMG_date = fieldName_EMG(index(2)+1:end);

    fieldName_CS = fields_CS{i};
    cs_times = cs_times_struct.(fieldName_CS);

    index = strfind(fieldName_CS, '_');
    CS_date = fieldName_CS(index(2)+1:end);

    if strcmp(TS_date, EMG_date) && strcmp(TS_date, CS_date)
        num_trials = length(cs_times);
        cr_detected = false(num_trials, 1);
        cr_onset_times = zeros(num_trials, 1);

        for trial = 1:num_trials
            % Find indices for baseline period (500ms before CS)
            baseline_idx = timestamps >= (cs_times(trial) - BASELINE_DURATION) & ...
                           timestamps < cs_times(trial);
            baseline = emg_data(baseline_idx);


            % Calculate baseline statistics
            baseline_mean = mean(baseline);
            baseline_std = std(baseline);
            threshold = baseline_mean + (STD_THRESHOLD * baseline_std);

            % Define response window (from CS + 50ms to CS + trace period end)
            response_idx = timestamps >= (cs_times(trial) + MIN_LATENCY ) & ...
                           timestamps < (cs_times(trial) + CS_DURATION + TRACE_DURATION);

            %late CRs
            %response_idx2 = timestamps >= (cs_times(trial) + CS_DURATION + TRACE_DURATION - adaptive) & ...
            %               timestamps < (cs_times(trial) + CS_DURATION + TRACE_DURATION);

            response_window = emg_data(response_idx);
            response_times = timestamps(response_idx);

            % Find periods where EMG exceeds threshold
            above_threshold = response_window > threshold;

            % Use runs test to find consecutive samples above threshold
            [runs_start, runs_length, run_times] = findRuns(above_threshold, response_times);

            % Convert run lengths from samples to duration in seconds
            run_durations = runs_length .* nanmean(diff(timestamps))

            % Check if any run meets the duration criterion
            valid_runs = run_durations >= MIN_CR_DURATION;

            if any(valid_runs)
                % CR detected - use first valid run
                first_valid_run = find(valid_runs, 1);
                cr_detected(trial) = true;
                cr_onset_times(trial) = run_times(first_valid_run);
            end
        end

        EMG_cr_detected.(fieldName_TS) = cr_detected;
        EMG_cr_onset_times.(fieldName_EMG) = cr_onset_times;
    else
        error('Dates do not match: TS_date = %s, EMG_date = %s, CS_date = %s', TS_date, EMG_date, CS_date);
    end
end


function [runs_start, runs_length, run_times] = findRuns(binary_signal, times)
% Helper function to find runs of ones in a binary signal
% Returns the start indices, lengths, and start times of runs

% Find transitions
transitions = diff([0; binary_signal; 0]);
run_starts = find(transitions == 1);
run_ends = find(transitions == -1) - 1;

runs_start = run_starts;
runs_length = run_ends - run_starts + 1;
run_times = times(run_starts);
