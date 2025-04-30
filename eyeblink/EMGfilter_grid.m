% Define parameter ranges for grid search
low_cutoffs = [30, 50, 75, 100];
high_cutoffs = [400, 500, 750, 990];
envelope_cutoffs = [65,70,80];

Fs = 2000;

% Initialize indexed storage arrays
param_list = {};
perc_correct_list = [];

% Original EMG data

raw_emg = test.EMG;
timestamps = test.EMG_ts;
CS_times = test.CS;

count = 1; % Explicit indexing

for i = 1:length(low_cutoffs)
    for j = 1:length(high_cutoffs)
        if high_cutoffs(j) <= low_cutoffs(i)
            continue; % Skip invalid frequencies
        end
        for k = 1:length(envelope_cutoffs)

            % Filtering parameters
            params.low = low_cutoffs(i);
            params.high = high_cutoffs(j);
            params.env = envelope_cutoffs(k);
            params.Fs = Fs;

            % Apply filtering
            nm = fieldnames(raw_emg);
            remg = raw_emg.(nm{1});

            filtered_emg = filter_emg2_params(remg, params);

            % Assign for CR detection compatibility
            test.temp_filtered.(nm{1}) = filtered_emg;

            % CR detection
            [EMG_cr_detected, ~] = detectCRs(test.temp_filtered, timestamps, CS_times);
            nm = fieldnames(EMG_cr_detected);
            cr = EMG_cr_detected.(nm{1});

            % Calculate percentage correct
            perc_correct = length(find(cr > 0)) / 50;

            % Store results explicitly indexed
            param_list{count, 1} = sprintf('low%d_high%d_env%d', params.low, params.high, params.env);
            perc_correct_list(count, 1) = perc_correct;

            count = count + 1; % Increment counter
        end
    end
end

% Create aligned result table
result_table = table(param_list, perc_correct_list, ...
    'VariableNames', {'Parameters', 'Perc_Correct'});

disp(result_table);
