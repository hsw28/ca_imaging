function [EMG_cr_detected, percentCRs, EMG_cr_onset_times, EMG_cr_types] = detectCRs(EMG_struct_unfiltered, EMG_ts_struct, us_times_struct, debug)

%Detects conditioned responses in EMG data
% filters EMG data using filter_emg
% tags noisy trials using flagNoisyBaselines
%determines if there was a CR response

%A CR was defined as an increase in integrated EMG activity that was greater than the mean baseline amplitude plus
%4 times the standard deviation of the baseline activity, for a minimum of 10 ms. The response also had to begin at
%least 35 ms after the onset of the CS, but before US onset. A response in the 200 ms before the onset of the US
%was defined as a late or adaptive CR. Finally, a short-latency or alpha response was defined as one that began
%within 35 ms of CS onset.
%
%CRs that returned to baseline before the onset of the US were considered poorly timed or “non-adaptive” CRs and
%were not counted for analysis. Any CR that began <35 ms after CS onset was considered an “alpha” response and was not counted.

%   EMG_struct: Struct containing raw EMG signals -- SIGNALS ARE FILTERED IN THIS FUNCTION
%   EMG_ts_struct: Struct containing timestamps corresponding to each EMG sample (in seconds)
%   cs_times_struct: Struct containing US onset times (in seconds) -- note this corrects for the 0.0083 offset
%
% Returns:
%   EMG_cr_detected: Struct indicating presence of CR for each trial
%   EMG_cr_onset_times: Struct of CR onset times (in seconds) for each trial



EMG_struct = EMG_struct_unfiltered;

if nargin < 4
    debug = false;
end




% Constants
BASELINE_DURATION = .5;
US_DELAY = 0.750;                  % Time from CS to US
MIN_LATENCY = 0.035;
MAX_LATENCY = 0.740;
ADAPTIVE_WINDOW_START = 0.73;
MIN_CR_DURATION = 0.010;
STD_THRESHOLD = 4;

fields_EMG = fieldnames(EMG_struct);
fields_TS = fieldnames(EMG_ts_struct);
fields_US = fieldnames(us_times_struct);

if length(fields_EMG) ~= length(fields_TS) | length(fields_EMG) ~= length(fields_US)
  error('your structure lengths do night align')
end

% Output structs
EMG_cr_detected = struct();
EMG_cr_onset_times = struct();
EMG_cr_types = struct();

for i = 1:numel(fields_TS)
    ts_field = fields_TS{i};
    timestamps = EMG_ts_struct.(ts_field);
    emg_field = fields_EMG{i};
    emg_data = EMG_struct.(emg_field);
    us_field = fields_US{i};
    us_times = us_times_struct.(us_field);


    num_trials = numel(us_times);
    cr_detected = zeros(num_trials, 1);
    cr_onsets = cell(num_trials, 1);
    cr_types = cell(num_trials, 1);

    if length(emg_data)>1 & length(us_times)>1
    emg_data = filter_emg(emg_data);

    [US_times_noise noisy_trials] = flagNoisyBaselines(emg_data, timestamps, us_times);

    us_times = US_times_noise;
    cs_times = us_times - US_DELAY;

      for t = 1:num_trials
          cs = cs_times(t);
          us = us_times(t);

          if isnan(cs)==1 %this means noisy trial
            cr_detected(t) = NaN;
            cr_onsets{t} = NaN;
            cr_types{t} = NaN;
            continue;
          end


          % Baseline
          baseline_idx = timestamps >= (cs - BASELINE_DURATION) & timestamps < cs;
          baseline = emg_data(baseline_idx);
          threshold = nanmean(baseline) + STD_THRESHOLD * nanstd(baseline);

          % Response window
          response_idx = timestamps >= cs & timestamps < us;
          response_window = emg_data(response_idx);
          response_times = timestamps(response_idx);
          above_threshold = response_window > threshold;

          [starts, lengths, times] = findRuns(above_threshold, response_times);
          durations = lengths .* nanmean(diff(timestamps));
          valid_idx = durations >= MIN_CR_DURATION;

          trial_onsets = [];
          trial_types = {};

          for j = find(valid_idx)'
              onset = times(j);
              rel_onset = onset - cs;

              if rel_onset < MIN_LATENCY
                  continue;  % alpha
              end

              % Check adaptive: still above threshold in last 20 ms before US
              %adaptive_idx = timestamps >= (cs) & timestamps < us;
              adaptive_idx = timestamps >= (us - 0.200) & timestamps < us;
              if any(emg_data(adaptive_idx) > threshold)
                  trial_onsets(end+1) = onset;
                  trial_types{end+1} = "adaptive";
              else
                  % Non-adaptive: skip
                  continue;
              end
          end

          if ~isempty(trial_onsets)
              cr_detected(t) = 1;
          end

          cr_onsets{t} = trial_onsets;
          cr_types{t} = trial_types;

          % Optional plotting
          if debug
              figure;
              plot(timestamps, emg_data); hold on;
              yline(threshold, 'r--');
              xline(cs, 'k--', 'CS');
              xline(us, 'm--', 'US');
              for m = 1:length(trial_onsets)
                  xline(trial_onsets(m), '--g', trial_types{m});
              end
              xline(us - 0.020, 'c--', 'US - 20 ms');
              xlabel('Time (s)'); ylabel('EMG');
          end
      end

      EMG_cr_detected.(ts_field) = cr_detected;
      EMG_cr_onset_times.(emg_field) = cr_onsets;
      EMG_cr_types.(ts_field) = cr_types;

      percentCRs.(ts_field) = length(find(cr_detected==1))./length(find(isnan(cr_detected)==0));

    else
      EMG_cr_detected.(ts_field) = NaN;
      EMG_cr_onset_times.(emg_field) = NaN;
      EMG_cr_types.(ts_field) = NaN;
      percentCRs.(ts_field) = NaN;
    end

end


function [starts, lengths, times] = findRuns(binary, t)
d = diff([0; binary(:); 0]);
starts = find(d == 1);
ends = find(d == -1) - 1;
lengths = ends - starts + 1;
times = t(starts);
