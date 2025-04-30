function emg_processed = filter_emg3(emg_raw)
% process_emg_signal - Full EMG processing pipeline
%   - Bandpass 100 Hz to 1000 Hz (adjusted for fs=2000Hz)
%   - Full-wave rectification
%   - Integration with a 10 ms time constant

% Inputs:
%   emg_raw - Raw EMG signal (vector)
%   fs - Sampling frequency (Hz)

% Output:
%   emg_processed - Final processed EMG signal

fs = 2000;

% --- Step 1: Bandpass filter (100 Hz to 1000 Hz or fs/2)
low_cutoff = 100;            % Hz
high_cutoff = min(1000, fs/2-1);  % Hz (cannot exceed fs/2)

[b_bp, a_bp] = butter(4, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
emg_bandpassed = filtfilt(b_bp, a_bp, emg_raw);

% --- Step 2: Full-wave rectification
emg_rectified = abs(emg_bandpassed);

% --- Step 3: Integration with 10 ms time constant
time_constant = 0.010; % 10 ms = 0.010 seconds
fc_integrator = 1 / (2*pi*time_constant);  % Corresponding cutoff frequency

[b_lp, a_lp] = butter(4, fc_integrator/(fs/2), 'low');
emg_processed = filtfilt(b_lp, a_lp, emg_rectified);
