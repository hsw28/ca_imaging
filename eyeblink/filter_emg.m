function emg_processed = filter_emg(emg_raw)
% process_emg_signal - Full EMG processing pipeline
%   - Bandpass 100 Hz to 1000 Hz (adjusted for fs=2000Hz)
%   - Full-wave rectification
%   - Integration with a 10 ms time constant

% Inputs:
%   emg_raw - Raw EMG signal (vector)
%   fs - Sampling frequency (Hz)

% Output:
%   emg_processed - Final processed EMG signal

%fs = 3.2000e+04;
fs=2000;

% --- Step 1: Bandpass filter (100 Hz to 1000 Hz or fs/2)
low_cutoff = 100;            % Hz
high_cutoff = min(1000, fs/2-1);  % Hz (cannot exceed fs/2)
%high_cutoff = 5000;  % Hz (cannot exceed fs/2)

[b_bp, a_bp] = butter(4, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
emg_bandpassed = filtfilt(b_bp, a_bp, emg_raw);

% --- Step 2: Full-wave rectification
emg_rectified = abs(emg_bandpassed);

% --- Step 3: Integration with 10 ms time constant
window_size = max(1, round(0.010 * fs));  % 10 ms
emg_processed = movmean(emg_rectified, window_size);
