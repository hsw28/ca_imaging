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


if nargin < 2
    fs = 2000;
end

% --- Step 1: Low-order causal Butterworth bandpass
%[b_bp, a_bp] = butter(2, [100 999]/(fs/2), 'bandpass');  % Low-order = shorter impulse
%emg_bandpassed = filter(b_bp, a_bp, emg_raw);

% --- Step 2: Rectify
emg_rectified = abs(emg_raw);

% --- Step 3: Causal 10 ms smoothing
window_size = round(0.010 * fs);
kernel = ones(window_size,1) / window_size;
emg_processed = filter(kernel, 1, emg_rectified);  % Causal only


end
