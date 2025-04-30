function EMG_final = filter_emg2(EMG_signal)

Fs = 2000;                 % Sampling frequency in Hz (adjust to your data)

% Step 1: Bandpass Filter Design (Butterworth)
low_cutoff = 100;   % Hz, adjust as needed (typically 30-100 Hz lower bound)
high_cutoff = 400; % Hz upper bound 1000
filter_order = 4;   % Standard Butterworth order for EMG

% Normalize frequencies by Nyquist frequency (Fs/2)
[b_bp, a_bp] = butter(filter_order, [low_cutoff, high_cutoff]/(Fs/2), 'bandpass');

% Apply bandpass filter using zero-phase filtering (filtfilt)
EMG_bandpassed = filtfilt(b_bp, a_bp, EMG_signal);

% Step 2: Rectification
EMG_rectified = abs(EMG_bandpassed);

% Step 3: Low-pass filter to extract envelope (e.g., 30 Hz)
envelope_cutoff = 30; % Hz, adjust as needed (~20-50 Hz typical, start w 30)
[b_lp, a_lp] = butter(filter_order, envelope_cutoff/(Fs/2), 'low');

% Apply low-pass filter (zero-phase filtering)
EMG_envelope = filtfilt(b_lp, a_lp, EMG_rectified);

% Optional: normalization (z-score or baseline normalization)
EMG_final = (EMG_envelope - mean(EMG_envelope)) / std(EMG_envelope);
