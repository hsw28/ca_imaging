function emg_filtered = filter_emg(emg_raw)

%% Continuous EMG data were initially filtered using a 4th order Butterworth band‐pass filter with
%cutoff frequency of 28 Hz and 250 Hz, following previous literature (Barker, Reeb‐Sutherland, & Fox, 2014).
%Mains noise was removed with a 50 Hz notch filter. Filtered continuous EMG signals were rectified and
%smoothed using a 4th order Butterworth low‐pass filter with a time constant of 3 ms corresponding to a
%cutoff frequency of 53.05 Hz (Blumenthal et al., 2005).
%https://pmc.ncbi.nlm.nih.gov/articles/PMC5298047/

% preprocess_emg - Preprocess EMG signal
%
% Inputs:
%   emg_raw - Raw EMG signal vector
%   fs      - Sampling frequency in Hz
%
% Output:
%   emg_filtered - Preprocessed EMG signal

fs = 2000;

% Design a 4th-order Butterworth band-pass filter (28–250 Hz)
[b_bp, a_bp] = butter(4, [28 250]/(fs/2), 'bandpass');
emg_bp = filtfilt(b_bp, a_bp, emg_raw);

% Design a 50 Hz notch filter
wo = 50/(fs/2);  % Normalized frequency
bw = wo/35;      % Bandwidth
[b_notch, a_notch] = iirnotch(wo, bw);
emg_notched = filtfilt(b_notch, a_notch, emg_bp);

% Rectify the signal
emg_rectified = abs(emg_notched);

% Design a 4th-order Butterworth low-pass filter with 53.05 Hz cutoff
[b_lp, a_lp] = butter(4, 53.05/(fs/2), 'low');
emg_filtered = filtfilt(b_lp, a_lp, emg_rectified);

end
