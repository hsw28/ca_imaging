function filtered_emg = filter_emg2_params(raw_emg, params)
    % params should include:
    % params.low, params.high (bandpass)
    % params.env (low-pass envelope)
    % params.Fs (sampling frequency)

    % Bandpass filtering
    [b_bp, a_bp] = butter(4, [params.low, params.high]/(params.Fs/2), 'bandpass');
    emg_bandpassed = filtfilt(b_bp, a_bp, raw_emg);

    % Rectification
    emg_rectified = abs(emg_bandpassed);

    % Envelope extraction (low-pass)
    [b_env, a_env] = butter(4, params.env/(params.Fs/2), 'low');
    filtered_emg = filtfilt(b_env, a_env, emg_rectified);
end
