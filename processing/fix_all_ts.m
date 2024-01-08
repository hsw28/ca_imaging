function [fixed_events fixed_EMG_ts fixed_EMG fixed_US_times fixed_CS_times fixed_pos] = fix_all_ts(event_files, EMG_ts, EMG, CS_times, US_times, pos)
%can use this whole function at once or individual functions within
%aligns everything to the same timestamps and starts the timestamps at time 0

%add CA_TS and CA_TS_fixed
%add CA_peak_data

%first find the time the miniscope starts recording
[starttime_raw fixed_events] = findMSstart(event_files);

%delete any EMG timestamps or corresponding EMGs that come before that time
%downsamples EMG
[fixed_EMG_ts fixed_EMG] = fixEMG_times(starttime_raw, EMG_ts, EMG);

fixed_US_times = fixUS(starttime_raw, US_times)
fixed_CS_times = fixUS(starttime_raw, CS_times)

%position time stamps start at zero so just change them to seconds.
%also changes px to cm
fixed_pos = fixPos_times(pos);

%CA time stamps start at zero so just change them to seconds
%do the same with peak data
%CA_TS_fixed = fixCA_times(CA_TS);
