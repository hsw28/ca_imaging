function [fixed_events fixed_EMG_ts fixed_EMG fixed_US_times fixed_pos] = fix_all_ts(event_files, EMG_ts, EMG, US_times, pos)
%can use this whole function at once or individual functions within
%aligns everything to the same timestamps and starts the timestamps at time 0

%add CA_TS and CA_TS_fixed

%first find the time the miniscope starts recording
[starttime_raw fixed_events] = findMSstart(event_files);

%delete any EMG timestamps or corresponding EMGs that come before that time
%downsamples EMG
[fixed_EMG_ts fixed_EMG] = fixEMG_times(starttime_raw, EMG_ts, EMG);

fixed_US_times = fixUS(starttime_raw, US_times)

%position time stamps start at zero so just change them to seconds
fixed_pos = fixPos_times(pos);

%CA time stamps start at zero so just change them to seconds
%CA_TS_fixed = fixCA_times(CA_TS);
