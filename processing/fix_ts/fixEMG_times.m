function [fixed_EMG_ts fixed_EMG] = fixEMG_times(starttime_raw_struct, EMG_ts_struct, EMG_struct)
fields_TS = fieldnames(EMG_ts_struct);
fields_EMG = fieldnames(EMG_struct);
fields_starts = fieldnames(starttime_raw_struct);

for i = 1:numel(fields_TS)
  fieldName_TS = fields_TS{i}
  fieldValue_TS = EMG_ts_struct.(fieldName_TS);
  EMG_ts = fieldValue_TS;

  index = strfind(fieldName_TS, '_');
  TS_date = fieldName_TS(index(2)+1:end);


  fieldName_EMG = fields_EMG{i};
  fieldValue_EMG = EMG_struct.(fieldName_EMG);
  EMG = fieldValue_EMG;

  index = strfind(fieldName_EMG, '_');
  EMG_date = fieldName_EMG(index(2)+1:end);

  fieldName_starts = fields_starts{i};
  fieldValue_starts = starttime_raw_struct.(fieldName_starts);
  starttime_raw = fieldValue_starts;

  index = strfind(fieldName_starts, '_');
  start_date = fieldName_starts(index(2)+1:end);

  if strcmp(TS_date, EMG_date)==1 && strcmp(TS_date, start_date)==1

    if EMG_ts(1,1)> 152985146928
      EMG_ts(:,1) = EMG_ts(:,1)./1000000;
    end

    if length(EMG)./length(EMG_ts)>500 %this means time to filter and resample and subtract start time
      %original sampling is 60hz for TS and 512x that for EMG
      %i want them both sampled at 1920, so ./16 for EMG and x32 for TS

      %first lets filter EMG
      EMG = EMGfilter150_960(EMG); %150-960hz using a bandpass filter with a blackman window order 100

      %then lets downsample EMG
      % Downsample the filtered signal to the target sample rate
      downsample_factor = 30720 / 1920;
      EMG = downsample(EMG, downsample_factor);


      %then upsample EMG_ts
      Fs_original = 60;
      Fs_target = 1920;
      t_original = (0:numel(EMG_ts)-1) / Fs_original;

      % Resample the signal using the 'resample' function
      t_target = (0:(numel(EMG_ts)*Fs_target/Fs_original)-1) / Fs_target;
      EMG_ts = interp1(t_original, EMG_ts, t_target, 'linear');


    %and sync start times
    %dontwant = find(EMG_ts<starttime_raw);
    %EMG_ts = EMG_ts(max(dontwant)+1:end);
    %EMG_ts = EMG_ts-starttime_raw;
    %EMG = EMG(max(dontwant+1):end);
    %EMG_ts = EMG_ts(1:length(EMG));
    EMG_ts = EMG_ts-starttime_raw;
    EMG_ts = EMG_ts(1:length(EMG));
    EMG_ts_struct.(fieldName_TS) = EMG_ts';
    EMG_struct.(fieldName_EMG) = EMG;
  end

  else
    TS_date
    EMG_date
    start_date
    error('dates do not match' )
  end

end
fixed_EMG_ts = EMG_ts_struct;
fixed_EMG = EMG_struct;

fprintf('EMGs fixed')
end
