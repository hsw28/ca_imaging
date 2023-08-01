function [fixed_events fixed_EMG_ts fixed_EMG pos_fixed] = fix_all_ts(event_files, EMG_ts, EMG, pos)
%can use this whole function at once or individual functions within
%aligns everything to the same timestamps and starts the timestamps at time 0

%add CA_TS and CA_TS_fixed

%%%%%%%%%%
%%%%%%%%%%MAKE IT CHECK THAT THE NAMES ARE THE SAME
%%%%%%%%%%

%first find the time the miniscope starts recording
[starttime_raw fixed_events] = findMSstart(event_files);

%delete any EMG timestamps or corresponding EMGs that come before that time
%downsamples EMG
[fixed_EMG_ts fixed_EMG] = fixEMG_times(starttime_raw, EMG_ts, EMG);

%position time stamps start at zero so just change them to seconds
pos_fixed = fixPos_times(pos);

%CA time stamps start at zero so just change them to seconds
%CA_TS_fixed = fixCA_times(CA_TS);

end

%- I believe for the new way I do it the recording starts on either SECOND ZERO or first FOUR, whichever is first
function [starttime_raw fixed_events] = findMSstart(event_files)
fields = fieldnames(event_files);
starttime_raw = struct();

for i = 1:numel(fields)
    fieldName = fields{i};
    ttl = event_files.(fieldName);

    if ttl(1,1)> 152985146928
      ttl(:,1) = ttl(:,1)./1000000;
    end
    ts0 = sort(find(ttl(:,2)==0));
    ts4 = sort(find(ttl(:,2)==4));
    ts0 = ts0(2);
    ts4 = ts4(1);
    if ts0<ts4
      start = ttl(ts0,1);
      starttime_raw.(fieldName) = ttl(ts0,1);
    else
      start = ttl(ts4,1);
      starttime_raw.(fieldName) = ttl(ts4,1);
    end
    fixedevents = ttl;
    fixedevents(:,1) = ttl(:,1)-start;
    event_files.(fieldName) = fixedevents;
end
fixed_events = event_files;
fprintf('events fixed')
end



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




function pos_fixed = fixPos_times(pos)
fields = fieldnames(pos);
for i = 1:numel(fields)
    fieldName = fields{i};
    current_pos = pos.(fieldName);
  if current_pos(1,1)>1 || current_pos(2,1)>1
    current_pos(:,1) = current_pos(:,1)/ 1000;
  end
pos.(fieldName) = current_pos;
end
pos_fixed = pos;
fprintf('positions fixed')
end

function CA_TS_fixed = fixCA_times(CA_TS)
fields = fieldnames(CA_TS);
  for i = 1:numel(fields)
    fieldName = fields{i};
    current_CA = CA_TS.(fieldName);
  if current_CA(1)>1 || current_CA(2)>1
    current_CA(1,:) = current_CA(1,:)/ 1000;
  end
CA_TS.(fieldName) = current_CA;
end
CA_TS_fixed = CA_TS;
fprintf('calcium times fixed')
end
