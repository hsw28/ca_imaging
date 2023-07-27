function [fixed_events fixed_EMG_ts fixed_EMG pos_fixed CA_TS_fixed] = fix_all_ts(event_files, EMG_ts, EMG, pos, CA_TS)
%can use this whole function at once or individual functions within
%aligns everything to the same timestamps and starts the timestamps at time 0

%first find the time the miniscope starts recording
[starttime_raw fixed_events] = findMSstart(event_files)

%delete any EMG timestamps or corresponding EMGs that come before that time
[fixed_EMG_ts fixed_EMG] = fixEMGs(starttime_raw, EMG_ts, EMG);

%position time stamps start at zero so just change them to seconds
pos_fixed = fixPos_times(pos);

%CA time stamps start at zero so just change them to seconds
CA_TS_fixed = fixPos_times(CA_TS);

%%%%%%%% CHECK DIMENSIONS ON THE POS FILE
%%%%%%%% CHECK TO MAKE SURE SIZES FOR EMG/EEG/EVENTS IS CORRECT FOR GOING DOWN TO SECONDS

%- I believe for the new way I do it the recording starts on either SECOND ZERO or first FOUR, whichever is first
function [starttime_raw fixed_events] = findMSstart(event_files)
fields = fieldnames(event_files);
starttime_raw = struct;

for i = 1:numel(fields)
    fieldName = fields{i};
    ttl = myStruct.(fieldName);

    if ttl(1,1)> 1000000000
      eventsfile(:,1) = eventsfile(:,1)./1000000;
    end
    ttl = eventsfile;
    ts0 = sort(find(ttl(:,2)==0));
    ts4 = sort(find(ttl(:,2)==4));
    ts0 = ts0(2);
    ts4 = ts4(1);
    if ts0<ts04
      starttime_raw = ttl(ts0,1);
    else
      starttime_raw = ttl(ts1,1);
    end
    fixedevents = eventsfile(:,1)-starttime_raw;

    event_files.(fieldName) = fixedevents;
    starttime_raw.(fieldName) = starttime_raw;
end
fixed_events = event_files;
end



function [fixed_EMG_ts fixed_EMG] = fixEMG_times(starttime_raw, EMG_ts, EMG)
fields_TS = fieldnames(EMG_ts);
fields_EMG = fieldnames(EMG);
fields_starts = fieldnames(starttime_raw);

for i = 1:numel(fields)
  fieldName_TS = fields_TS{i};
  fieldValue_TS = myStruct.(fieldName_TS);
  EMG_ts = fieldValue_EMG_TS;

  fieldName_EMG = fields_EMG{i};
  fieldValue_EMG = myStruct.(fieldName_EMG);
  EMG = fieldValue_EMG;

  fieldName_starts = fields_starts{i};
  fieldValue_starts = myStruct.(fieldName_starts);
  starttime_raw = fieldValue_starts;

  if EMG_ts(1,1)> 1000000000
    EMG_ts(:,1) = EMG_ts(:,1)./1000000;
  end
  dontwant = find(EMG_ts<starttime_raw);
  EMG_ts = EMG_ts(max(dontwant+1):end);
  EMG_ts = EMG_ts-starttime_raw;
  EMG = EMG(max(dontwant+1):end);

  fixed_EMG_ts.(fieldName_TS) = EMG_ts;
  fixed_EMG.(fieldName_EMG) = EMG;

end
end



function pos_fixed = fixPos_times(pos)
fields = fieldnames(pos);
for i = 1:numel(fields)
    fieldName = fields{i};
    current_pos = myStruct.(fieldName);
  if current_pos(1)>1 || current_pos(2)>1
    current_pos(1,:) = current_pos(1,:)/ 1000;
  end
pos.(fieldName) = current_pos;
end
pos_fixed = pos;
end

function CA_TS_fixed = fixPos_times(CA_TS)
fields = fieldnames(CA_TS);
  for i = 1:numel(fields)
  current_CA = fields{i};
  if current_CA(1)>1 || current_CA(2)>1
    current_CA(1,:) = current_CA(1,:)/ 1000;
  end
CA_TS.(fieldName) = current_CA;
end
CA_TS_fixed = CA_TS
end
