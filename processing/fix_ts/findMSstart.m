function [starttime_raw fixed_events] = findMSstart(event_files)
  %- I believe for the new way I do it the recording starts on either SECOND ZERO or first FOUR, whichever is first


fields = fieldnames(event_files);
starttime_raw = struct();

for i = 1:numel(fields)
    fieldName = fields{i}
    ttl = event_files.(fieldName);
    if length(ttl)<=1
      starttime_raw.(fieldName) = 0;
      event_files.(fieldName) = 0;
      fprintf('empty event')
    else

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
        ts0 = sort(find(ttl(:,2)>0));
        starttime_raw.(fieldName) = ttl(ts0(1),1);

        fixedevents(:,1) = ttl(:,1)-start;
        event_files.(fieldName) = fixedevents;
end
fixed_events = event_files;
%fprintf('events fixed')
end
