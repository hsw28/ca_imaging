function fixed_US_times = fixUS(starttime_raw_struct, US_times_struct)
  %starttimes is from findMSstart

  fields_US = fieldnames(US_times_struct);
  fields_starts = fieldnames(starttime_raw_struct);

  for i = 1:numel(fields_US)
    fieldName_US = fields_US{i}
    fieldValue_US = US_times_struct.(fieldName_US);
    US = fieldValue_US;

    index = strfind(fieldName_US, '_');
    US_date = fieldName_US(index(2)+1:end);

    fieldName_starts = fields_starts{i};
    fieldValue_starts = starttime_raw_struct.(fieldName_starts);
    starttime_raw = fieldValue_starts;

    index = strfind(fieldName_starts, '_');
    start_date = fieldName_starts(index(2)+1:end);

    if US(1)> 152985146928
      US = US./1000000;
    end

    if strcmp(US_date, start_date)==1

        if starttime_raw ==0
          warning('your start time is zero')
        %  starttime_raw = US(1)+.0024; %If CS
          starttime_raw = US(1)+(1-0024)
          US = US-starttime_raw;
        else
          US = US-starttime_raw;
        end
        if length(US)<5
          US = NaN;
        end
    else
      error('dates do not match' )
    end

    US_times_struct.(fieldName_US) = US;
  end
  fixed_US_times =  US_times_struct;
  fprintf('USs fixed')
end
