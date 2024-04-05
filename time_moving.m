function moving_percent = time_moving(pos_struct, CSUS_id, velthreshold)
  %gives you the percent of time the animal is moving

fields_pos = fieldnames(pos_struct);
fields_id =  fieldnames(CSUS_id);



for i = 1:numel(fields_pos)

    fieldName_pos = fields_pos{i};
    fieldValue_pos = pos_struct.(fieldName_pos);
    pos = fieldValue_pos;

    fieldName_id = fields_id{i};
    fieldValue_id = CSUS_id.(fieldName_id);
    ID = fieldValue_id;

    index = strfind(fieldName_pos, '_');
    spikes_date = fieldName_pos(index(2)+1:end)

    pos = smoothpos(pos);




    trials = find(ID(1,:)~=0);
    if length(trials)>10
      start = 1;
      stop = max(trials);
      stop = length(pos);
    else
      start = 1;
      stop = length(pos);
    end

    if stop>length(pos)
      stop = length(pos)
      warning('your stop time is after your position')
    end
    pos = pos([start:stop],:);
    vel = ca_velocity(pos);
    goodvel = find(vel(1,:)>=velthreshold);
    percent = length(goodvel)./length(vel);
    av_speed = nanmean(vel(1,:));

    moving_percent.(sprintf('moving_%s', spikes_date)) = [av_speed, percent];
end
