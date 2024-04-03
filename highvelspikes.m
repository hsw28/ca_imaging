function f = highvelspikes(spike_structure, pos_structure, dim, velthreshold, CA_timestamps)

if dim>velthreshold
  warning('did you reverse your dimension and vel threshold?')
end

  fields_spikes = fieldnames(spike_structure);
  fields_pos = fieldnames(pos_structure);
  fields_cats = fieldnames(CA_timestamps);

  if numel(fields_spikes) ~= numel(fields_pos)
    error('your spike and pos structures do not have the same number of values. you may need to pad your US structure for exploration days')
  end

  if numel(fields_pos) ~= numel(fields_cats)
    error('your pos and timestamp structures do not have the same number of values. you may need to pad your US structure for exploration days')
  end


  for i = 1:numel(fields_spikes)

        fieldName_spikes = fields_spikes{i};
        fieldValue_spikes = spike_structure.(fieldName_spikes);
        peaks_time = fieldValue_spikes;


        index = strfind(fieldName_spikes, '_');
        spikes_date = fieldName_spikes(index(2)+1:end)

        fieldName_pos = fields_pos{i};
        fieldValue_pos = pos_structure.(fieldName_pos);
        pos = fieldValue_pos;

        if size(pos,2)>3
          error('you are not using a fixed position')
        end

        if size(pos,2)>size(pos,1)
          pos = pos';
        end

        fieldName_cats = fields_cats{i};
        curr_CA_timestamps = CA_timestamps.(fieldName_cats);

        index = strfind(fieldName_spikes, '_');
        pos_date = fieldName_spikes(index(2)+1:end);

        if length(peaks_time) <5
          highspeedspikes_struct.(sprintf('MI_%s', spikes_date)) = NaN;
          continue
        end

        if (pos(1,1)-pos(end,1))./length(pos) < 1
          pos = convertpostoframe(pos, curr_CA_timestamps);
        end

        mutinfo = NaN(size(peaks_time,1),1);

        tm = pos(:, 1);
        biggest = max(peaks_time(:));
        [minValue,closestIndex] = min(abs(biggest-tm));
        pos = pos(1:closestIndex, :);


  velthreshold = 2;
  vel = ca_velocity(pos);
  goodvel = find(vel(1,:)>=velthreshold);
  goodtime = pos(goodvel, 1);
  goodpos = pos(goodvel,:);


  mintime = vel(2,1);
  maxtime = vel(2,end);
  tm = vel(2,:);

  numunits = size(peaks_time,1);

  HSS = NaN(numunits, 1);
  if numunits<=1
    highspeedspikes_struct.(sprintf('MI_%s', spikes_date)) = NaN;
    warning('you have no spikes')
  else
  for k=1:numunits
    highspeedspikes = [];

    [c indexmin] = (min(abs(peaks_time(k,:)-mintime))); %
    [c indexmax] = (min(abs(peaks_time(k,:)-maxtime))); %
    currspikes = peaks_time(k,indexmin:indexmax);


    for ii=1:length(currspikes) %finding if in good vel
      [minValue,closestIndex] = min(abs(currspikes(ii)-goodtime));

      if minValue <= 1/15 %if spike is within 1 second of moving. no idea if good time
        highspeedspikes(end+1) = currspikes(ii);
      end;
    end

  %want highspeedspikes


    HSS(k) = length(highspeedspikes);

  end

  highspeedspikes_struct.(sprintf('MI_%s', spikes_date)) = HSS';
  end
  end

  f = highspeedspikes_struct;
