function f = mutualinfo_openfield_noCSUS(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, CSUS_id_struct)
%finds mutual info for a bunch of cells
%little did I know i already had code for this: ca_mutualinfo_openfield.m

tic



fields_spikes = fieldnames(spike_structure);
fields_pos = fieldnames(pos_structure);
fields_cats = fieldnames(CA_timestamps);
fields_CSUS = fieldnames(CSUS_id_struct);


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

      fieldName_CSUS = fields_CSUS{i};
      CSUS_id = CSUS_id_struct.(fieldName_CSUS);

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
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        continue
      end

      if (pos(1,1)-pos(end,1))./length(pos) < 1
        pos = convertpostoframe(pos, curr_CA_timestamps);
      end

      mutinfo = NaN(size(peaks_time,1),1);

      tm = pos(:, 1);
      biggest = max(peaks_time(:));
      [minValue,closestIndex] = min(abs(biggest-tm));

      pos = smoothpos(pos);


      goodCSUS = find(CSUS_id(1,:)>0);

      good_CSUStime = pos(goodCSUS,1);
      good_CSUSpos = pos(goodCSUS,:);

      vel = ca_velocity(pos);
      goodvel = find(vel(1,:)>=velthreshold);
      goodtime = pos(goodvel, 1);
      goodpos = pos(goodvel,:);
      goodvel = setdiff(goodvel, goodCSUS);


      mintime = vel(2,1);
      maxtime = vel(2,end);
      tm = vel(2,:);

      numunits = size(peaks_time,1);

      if numunits<=1
        mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = NaN;
        warning('you have no spikes')
      else
      for k=1:numunits
        highspeedspikes = [];

        [c indexmin] = (min(abs(peaks_time(k,:)-mintime))); %
        [c indexmax] = (min(abs(peaks_time(k,:)-maxtime))); %
        currspikes = peaks_time(k,indexmin:indexmax);


        for ii=1:length(currspikes) %finding if in good vel
          [minValue_CSUS,closestIndex] = min(abs(currspikes(ii)-good_CSUStime));
          [minValue_vel,closestIndex] = min(abs(currspikes(ii)-goodtime));
          if minValue_CSUS <= 1/15 & isnan(currspikes(ii))==0 %being CSUS takes precedence
            continue;
          elseif minValue_vel <= 1/15 & isnan(currspikes(ii))==0
            highspeedspikes(end+1) = currspikes(ii);
          end
        end



%want highspeedspikes



  set(0,'DefaultFigureVisible', 'off');
  if length(highspeedspikes)>0
  [rate totspikes totstime colorbar spikeprob occprob] = CA_normalizePosData(highspeedspikes, goodpos, dim, 1.000);
          if (size(spikeprob,1)) < (size(spikeprob,2))
            spikeprob = spikeprob';
          end
          if (size(occprob,1)) < (size(occprob,2))
            occprob = occprob';
          end
  mutinfo(k) = mutualinfo([spikeprob, occprob]);
  else
    mutinfo(k) = NaN;
  end

end

mutualinfo_struct.(sprintf('MI_%s', spikes_date)) = mutinfo';
end
end

f = mutualinfo_struct;
