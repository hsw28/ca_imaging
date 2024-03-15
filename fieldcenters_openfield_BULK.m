function field_centers  = fieldcenters_openfield_BULK(peaks_time_struct, pos_struct, dim, velthreshold)
  %good cells IS AN OPTIONAL INPUT and are indices of the cells you know have fields
  % field ceenters are the highest spiking point, not the geometric center
  %rates returns max rate, av rate, min rate



  fields_peaks_time = fieldnames(peaks_time_struct);
  fields_pos = fieldnames(pos_struct);

  if numel(fields_peaks_time) ~= numel(fields_pos)
    error('your spike and pos structures do not have the same number of values. you may need to pad your US structure for exploration days')
  end


  for i = 1:numel(fields_pos)
        fieldName_spikes = fields_peaks_time{i};
        fieldValue_spikes = peaks_time_struct.(fieldName_spikes);
        peaks_time = fieldValue_spikes;


        index = strfind(fieldName_spikes, '_');
        spikes_date = fieldName_spikes(index(2)+1:end)

        fieldName_pos = fields_pos{i};
        fieldValue_pos = pos_struct.(fieldName_pos);
        pos = fieldValue_pos;

        vel = ca_velocity(pos);
        goodvel = find(vel(1,:)>=velthreshold);
        goodtime = pos(goodvel, 1);
        goodpos = pos(goodvel,:);

        mintime = vel(2,1);
        maxtime = vel(2,end);

        numunits = size(peaks_time,1);
        maxrate = NaN(2,numunits);

        for k=1:numunits

          highspeedspikes = [];

          [c indexmin] = (min(abs(peaks_time(k,:)-mintime))); %
          [c indexmax] = (min(abs(peaks_time(k,:)-maxtime))); %
          currspikes = peaks_time(k,indexmin:indexmax);

          for i=1:length(currspikes) %finding if in good vel
            [minValue,closestIndex] = min(abs(currspikes(i)-goodtime));
            if minValue <= 1 %if spike is within 1 second of moving. no idea if good time
              highspeedspikes(end+1) = currspikes(i);
            end;
          end


          %subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)
          set(0,'DefaultFigureVisible', 'off');



          fr = ca_firingrate(currspikes, pos);

          if fr > .0000000001 && length(highspeedspikes)>0

            [rate totspikes totstime colorbar spikeprob occprob] = CA_normalizePosData(highspeedspikes,goodpos,dim, 1.000);
            rate;
            [maxval, maxindex] = max(rate(:));
            [x,y] = ind2sub(size(rate), maxindex);
            maxrate(1, k) = x*dim;
            maxrate(2, k) = y*dim;
          else
            maxrate(1, k) = NaN;
            maxrate(2, k) = NaN;
          end

        end

        field_centers.(sprintf('centers_%s', spikes_date)) = maxrate';
  end
