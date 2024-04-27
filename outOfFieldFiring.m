function f = outOfFieldFiring(place_MI, spikes, pos, field_center, CSUS_id, curr_CA_timestamps)
%ex: outOfFieldFiring(rat0313.MI_shuff.MI_2023_05_31, rat0313.Ca_peaks.CA_peaks_2023_05_31, rat0313.pos.pos_2023_05_31, rat0313.field_centers.centers_2023_05_31, rat0313.CSUS_id.CSUS_id_2023_05_31, rat0313.Ca_ts.CA_time_2023_05_31)

%looking at mutual info to see if it is a place cell
place_MI = place_MI(:,3);
pc = find(place_MI>=.95);

%only looking at spikes of cells that are place cells
spikes = spikes(pc,:);
field_center = field_center(pc,:);

%only looking at spikes that occur when the animal is >4cm/s OR during a CSUS period

%first getting pos in the right format and finding velocities
if (pos(1,1)-pos(end,1))./length(pos) < 1
  pos = convertpostoframe(pos, curr_CA_timestamps);
end


pos = smoothpos(pos);
vel = ca_velocity(pos);


goodvel = find(vel(1,:)>=4);
goodCSUS = find(CSUS_id(1,:)>0);
good_veltime = pos(goodvel, 1);
good_CSUStime = pos(goodCSUS,1);

good_velpos = pos(goodvel, :);
good_CSUSpos = pos(goodCSUS,:);


%for each cell, go through the spikes and find if in good velocity or CSUS
highspeedspikes = [];
CSUSspikes = [];
dist = NaN(size(spikes,1),4);
for i=1:size(spikes,1) %picking cells
  currspikes = spikes(i,:);
  CSUSspikes = [];
  highspeedspikes = [];
    for ii=1:length(currspikes) %finding if in good vel
      [minValue_CSUS,closestIndex] = min(abs(currspikes(ii)-good_CSUStime));
      [minValue_vel,closestIndex] = min(abs(currspikes(ii)-good_veltime));

      if minValue_CSUS <= 1/15 && isnan(currspikes(ii))==0 %being CSUS takes precedence
        CSUSspikes(end+1) = currspikes(ii);
      elseif minValue_vel <= 1/15 && isnan(currspikes(ii))==0
        highspeedspikes(end+1) = currspikes(ii);
      end
    end %now we have if the spike was during CSUS period or during running


    %time to see where the spikes for this cell are

    cellcenter = field_center(i,:);
    d_CSUS = NaN;
    d_fast = NaN;
    p_CSUS = [NaN,NaN];
    p_fast = [NaN,NaN];
    dim = 2.5;
    if length(CSUSspikes>1)
      [rate totspikes totstime colorbar spikeprob occprob] = CA_normalizePosData(CSUSspikes,good_CSUSpos,2.5, 1.000);
      [maxval, maxindex] = max(rate(:));
      [x,y] = ind2sub(size(rate), maxindex);
      maxrate_csus = [x*dim,y*dim];

      %  CSUSspikes_place = ca_placeevent(CSUSspikes, pos);
      %  CSUSspikes_place = CSUSspikes_place(2:3,:);
      %  for z = 1:length(CSUSspikes)
      %    curr_center = CSUSspikes_place(:,z)';
      %    points = [cellcenter; curr_center];
      %    p_CSUS(end+1,:) = curr_center;
      %    d_CSUS(end+1) = pdist(points, 'euclidean');
      %  end
    else
        d_CSUS(end+1) = NaN;
    end


    if length(highspeedspikes>1)
      [rate totspikes totstime colorbar spikeprob occprob] = CA_normalizePosData(highspeedspikes,good_velpos,2.5, 1.000);
      [maxval, maxindex] = max(rate(:));
      [x,y] = ind2sub(size(rate), maxindex);
      maxrate_fast = [x*dim,y*dim];

    %  fast_spikes_place = ca_placeevent(highspeedspikes, pos);
    %  fast_spikes_place = fast_spikes_place(2:3,:);
    %  for z = 1:length(highspeedspikes)
    %    curr_center = fast_spikes_place(:,z)';
    %    points = [cellcenter; curr_center];
    %    p_fast(end+1,:) = curr_center;
    %    d_fast(end+1) = pdist(points, 'euclidean');
    %  end
    else
        d_fast(end+1) = NaN;
    end



    dist(i,1:2) = maxrate_csus;
    dist(i,3:4) = maxrate_fast;




  end
   f = dist;


%{
    %below is for finding firing locations, not sure if I will use
    if length(CSUSspikes>1)
      CSUS_rate = ca_firingPerPos(pos(good_CSUStime, :), CSUSspikes, 2.5, 0);
    else
      CSUS_rate = NaN;
    end
    if length(highspeedspikes>1)
      fast_rate = ca_firingPerPos(pos(goodvel, :), highspeedspikes, 2.5, 0);
    else
      fast_rate = NaN;
    end
%}
