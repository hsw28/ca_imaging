function f = pos_maps_directional(peaks_time, pos, dim)
%plots place cell maps for a bunch of 'cells'



velthreshold = 12;
vel = cavelocity(pos);
vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);

figure

numunits = size(peaks_time,1);
ha = tight_subplot(2,ceil((numunits)),[.01 .03],[.1 .01],[.01 .01]);

for k=1:numunits
  highspeedspikes = [];
  for i=1:length(peaks_time(k,:)) %finding if in good vel
    [minValue,closestIndex] = min(abs(peaks_time(k,i)-goodtime));

    if minValue <= 1 %if spike is within 1 second of moving. no idea if good time
      highspeedspikes(end+1) = peaks_time(k,i);
    end;
  end

  %gets direction
  fwd = [];
  bwd = [];
  for z = 1:length(highspeedspikes)
    [minValue,closestIndex] = min(abs(pos(:,1)-highspeedspikes(z)));
    if pos(closestIndex-15,2)-pos(closestIndex+15,2)>0
      fwd(end+1) = highspeedspikes(z);
    else
      bwd(end+1) = highspeedspikes(z);
    end
  end

  %subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)

  %plot fwd

  axes(ha(k));
  pos(:,3) = 1;
  [rate totspikes totstime colorbarf] = normalizePosData(fwd,pos,dim, 2.5);
  camutualinfo(rate)

  q = ceil((numunits))*2;
  axes(ha(q./2+k));
  pos(:,3) = 1;
  [rate totspikes totstime colorbarb] = normalizePosData(bwd,pos,dim, 2.5);
  camutualinfo(rate)

  if max(colorbarf) > max(colorbarb)
    axes(ha(q./2+k));
    set(gca, 'clim', colorbarf);
    colorbar
    axes(ha(k));
    set(gca, 'clim', colorbarf);
    colorbar
  else
    axes(ha(k));
    set(gca, 'clim', colorbarb);
    colorbar
    axes(ha(q./2+k));
    set(gca, 'clim', colorbarb);
    colorbar
  end


  %normalizePosData(peaks_time(k,:),pos,dim);
end
