function f = ca_mutualinfo(peaks_time, pos, dim)
%plots place cell maps for a bunch of 'cells'


mutinfo = NaN(2, size(peaks_time,1));

velthreshold = 8;
vel = ca_velocity(pos);
vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);
goodpos = pos(goodvel,:);

numunits = size(peaks_time,1);

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


  set(0,'DefaultFigureVisible', 'off');

  if length(fwd) >5
  goodpos(:,3) = 1;
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(fwd,goodpos,dim, 2.5);

  mutinfo(1, k) = mutualinfo([spikeprob', occprob']);
  else
    mutinfo(1, k) = NaN;
  end


  if length(bwd) >5
  goodpos(:,3) = 1;
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(bwd,goodpos,dim, 2.5);
  mutinfo(2, k) = mutualinfo([spikeprob', occprob']);
  else
    mutinfo(2, k) = NaN;
  end


end

f = mutinfo';
