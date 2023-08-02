function f = ca_mutualinfo_openfield(peaks_time, pos, dim, velthreshold)
%finds mutual info for a bunch of cells

%%%%%%%%%%%%%%NEED TO FIGURE OUT SOME SMOOTHING
tic
mutinfo = NaN(1, size(peaks_time,1));




vel = ca_velocity(pos);
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);
goodpos = pos(goodvel,:);

mintime = vel(2,1);
maxtime = vel(2,end);

numunits = size(peaks_time,1);

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

  if fr > .01 && length(highspeedspikes)>0

    [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(highspeedspikes,goodpos,dim, 6.85);
    totspikes
    mutinfo(1, k) = mutualinfo([spikeprob', occprob']);
  else
    mutinfo(1, k) = NaN;
  end

end

f = mutinfo';
toc
