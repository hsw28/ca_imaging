function f = ca_mutualinfo(peaks_time, pos, dim)
%finds mutual info for a bunch of cells

tic
mutinfo = NaN(2, size(peaks_time,1));



velthreshold = 12;
vel = ca_velocity(pos)
%vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30.0005); %originally had this at 30, trying with 15 now
goodvel = find(vel(1,:)>=velthreshold)
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

  %gets direction
  fwd = [];
  bwd = [];
  for z = 1:length(highspeedspikes)
    [minValue,closestIndex] = min(abs(pos(:,1)-highspeedspikes(z)));
    if pos(max(closestIndex-15, 1),2)-pos(min(closestIndex+15,length(pos)),2)>0
      fwd(end+1) = highspeedspikes(z);
    else
      bwd(end+1) = highspeedspikes(z);
    end
  end

length(fwd);
length(bwd);
length(highspeedspikes);


  %subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)


  set(0,'DefaultFigureVisible', 'off');

 fr = ca_firingrate(currspikes, pos);

  if fr > .01 && length(fwd)>0
  goodpos(:,3) = 1;
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(fwd,goodpos,dim, 2.5);

  mutinfo(1, k) = mutualinfo([spikeprob', occprob']);
  else
    mutinfo(1, k) = NaN;
  end

  if fr > .01 && length(bwd)>0
  goodpos(:,3) = 1;
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(bwd,goodpos,dim, 2.5);
  mutinfo(2, k) = mutualinfo([spikeprob', occprob']);
  else
    mutinfo(2, k) = NaN;
  end


end

f = mutinfo';
toc
