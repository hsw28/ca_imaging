function f = ca_mutualinfo(peaks_time, pos, dim, varargin)
%finds mutual info for a bunch of cells

tic
mutinfo = NaN(2, size(peaks_time,1));

if size(varargin)==0
psize = 1.000; %some REAL ratio of pixels to cm -- 3.5 for wilson, 2.5 for disterhoft linear, 1.000 for eyeblink
else
psize = varargin;
end



velthreshold = 2;
vel = ca_velocity(pos)
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

    if minValue <= 1/15 %if spike is within 1 second of moving. no idea if good time
      highspeedspikes(end+1) = currspikes(i);
    end;
  end

  %gets direction
  fwd = [];
  bwd = [];
  for z = 1:length(highspeedspikes)
    [minValue,closestIndex] = min(abs(pos(:,1)-highspeedspikes(z)));
    if pos(max(closestIndex-7, 1),2)-pos(min(closestIndex+7,length(pos)),2)>0
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
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(fwd,goodpos,dim, psize);

  mutinfo(1, k) = mutualinfo([spikeprob', occprob']);
  else
    mutinfo(1, k) = NaN;
  end

  if fr > .01 && length(bwd)>0
  goodpos(:,3) = 1;
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(bwd,goodpos,dim, psize);
  mutinfo(2, k) = mutualinfo([spikeprob', occprob']);
  else
    mutinfo(2, k) = NaN;
  end


end

f = mutinfo';
toc
