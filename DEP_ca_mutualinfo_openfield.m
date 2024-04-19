function f = DEPca_mutualinfo_openfield(peaks_time, pos, dim, velthreshold, varargin)
%finds mutual info for a bunch of cells

if size(varargin)==0
psize = 6.85; %some REAL ratio of pixels to cm -- 3.5 for wilson, 2.5 for disterhoft linear, 6.85 for eyeblink
else
psize = varargin;
end


%%%%%%%%%%%%%%NEED TO FIGURE OUT SOME SMOOTHING
tic
numunits = size(peaks_time,1);
mutinfo = NaN(numunits,1);




vel = ca_velocity(pos);
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);
goodpos = pos(goodvel,:);

mintime = vel(2,1);
maxtime = vel(2,end);



for k=1:numunits
  highspeedspikes = [];

  [c indexmin] = (min(abs(peaks_time(k,:)-mintime))); %
  [c indexmax] = (min(abs(peaks_time(k,:)-maxtime))); %
  currspikes = peaks_time(k,indexmin:indexmax);

  for i=1:length(currspikes) %finding if in good vel
    [minValue,closestIndex] = min(abs(currspikes(i)-goodtime));
    if minValue <= 1/7.5 %if spike is within 1 second of moving. no idea if good time
      highspeedspikes(end+1) = currspikes(i);
    end;
  end


  %subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)
  set(0,'DefaultFigureVisible', 'off');

  fr = ca_firingrate(currspikes, pos);

  if fr > .01 && length(highspeedspikes)>0

    [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(highspeedspikes,goodpos,dim, psize);
    totspikes
    mutinfo(k, 1) = mutualinfo([spikeprob', occprob']);
  else
    mutinfo(k, 1) = NaN;
  end

end

f = mutinfo;
toc
