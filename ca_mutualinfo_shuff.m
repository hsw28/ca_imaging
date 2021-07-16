function f = ca_mutualinfo_shuff(peaks_time, pos, dim, num_times_to_run, ca_MI)
%finds 95% top shuffled mutual info for X number of runs
%put in ca_mutualinfo so you know what to skip
 tic

mutinfo = NaN(2, size(peaks_time,1));

velthreshold = 12;
vel = ca_velocity(pos);
mintime = vel(2,1);
maxtime = vel(2,end);
vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);
goodpos = pos(goodvel,:);

figure

numunits = size(peaks_time,1);

for k=1:numunits
  [c indexmin] = (min(abs(peaks_time(k,:)-mintime))); %how close the REM time is to velocity-- index is for REM time
  [c indexmax] = (min(abs(peaks_time(k,:)-maxtime))); %how close the REM time is to velocity
  currspikes = peaks_time(k,indexmin:indexmax);

  if isnan(ca_MI(k,1))==1 && isnan(ca_MI(k,2))==1
    mutinfo(1, k) = NaN;
    mutinfo(2, k) = NaN;

  else
  highspeedspikes = [];
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

  %subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)


  set(0,'DefaultFigureVisible', 'off');

fwdshuf = [];
bwdshuf = [];

for l = 1:num_times_to_run

  if isnan(ca_MI(k,1))==0
  goodpos(:,3) = 1;
  shufffwd = randsample(goodtime, length(fwd));
  shufffwd = sort(shufffwd);
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(shufffwd,goodpos,dim, 2.5);

  fwdshuf(end+1) = mutualinfo([spikeprob', occprob']);
  else
    fwdshuf(end+1) = NaN;
  end


  if isnan(ca_MI(k,2))==0
  goodpos(:,3) = 1;
  shuffbwd = randsample(goodtime, length(bwd));
  shuffbwd = sort(shuffbwd);
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(shuffbwd,goodpos,dim, 2.5);
  bwdshuf(end+1) = mutualinfo([spikeprob', occprob']);
  else
    bwdshuf(end+1) = NaN;
  end
end

topMI = floor(num_times_to_run*.95);
fwdshuf = sort(fwdshuf);
bwdshuf = sort(bwdshuf);
mutinfo(1, k) = fwdshuf(topMI);
mutinfo(2, k) = bwdshuf(topMI);

end
end

toc
f = mutinfo';
