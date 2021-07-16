function f = ca_mutualinfo_shuff_all(peaks_time, pos, dim, num_times_to_run)
%finds 95% top shuffled mutual info for X number of runs
%shuffles spike trains for mutual info in all times (NOT only in run time)


mutinfo = NaN(2, size(peaks_time,1));

velthreshold = 8;
vel = ca_velocity(pos);
vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);
goodpos = pos(goodvel,:);

figure

numunits = size(peaks_time,1);

for k=1:numunits

  for l = 1:num_times_to_run

  shuff = randsample(pos(:, 1), length(numunits));
  shuff = sort(shuff);

  highspeedspikes = [];
  for i=1:length(shuff) %finding if in good vel
    [minValue,closestIndex] = min(abs(shuff(i)-goodtime));

    if minValue <= 1 %if spike is within 1 second of moving. no idea if good time
      highspeedspikes(end+1) = shuff(i);
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

fwdshuf = NaN(num_times_to_run,1);
bwdshuf = NaN(num_times_to_run,1);



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

topMI = floor(num_times_to_run*.95);
fwdshuf = sort(fwdshuf);
bwdshuf = sort(bwdshuf);
mutinfo(1, k) = fwdshuf(topMI);
mutinfo(2, k) = bwdshuf(topMI);

end

f = mutinfo';
