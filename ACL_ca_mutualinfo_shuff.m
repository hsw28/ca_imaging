function f = ACL_ca_mutualinfo_shuff(peaks_time, pos, dim, num_times_to_run, ca_MI, varargin)
%finds 95% top shuffled and 99% top shuffled mutual info for X number of runs
%put in ca_mutualinfo so you know what to skip
 tic

 if size(varargin)==0
 psize = 1.000; %some REAL ratio of pixels to cm -- 3.5 for wilson, 2.5 for disterhoft linear, 1.000 for eyeblink
 else
 psize = varargin;
 end

 if size(pos,2) > size(pos,1)
   pos = pos';
 end

 if size(peaks_time,2) > size(peaks_time,1)
   peaks_time = peaks_time';
 end


mutinfo = NaN(4, size(peaks_time,2));
stddev3 = NaN(2, size(peaks_time,2));

velthreshold = 8;
vel = ca_velocity(pos);
mintime = vel(2,1);
maxtime = vel(2,end);
%vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);
goodpos = pos(goodvel,:);

figure

numunits = size(peaks_time,2);

for k=1:numunits
  want = find(peaks_time(:,k)>0);
  currspikes = pos(want,1);

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

fwdshuf = NaN(num_times_to_run,1);
bwdshuf = NaN(num_times_to_run,1);
goodpos(:,3) = 1;

parfor l = 1:num_times_to_run

  if isnan(ca_MI(k,1))==0
  shufffwd = randsample(goodtime, length(fwd));
  shufffwd = sort(shufffwd);
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(shufffwd,goodpos,dim, 1);

  fwdshuf(l) = mutualinfo([spikeprob', occprob']);
  else
    fwdshuf(l) = NaN;
  end

  if isnan(ca_MI(k,2))==0
  shuffbwd = randsample(goodtime, length(bwd));
  shuffbwd = sort(shuffbwd);
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(shuffbwd,goodpos,dim, 1);
  bwdshuf(l) = mutualinfo([spikeprob', occprob']);
  else
    bwdshuf(l) = NaN;
  end
end


topMI5 = floor(num_times_to_run*.95);
topMI1 = floor(num_times_to_run*.99);
fwdshuf = sort(fwdshuf);
bwdshuf = sort(bwdshuf);

mutinfo(1, k) = fwdshuf(topMI5);

[c index] = (min(abs(ca_MI(k)-fwdshuf)));
if isnan(index)==0
  rank = index./num_times_to_run;
  mutinfo(2, k) = rank;
else
  mutinfo(2,k) = NaN;
end


mutinfo(3, k) = bwdshuf(topMI5);

[c index] = (min(abs(ca_MI(k)-bwdshuf)));
if isnan(index)==0
  rank = index./num_times_to_run;
  mutinfo(4, k) = rank;
else
  mutinfo(4,k) = NaN;
end



%stddev3(1,k) = nanmean(fwdshuf)+(3*nanstd(fwdshuf));
%stddev3(2,k) = nanmean(bwdshuf)+(3*nanstd(bwdshuf));

end
end

toc
f = mutinfo';
f% = stddev3';



%{
set(0,'DefaultFigureVisible', 'on');
figure
subplot(2,1,1)
histogram(fwdshuf, 'BinWidth', .01, 'Normalization','probability')
vline(mutinfo(1))
vline(ca_MI(1), 'g')
%xlabel('Mutual Information')
%ylabel('Occurance (%)')
subplot(2,1,2)
histogram(bwdshuf, 'BinWidth', .01, 'Normalization','probability')
vline(mutinfo(2))
vline(ca_MI(2), 'g')
%}
