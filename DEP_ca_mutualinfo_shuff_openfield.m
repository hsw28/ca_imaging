function f = DEPca_mutualinfo_shuff_openfield(peaks_time, pos, dim, num_times_to_run, ca_MI, velthreshold)
%finds 95% top shuffled and 99% top shuffled mutual info for X number of runs
%put in ca_mutualinfo so you know what to skip

%%%%%%%%%%%%%%NEED TO FIGURE OUT SOME SMOOTHING
tic

mutinfo = NaN(2, size(peaks_time,1));
stddev3 = NaN(2, size(peaks_time,1));


%velthreshold = 6;
vel = ca_velocity(pos);
mintime = vel(2,1);
maxtime = vel(2,end);
%vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);
goodpos = pos(goodvel,:);



numunits = size(peaks_time,1);

for k=1:numunits

  currspikes = peaks_time(k,:);
  if isnan(ca_MI(k,1))==1
    mutinfo(1, k) = NaN;
    mutinfo(2, k) = NaN;

  else
  highspeedspikes = [];
  for i=1:length(currspikes) %finding if in good vel
    [minValue,closestIndex] = min(abs(currspikes(i)-goodtime));

    if minValue <= 1/7.5 %if spike is within 1 second of moving. no idea if good time
      highspeedspikes(end+1) = currspikes(i);
    end;
  end


shuf = NaN(num_times_to_run,1);

parfor l = 1:num_times_to_run

  if isnan(ca_MI(k))==0 && length(highspeedspikes)>1
  shufff = randsample(goodtime, length(highspeedspikes));
  shufff = sort(shufff);

  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(shufff,goodpos,dim, 6.85);

  shuf(l) = mutualinfo([spikeprob', occprob']);
  else
    shuf(l) = NaN;
  end


end


topMI5 = floor(num_times_to_run*.95);
topMI1 = floor(num_times_to_run*.99);
shuf = sort(shuf);
mutinfo(1, k) = shuf(topMI5);
mutinfo(2, k) = shuf(topMI1);


%stddev3(1,k) = nanmean(fwdshuf)+(3*nanstd(fwdshuf));

end
end

toc
f = mutinfo';
%f = stddev3';



%{
(0,'DefaultFigureVisible', 'on');
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
