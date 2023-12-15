function f = ACL_ca_firingrate(spiketimes, pos)
  %finds firing rate for not movement and for moveement. just change what you want in the code bc too lazy

timelength = pos(end,1)-pos(1,1);

%rate for not movement %%%%%%%%%%%%%%
rate= [];
for k = 1:size(spiketimes,1)
  want = find(isnan(spiketimes(k,:))==0);
  rate(end+1) = length(want)./timelength;
end

length(find(rate<=.01));
f = rate; %ratee for all

%{
%rate for movement %%%%%%%%%%%%%%
ratemov= [];
peaks_time= spiketimes;
velthreshold = 12;
vel = ca_velocity(pos);
vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30.0005); %originally had this at 30, trying with 15 now
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


  ratemov(end+1) = length(highspeedspikes)./(length(goodtime)*30);

end

length(find(ratemov<.01))

f = ratemov;
%}
