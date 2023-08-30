function allframes = convertpostoframe(pos, CA_timestamps)
%converts pos to same frames/times as CA imaging


timestamps = CA_timestamps;
if isa(timestamps,'table')
  timestamps = table2array(timestamps);
  timestamps = timestamps(:,2);
end

if size(timestamps,2)==3
  timestamps = timestamps(:,2);
end

if timestamps(5)>2
timestamps = timestamps./1000;
end


pos_time = pos(:,1);
allframes = NaN(floor(length(timestamps)./2), 3);
for k=1:length(allframes)

  currconv = timestamps(k);
  [c index] = min(abs(pos_time-currconv));

  allframes(k, :) =  pos(index, :);

end
