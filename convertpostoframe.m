function newpos = convertpostoframe(pos, CA_timestamps)
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
j=2;
k=1;
while j <= length(timestamps) && k<=length(allframes)


  currconv = timestamps(j);
  [c index] = min(abs(pos_time-currconv));


  allframes(k, :) =  pos(index, :);
    j=j+2;
    k = k+1;
end

newpos = allframes;

%{
xpos = newpos(:,2);
ypos = newpos(:,3);
if for_rec_1_for_oval_2 == 2 %oval
  xpos = xpos*.14;
  xpos = xpos-min(xpos);
  ypos = ypos*.15;
  ypos = ypos-min(ypos);
else
  %sq/rectangle. this is to set the bounds but the tracking seems good without so?
    xpos = xpos*.19;
    xpos = xpos-min(xpos);
    ypos = ypos*.15;
    ypos = ypos-min(ypos);
end
newpos(:,2) = xpos;
newpos(:,3) = ypos;
end
%}
