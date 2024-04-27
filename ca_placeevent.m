function f = ca_placeevent(events, pos)

postime = pos(:,1);
event_pos = NaN(3,length(events));
for i=1:length(events)
  [minValue_vel,closestIndex] = min(abs(events(i)-postime));

  event_pos(2:3,i) = pos(closestIndex,2:3);
  event_pos(1,i) = events(i);
end

f = event_pos;
