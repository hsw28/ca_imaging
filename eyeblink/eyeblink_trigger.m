function f = eyeblink_trigger(trigger, event)
  %put in trigger times and matrix of event times, and gives you number of spikes in each PSTH window for each cell

event =  event';

summer = NaN(1,size(event, 2));
for k=1:size(event, 2)
  good = ~isnan(event(:,k));
  good = event(good,k);
  summer(k) = sum(psth(trigger, good));
end

f = summer';
