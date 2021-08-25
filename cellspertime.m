function f = cellspertime(spiketime)
%compares number of spikes for 5 min of run at the beginning of trial vs 5 min run at end

last = max(spiketime(:));

startnum = [];
endnum = [];
for k=1:size(spiketime,1)
  temp1 = find(spiketime(k,:)>300);
  temp2 = find(spiketime(k,:)<600);
  startnum(end+1) = length(intersect(temp1,temp2));
  if startnum(end)==0
    startnum(end) = NaN;
  end
  temp1 = find(spiketime(k,:)>(last-600));
  temp2 = find(spiketime(k,:)<(last-300));
  endnum(end+1) = length(intersect(temp1, temp2));
  if endnum(end)==0
  endnum(end) = NaN;
end
end


[h,p,ci,stats] = ttest(startnum,endnum, 'Tail','left');


f = [startnum; endnum];
