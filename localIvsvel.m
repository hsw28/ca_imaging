function f = localIvsvel(localI, pos)
  %bbins time on track and plots local Is of place cells against velocity in field location


binsize = 5;

vel = ca_velocity(pos);
vel(1,:) = smooth(vel(1,:), 15);
vel = vel(1,:);
fast = find(vel>12);

pos = pos(:,2)./2.5; %KEEP
pos = pos(1:length(vel));
pos = pos(fast);



numbins = floor((max(pos)-min(pos))./binsize);
binsize = (max(pos)-min(pos))./numbins;
edges = floor(min(pos)):binsize:ceil(max(pos));

[N,EDGES, BIN] = histcounts(pos, edges);
avspeed = [];
for k=1:length(N)
  want = find(BIN==k);
  avspeed(end+1) = nanmean(vel(want));
end


pairs = [];
for j = 1:2
for k=1:size(localI,1)
    field = localI(k,j*2).*4;
        I = localI(k,(j*2)-1);
      if isnan(field)==0 && isnan(I)==0
        dif = (EDGES)-field;
        dif = min(find(dif>0));
    if length(dif) == 0
      dif = length(EDGES);
    end
    both = [avspeed(dif-1), I];
    pairs = [pairs;both];
  else
    both = [NaN, I];
    pairs = [pairs;both];
  end
  end
end

f = pairs;
figure
scatter(pairs(:,1), pairs(:,2))

fitline(pairs(:,1), pairs(:,2));
