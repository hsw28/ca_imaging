function f = boxplot_grouped(x, y, x_divisions)
%makes boxplot of scatter plot data but groups your x variable in divisions of the size specified

x_edges = min(x):x_divisions:max(x);


[N,EDGES,BIN] = histcounts(x,x_edges);
%BIN is an array of the same size as X whose elements are the bin indices for the corresponding elements in X.


space = floor(length(x_edges)./3);
categ = NaN(1,length(x));
meanz = [];
stdz = [];
firstthird =[];
lastthird = [];
for k = 1:length(x_edges)
  curbin = find(BIN==k);
  categ(curbin) = EDGES(k);
  meanz(end+1) = nanmean(y(curbin));
  stdz(end+1) = std(y(curbin))./sqrt(length(curbin));
  if k<space
      firstthird(end+1:end+length(curbin)) =y(curbin);
  elseif k>=length(x_edges)-space
    lastthird(end+1:end+length(curbin)) = y(curbin);
  end
end


%figure
%boxplot(y, categ)


figure
bar(x_edges, meanz)
hold on
er = errorbar(x_edges, meanz, stdz);
er.Color = [0 0 0];
er.LineStyle = 'none';
xlabel('Binned Place Field Locations (cm)')
ylabel('Average Local Morans I')

f = meanz;


%firstthird
%lastthird
%mean(firstthird)
%mean(lastthird)
[ a b c d] = ttest2(firstthird,lastthird);
