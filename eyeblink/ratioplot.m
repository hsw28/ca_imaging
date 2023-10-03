function f = ratioplot(ratio)

%{
figure
hold on

for i = 1:size(ratio,1)
  want = ~isnan(ratio(i,:));
  if length(find(want>0))>1
%    want(16) = 1;
    len = [1:1:size(ratio,2)];
    xplot = len(want);
    yplot = ratio(i,xplot);
    plot(xplot,yplot,'o-')
  end
end

plot(nanmean(ratio), 'LineWidth', 3, 'Color', 'k');

xlabel('Training day')
ylabel('US Response (total calcium events)')


figure
hold on
for i = 1:size(ratio,1)
  av3 = [];
  j = 1;
  while j <= (size(ratio,2))
    if j+2<=(size(ratio,2))
    av3(end+1) = nanmean(ratio(i,j:j+2));
    j = j+3;
    else
    av3(end+1) = nanmean(ratio(i,j:end));
    j = j+3;
    end
  end
  want = ~isnan(av3);
  if length(find(want>0))>1
    len = [1:1:length(av3)];
    xplot = len(want);
    yplot = av3(xplot);
    plot(xplot,yplot,'o-')
  end
end
%}


figure
hold on
allav = [];
for i = 1:size(ratio,1)
  av3 = [];
  j = 1;
  while j <= (size(ratio,2))
    if j+2<=(size(ratio,2))
    av3(end+1) = nanmean(ratio(i,j:j+2));
    j = j+3;
    else
    av3(end+1) = nanmean(ratio(i,j:end));
    j = j+3;
    end
  end
    xplot = [1:1:length(av3)];
    yplot = av3(xplot);
    xplot = xplot.*3;
    plot(xplot,yplot,'o-')
    allav = [allav; av3];
end

plot(xplot,nanmean(allav), 'LineWidth', 3, 'Color', 'k');

xlabel('Training day')
ylabel('US Response (total calcium events)')
title('Responses averaged across 3 days')
