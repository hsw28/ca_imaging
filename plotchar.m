function f = plotchar(cellcenters, characteristic)
%plots cell centers in a color corresponding to a (numerical) characteristic


figure
colors = (characteristic - min(characteristic))./max(characteristic); %normalize to 0:1
colors = colors*100;
%colors = colors(indexes3,:);
c = colorbar;
%set(gca, 'clim', [min(indexes2),100]);
sizes = 100;

length(cellcenters(3,:))
size(colors)
scatter(cellcenters(3,:), cellcenters(4,:), sizes, colors, 'filled')
