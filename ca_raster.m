function f = ca_raster(peak_matrix)
%plots input times

frames = size(peak_matrix,2);

raster((peak_matrix'), 1, max(peak_matrix(:)));
title('Ca2+ Imaging Raster Plot')
xlabel('Seconds')
ylabel('Cells')

%set(gca,'xtick',0:frames/length_in_seconds*60:frames);
%xt = get(gca, 'XTick');
%set(gca, 'XTick', xt, 'XTickLabel', xt/frames*length_in_seconds);
