function f = ca_raster(peak_matrix, length_in_seconds)

frames = size(peak_matrix,2);

spike_val = NaN(size(peak_matrix));
[x, y] = find(peak_matrix>0)
for k=1:max(x)
  temp = find(x==k);
  spike_val(k, 1:length(temp)) = y(temp);
end

raster(spike_val', min(y), max(y));
title('Ca2+ Imaging Raster Plot')
xlabel('Seconds')
ylabel('Cells')
set(gca,'xtick',0:frames/length_in_seconds*60:frames);
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt/frames*length_in_seconds);
