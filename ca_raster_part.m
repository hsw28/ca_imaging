function f = ca_raster_part(peak_matrix, mintime, maxtime)
%plots input times

[tempx tempy] = find(peak_matrix>=mintime & peak_matrix<maxtime); %x is cell number, y is time

[tempx,I] = sort(tempx);
tempy = tempy(I);


new_matrix = NaN(size(peak_matrix));


k = 1;
z = 1;
val = 1;
uniquex = unique(tempx);
while z<length(tempx)
    if tempx(z) == uniquex(k)
      new_matrix(k, val) = peak_matrix(tempx(z),tempy(z));
      z=z+1;
      val = val+1;
    else
      k = k+1;
      new_matrix(k, val) = peak_matrix(tempx(z),tempy(z));
      z=z+1;
      val = 1;
    end

end

f = new_matrix;



raster((new_matrix'), 1, max(new_matrix(:)));
title('Ca2+ Imaging Raster Plot')
xlabel('Seconds')
ylabel('Cells')
axis([mintime maxtime, 0 length(unique((tempy)))])

%set(gca,'xtick',0:frames/length_in_seconds*60:frames);
%xt = get(gca, 'XTick');
%set(gca, 'XTick', xt, 'XTickLabel', xt/frames*length_in_seconds);
