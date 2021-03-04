function f = pos_maps(peaks_time, pos, dim)
%plots place cell maps for a bunch of 'cells'

figure

numunits = size(peaks_time,1);
for k=1:numunits
  subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)
  want = peaks_time(k,:)>200;
  cur = peaks_time(k,want);

  normalizePosData(cur,pos,dim);

  %normalizePosData(peaks_time(k,:),pos,dim);
end
