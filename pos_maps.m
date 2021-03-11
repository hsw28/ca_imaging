function f = pos_maps(peaks_time, pos, dim)
%plots place cell maps for a bunch of 'cells'
%need to eliminate low velocity periods

velthreshold = 8;
vel = velocity(pos);
vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);

figure

numunits = size(peaks_time,1);
for k=1:numunits
  highspeedspikes = [];
  for i=1:length(peaks_time(k,:)) %finding if in good vel
    [minValue,closestIndex] = min(abs(peaks_time(k,i)-goodtime));

    if minValue <= 1 %if spike is within 1 second of moving. no idea if good time
      highspeedspikes(end+1) = peaks_time(k,i);
    end;
  end

  subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)


  normalizePosData(highspeedspikes,pos,dim);

  %normalizePosData(peaks_time(k,:),pos,dim);
end
