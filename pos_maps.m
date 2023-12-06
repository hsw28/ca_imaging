function f = pos_maps(peaks_time, pos, dim, velthreshold)
%plots place cell maps for a bunch of 'cells'




vel = ca_velocity(pos);
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);



numunits = size(peaks_time,1);
ha = tight_subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)),[.01 .03],[.1 .01],[.01 .01]);

for k=1:numunits
  highspeedspikes = [];
  for i=1:length(peaks_time(k,:)) %finding if in good vel
    [minValue,closestIndex] = min(abs(peaks_time(k,i)-goodtime));

    if minValue <= 1 %if spike is within 1 second of moving. no idea if good time
      highspeedspikes(end+1) = peaks_time(k,i);
    end;
  end

  %subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)

  axes(ha(k));
  [rate totspikes totstime colorbar spikeprob occprob] = normalizePosData(highspeedspikes,pos,dim, 1.000);

  sigma = 1; % set sigma to the value you need
  sz = 2*ceil(2.6 * sigma) + 1; % See note below
  mask = fspecial('gauss', sz, sigma);
 rate = nanconv(rate, mask, 'same');

 %%%%%%%

 [nr,nc] = size(rate);
 colormap('parula');
 %lower and higher three percent of firing sets bounds
 numrate = rate(~isnan(rate));
 numrate = sort(numrate(:),'descend');
 maxratefive = min(numrate(1:ceil(length(numrate)*0.03)));
 numrate = sort(numrate(:),'ascend');
 minratefive = max(numrate(1:ceil(length(numrate)*0.03)));


 pcolor([rate nan(nr,1); nan(1,nc+1)]);
 shading flat;
 set(gca, 'ydir', 'reverse');
 if minratefive ~= maxratefive
   set(gca, 'clim', [minratefive*2, maxratefive]);
end

colorbar = [minratefive*1.5, maxratefive];

end
