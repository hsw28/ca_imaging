function f = pos_maps_directional(peaks_time, pos, dim)
%plots place cell maps for a bunch of 'cells'



velthreshold = 8;
vel = ca_velocity(pos);
vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
goodvel = find(vel(1,:)>=velthreshold);
goodtime = pos(goodvel, 1);
goodpos = pos(goodvel,:);

figure

numunits = size(peaks_time,1);
ha = tight_subplot(2,ceil((numunits)),[.01 .03],[.1 .01],[.01 .01]);

for k=1:numunits
  highspeedspikes = [];
  for i=1:length(peaks_time(k,:)) %finding if in good vel
    [minValue,closestIndex] = min(abs(peaks_time(k,i)-goodtime));

    if minValue <= 1/7.5 %if spike is within 1 second of moving. no idea if good time
      highspeedspikes(end+1) = peaks_time(k,i);
    end;
  end

  %gets direction
  fwd = [];
  bwd = [];
  for z = 1:length(highspeedspikes)
    [minValue,closestIndex] = min(abs(pos(:,1)-highspeedspikes(z)));
    if pos(max(closestIndex-15, 1),2)-pos(min(closestIndex+15,length(pos)),2)>0
      fwd(end+1) = highspeedspikes(z);
    else
      bwd(end+1) = highspeedspikes(z);
    end
  end

  %subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)

  %plot fwd

  axes(ha(k));
  goodpos(:,3) = 1;
  [rate totspikes totstime colorbarf] = normalizePosData(fwd,goodpos,dim, 2.5);
  colorbar('southoutside');
  set(gca,'ytick',[])
ax = gca;
  axis([0 264/dim, 1 2])
tick_scale_factor = dim;
ax.XTickLabel = ax.XTick * tick_scale_factor;



  q = ceil((numunits))*2;
  axes(ha(q./2+k));
  goodpos(:,3) = 1;
  [rate totspikes totstime colorbarb] = normalizePosData(bwd,goodpos,dim, 2.5);
colorbar('southoutside');
set(gca,'ytick',[])
ax = gca;
  axis([0 264/dim, 1 2])
tick_scale_factor = dim;
ax.XTickLabel = (ax.XTick) * tick_scale_factor;



if max(colorbarf) < .1
        axes(ha(k));
  set(gca, 'clim', [0,.1]);
  c = colorbar('southoutside');
  %c.Label.String = 'spikes per second per bin';
  %xlabel('cm')
end
if max(colorbarb) < .1
    axes(ha(q./2+k));
  set(gca, 'clim', [0,.1]);
colorbar('southoutside');
end

%{
  mn = mean([max(colorbarb), max(colorbarf)])
  axes(ha(q./2+k));
  set(gca, 'clim', [0 mn]);
  colorbar
  axes(ha(k));
  set(gca, 'clim', [0 mn]);
  colorbar

  if max(colorbarf) > max(colorbarb)
    axes(ha(q./2+k));
    set(gca, 'clim', colorbarf);
    colorbar
    axes(ha(k));
    set(gca, 'clim', colorbarf);
    colorbar
  else
    axes(ha(k));
    set(gca, 'clim', colorbarb);
    colorbar
    axes(ha(q./2+k));
    set(gca, 'clim', colorbarb);
    colorbar
  end

%}

end
