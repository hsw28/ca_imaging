function fieldcent  = fieldcenters_openfield(peaks_time, pos, dim, velthreshold, goodcells)
  %good cells IS AN OPTIONAL INPUT and are indices of the cells you know have fields
  % field ceenters are the highest spiking point, not the geometric center
  %rates returns max rate, av rate, min rate


  if nargin < 5
      goodcells = (1:size(peaks_time,1)); % Example: considering all cells as good
  end

  vel = ca_velocity(pos);
  goodvel = find(vel(1,:)>=velthreshold);
  goodtime = pos(goodvel, 1);
  goodpos = pos(goodvel,:);

  mintime = vel(2,1);
  maxtime = vel(2,end);

  numunits = size(peaks_time,1);
  maxrate = NaN(2,numunits);

  for k=1:numunits
    if length(find(k==goodcells))<1
      maxrate(1, k) = NaN;
      maxrate(2, k) = NaN;
    else

    highspeedspikes = [];

    [c indexmin] = (min(abs(peaks_time(k,:)-mintime))); %
    [c indexmax] = (min(abs(peaks_time(k,:)-maxtime))); %
    currspikes = peaks_time(k,indexmin:indexmax);

    for i=1:length(currspikes) %finding if in good vel
      [minValue,closestIndex] = min(abs(currspikes(i)-goodtime));
      if minValue <= 1/7.5 %if spike is within 1 second of moving. no idea if good time
        highspeedspikes(end+1) = currspikes(i);
      end;
    end


    %subplot(ceil(sqrt(numunits)),ceil(sqrt(numunits)), k)
    set(0,'DefaultFigureVisible', 'off');



    fr = ca_firingrate(currspikes, pos);

    if fr > .000001 && length(highspeedspikes)>0

      [rate totspikes totstime colorbar spikeprob occprob] = CA_normalizePosData(highspeedspikes,goodpos,dim, 1.000);
      rate;
      rate = smoothdata(rate, 'gaussian', dim);
      [maxval, maxindex] = max(rate(:));
      [x,y] = ind2sub(size(rate), maxindex);
      maxrate(1, k) = x*dim;
      maxrate(2, k) = y*dim;
    else
      maxrate(1, k) = NaN;
      maxrate(2, k) = NaN;
    end

  end

  end

  fieldcent = maxrate';
