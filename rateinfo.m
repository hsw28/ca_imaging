function rates  = rateinfo(peaks_time, goodcells, pos, dim)
%rates returns max rate, min rate,  av rate, for fwd and backwrd


  velthreshold = 12;
  vel = ca_velocity(pos);
  vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30.0005); %originally had this at 30, trying with 15 now
  goodvel = find(vel(1,:)>=velthreshold);
    goodvel = find(vel(1,:)<=5);
  goodtime = pos(goodvel, 1);
  goodpos = pos(goodvel,:);

  mintime = vel(2,1);
  maxtime = vel(2,end);

  numunits = size(peaks_time,1);

  maxrate = [];
  fwdmax = [];
  bkmax = [];
  fwdmin = [];
  bkmin = [];
  fwdav = [];
  bkav = [];
  av = [];
  allmx = [];
  for k=1:numunits
    if length(find(goodcells(:,1)==k))==0
      maxrate(1, k) = NaN;
      maxrate(2, k) = NaN;
    else
    goodcelldir = find(goodcells(:,1)==k);
    goodcelldir = goodcells(goodcelldir,:);

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


    set(0,'DefaultFigureVisible', 'off');


    if length(find(goodcelldir==1))>0 && length(fwd)>0
    goodpos(:,3) = 1;
    [rate1 totspikes totstime colorbar spikeprob occprob] = normalizePosData(fwd,goodpos,dim, 2.5);
    rate1 = smoothdata(rate1, 'gaussian', dim);
    [maxval,maxindex] = max(rate1);
    [minval,minindex] = min(rate1);
    avval = nanmean(rate1);
    fwdmax(end+1) = maxval;
    fwdmin(end+1) = minval;
    fwdav(end+1) = avval;
    else
      rate1 = NaN;
      maxval = NaN;
      fwdmax(end+1) = NaN;
      fwdmin(end+1) = NaN;
      fwdav(end+1) = NaN;
    end


   if length(find(goodcelldir==2))>0 && length(bwd)>0
    goodpos(:,3) = 1;
    [rate2 totspikes totstime colorbar spikeprob occprob] = normalizePosData(bwd,goodpos,dim, 2.5);
    rate2 = smoothdata(rate2, 'gaussian', dim);
    [maxval2,maxindex] = max(rate2);
    [minval,minindex] = min(rate2);
    avval = nanmean(rate2);
    bkmax(end+1) = maxval;
    bkmin(end+1) = minval;
    bkav(end+1) = avval;
    else
      rate2 = NaN;
      maxval2 = NaN;
      bkmax(end+1) = NaN;
      bkmin(end+1) = NaN;
      bkav(end+1) = NaN;
    end

av(end+1) = nanmean([rate1, rate2]);
allmx(end+1) = max(maxval, maxval2);

  end

  end

rates = [fwdmax; fwdmin; fwdav; bkmax; bkmin; bkav]';
rates = [av; allmx]';
