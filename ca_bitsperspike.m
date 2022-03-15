function f = ca_bitsperspike(posData, clusters, dim)
  %IN PREP
  %need to change velocity to own function with dimensions correct
  %right now firing/per takes a structure

%outputs name, bits per spike, and mean firing rate normalized by position

%determine how many spikes & pos files

spikenum = size(clusters,2);

output = {'cluster name'; 'bits/spike'; 'mean rate'};




  %Sum of (occprobs * mean firing rate per bin / overall mean rate) * log2 (mean firing rate per bin / overall mean rate)
  psize = 2.5 * dim;
  xvals = posData(:,2);
  yvals = posData(:,3);
  xmin = min(posData(:,2));
  ymin = min(posData(:,3));
  xmax = max(posData(:,2));
  ymax = max(posData(:,3));


  xbins = ceil((xmax-xmin)/psize); %number of x
  ybins = ceil((ymax-ymin)/psize); %number of y


  xinc = xmin +(0:xbins)*psize; %makes a vectors of all the x values at each increment
  yinc = ymin +(0:ybins)*psize; %makes a vector of all the y values at each increment

  velthreshold = 6;
  vel = velocity(posData);
  size(vel)
  vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
  fastvel = find(vel(1,:) > velthreshold);
  totaltime = length(fastvel)./7.5;
  posDataFast = posData(fastvel, :);
  xvalsFast = posDataFast(:,2);
  yvalsFast = posDataFast(:,3);

  %occupancy
  occ = zeros(xbins, ybins);
  testing = 0;
  for x = (1:xbins)
    for y = (1:ybins)
      if x<xbins & y<ybins
      occx = find(xvalsFast>=xinc(x) & xvalsFast<xinc(x+1));
      occy = find(yvalsFast>=yinc(y) & yvalsFast<yinc(y+1));
      elseif x==xbins & y<ybins
      occx = find(xvalsFast>=xinc(x));
      occy = find(yvalsFast>=yinc(y) & yvalsFast<yinc(y+1));
      elseif x<xbins & y==ybins
      occx = find(xvalsFast>=xinc(x) & xvalsFast<xinc(x+1));
      occy = find(yvalsFast>=yinc(y));
      elseif x==xbins & y==ybins
      occx = find(xvalsFast>=xinc(x));
      occy = find(yvalsFast>=yinc(y));
      end

      if length(intersect(occx, occy)) == 0
      occ(x,y) = NaN;
      else
      occ(x,y) = length(intersect(occx, occy));
      end
  end
  end

numocc = occ(~isnan(occ));
occtotal = sum(((numocc)), 'all');
occprobs = occ./(occtotal)
occprobs = chartinterp(occprobs);



%if length(~isnan(chart(:))==1)>12000

%Sum of (occprobs * mean firing rate per bin / meanrate) * log2 (mean firing rate per bin / meanrate)

%spike rates
cnames = {};

chart = normalizePosData([1], posDataFast, 5);
length(~isnan(chart(:))==1)





for c = 1:(spikenum)
  clust = clusters(c,:);
    [clustmin indexmin] = min(abs(posData(1,1)-clust));
    [clustmax indexmax] = min(abs(posData(end,1)-clust));


    fxmatrix = firingPerPos(posData, clust, dim, 1, 30, occ);


    assvel = assignvelOLD(clust, vel);
    fastspikeindex = find(assvel > velthreshold);
    %meanrate = length(fastspikeindex)./(totaltime); %WANT ONLY AT HIGH VEL


    fxclust = fxmatrix;
    length(~isnan(fxclust(:))==1);
    fxclust = chartinterp(fxclust);
    meanrate = nanmean(fxclust(:));
    fxclust = ndnanfilter(fxclust, 'gausswin', [10/dim, 10/dim], 2, {}, {'symmetric'}, 1);

    neg = find(fxclust(:)<0);
    fxclust(neg) = eps;


    oldbits = 0;
    newbits = 0;
    bitsper = 0;
    bigX = 0;
    bigY = 0;
    for x = (1:xbins) %WANT TO PERMUTE THROUGH EACH SQUARE OF SPACE SKIPPING NON OCCUPIED SQUARES. SO EACH BIN SHOULD HAVE TWO COORDINATES
      for y = (1:ybins)
        if occprobs(x,y)>0 & ~isnan(fxclust(x,y))==1 & ~isnan(occprobs(x,y))==1


        newbits = (occprobs(x,y) .* (fxclust(x,y) ./ meanrate) * log2((fxclust(x,y) ./ meanrate)));
      %  if newbits>5
      %    x
      %    y
      %  end

        bitsper = bitsper + newbits; %if you want per location, assign this to a matrix


        if newbits > oldbits
          oldbits = newbits;

          xbins = ceil((xmax-xmin)/psize);
          ybins = ceil((ymax-ymin)/psize);

          bigX = (xmax-xmin)./xbins * x + xmin;
          bigY = (ymax-ymin)./ybins * y + ymin;
        end



        end
      end
    end


    bigXall = [bigXall, bigX];
    bigYall = [bigYall, bigY];
    if meanrate <.05
      bitsper = NaN;
    end

    newdata = {name; bitsper; meanrate};

    output = horzcat(output, newdata);


end




f = output';
