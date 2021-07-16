function [f allsizescenters shuffledskew] = CAplacefieldnumDIRECTIONAL(clusters,posData, dim, ca_MI_good)

      %set(0,'DefaultFigureVisible', 'off');

%determine how many spikes & pos files

posData(:,3) = 1;
allsizes = [];
allcenterXmean = [];
allcenterYmean = [];
allcenterXmax = [];
allcenterYmax = [];
allskew = [];
alldir = [];
alldirskew = [];

allavpfrate = [];
allmaxpfrate = [];
shuffledskew = NaN(500, 500);
totsnum = 0;
FRAI = [];



%output = {'cluster name'; 'cluster size'; 'direction'; 'num of fields'; 'field size in cm'; 'centermax'; 'centermean'; 'skewness'};
output = {'cluster name'; 'cluster size'; 'direction'; '1=to, 2=away'; 'field size in cm'; 'centermax X'; 'centermax Y';'av field rate'; 'max field rate'};


  velthreshold = 8;
  vel = ca_velocity(posData);
  vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30); %originally had this at 30, trying with 15 now
  fastvel = find(vel(1,:) > velthreshold);
  totaltime = length(fastvel)./30;
  posDataFast = posData(fastvel, :);
  xvalsFast = posDataFast(:,2);
  yvalsFast = posDataFast(:,3);

  psize = 2.5 * dim;
  xvals = posDataFast(:,2);
  yvals = posDataFast(:,3);
  xmin = min(posDataFast(:,2));
  ymin = min(posDataFast(:,3));
  xmax = max(posDataFast(:,2));
  ymax = max(posDataFast(:,3));
  xbins = ceil((xmax)/psize); %number of x
  xinc = (0:xbins)*psize; %makes a vectors of all the x values at each increment

  %occupancy
  occ = zeros(xbins);
  testing = 0;
  for x = (1:xbins)
      if x<xbins
      occx = find(xvalsFast>=xinc(x) & xvalsFast<xinc(x+1));
      elseif x==xbins
      occx = find(xvalsFast>=xinc(x));
      end

      if length((occx)) == 0
      occ(x) = NaN;
      else
      occ(x) = length(occx);
      end
  end


numocc = occ(~isnan(occ));
occtotal = sum(((numocc)), 'all');
occprobs = occ./(occtotal);

%Sum of (occprobs * mean firing rate per bin / meanrate) * log2 (mean firing rate per bin / meanrate)

%spike rates


for z = 1:length(ca_MI_good)

    c = ca_MI_good(z,1);
    clust = clusters(c,:);
    clustsize = length(find(isnan(clust)==0));
    [clustmin indexmin] = min(abs(posData(1,1)-clust));
    [clustmax indexmax] = min(abs(posData(end,1)-clust));
    clust = clust(indexmin:indexmax);

    assvel = assignvelOLD(clust, vel);
    fastspikeindex = find(assvel > velthreshold);
    %meanrate = length(fastspikeindex)./(totaltime); %WANT ONLY AT HIGH VEL

    torewardspikes =[];
    awayrewardspikes=[];
    for z = 1:length(fastspikeindex)
      [minValue,closestIndex] = min(abs(posData(:,1)-fastspikeindex(z)));
      if posData(max(closestIndex-15, 1),2)-posData(min(closestIndex+15,length(posData)),2)>0
        torewardspikes(end+1) = fastspikeindex(z);
      else
        awayrewardspikes(end+1) = fastspikeindex(z);
      end
    end


      %%%%%%%%%%%%

    %NOW DO CHARTS AND EVERYTHING FOR BOTH
    %spiking normalization chart
  for zdir=1:2;
      if zdir ==1;
        spikestochart = torewardspikes;
        dir = 'to';
        currentdir = 1;
      else
        spikestochart = awayrewardspikes;
        dir = 'away';
        currentdir = 2;
      end


    spikestochart = cutclosest(posDataFast(1,1), posDataFast(end,1), spikestochart, spikestochart);
    if length(spikestochart)<2
      continue
    end


    [chart spikegraph timegraph] = normalizePosData(spikestochart, posDataFast, dim, 2.5);
    chart(find(chart==0)) = eps;

    spikegraph = (spikegraph);
    timegraph = (timegraph);


    %smoothing spiking normalization
    chart = smoothdata(chart, 'gaussian', 2.5/dim);
    spikegraph = smoothdata(spikegraph, 'gaussian', 2.5/dim);
    timegraph = smoothdata(timegraph, 'gaussian', 2.5/dim);




    chartlin = sort(chart(:));
    chartlin = chartlin(~isnan(chartlin));
    chartnozero = mean(chartlin(find(chartlin>0)));
    chartlinstd = std(chartlin);
    chartlinmedian = median(chartlin);
    chartlinmean = mean(chartlin);



    %divided into 5cm by 5cm bins
    %place fields are 5 or more adjacent picels with a firing rate >3xmean unit rate
    %finds mean rate


    meanrate = nanmean(chart(:));
    maxrate = max(chart(:));
    chart(isnan(chart)) = 0;
    spikegraph(isnan(spikegraph)) = 0;
    timegraph(isnan(timegraph)) = 0;



    %finds maxes
    actualmax = imregionalmax(chart); %maxes are where there is a 1 on this chart
    [I,J] = find(actualmax==1);
    linearmax = sub2ind(size(actualmax), I, J); %linear indices

    %finds areas where firing > 3x mean. those are marked with a 1 on the chart
      %[I,J] = find(chart>=chartlinmean+(1*(chartlinstd)));
      %[I,J] = find(chart>=3*meanrate);


%We defined the extent of a place field as all connected occupancy bins
%whose smoothed event rate exceeded 1/3 of the peak event rate occupancy bin.
    [I,J] = find(chart>=(max(chart)*.33));


    chartmax = zeros(size(chart));


    for p=1:length(I)
      chartmax(I(p),J(p)) = 1;
    end

    CC = bwconncomp(chartmax,4); %finds pixels connected by sides
    numfields = 0;
    fieldsize = [];
    centermax = [];
    centermean =[];


    for z=1:length(CC.PixelIdxList)


      [Yindex, Xindex] = ind2sub(size(chartmax),cell2mat(CC.PixelIdxList(z)));
      XM = length(unique(Xindex));

        fsize = XM*dim;


        curr = (chart(Yindex, Xindex));
        currspikegraph = (spikegraph(Yindex, Xindex));
        currtimegraph= (timegraph(Yindex, Xindex));


        fieldsize(end+1) = fsize;

        %find all instances where animal goes through place field
        %finds indices of place fields
        [centerY, centerX] = ind2sub(size(chartmax),cell2mat(CC.PixelIdxList(z)));
        %converts indices to actual values
        realX = (xmax)./xbins .* centerX;

        currentrates = chart(centerY, centerX);


        avpfrate = nanmean(currentrates(:));
        maxpfrate = max(currentrates(:));


        numfields = numfields+1;
        [centerY, centerX] = ind2sub(size(chartmax),cell2mat(CC.PixelIdxList(z)));

        xbins = ceil((xmax)/psize);
        ybins = ceil((ymax)/psize);



        newchart = zeros(size(chart));
        newchart(centerY, centerX) = chart(centerY, centerX);
        M = max(newchart, [], 'all');
        [centerYmax, centerXmax] = find(newchart==M);
        centerYmax = ybins-centerYmax; %here?
        newchart(find(newchart==0))=NaN;

        newspikegraph  = zeros(size(chart));
        newtimegraph  = zeros(size(chart));
        newspikegraph(centerY, centerX) = spikegraph(centerY, centerX);
        newtimegraph(centerY, centerX) = timegraph(centerY, centerX);
        newspikegraph(find(newspikegraph==0))=NaN;
        newtimegraph(find(newtimegraph==0))=NaN;






        centerXmean = nanmean(centerX);
        centerYmean = nanmean(centerY);
        centerYmean = ybins-centerYmean;
        centerXmean = (xmax)./xbins * centerXmean;
        centerYmean = (ymax)./ybins * centerYmean;





        newflattened = nanmean(chart,1);
        flatstart = (find(newflattened>0));
        flattenednewspikegraph = nanmean(newspikegraph,1);
        flattenednewtimegraph = nanmean(newtimegraph,1);

        newflattened = newflattened(flatstart(1):end);
        flattenednewspikegraph = flattenednewspikegraph(flatstart(1):end);
        flattenednewtimegraph = flattenednewtimegraph(flatstart(1):end);


        %SHUFFLE
        %SHUFFLE
        totsnum = totsnum+1;
        %for sss=1:500 %SHUFFLE
        %  newflattened = newflattened(~isnan(newflattened)); %SHUFFLE
        %  newflattened = newflattened(randperm(length(newflattened))); %SHUFFLE






        counter = 0;
        flatmean = 0;
        countersum = 0;


        for kk = 1:length(newflattened)
          if newflattened(kk)>0
          flatmean = flatmean+(kk*newflattened(kk));
          counter = counter+1;
          countersum = countersum+newflattened(kk);
          end
        end
        flatmean = flatmean./countersum;


        flatmom = 0;
        temp = [];
        for kk = 1:length(newflattened)
          if newflattened(kk)>0
            kk-flatmean;
            temp(end+1)= ((kk-flatmean)^3)*newflattened(kk);
          flatmom = flatmom+((kk-flatmean)^3)*newflattened(kk);
          end
        end








        centerXmax = (xmax)./xbins * centerXmax;
        centermax(end+1) = centerXmax;

        centerXmean = nanmean(centerX);
        centerXmean = (xmax)./xbins * centerXmean;
        centermean(end+1) = centerXmean;


        allsizes(end+1) = fieldsize(end);
        allcenterXmean(end+1) = centerXmean;
        allcenterYmean(end+1) = centerYmean;
        allcenterXmax(end+1) = centerXmax(1);
        allcenterYmax(end+1) = centerYmax(1);
        alldir(end+1) = currentdir;
        allavpfrate(end+1) = avpfrate;
        allmaxpfrate(end+1) = maxpfrate;

        newdata = {c; clustsize; dir; alldir(end); allsizes(end); allcenterXmax(end); allcenterYmax(end); avpfrate(end); maxpfrate(end)};
        output = horzcat(output, newdata);




    end

    %maximum firing rate >= .5hz

    %newdata = {name; clustsize; dir; numfields; fieldsize; centermean; centermax; skewness};
    %output = horzcat(output, newdata);
  end




end

f = output';
