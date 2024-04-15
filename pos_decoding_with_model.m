function [values dec_error] = pos_decoding_with_model(model, newClusters, newPos)



  pos = newPos;
  posData = smoothpos(pos);

  spike_clusters = newClusters;

  if isa(spike_clusters,'double')==1
      for k = 1:size(spike_clusters,1)
          name = num2words(k);
          for z = 1:length(name)
              if strfind(name(z),'-')
                  name(z) = ' ';
              end
          end
          name = name(find(~isspace(name)));
          curclus = spike_clusters(k,:);
          curclus = curclus(find(~isnan(curclus)));
          cluststruct.(name) = curclus;
      end
      clusters = cluststruct;
  end  % <-- This is the missing end statement


  %find number of clusters
  clustname = (fieldnames(clusters));
  numclust = length(clustname);



  fxmatrix = model.fxmatrix;
  tdecode = model.tdecode;
  dim = model.dim;
  velthreshold = model.velthreshold;

  ogspikes = length(fieldnames(fxmatrix));


  if ogspikes ~= numclust
    error('did you only use spikes common between training and running the model?')
  end




maxprob = [];
spikenum = 1;
times = [];
percents = [];
numcel = [];
maxx = [];
maxy = [];
same = 0;

%BIN
psize = 1.000 * dim;
xvals = posData(:,2);
yvals = posData(:,3);
xmin = min(posData(:,2));
ymin = min(posData(:,3));
xmax = max(posData(:,2));
ymax = max(posData(:,3));
xbins = ceil((xmax-xmin)/psize); %number of x
ybins = ceil((ymax-ymin)/psize); %number of y
if ybins ==0
  ybins = 1;
end

xinc = xmin +(0:xbins)*psize; %makes a vectors of all the x values at each increment
yinc = ymin +(0:ybins)*psize; %makes a vector of all the y values at each increment

%xmin = 25;
%ymin = 50;
%xmax = 175;
%ymax = 150;

%{
if xmin>min(posData(:,2))
  error('your xmin value is too large')
end
if ymin>min(posData(:,3))
  error('your ymin value is too large')
end
if xmax<max(posData(:,2))
  max(posData(:,2))
    error('your xmax value is too small')
end
if ymax<max(posData(:,3))
    error('your ymax value is too small')
end
%}


occ = zeros(xbins, ybins);
testing = 0;
for x = (1:xbins)
  for y = (1:ybins)
    if x<xbins & y<ybins
      occx = find(xvals>=xinc(x) & xvals<xinc(x+1));
      occy = find(yvals>=yinc(y) & yvals<yinc(y+1));
    elseif x==xbins & y<ybins
      occx = find(xvals>=xinc(x));
      occy = find(yvals>=yinc(y) & yvals<yinc(y+1));
    elseif x<xbins & y==ybins
      occx = find(xvals>=xinc(x) & xvals<xinc(x+1));
      occy = find(yvals>=yinc(y));
    elseif x==xbins & y==ybins
      occx = find(xvals>=xinc(x));
      occy = find(yvals>=yinc(y));
    end
    if length(intersect(occx, occy)) ==0
      occ(x,y) = 0;
    else
      occ(x,y) = 1;
    end
end
end



vel = ca_velocity(posData);

timee = pos(:,1);
[cc indexmin] = min(abs(posData(1,1)-timee));
[cc indexmax] = min(abs(posData(end,1)-timee));
decodetimevector = timee(indexmin:indexmax);
if length(decodetimevector)<10
  timevector = timee;
else
  timevector = decodetimevector;
end


pos_samp_per_sec = length(posData(:,1))./(posData(end,1)-posData(1,1));
%time decoding
tdecodesec = tdecode;
t = round(pos_samp_per_sec*tdecode);


n =0;
nivector = zeros((numclust),1);
tm = 1;

while tm < (length(timevector)-t)

  goodvel = find(vel(2,:)>=timevector(tm) & vel(2,:)<timevector(tm+t));

%  if nanmean((vel(1,goodvel)))>velthreshold & nanmedian((vel(1,goodvel)))>velthreshold
  if length(find(vel(1,goodvel)>velthreshold)) >= length(goodvel)*.75
   %find spikes in each cluster for time
   nivector = zeros((numclust),1);
   for c=1:numclust   %permute through cluster
     name = char(clustname(c));
     %find number of cells that fire during each period
     nivector(c) = length(find(clusters.(name)>=timevector(tm) & clusters.(name)<timevector(tm+t)));
   end

      %for the cluster, permute through the different positions
      endprob = zeros(xbins, ybins);
    endprob = zeros(size(fxmatrix.one));
      %for x = (1:xbins) %WANT TO PERMUTE THROUGH EACH SQUARE OF SPACE SKIPPING NON OCCUPIED SQUARES. SO EACH BIN SHOULD HAVE TWO COORDINATES
        %for y = (1:ybins)
      for x = 1:size(fxmatrix.one,1)
        for y = 1:size(fxmatrix.one,2)
          productme =0;
          expme = 0;
          c = 1;

        %  if occ(x,y) == 0 %means never went there, dont consider
        %    endprob(x,y) = NaN;
        %    break
        %  end

          for c=1:numclust  %permute through cluster
              ni = nivector(c);
              name = char(clustname(c));
              fx = fxmatrix.(name);

              if length(fx)<2
                continue
              end


              if x>size(fx,1)
                x = size(fx,1);
              end
              if y>size(fx,2)
                y = size(fx,2);
              end

              fx = (fx(x, y));
              productme = productme + (ni)*log(fx);  %IN
              expme = (expme) + (fx);
               % goes to next cell, same location

          end
          numcel(end+1) = (ni);
          % now have all cells at that location
          tmm = tdecodesec;

        endprob(x, y) = (productme) + (-tmm.*expme); %IN
        end
        end



        [maxvalx, maxvaly] = find(endprob == max(endprob(:)));

        mp = max(endprob(:))-12;

      endprob = exp(endprob-mp);


         %finds indices
        conv = 1./sum(endprob(~isnan(endprob)), 'all');
        endprob = endprob.*conv; %matrix of percents
        %percents = vertcat(percents, endprob);

        percents(end+1) = max(endprob(:)); %finds confidence
        if length(maxvalx) > 1 %if probs are the sample, randomly pick one and print warning
            same = same+1;
            maxvalx = datasample(maxvalx, 1);
            maxvaly = datasample(maxvaly, 1);

        end

            if length(maxvalx)<1 | length(maxvaly) <1
              maxx(end+1) = NaN;
              maxy(end+1) = NaN;
            else
              maxx(end+1) = (xinc(maxvalx)); %translates to x and y coordinates
              maxy(end+1) = (yinc(maxvaly));
            end


    else
      %means vel is too low
      maxx(end+1) = NaN;
      maxy(end+1) = NaN;
      percents(end+1) = NaN;
    end

        times(end+1) = timevector(tm);

    %if want overlap
  %  if tdecodesec>.5
  %    tm = round(tm+(t/2));
  %  else
      tm = tm+t;
  %  end

    n = n+1;
    if rem(n,500)==0
      n;
    end
end

warning('your probabilities were the same')
same = same
maxx = maxx+(psize/2);
maxy = maxy+(psize/2);
values = [maxx; maxy; percents; times];



values;

dec_error = ca_decodederror(values, posData, tdecode);
error_av = nanmean(dec_error(:,6));
error_med = nanmedian(dec_error(:,6));
