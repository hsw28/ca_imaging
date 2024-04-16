function model = pos_decoding_model(pos, spike_clusters, tdecode, dim, velthreshold)
%decodes position and outputs decoded x, y, confidence(in percents), and time.
%dim is bin size in cm
%tdecode is decoding in seconds


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
end


posData = pos;
%posData = fixpos(posData);


timee = pos(:,1);


[cc indexmin] = min(abs(posData(1,1)-timee));
[cc indexmax] = min(abs(posData(end,1)-timee));
decodetimevector = timee(indexmin:indexmax);
if length(decodetimevector)<10
  timevector = timee;
else
  timevector = decodetimevector;
end


posData = smoothpos(posData);

tdecodesec = tdecode;


%find number of clusters
clustname = (fieldnames(clusters));
numclust = length(clustname);

%BIN
psize = 1.000 * dim;

xvals = posData(:,2);
yvals = posData(:,3);
xmin = min(posData(:,2));
ymin = min(posData(:,3));
xmax = max(posData(:,2));
ymax = max(posData(:,3));

%xmin = 25;
%ymin = 50;
%xmax = 175;
%ymax = 150;

if xmin>min(posData(:,2))
  error('your xmin value is too large')
end
if ymin>min(posData(:,3))
  error('your ymin value is too large')
end
if xmax<max(posData(:,2))
    error('your xmax value is too small')
end
if ymax<max(posData(:,3))
    error('your ymax value is too small')
end




xbins = ceil((xmax-xmin)/psize); %number of x
ybins = ceil((ymax-ymin)/psize); %number of y
if ybins ==0
  ybins = 1;
end


xinc = xmin +(0:xbins)*psize; %makes a vectors of all the x values at each increment
yinc = ymin +(0:ybins)*psize; %makes a vector of all the y values at each increment



% for each cluster,find the firing rate at esch location
fxmatrix = ca_firingPerPos(posData, clusters, dim, tdecodesec, 7.5, velthreshold);
names = (fieldnames(fxmatrix));
for k=1:length(names)
    curname = char(names(k));
    now = fxmatrix.(curname);
    if size(fxmatrix.(curname),2)>1 & max(now)>0

      fxmatrix.(curname) = ndnanfilter(fxmatrix.(curname), 'gausswin', [dim*2/dim, dim*2/dim], 2, {}, {'replicate'}, 1);

    end

  current = fxmatrix.(curname);
  current(isnan(current)) = eps;
  fxmatrix.(curname) = current;


end




model = struct('fxmatrix', fxmatrix, 'tdecode', tdecode, 'dim', dim, 'velthreshold', velthreshold)
