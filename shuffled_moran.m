function percents = shuffled_moran(cellcenter, fieldcenter, num_times_to_shuffle)

%average soma size is 21um https://synapseweb.clm.utexas.edu/dimensions-dendrites
%so radius is aout 11um or 0.011mm
%thats 2 pixels out of 200 if view field is 1.1mm

if isstruct(cellcenter)==1
  cellnames = fieldnames(cellcenter);
  fieldcenternames = fieldnames(fieldcenter);
  daynum = length(cellnames);
else
  daynum =1;
end

neighper = [];
disper = [];
morall = [];

for z=1:(daynum)
  if isstruct(cellcenter)==1
    name = char(cellnames(z));
    currentcellcenter = cellcenter.(name);
    name = char(fieldcenternames(z));
    currentfieldcenter = fieldcenter.(name);
  else
    currentcellcenter = cellcenter;
    currentfieldcenter = fieldcenter;
  end

  currentbicenter = NaN(size(currentcellcenter,2),1);
  for k=1:size(currentcellcenter,2)
    f1f = currentfieldcenter(k, 1); %field center forward
    f1b = currentfieldcenter(k, 2); %field center backward
    if isnan(f1f)==1 && isnan(f1b)==1
      currentbicenter(k) = NaN;
    elseif isnan(f1f)==1
      currentbicenter(k) = f1b;
    elseif isnan(f1b)==1
      currentbicenter(k) = f1f;
    elseif abs(f1f-f1b) <= 15 %same field
      currentbicenter(k) = mean(f1f, f1b);
    elseif abs(f1f-f1b) > 15
      currentbicenter(k) = NaN;
    end
  end


  nans = isnan(currentbicenter);
  currentbicenter = currentbicenter(~nans);

  %this divides track into three segments <-- gets 12/24 with all inputs same
  %tracklength = max(currentbicenter)-min(currentbicenter);
  %[N,EDGES,BIN] = histcounts(currentbicenter, 'BinEdges', [min(currentbicenter), min(currentbicenter)+(tracklength./3), min(currentbicenter)+(2*tracklength./3), max(currentbicenter)]);
  %currentbicenter = BIN;

  currentcellcenter1 = currentcellcenter(1,~nans);
  currentcellcenter1 = (currentcellcenter1);
  currentcellcenter2 = currentcellcenter(2,~nans);
  currentcellcenter2 = (currentcellcenter2);
  currentcellcenter = [currentcellcenter1;currentcellcenter2];



  [W1 W W3] = xy2cont(currentcellcenter(1,:), currentcellcenter(2,:));
  neighweights = normw(W);



  disweights = dist(currentcellcenter);

  %disweights = make_nnw(currentcellcenter(1,:)', currentcellcenter(2,:)', 6, 4);
  %disweights = make_neighborsw(currentcellcenter(1,:)', currentcellcenter(2,:)', 6);
  %disweights = pdweight(currentcellcenter(1,:)', currentcellcenter(2,:)',0,50,1);


  %zzz = find(disweights(:)<3); %centers less than 4 pixels apart
  %disweights(zzz) = 0;

  %qqq = abs(currentbicenter'-currentbicenter);
  %qqq = find(qqq(:)<10); %fields less than 10cm apart
  %zzz = intersect(zzz,qqq); %cells that are close in cell location and field location




  notzero = find(disweights(:)>0);
  %disweights(notzero) = disweights(notzero).^3;
  disweights(notzero) = 1./disweights(notzero);
  disweights = normw(disweights);

  subplot(ceil(daynum./5), ceil(daynum./5), z);
  spatiallag(disweights, currentbicenter);



  neighmor = (moran(currentcellcenter(1:2,:), currentbicenter, neighweights));
  dismore = (moran(currentcellcenter(1:2,:), currentbicenter, disweights));
  morall(end+1) = dismore;


  mi_neigh =[];
  mi_dist = [];
  for k=1:num_times_to_shuffle
      fieldshuff = randperm(length(currentbicenter));
      newfield = currentbicenter(fieldshuff);
      a = (moran(currentcellcenter(1:2,:), newfield, neighweights));
      b = (moran(currentcellcenter(1:2,:), newfield, disweights));

      mi_neigh(end+1) = a;
      mi_dist(end+1) = b;

    end

    mi_neigh =  mi_neigh';
    mi_dist = mi_dist';

    neighper(end+1) = length(find(neighmor>mi_neigh))./num_times_to_shuffle;
    disper(end+1) = length(find(dismore>mi_dist))./num_times_to_shuffle;
end

percents = [neighper; disper; morall]';

length(find(percents(:,1)>=.95))
length(find(percents(:,2)>=.95))
