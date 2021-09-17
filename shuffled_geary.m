function percents = shuffled_geary(cellcenter, fieldcenter, num_times_to_shuffle)

if isstruct(cellcenter)==1
  cellnames = fieldnames(cellcenter);
  fieldcenternames = fieldnames(fieldcenter);
  daynum = length(cellnames);
else
  daynum =1;
end

neighper = [];
disper = [];
gearallneigh = [];
gearalldist = [];

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
  currentcellcenter1 = currentcellcenter(1,~nans);
  currentcellcenter2 = currentcellcenter(2,~nans);
  currentcellcenter = [currentcellcenter1;currentcellcenter2];



  [W1 W W3] = xy2cont(currentcellcenter(1,:), currentcellcenter(2,:));
  neighweights = normw(W);



  %disweights = make_nnw(currentcellcenter(1,:)', currentcellcenter(2,:)', 4, 4);
  %disweights = make_neighborsw(currentcellcenter(1,:)', currentcellcenter(2,:)', 2);
  %disweights = pdweight(currentcellcenter(1,:)', currentcellcenter(2,:)',0,3,1);




  disweights = dist(currentcellcenter);
  notzero = find(disweights(:)>0);
  %disweights(notzero) = 2.^disweights(notzero);
  disweights(notzero) = 1./disweights(notzero);
  disweights = normw(disweights);


  neighmor = (geary(currentcellcenter(1:2,:), currentbicenter, neighweights));
    gearallneigh(end+1) = neighmor;
  dismore = (geary(currentcellcenter(1:2,:), currentbicenter, disweights));
  gearalldist(end+1) = dismore;




  mi_neigh =[];
  mi_dist = [];
  for k=1:num_times_to_shuffle
      fieldshuff = randperm(length(currentbicenter));
      newfield = currentbicenter(fieldshuff);
      a = (geary(currentcellcenter(1:2,:), newfield, neighweights));
      b = (geary(currentcellcenter(1:2,:), newfield, disweights));

      mi_neigh(end+1) = a;
      mi_dist(end+1) = b;

    end

    mi_neigh =  mi_neigh';
    mi_dist = mi_dist';

    neighper(end+1) = length(find(neighmor>mi_neigh))./num_times_to_shuffle;
    disper(end+1) = length(find(dismore>mi_dist))./num_times_to_shuffle;
end


percents = [neighper; disper; gearallneigh; gearalldist]';


length(find(percents(:,1)>=.95))
length(find(percents(:,2)>=.95))
