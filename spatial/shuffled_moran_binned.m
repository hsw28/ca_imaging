function percents = shuffled_moran_binned(cellcenter, fieldcenter, num_times_to_shuffle)
%does moran but bins track into reward area and not reward area
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
meandist = [];
cellnum = [];

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



  %this divides track into 5 segments <-- gets 14/24 with all inputs same
  tracklength = max(currentbicenter)-min(currentbicenter);
  [N,EDGES,BIN] = histcounts(currentbicenter, 'BinEdges', [min(currentbicenter), min(currentbicenter)+(tracklength./5), min(currentbicenter)+(2*tracklength./5), min(currentbicenter)+(3*tracklength./5), min(currentbicenter)+(4*tracklength./5), max(currentbicenter)]);
  %BIN(find(BIN==5)) = 1; %reward
  BIN(find(BIN==2)) = 3; %not reward
  BIN(find(BIN==4)) = 3;
  currentbicenter = BIN;


  currentcellcenter1 = currentcellcenter(1,~nans);
  currentcellcenter1 = (currentcellcenter1);
  currentcellcenter2 = currentcellcenter(2,~nans);
  currentcellcenter2 = (currentcellcenter2);
  currentcellcenter = [currentcellcenter1;currentcellcenter2];

  %for plotting moran binned
  figure
  sizes = 100;
    colormatrix = NaN(length(BIN), 3);
  for c=1:length(BIN)
    if BIN(c) ==1
      colormatrix(c,:)=[0 0 1]; %reward is blue
    elseif BIN(c) ==2
      colormatrix(c,:)=[1 0 0]; %not reward os red
    end
  end
  %subplot(ceil(daynum./6), ceil(daynum./4), z);
  scatter(jitter(currentcellcenter(1,:)), currentcellcenter(2,:), sizes, colormatrix, 'filled')




  [W1 W W3] = xy2cont(currentcellcenter(1,:), currentcellcenter(2,:));
  neighweights = normw(W);



  disweights = dist(currentcellcenter);
  cellnum(end+1) = length(disweights);

  notzero = find(disweights(:)>0);
  disweights(notzero) = 1./disweights(notzero);
  disweights = normw(disweights);

  %PLOTTING
  %figure
  %subplot(ceil(daynum./6), ceil(daynum./4), z);
  plotval = spatiallag(disweights, currentbicenter);

  figure
  %plotval(lh) = 1; %light grey [0.85 0.85 0.85]
  %plotval(hh) = 2; % red [1 0 0]
%  plotval(ll) = 3; %blue [0 0 1]
%  plotval(hl) = 4; %dark grey [0.4 0.4 0.4]
  lh = [0.85 0.85 0.85];
  hh= [1 0 0];
  ll = [0 0 1];
  hl = [0.4 0.4 0.4];
  colormatrix = NaN(length(plotval), 3);
  for c=1:length(plotval)
    if plotval(c) ==1
      colormatrix(c,:)=lh;
    elseif plotval(c) ==2
      colormatrix(c,:)=hh;
    elseif plotval(c) ==3
      colormatrix(c,:)=ll;
    elseif plotval(c) ==4
      colormatrix(c,:)=hl;
    end
  end
  sizes = 100;
  scatter(jitter(currentcellcenter(1,:)), currentcellcenter(2,:), sizes, colormatrix, 'filled')
  %title('Both Directions')
  axis([0 200,0 200])




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
    meandist(end+1) = nanmean(mi_dist);

    neighper(end+1) = length(find(neighmor>mi_neigh))./num_times_to_shuffle;
    disper(end+1) = length(find(dismore>mi_dist))./num_times_to_shuffle;
end

percents = [neighper; disper; morall; meandist; cellnum]';

length(find(percents(:,1)>=.95))
length(find(percents(:,2)>=.95))
