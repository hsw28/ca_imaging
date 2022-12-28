function percents = shuffled_moran(cellcenter, fieldcenter, num_times_to_shuffle, varargin)

%average soma size is 21um https://synapseweb.clm.utexas.edu/dimensions-dendrites
%so radius is aout 11um or 0.011mm
%thats 2 pixels out of 200 if view field is 1.1mm
%We thus considered a candidate set of cells to be the same neuron if all pair-wise separations were ≤6 μm.
%VARARGIN IS WEIGHTS MATRIX

if isstruct(cellcenter)==1
  cellnames = fieldnames(cellcenter);
  fieldcenternames = fieldnames(fieldcenter);
  daynum = length(cellnames);
else
  daynum =1;
end

neighper = [];
disper = [];
disper_less = [];
morall = [];
meandist = [];
cellnum = [];
perc_below = [];
shuff = [];




for z=1:(daynum)
if length(varargin)<1
  if isstruct(cellcenter)==1
    name = char(cellnames(z));
    currentcellcenter = cellcenter.(name);
    name = char(fieldcenternames(z));
    currentfieldcenter = fieldcenter.(name);
    shuff = [];
  else
    currentcellcenter = cellcenter;
    currentfieldcenter = fieldcenter;
  end

%{
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
    elseif abs(f1f-f1b) <= 20 %same field
      currentbicenter(k) = mean(f1f, f1b);
    elseif abs(f1f-f1b) > 20
      currentbicenter(k) = NaN;
    end
  end
%}

%comment this in for forward
currentbicenter = currentfieldcenter(:, 1);
%comment tis in for backward
%currentbicenter = currentfieldcenter(:, 2);

  nans = isnan(currentbicenter);
  currentbicenter = currentbicenter(~nans);



  currentcellcenter1 = currentcellcenter(1,~nans);
  currentcellcenter1 = (currentcellcenter1);
  currentcellcenter2 = currentcellcenter(2,~nans);
  currentcellcenter2 = (currentcellcenter2);
  currentcellcenter = [currentcellcenter1;currentcellcenter2];

  %for plotting moran binned
%{
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
  scatter(jitter(currentcellcenter(1,:)), currentcellcenter(2,:), sizes, colormatrix, 'filled')

%}


  [W1 W W3] = xy2cont(currentcellcenter(1,:), currentcellcenter(2,:));
  neighweights = normw(W);




  disweights = dist(currentcellcenter);
  cellnum(end+1) = length(disweights);

%FOR ONLY COUNTING CELLS MORE THAN AMOUNT
  %zero1 = find(disweights(:)<discutoff);
  %zero2 = find(disweights(:)>0);
  %zero = intersect(zero1, zero2);
  %perc_below(end+1) = length(zero)./length(zero2);

  zero = find(disweights(:)==0);


  if length(zero)>(0)
  disweights(zero) = 0;
  end
  notzero = find(disweights(:)>0);
  disweights(notzero) = 1./disweights(notzero);
  disweights = normw(disweights);




  %disweights = make_nnw(currentcellcenter(1,:)', currentcellcenter(2,:)', 6, 4);
  %disweights = make_neighborsw(currentcellcenter(1,:)', currentcellcenter(2,:)', 6);
  %disweights = pdweight(currentcellcenter(1,:)', currentcellcenter(2,:)',0,50,1);



  %PLOTTING
  figure
  subplot(ceil(daynum./6), ceil(daynum./4), z);
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
      %colormatrix(c,:)=ll;
      colormatrix(c,:)=ll;
    elseif plotval(c) ==4
      %colormatrix(c,:)=hl;
      colormatrix(c,:)=hl;
    end
  end
  sizes = 100;
  scatter((currentcellcenter(1,:)), currentcellcenter(2,:), sizes, colormatrix, 'filled')
  title('Both Directions colors')
  axis([0 200,0 200])

  f = [currentcellcenter(1,:), currentcellcenter(2,:)];

  figure
  colors = ((currentbicenter) - min(currentbicenter))./max(currentbicenter); %normalize to 0:1
  colors = colors*100;
  %colors = colors(indexes3,:);
  c = colorbar;
  %set(gca, 'clim', [min(indexes2),100]);
  sizes = 100;
  scatter((currentcellcenter(1,:)), currentcellcenter(2,:), sizes, colors, 'filled')
  set(gcf,'renderer','painters');
  colorbar
  axis([0 200,0 200])



  neighmor = (moran(currentcellcenter(1:2,:), currentbicenter, neighweights));
  dismore = (moran(currentcellcenter(1:2,:), currentbicenter, disweights));
  morall(end+1) = dismore;







elseif length(varargin)==1
  %neighweights = cell2mat(varargin);
  disweights = cell2mat(varargin);
  currentbicenter = fieldcenter;
  %neighmor = (moran(currentcellcenter(1:2,:), currentbicenter, neighweights));
  currentcellcenter = cellcenter;
  dismore = (moran(currentcellcenter(1:2,:), currentbicenter, disweights));
  morall(end+1) = dismore;
end


  mi_neigh =[];
  mi_dist = [];
  for k=1:num_times_to_shuffle
      fieldshuff = randperm(length(currentbicenter));
      newfield = currentbicenter(fieldshuff);
      %a = (moran(currentcellcenter(1:2,:), newfield, neighweights));
      b = (moran(currentcellcenter(1:2,:), newfield, disweights));

      %mi_neigh(end+1) = a;
      mi_dist(end+1) = b;

      figure
      subplot(2,1,1)
      colors = ((newfield) - min(newfield))./max(newfield); %normalize to 0:1
      colors = colors*100;
      %colors((colors<20))==0;
      %colors((colors>80))==100;
      %colors = colors(indexes3,:);
      c = colorbar;
      %set(gca, 'clim', [min(indexes2),100]);
      sizes = 100;
      scatter((currentcellcenter(1,:)), currentcellcenter(2,:), sizes, colors, 'filled')
      set(gcf,'renderer','painters');
      colorbar
      axis([0 200,0 200])

      subplot(2,1,2)
      length(find(~isnan(colors)))
      fullmat = NaN(200,200);
      for k=1:length(currentcellcenter(1,:))
        fullmat(currentcellcenter(1,k), currentcellcenter(2,k)) = colors(k);
      end
      [tempmat temp] = ndnanfilter(fullmat, 'gausswin', [10 10]);
      win = [[temp; flip(temp)], flip([temp; flip(temp)], 2)];
      [fullmat a] = ndnanfilter(fullmat, win, [10 10]);
      h = pcolor(fullmat')
      set(h, 'EdgeColor', 'none');
      axis([0 200,0 200])
      set(gcf,'renderer','painters');


      plotval = spatiallag(disweights, newfield);
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
          %colormatrix(c,:)=ll;
          colormatrix(c,:)=hh;
        elseif plotval(c) ==4
          %colormatrix(c,:)=hl;
          colormatrix(c,:)=lh;
        end
      end
      sizes = 100;
      scatter((currentcellcenter(1,:)), currentcellcenter(2,:), sizes, colormatrix, 'filled')
      title('Both Directions')
      axis([0 200,0 200])
    end

    mi_neigh =  mi_neigh';
    mi_dist = mi_dist'
    meandist(end+1) = nanmean(mi_dist);

    %neighper(end+1) = length(find(neighmor>mi_neigh))./num_times_to_shuffle;
    disper(end+1) = length(find(dismore>mi_dist))./num_times_to_shuffle;
    disper_less(end+1) = length(find(dismore<mi_dist))./num_times_to_shuffle;
end

%percents = [neighper; disper; morall; meandist; cellnum]';
percents = [disper; morall; meandist; disper_less]'; %percent greater than shuffles, actual moran, mean shuffled moran, % less than shuffled



%{
length(find(percents(:,1)>=.95))
length(find(percents(:,2)>=.95))


figure
histogram(mi_dist(:),'BinWidth', .001)
vline(morall)

nanmean(perc_below)
nanstd(perc_below)
%}
