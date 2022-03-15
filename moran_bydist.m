function [f num] = moran_bydist(cellcenter, fieldcenter, varargin)
%finds local moran's by increments of distance in 0.05mm

%varargin is yyour wieight matrrix if you like
%200pixels,  view field is 1.1mm
%each pixel is 0.0055mm
% 9.0909 pixels is .05mm


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


%200pixels,  view field is 1.1mm
%each pixel is 0.0055mm
% 9.0909 pixels is .05mm
%1.8182 is 0.01
%0.9091 is .005



dismoran_shuff = NaN(50, 4, (daynum));
dismoran = NaN((daynum),50);
num = NaN((daynum),50);

for z=1:(daynum)
  discutoff=0;
  if isstruct(cellcenter)==1
    name = char(cellnames(z))
    currentcellcenter = cellcenter.(name);
    name = char(fieldcenternames(z))
    currentfieldcenter = fieldcenter.(name);
  else
    currentcellcenter = cellcenter;
    currentfieldcenter = fieldcenter;
  end



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







for dc=1:50 %should be 50 times




  disweights = dist(currentcellcenter);



  %zero1 = find(disweights(:)>discutoff); %FOR ONLY COUNTING CELLS LESS THAN AMOUNT
  zero1 = find(disweights(:)<discutoff); %FOR ONLY COUNTING CELLS MORE THAN AMOUNT
  %zero1 = [find(disweights(:)<discutoff-9.0909); find(disweights(:)>(discutoff))];  %window


  zero2 = find(disweights(:)>0);
  zero = intersect(zero1, zero2);

  if length(zero)>(0)
    disweights(zero) = 0;
  end

      notzero = find(disweights(:)>0);
      disweights(notzero) = 1./disweights(notzero);
      disweights = normw(disweights);

      val = (moran(currentcellcenter, currentbicenter, disweights)); %for just moran
      val_shuff = shuffled_moran(currentcellcenter, currentbicenter, 1000, disweights); %for shuffled moran

      dismoran(z,dc) = val;
      num(z,dc) = sqrt(length(notzero));


      dismoran_shuff(dc, :, (z))  = val_shuff;

      discutoff = discutoff+0.9091;

end

end




  f =   dismoran;
  f = dismoran_shuff;

%{
  histogram(f, 'BinWidth', .1, 'Normalization','probability')
  ylabel('Proportion of Cells')
  xlabel('Local Morans I')
  vline(mean(f))


  figure
  colors = f; %normalize to 0:1
  %colors = colors(indexes3,:);
  c = colorbar;
  %set(gca, 'clim', [min(indexes2),100]);
  sizes = 100;
  scatter(jitter(currentcellcenter(1,:)), currentcellcenter(2,:), sizes, colors, 'filled')
  title('Both Directions')
  axis([0 200,0 200])
  colorbar



%err = std(f) ./ sqrt(length(f))
%errorbar([.05:.05:1.1], mean(f), err)

figure

plot([0:.005:.25-.005], nanmean(f))
xlabel('mm apart')
ylabel('mean local Morans I')

set(0,'DefaultFigureVisible', 'on');

%}
