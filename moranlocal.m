function f = moranlocal(cellcenter, fieldcenter, varargin)
%finds local moran's i
%varargin is yyour wieight matrrix if you like


%{
However, Moran’s 𝐼 is does not indicate areas within the map where specific types
of values (e.g. high, low) are located. In other words, Moran’s I can tell us whether
values in our map cluster together (or disperse) overall, but it will not inform us
about where specific clusters (or outliers) are.

The core idea of a local Moran’s 𝐼𝑖 is to identify cases in which the value of an
observation and the average of its surroundings is either more similar (HH or LL in
the scatterplot above) or dissimilar (HL, LH) than we would expect from pure chance.
The mechanism to do this is similar to the one in the global Moran’s I, but applied
in this case to each observation.

it is important to keep in mind that the high positive values arise from value
similarity in space, and this can be due to either high values being next to high
values or low values next to low values. The local 𝐼𝑖 values alone cannot distinguish
between these two.

The values in the left tail of the density represent locations displaying negative
spatial association. There are also two forms, a high value surrounded by low values,
or a low value surrounded by high valued neighboring observations. And, again,
the 𝐼𝑖 value cannot distinguish between the two cases.
%}

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
  for k=1:size(currentfieldcenter,1)
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
  currentcellcenter1 = (currentcellcenter1);
  currentcellcenter2 = currentcellcenter(2,~nans);
  currentcellcenter2 = (currentcellcenter2);
  currentcellcenter = [currentcellcenter1;currentcellcenter2];


  disweights = dist(currentcellcenter);
  cellnum(end+1) = length(disweights);


  notzero = find(disweights(:)>0);
  disweights(notzero) = 1./disweights(notzero);
  disweights = normw(disweights);


  y = currentcellcenter;
  x = currentbicenter;
  i = y(1,:);
  j = y(2,:);
  n = length(x);

  zi = x-mean(x);
  secondmoment = sum(zi.^2)./n;

  localIs = [];
  for k = 1:length(x)
      zi = x(k)-mean(x);
      mul = zi./secondmoment;
      summer = sum(disweights(k,:).*x');
      localm = mul.*summer;
      localIs(end+1) = localm;
    end
  end

  f =   localIs';
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

f = [localIs; currentcellcenter]';'
