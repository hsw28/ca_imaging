function f = moranlocal_bydist(cellcenter, fieldcenter, varargin)
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
discutoff=9.0909;
dismoran = NaN(22, length(daynum));
for dc=1:22 %should be 21 times

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

  zero = find(disweights(:)>discutoff);
  disweights(zero) = 0;

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
      %localIs(end+1) = localm;

      %dismoran()
      dismoran(dc, k) = localm;
    end


end
    discutoff = discutoff+9.0909;
end

  f =   dismoran';


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
%}

figure
%err = std(f) ./ sqrt(length(f))
%errorbar([.05:.05:1.1], mean(f), err)
plot([.05:.05:1.1], mean(f))
xlabel('mm apart')
ylabel('mean local Morans I')
