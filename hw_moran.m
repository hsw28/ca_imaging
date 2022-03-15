function [weights cellcenterf] = hw_moran(cellcenter, fieldcenter)
  %DEPRICATED
  %runs a moran test for shared neighbors and distance matrix
%see https://cran.r-project.org/web/packages/Irescale/vignettes/rectifiedI.html

%forward
nans = isnan(fieldcenter(:,1));
fieldcenterf = fieldcenter(~nans,1);
cellcenter1 = cellcenter(1,~nans);
cellcenter2 = cellcenter(2,~nans);
cellcenterf = [cellcenter1;cellcenter2];


[W1 W W3] = xy2cont(cellcenterf(1,:), cellcenterf(2,:));
neighweights = W;
disweights = dist(cellcenterf);
notzero = find(disweights(:)>0);
disweights(notzero) = 1./disweights(notzero);


prt(moran(cellcenterf(1:2,:), fieldcenterf, neighweights));
prt(moran(cellcenterf(1:2,:), fieldcenterf, disweights));

%backward
nans = isnan(fieldcenter(:,2));
fieldcenterb = fieldcenter(~nans,2);
cellcenter1 = cellcenter(1,~nans);
cellcenter2 = cellcenter(2,~nans);
cellcenterb = [cellcenter1;cellcenter2];


[W1 W W3] = xy2cont(cellcenterb(1,:), cellcenterb(2,:));
neighweights = W;
disweights = dist(cellcenterb);
notzero = find(disweights(:)>0);
disweights(notzero) = 1./disweights(notzero);



prt(moran(cellcenterb(1:2,:), fieldcenterb, neighweights))
prt(moran(cellcenterb(1:2,:), fieldcenterb, disweights))


%both
bicenter = NaN(size(cellcenter,2),1);
for k=1:size(cellcenter,2)
  f1f = fieldcenter(k, 1); %field center forward
  f1b = fieldcenter(k, 2); %field center backward
  if isnan(f1f)==1 && isnan(f1b)==1
    bicenter(k) = NaN;
  elseif isnan(f1f)==1
    bicenter(k) = f1b;
  elseif isnan(f1b)==1
    bicenter(k) = f1f;
  elseif abs(f1f-f1b) <= 15 %same field
    bicenter(k) = mean(f1f, f1b);
  elseif abs(f1f-f1b) > 15
    bicenter(k) = NaN;
  end
end

nans = isnan(bicenter);
bicenter = bicenter(~nans);
cellcenter1 = cellcenter(1,~nans);
cellcenter2 = cellcenter(2,~nans);
cellcenter = [cellcenter1;cellcenter2];


[W1 W W3] = xy2cont(cellcenter(1,:), cellcenter(2,:));
neighweights = normw(W);
disweights = dist(cellcenter);
notzero = find(disweights(:)>0);
disweights(notzero) = 1./disweights(notzero);
disweights = normw(disweights);


(moran(cellcenter(1:2,:), bicenter, neighweights));
(moran(cellcenter(1:2,:), bicenter, disweights));

figure
axis([min(cellcenter1)-1, max(cellcenter1)+1, min(cellcenter2)-1, max(cellcenter2)+1])
for ii = 1:length(cellcenter1)
    t = text(cellcenter1(ii),cellcenter2(ii),num2str(ii));
    hold on
end

figure
for ii = 1:length(cellcenter1)
  jj = ii+1;
  while jj <=length(cellcenter1)
    if neighweights(jj,ii)>0
      hold on
      x = [cellcenter1(ii), cellcenter1(jj)];
      y = [cellcenter2(ii), cellcenter2(jj)];
      plt = plot(x,y);
      color = plt.Color;
      x2 = [(cellcenter1(jj)+cellcenter1(ii))./2];
      y2 = [(cellcenter2(jj)+cellcenter2(ii))./2];
      plt = plot(x2,y2, 'LineStyle', 'none');
      num = round(neighweights(jj,ii)*100)/100;
      label(plt,num2str(num),'location','center', 'color', color)
      %plot(x,y,'Color', color);

    end
    jj = jj+1;
  end
end
