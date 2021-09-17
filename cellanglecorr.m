function [distances distancebyangle temp] = cellanglecorr(cellcenters, fieldcenters, angleincrement)
  %finds the distance between field centers and the ANGLE beetween cell ceenters
  %determines if there is a correlation
  %varargin is the cell number if you only want to do one cell. put in zero if you want all cells non repeating, and no input for all cells repeating




bicenter = NaN(size(cellcenters,2),1);
for k=1:size(cellcenters,2)
  f1f = fieldcenters(k, 1); %field center forward
  f1b = fieldcenters(k, 2); %field center backward
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

celldistance = [];
forwarddistance = [];
backwarddistance = [];
bothdistance = [];
slope = [];

  for k=1:size(cellcenters,2)
    p1x = cellcenters(1,k);
    p1y = cellcenters(2,k); %cell center
    f1f = fieldcenters(k, 1); %field center forward
    f1b = fieldcenters(k, 2); %field center backward
    f1both = bicenter(k);
    for j=k+1:size(cellcenters,2) %goes through other cells

        %finding distance between cell centers
         p2x = cellcenters(1,j);
         p2y = cellcenters(2,j);

         if p2x > p1x && isnan(p2x)==0 && isnan(p1x)==0
           slope(end+1) = atand((p2y-p1y)./(p2x-p1x));
         elseif isnan(p2x)==0 && isnan(p1x)==0
           slope(end+1) = atand((p1y-p2y)./(p1x-p2x));
         else
           slope(end+1) = NaN;
         end

         celldistance(end+1) = abs(norm([p1x,p1y]-[p2x,p2y]));

         %distance between field centers
         f2f = fieldcenters(j, 1); %field center forward
         f2b = fieldcenters(j, 2); %field center backward
         f2both = bicenter(j);
         forwarddistance(end+1) = abs(f2f-f1f);
         backwarddistance(end+1) = abs(f2b-f1b);
         bothdistance(end+1) = abs(f2both-f1both);



    end
  end


  set(0,'DefaultFigureVisible', 'on');


k = min(slope(~isnan(slope)));
formean = [];
backmean = [];
bothmean = [];
binedge = [];
count = 1;
allnum = NaN(length(bothdistance),1);
while k< (max(slope)+angleincrement)
  wanted = find(slope>=k & slope<k+angleincrement);

  if  length(wanted)==0 & abs(k-max(slope))<=angleincrement
  break
  else
    formean(end+1) = nanmean(forwarddistance(wanted));
    backmean(end+1) = nanmean(backwarddistance(wanted));
    bothmean(end+1) = nanmean(bothdistance(wanted))
    %allnum(count, 1:length(wanted)) = bothdistance(wanted(:));
    %nanmean(allnum(count,1:length(wanted)))
    %allnum(count, length(wanted):end) = NaN;
    %nanmean(allnum(count,:))
    allnum(wanted) = count;
  end


  binedge(end+1) = k;
  count = count+1;
  k = k+angleincrement;
end

figure
subplot(1,3,1)
plot(binedge, formean);
xlabel('Bin Edge')
ylabel('Mean Distance between fields (Forward)')
subplot(1,3,2)
plot(binedge, backmean);
xlabel('Bin Edge')
ylabel('Mean Distance between fields (Backwards)')
distancebyangle = [binedge; formean; backmean; bothmean];
subplot(1,3,3)
plot(binedge, bothmean);
xlabel('Bin Edge')
ylabel('Mean Distance between fields (Both directions)')
[p,tbl,stats]  = anovan(bothdistance, allnum)
temp = [bothdistance; allnum'];

distances = [celldistance; slope; forwarddistance; backwarddistance; bothdistance];


%plotting for all no repeating or single cell

figure
subplot(1,3,1)
hold on
idx = (isnan(slope));
x = slope(~idx);
y = forwarddistance(~idx);
idx = (isnan(y));
x = x(~idx);
y = y(~idx);
scatter(x, y)
coeffs = polyfit(x, y, 1);




polydata = polyval(coeffs,x);
sstot = sum((y - mean(y)).^2);
ssres = sum((y - polydata).^2);
rsquared = 1 - (ssres / sstot); % get r^2 value
stats = fitlm(x,y);
pvalforward = stats.Coefficients.pValue(2);
y = polyval(coeffs,x);
if pvalforward <= .05
  plot(x, y, 'r', 'LineWidth', 5) % best fit line
else
  plot(x, y, 'black', 'LineWidth', 5) % best fit line
end
str1 = {'p value' pvalforward, 'r2 value' rsquared};
text(25,max(forwarddistance)*.8,str1,'FontSize',12);
title('Slope between cells versus forward field distance')
xlabel('Slope between cells')
ylabel('Distance between fields')



subplot(1,3,2)
hold on
idx = (isnan(slope));
x = slope(~idx);
y = backwarddistance(~idx);
idx = (isnan(y));
x = x(~idx);
y = y(~idx);

scatter(x, y);
coeffs = polyfit(x, y, 1);

polydata = polyval(coeffs,x);
sstot = sum((y - mean(y)).^2);
ssres = sum((y - polydata).^2);
rsquared = 1 - (ssres / sstot); % get r^2 value
stats = fitlm(x,y);
pvalbackwards = stats.Coefficients.pValue(2);
y = polyval(coeffs,x);
if pvalbackwards <= .05
  plot(x, y, 'r', 'LineWidth', 5) % best fit line
else
  plot(x, y, 'black', 'LineWidth', 5) % best fit line
end
str1 = {'p value' pvalbackwards, 'r2 value' rsquared};
text(25,max(backwarddistance)*.8,str1,'FontSize',12);
title('Slope between cells versus backward field distance')
xlabel('Slope between cells')
ylabel('Distance between fields')


subplot(1,3,3)
hold on
idx = (isnan(slope));
x = slope(~idx);
y = bothdistance(~idx);
idx = (isnan(y));
x = x(~idx);
y = y(~idx);

scatter(x, y);
coeffs = polyfit(x, y, 1);

polydata = polyval(coeffs,x);
sstot = sum((y - mean(y)).^2);
ssres = sum((y - polydata).^2);
rsquared = 1 - (ssres / sstot); % get r^2 value
stats = fitlm(x,y);
pvalbackwards = stats.Coefficients.pValue(2);
y = polyval(coeffs,x);
if pvalbackwards <= .05
  plot(x, y, 'r', 'LineWidth', 5) % best fit line
else
  plot(x, y, 'black', 'LineWidth', 5) % best fit line
end
str1 = {'p value' pvalbackwards, 'r2 value' rsquared};
text(25,max(backwarddistance)*.8,str1,'FontSize',12);
title('Slope between cells versus field distance (both directions)')
xlabel('Slope between cells')
ylabel('Distance between fields')
