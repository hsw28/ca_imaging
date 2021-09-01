function [distances distancebyangle] = cellanglecorr(cellcenters, fieldcenters, angleincrement)
  %finds the distance between field centers and the ANGLE beetween cell ceenters
  %determines if there is a correlation
  %varargin is the cell number if you only want to do one cell. put in zero if you want all cells non repeating, and no input for all cells repeating


celldistance = [];
forwarddistance = [];
backwarddistance = [];
slope = [];

  for k=1:size(cellcenters,2)
    p1x = cellcenters(1,k);
    p1y = cellcenters(2,k); %cell center
    f1f = fieldcenters(k, 1); %field center forward
    f1b = fieldcenters(k, 2); %field center backward
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
         forwarddistance(end+1) = abs(f2f-f1f);
         backwarddistance(end+1) = abs(f2b-f1b);



    end
  end


  set(0,'DefaultFigureVisible', 'on');


k = min(slope(~isnan(slope)));
formean = [];
backmean = [];
binedge = [];
while k< (max(slope)+angleincrement)
  wanted = find(slope>=k & slope<k+angleincrement);
  formean(end+1) = nanmean(forwarddistance(wanted));
  backmean(end+1) = nanmean(backwarddistance(wanted));
  binedge(end+1) = k;
  k = k+angleincrement;
end
figure
subplot(1,2,1)
plot(binedge, formean);
xlabel('Bin Edge')
ylabel('Mean Distance between fields (Forward)')
subplot(1,2,2)
plot(binedge, backmean);
xlabel('Bin Edge')
ylabel('Mean Distance between fields (Forward)')
distancebyangle = [binedge; formean; backmean];

distances = [celldistance; slope; forwarddistance; backwarddistance];


%plotting for all no repeating or single cell

figure
subplot(1,2,1)
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



subplot(1,2,2)
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
