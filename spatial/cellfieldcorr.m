function [distances pvals] = cellfieldcorr(cellcenters, fieldcenters, varargin)
  %finds the distance between field centers and the distance beetween cell ceenters
  %determines if there is a correlation
  %varargin is the cell number if you only want to do one cell. put in zero if you want all cells non repeating, and no input for all cells repeating


celldistance = [];
forwarddistance = [];
backwarddistance = [];
bothdistance = [];

bicenter = fieldcenters;
%{
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
%}

%for one cell
if length(cell2mat(varargin))>0 && cell2mat(varargin)~=0
  i = cell2mat(varargin);
  p1 = [cellcenters(1,i); cellcenters(2,i)]; %cell center
  f1f = fieldcenters(i, 1); %field center forward
  f1b = fieldcenters(i, 2); %field center backward
  f1a = bicenter(k);
  for j=1:size(cellcenters,2) %goes through other cells
      if j == i
        forwarddistance(end+1) = NaN;
        backwarddistance(end+1) = NaN;
        celldistance(end+1) = NaN;
        bothdistance(end+1) = NaN;

      else
      %finding distance between cell centers
       p2 = [cellcenters(1,j); cellcenters(2,j)];
       celldistance(end+1) = abs(norm(p1 - p2));

       %distance between field centers
       f2f = fieldcenters(j, 1); %field center forward
       f2b = fieldcenters(j, 2); %field center backward
       f2a = bicenter(j);
       bothdistance(end+1) = abs(f2a-f1a);
       forwarddistance(end+1) = abs(f2f-f1f);
       backwarddistance(end+1) = abs(f2b-f1b);
     end
  end

%for all cells, not repeating
elseif length(cell2mat(varargin))>0 && cell2mat(varargin)==0
  for k=1:size(cellcenters,2)
    p1 = [cellcenters(1,k); cellcenters(2,k)]; %cell center
    f1f = fieldcenters(k, 1); %field center forward
    f1b = fieldcenters(k, 2); %field center backward
    f1a = bicenter(k);
    for j=k+1:size(cellcenters,2) %goes through other cells

        %finding distance between cell centers
         p2 = [cellcenters(1,j); cellcenters(2,j)];
         celldistance(end+1) = abs(norm(p1 - p2));

         %distance between field centers
         f2f = fieldcenters(j, 1); %field center forward
         f2b = fieldcenters(j, 2); %field center backward
         f2a = bicenter(j);
         bothdistance(end+1) = abs(f2a-f1a);
         forwarddistance(end+1) = abs(f2f-f1f);
         backwarddistance(end+1) = abs(f2b-f1b);


    end
  end

%for all cells, repeating
else

  pvalforward = NaN(size(cellcenters,2),1);
  pvalbackwards = NaN(size(cellcenters,2),1);
  rsquaredforward = NaN(size(cellcenters,2),1);
  rsquaredbackward = NaN(size(cellcenters,2),1);
  for k=1:size(cellcenters,2)
    celldistance = [];
    forwarddistance = [];
    backwarddistance = [];
    p1 = [cellcenters(1,k); cellcenters(2,k)]; %cell center
    f1f = fieldcenters(k, 1); %field center forward
    f1b = fieldcenters(k, 2); %field center backward
    f1a = bicenter(k);
    for j=1:size(cellcenters,2) %goes through other cells
      if j == k
        forwarddistance(end+1) = NaN;
        backwarddistance(end+1) = NaN;
        celldistance(end+1) = NaN;
        bothdistance(end+1) = NaN;
       else
        %finding distance between cell centers
         p2 = [cellcenters(1,j); cellcenters(2,j)];
         celldistance(end+1) = abs(norm(p1 - p2));

         %distance between field centers
         f2f = fieldcenters(j, 1); %field center forward
         f2b = fieldcenters(j, 2); %field center backward
         f2a = bicenter(j);
        bothdistance(end+1) = abs(f2a-f1a);
         forwarddistance(end+1) = abs(f2f-f1f);
         backwarddistance(end+1) = abs(f2b-f1b);
        end
      end


    idx = (isnan(forwarddistance));
    x = celldistance(~idx);
    y = forwarddistance(~idx);
    if length(x)>0
    stats = fitlm(x,y);
    pval = stats.Coefficients.pValue(2);
    pvalforward(k) = pval;
    coeffs = polyfit(x, y, 1);
    polydata = polyval(coeffs,x);
    sstot = sum((y - mean(y)).^2);
    ssres = sum((y - polydata).^2);
    rsquaredforward(k) = 1 - (ssres / sstot); % get r^2 value
    end

    idx = (isnan(forwarddistance));
    x = celldistance(~idx);
    y = backwarddistance(~idx);
    if length(x)>0
   stats = fitlm(x,y);
   pval = stats.Coefficients.pValue(2);
   pvalbackwards(k) = pval;
   coeffs = polyfit(x, y, 1);
   polydata = polyval(coeffs,x);
   sstot = sum((y - mean(y)).^2);
   ssres = sum((y - polydata).^2);
   rsquaredbackward(k) = 1 - (ssres / sstot); % get r^2 value
  end
  end
  pvalr2forward = [pvalforward, rsquaredforward];
  pvalr2backwards = [pvalbackwards, rsquaredbackward];
end


distances = [celldistance;forwarddistance;backwarddistance];


%pvals = [pvalr2forward; pvalr2backwards]



%plotting for all no repeating or single cell
if length(cell2mat(varargin))>0
figure
subplot(1,3,1)
hold on
idx = (isnan(forwarddistance));
x = celldistance(~idx);
y = forwarddistance(~idx);
scatter(x, y)
coeffs = polyfit(x, y, 1);
slope = coeffs(1);
polydata = polyval(coeffs,x);
sstot = sum((y - mean(y)).^2);
ssres = sum((y - polydata).^2);
rsquaredforward = 1 - (ssres / sstot); % get r^2 value
stats = fitlm(x,y);
pvalforward = stats.Coefficients.pValue(2);
y = polyval(coeffs,x);
if pvalforward <= .05
  plot(x, y, 'r', 'LineWidth', 5) % best fit line
else
  plot(x, y, 'black', 'LineWidth', 5) % best fit line
end
str1 = {'p value' pvalforward, 'r2 value' rsquaredforward};
text(25,max(forwarddistance)*.8,str1,'FontSize',12);
title('Distance between cells versus forward field distance')
xlabel('Distance between cells')
ylabel('Distance between fields')

set(0,'DefaultFigureVisible', 'on');


subplot(1,3,2)
hold on
idx = (isnan(backwarddistance));
x = celldistance(~idx);
y = backwarddistance(~idx);
scatter(x, y);
coeffs = polyfit(x, y, 1);
slope = coeffs(1);
polydata = polyval(coeffs,x);
sstot = sum((y - mean(y)).^2);
ssres = sum((y - polydata).^2);
rsquaredbackward = 1 - (ssres / sstot); % get r^2 value
stats = fitlm(x,y);
pvalbackwards = stats.Coefficients.pValue(2);
y = polyval(coeffs,x);
if pvalbackwards <= .05
  plot(x, y, 'r', 'LineWidth', 5) % best fit line
else
  plot(x, y, 'black', 'LineWidth', 5) % best fit line
end
str1 = {'p value' pvalbackwards, 'r2 value' rsquaredbackward};
text(25,max(backwarddistance)*.8,str1,'FontSize',12);
title('Distance between cells versus backward field distance')
xlabel('Distance between cells')
ylabel('Distance between fields')



subplot(1,3,3)
hold on
idx = (isnan(bothdistance));
x = celldistance(~idx);
y = bothdistance(~idx);
scatter(x, y);
coeffs = polyfit(x, y, 1);
slope = coeffs(1);
polydata = polyval(coeffs,x);
sstot = sum((y - mean(y)).^2);
ssres = sum((y - polydata).^2);
rsquaredboth = 1 - (ssres / sstot); % get r^2 value
stats = fitlm(x,y);
pvalboth = stats.Coefficients.pValue(2);
y = polyval(coeffs,x);
if pvalboth <= .05
  plot(x, y, 'r', 'LineWidth', 5) % best fit line
else
  plot(x, y, 'black', 'LineWidth', 5) % best fit line
end
str1 = {'p value' pvalboth, 'r2 value' rsquaredboth};
text(25,max(backwarddistance)*.8,str1,'FontSize',12);
title('Distance between cells versus both field distance')
xlabel('Distance between cells')
ylabel('Distance between fields')

pvalr2forward = [pvalforward, rsquaredforward];
pvalr2backwards = [pvalbackwards, rsquaredbackward];
pvalr2both = [pvalboth, rsquaredboth];

pvals = [pvalr2forward; pvalr2backwards; pvalr2both]
end
