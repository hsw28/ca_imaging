function plotval= spatiallag(weightmatrix, values)
  %finds spatial lags and makes a moran plot

weightmatrix = normw(weightmatrix);

if size(values,1)>size(values,2)
 values = values';
end

sx = sqrt(((1./(length(values))).*(sum(values.^2)))-(mean(values).^2));
x = (values-mean(values))./ sx;

lag = weightmatrix .* x;
lag = sum(lag')';


y = lag';


figure
scatter(x, lag);

hold on
coeffs = polyfit(x, y, 1);



polydata = polyval(coeffs,x);
sstot = sum((y - mean(y)).^2);
ssres = sum((y - polydata).^2);
rsquared = 1 - (ssres / sstot); % get r^2 value
stats = fitlm(x,y);
pvalforward = stats.Coefficients.pValue(2);
y = polyval(coeffs,x);
title('Spatial Lag')
xlabel('Standardized Values')
ylabel('Spatial Lag')
if pvalforward <= .05
  plot(x, y, 'r', 'LineWidth', 5) % best fit line
else
  plot(x, y, 'black', 'LineWidth', 5) % best fit line
end
str1 = {'p value' pvalforward, 'r2 value' rsquared, 'slope' coeffs(1)}
text(max(x)*.7,min(lag),str1,'FontSize',8);
vline(0)
hline(0)



smallx = find(x<0);
bigx = find(x>0);
smalllag = find(lag<0);
biglag = find(lag>0);
lh = intersect(smallx, biglag);
hh = intersect(bigx, biglag);
ll = intersect(smallx, smalllag);
hl = intersect(bigx, smalllag);
plotval = NaN(length(x),1);
plotval(lh) = 1;
plotval(hh) = 2;
plotval(ll) = 3;
plotval(hl) = 4;
