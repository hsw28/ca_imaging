function f = spatiallag(weightmatrix, values)
  %finds spatial lags and makes a moran plot

weightmatrix = normw(weightmatrix);

if size(values,1)>size(values,2)
 values = values';
end


lag = weightmatrix .* values;
lag = mean(lag')';

x = values;
y = lag';
scatter(values, lag);

hold on
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
