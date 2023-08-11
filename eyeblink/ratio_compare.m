function f = ratio_compare(ratios1, ratios2, days_to_average)
  %averages ratio response for ratios one and ratios 2
  %then scatter plots and finds best fit line

  ratios1 = nanmean(ratios1');
  ratios2 = nanmean(ratios2');

  x = ratios1;
  y = ratios2;

  ratios1_good = (find(~isnan(ratios1)));
  ratios2_good = (find(~isnan(ratios2)));

  all_good = intersect(ratios1_good, ratios2_good);


  x = ratios1(all_good);
  y = ratios2(all_good);



% Calculate the best-fit line coefficients using polyfit
coefficients = polyfit(x, y, 1);

% Calculate the corresponding y values for the best-fit line
yFit = polyval(coefficients, x);

% Create a scatter plot
scatter(x, y, 'b', 'filled'); % Original data points
hold on;
plot(x, yFit, 'r'); % Best-fit line

xlabel('X');
ylabel('Y');
title('Scatter Plot with Best-Fit Line');
legend('Data', 'Best-Fit Line');
grid on;

coeffs = polyfit(x, y, 1);
slope = coeffs(1); % get slope of best fit line
intercept = coeffs(2);
% Get fitted values
polydata = polyval(coeffs,x);
sstot = sum((y - mean(y)).^2);
ssres = sum((y - polydata).^2);
rsquared = 1 - (ssres / sstot) % get r^2 value

stats = fitlm(x,y);
pval = stats.Coefficients.pValue(2)

hold off;
