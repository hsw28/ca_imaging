function f = fitline(x, y)
%put in a structure or vector and get p val, r2, slope

if isstruct(x)==1
  Xnames = fieldnames(x);
  Xstructnum = length(Xnames);
  Ynames = fieldnames(y);
else
  Xstructnum =1;
end

pval = [];
rsquared = [];
slope = [];

for k=1:Xstructnum
  if isstruct(x)==1
    name = char(Xnames(k));
    currentX = x.(name);
    name = char(Ynames(k));
    currentY = y.(name);
  else
    currentX = x;
    currentY = y;
  end

%  currentX = currentX(:,3);
%  currentY = currentY(:,4);

  xnan = find(isnan(currentX)==0);
  currentX = currentX(xnan);
  currentY = currentY(xnan);
  ynan = find(isnan(currentY)==0);
  currentX = currentX(ynan);
  currentY = currentY(ynan);

  coeffs = polyfit(currentX, currentY, 1);
  polydata = polyval(coeffs,currentX);
  sstot = sum((currentY - mean(currentY)).^2);
  ssres = sum((currentY - polydata).^2);
  rsquared(end+1) = 1 - (ssres / sstot); % get r^2 value
  stats = fitlm(currentX,currentY);
  pval(end+1) = stats.Coefficients.pValue(2);
  slope(end+1) = coeffs(1);

end

f = [pval', rsquared', slope'];
