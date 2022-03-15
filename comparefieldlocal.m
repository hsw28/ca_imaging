function f = comparefieldlocal(alignmentdata, varargin)
  %makes a scattre plot showing placefieelds across sessions
  %can put in daysrecorded as a vector like [1, 2, 10, 14, 15, 16] to see how r2 values change

%figure
daysrecorded = cell2mat(varargin);

rs = [];
daysapart = [];

total=0;
pval = [];
for k=1:size(alignmentdata, 2)-1
  for z=k+1:size(alignmentdata, 2)
    subplot(5,5,total+z-1)
    x = alignmentdata(:,k);
    y = alignmentdata(:,z);
    xx = find(x>0);
    yy = find(y>0);
    xxyy = intersect(xx,yy);
    x = x(xxyy);
    y = y(xxyy);
  if length(x)>10
    scatter(x, y, 'filled', 'black')
    hold on

    coeffs = polyfit(x, y, 1);
    posslope = coeffs(1);
    polydata = polyval(coeffs,x);
    sstot = sum((y - mean(y)).^2);
    ssres = sum((y - polydata).^2);
    posrsquared = 1 - (ssres / sstot); % get r^2 value
    rs(end+1) = posrsquared;

    stats = fitlm(x,y);
    pospval = stats.Coefficients.pValue(2);
    y = polyval(coeffs,x);



    daysapart(end+1) = daysrecorded(z)-daysrecorded(k);


    if round(pospval)<=.05
      plot(x, y, 'r', 'LineWidth', 1.5) % best fit line
    else
      plot(x, y, 'black', 'LineWidth', 1) % best fit line
    end
    %str1 = {'p=' round(pospval,3), 'r2=' round(posrsquared,2)};
    %text(max(x)*.2,max(y)*.8,str1,'FontSize',12);

  end
    pval(end+1) = pospval;

  end
    total = total+5;
end

  f = [rs; daysapart]';
