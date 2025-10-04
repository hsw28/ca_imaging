function [f,y] = ksdensity_fallback(x)
try
    [f,y] = ksdensity(x);
catch
    nb = max(10, round(sqrt(numel(x))));
    [cnt,edges] = histcounts(x, nb, 'Normalization','pdf');
    y = movmean(edges,2,'Endpoints','discard');
    f = movmean(cnt,3,'Endpoints','shrink');
    yq = linspace(min(y), max(y), 200);
    f  = interp1(y, f, yq, 'pchip','extrap');
    y  = yq;
end
end
