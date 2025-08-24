function s = p2stars_scalar(p)
    if isnan(p) || p>=0.05
        s = 'n.s.';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    else
        s = '*';
    end
end
