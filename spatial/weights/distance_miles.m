function D = distance_miles(xc,yc)
% PURPOSE: Computes the list of pairwise distances (in miles) 
%          for a given set of latt,long points
% ----------------------------------------------------------
% Usage: D = distance_miles(latt,long)
% where: latt,long are vectors of latitude, longitude coordinates for each location
% ----------------------------------------------------------
% Returns: D = (n x n)-matrix of pairwise distances

% Written by: Jim LeSage 2/2011

% miles = sqrt(x * x + y * y)
% where x = 69.1 * (lat2 - lat1)
% and y = 69.1 * (lon2 - lon1) * cos(lat1/57.3)


n = length(xc) ;  %number of locations
for i=1:n;
    for j=1:n;
        X = 69.1*(xc(i,1) - xc(j,1));
        Y = 69.1*(yc(i,1) - yc(j,1))*cos(xc(j,1)/57.3);
        D(i,j) = sqrt(X*X + Y*Y);
    end;
end;

