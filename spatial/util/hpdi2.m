function bounds = hpdi2(adraw,hperc)
% PURPOSE: Computes an hperc-percent HPDI for a vector of MCMC draws
% --------------------------------------------------------------------
% Usage: bounds = hpdi(draws,hperc);
% where draws = an ndraw by nvar matrix
%       hperc = 0 to 1 value for hperc percentage point
% --------------------------------------------------------------------
% RETURNS:
%         bounds = a nvar x 2 matrix with 
%         bounds(:,1) = hperc percentage point
%         bounds(:,2) = 1-hperc percentage point
%          e.g. if hperc = 0.05
%          bounds(:,1) = 0.05 point for each vector in the matrix
%          bounds(:,2) = 0.95 point  for each vector in the matrix
% --------------------------------------------------------------------

% Written by Gary Koop
% documented by J.P. LeSage

% This function takes a vector of MCMC draws and calculates
%a hperc-percent HPDI
[ndraw,ncols]=size(adraw);

p = hperc;
qlow = quantile(adraw,p);
p = 1-hperc;
qhigh = quantile(adraw,p);

bounds = zeros(ncols,2);
for i=1:ncols;
bounds(i,1) = qlow(1,i);
bounds(i,2) = qhigh(1,i);
end;

