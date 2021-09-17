function [lbounds, means, ubounds] = summary_interval(adraw)
% PURPOSE: Computes an hperc-percent credible interval for a vector of MCMC draws
% --------------------------------------------------------------------
% Usage: bounds = cr_interval(draws,hperc);
% where draws = an ndraw by nvar matrix
%       hperc = 0 to 1 value for hperc percentage point
% --------------------------------------------------------------------
% RETURNS:
%         means   = a nvar x 1 vector with means
%         ubounds = a nvar x 2 vector with 
%         ubounds(nvar,2) = 0.95 and 0.99 percentage points
%         lbounds = a nvar x 2 vector with
%         lbounds(nvar,2) = 0.01 and 0.05 percentage points
% --------------------------------------------------------------------

% written by:
% James P. LeSage, 4/2018
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com


% This function takes a vector of MCMC draws and calculates
% an hperc-percent credible interval
[ndraw,ncols]=size(adraw);

means = mean(adraw)';

lower01 = quantile(adraw,0.01)';
lower05 = quantile(adraw,0.05)';

upper95 = quantile(adraw,0.95)';
upper99 = quantile(adraw,0.01)';

lbounds = [lower01 lower05];
ubounds = [upper95 upper01];


