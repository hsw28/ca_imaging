function x = beta_inv(p, a, b, crit)
% PURPOSE: inverse of the cdf (quantile) of the beta(a,b) distribution
%--------------------------------------------------------------
% USAGE: x = beta_inv(p,a,b)
% where:   p = vector of probabilities
%          a = beta distribution parameter, a = scalar
%          b = beta distribution parameter  b = scalar
% NOTE: mean [beta(a,b)] = a/(a+b), variance = ab/((a+b)*(a+b)*(a+b+1))
%--------------------------------------------------------------
% RETURNS: x at each element of p for the beta(a,b) distribution
%--------------------------------------------------------------
% SEE ALSO: beta_d, beta_pdf, beta_inv, beta_rnd
%--------------------------------------------------------------

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg
% documentation modified by LeSage to
% match the format of the econometrics toolbox

if (nargin ~= 3 & nargin ~= 4)
    error('Wrong # of arguments for beta_inv');
end
 
if any(any((a<=0)|(b<=0)))
   error('beta_inv parameter a or b is non-positive');
end
if any(any(abs(2*p-1)>1))
   error('beta_inv: A probability should be 0<=p<=1');
end
if(nargin == 3)
    crit.maxiters = 1000;
    crit.tol = 4096*eps;
end
x = a ./ (a+b);
dx = 1; iters = 1;
while any(any(abs(dx)>crit.tol*max(x,1))) & iters < crit.maxiters
   dx = (betainc(x,a,b) - p) ./ beta_pdf(x,a,b);
   x = x - dx;
   x = x + (dx - x) / 2 .* (x<0);
   iters = iters+1;
end
    

 
