function mi = ca_mutualinfo(p)
%MUTUALINFO computes mutual information
%
%  c=MUTUALINFO(p) computes the mutual information from joint probability
%  distribution p.
%
%  c=MUTUALINFO(p,base) computes mutual information using the logarithm
%  with the specified base (default=2).
%


base = 2;

%normalize
p = p./nansum(p(:));

%compute marginals
m1 = nansum( p, 1 );
m2 = nansum( p, 2 );

%compute mutual information
mi = p.*log(p./bsxfun(@times,m1,m2))./log(base);
mi = nansum( mi(:) );
