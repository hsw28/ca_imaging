% PURPOSE: An example of using sar_g() Gibbs sampling
%          spatial autoregressive model
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sar_gd (see also sar_gd2 for a large data set)
%---------------------------------------------------

clear all;

n = 200;
latt = rand(n,1);
long = rand(n,1);

W = make_neighborsw(latt,long,3);

k = 3;
x = randn(n,k);

beta = ones(k,1);

rho = 0.6;

e = randn(n,1);

y = (speye(n) - rho*W)\(x*beta + e);

% 1st order contiguity matrix for
% Anselin's Columbus crime dataset
% stored in sparse matrix format [i, j, s] = find(W);
% so that W = sparse(i,j,s); reconstructs the 49x49 matrix
% NOTE: already row-standardized



info.lflag = 0;
result0 = sar(y,x,W,info);
prt(result0);

ndraw = 2500;
nomit = 500;

prior2.novi = 1; % homoscedastic model
% uses default numerical integration for rho
prior2.lflag = 0;
results2 = sar_g(y,x,W,ndraw,nomit,prior2);
results2.tflag = 'tstat';
prt(results2,vnames);

% Gibbs sampling function heteroscedastic prior
% to maximum likelihood estimates
prior.rval = 4;
prior.lflag = 0;
results = sar_g(y,x,W,ndraw,nomit,prior);
results.tflag = 'tstat';
prt(results,vnames);


[h1,f1,y1] = pltdens(results.pdraw);
[h2,f2,y2] = pltdens(results2.pdraw);
plot(y1,f1,'.r',y2,f2,'.g');
legend('heteroscedastic','homoscedastic');
title('posterior distributions for rho');
xlabel('rho values');

probs = model_probs(results,results2);

fprintf(1,'posterior model probabilities \n');
in2.rnames = strvcat('Models','heteroscedastic model','homoscedastic model');
in2.cnames = strvcat('Posterior model probabilities');
in2.fmt = '%12.6f';

mprint(probs,in2);

