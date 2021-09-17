% simulate model using W, then calculate log-marginals

% y = rho*W*y + x*beta + W*x*gamma + e

% DGP: 
% yt = inv(I_n - rho*W) *(X*beta + W*X*gamma + varepsilon)

clear all;

load schools.dat;
% col 1 = school district ID
% col 2 = longitude centroid for the district
% col 3 = latitude centroid for the district

long = schools(:,2);
latt = schools(:,3);

W = make_neighborsw(latt,long,6);
  
[n,junk] = size(W);

N = n;
  
rng(86573);

sigx = 1;
x = randn(N,1)*sqrt(sigx);
x1 = x;

x = randn(N,1)*sqrt(sigx);
x2 = x;

x = randn(N,1)*sqrt(sigx);
x3 = x;

% simulate y


xo = [x1 x2 x3];

tmp = [ -0.5 1 0.5];
beta = tmp';

tmp = [-1 0.5 1]; % make sure these don't equal -rho*beta
gamm = tmp';

sige = 1;

alpha = 10;


% ============================================

xmat = [ones(N,1) xo W*xo]; % model includes W*x-variables

beta_gamma = [alpha
    beta
    gamm];


    
    rho = 0.5;
    lam = 0.4;

 
F = speye(N) - rho*W;
G = speye(N) - lam*W;

xmat = [ones(N,1) xo W*xo]; % model includes W*x-variables

beta_gamma = [alpha
    beta
    gamm];

eterm = randn(N,1)*sqrt(sige);

tmp = [xmat*beta_gamma];


y_slx = tmp + eterm;

result1 = lmarginal_cross_section(y_slx,xo,W); 

fprintf(1,'true model is SLX \n');
fprintf(1,'time taken is: %16.4f seconds \n',result1.time);

in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out = [result1.lmarginal result1.probs];

mprint(out,in);


y_sdm = F\(tmp + eterm);

result2 = lmarginal_cross_section(y_sdm,xo,W); 

fprintf(1,'true model is SDM \n');
fprintf(1,'time taken is: %16.4f seconds \n',result2.time);

fprintf('rho = %10.4f  \n',0.5);

in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out = [result2.lmarginal result2.probs];

mprint(out,in);


y_sdem = tmp + G\eterm;

result3 = lmarginal_cross_section(y_sdem,xo,W); 

fprintf(1,'true model is SDEM \n');
fprintf(1,'time taken is: %16.4f seconds \n',result3.time);

fprintf('lambda = %10.4f  \n',0.4);

in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out = [result3.lmarginal result3.probs];

mprint(out,in);

