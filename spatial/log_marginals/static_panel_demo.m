% simulate model using W, then calculate log-marginals

% yt = rho*W*yt + xt*beta +W*xt*gamma + e

% DGP: 
% yt = inv[I_n*t - rho*(I_t \otimes W) x (X*beta + W*X*gamma + varepsilon)
% no fixed effects used here


clear all;

load queenW.data;

Wc = normw(queenW);
W = Wc;
  
[n,junk] = size(W);

N = n;

T=30; % number of time periods
  
rng(86573); % sets random number generator seed value

sigx = 1;
x = randn(N,T)*sqrt(sigx);
x1 = x;

x = randn(N,T)*sqrt(sigx);
x2 = x;

x = randn(N,T)*sqrt(sigx);
x3 = x;

% simulate yt

xvec1 = (vec(x1));

xvec2 = (vec(x2));

xvec3 = (vec(x3));

xo = [xvec1 xvec2 xvec3];

tmp = [ -0.5 1 0.5];
beta = tmp';

tmp = [-1 0.5 1]; % make sure these don't equal -rho*beta
gamm = tmp';

sige = 1;

alpha = 10;


% ============================================
kron2= sparse(kron(speye(T),W));

xmat = [ones(N*T,1) xo kron2*xo]; % model includes W*x-variables

beta_gamma = [alpha
    beta
    gamm];


    
    rho = 0.5;
    lam = 0.4;

 
F = kron(speye(N),speye(T)) - rho*kron(speye(T),W);
G = kron(speye(N),speye(T)) - lam*kron(speye(T),W);

kron2= sparse(kron(speye(T),W));

xmat = [ones(N*T,1) xo kron2*xo]; % model includes W*x-variables

beta_gamma = [alpha
    beta
    gamm];

eterm = randn(N*T,1)*sqrt(sige);

% generate fixed effects
% tmp = eye(N);
% FE = tmp(:,1:N-1);
% bfe = unif_rnd(N-1,-1,1);

tmp = [xmat*beta_gamma];


y_slx = tmp + eterm;
% ted = 0;
% [y_slxs,xslx,No,To,Wf3]=demeanF(y_slx,xo,N,T,ted,W);

resulto = ols(y_slx,xmat);

result1 = lmarginal_static_panel(y_slx,xo,W,N,T); 

fprintf(1,'true model is SLX \n');
in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out = [result1.lmarginal result1.probs];

mprint(out,in);


y_sdm = F\(tmp + eterm);
% ted = 0;
% [y_sdms,xsdms,No,To,Wf]=demeanF(y_sdm,xo,N,T,ted,W);

result2 = lmarginal_static_panel(y_sdm,xo,W,N,T); 

fprintf(1,'true model is SDM \n');
fprintf('rho = %10.4f  \n',0.5);

in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out = [result2.lmarginal result2.probs];

mprint(out,in);


y_sdem = tmp + G\eterm;
% [y_sdems,xsdems,No,To,Wf2]=demeanF(y_sdem,xo,N,T,ted,W);


result3 = lmarginal_static_panel(y_sdem,xo,W,N,T); 

fprintf(1,'true model is SDEM \n');
fprintf('lambda = %10.4f  \n',0.4);

in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out = [result3.lmarginal result3.probs];

mprint(out,in);

