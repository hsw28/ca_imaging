% simulate model using W, then calculate log-marginals

% yt = rho*W*yt + xt*beta +W*xt*gamma + e

% DGP: 
% yt = inv[I_n*t - rho*(I_t \otimes W) x (X*beta + W*X*gamma + varepsilon)
% fixed effects used here


clear all;

load queenW.data;

Wc = normw(queenW);
  
[n,junk] = size(Wc);

N = n;

[a,b] = xlsread('states_borders.xls');
state_names = strvcat(b(2:end,1));

% creates a spatial weight matrix based on lengths of state borders in
% common with neighboring states
wraw = a(1:end,2:end);

n = 49;
wout = zeros(n,n);
for i=1:n;
    for j=1:n;
        wout(i,j) = wraw(i,j);
    end;
end;
 
for i=1:n;
    for j=1:n;
        wout(j,i) = wout(i,j);
    end;
end;
        
Wm = normw(wout); % normalizes the spatial weight matrix
W = Wm;

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

sige = 2;

alpha = 10;


% ============================================
kron2= sparse(kron(speye(T),W));

xmat = [ones(N*T,1) xo kron2*xo]; % model includes W*x-variables

beta_gamma = [alpha
    beta
    gamm];


    
    rho = 0.6;
    lam = 0.6;

 
F = kron(speye(N),speye(T)) - rho*kron(speye(T),W);
G = kron(speye(N),speye(T)) - lam*kron(speye(T),W);

kron2= sparse(kron(speye(T),W));

xmat = [ones(N*T,1) xo kron2*xo]; % model includes W*x-variables

beta_gamma = [alpha
    beta
    gamm];

eterm = randn(N*T,1)*sqrt(sige);

% generate fixed effects
tmp = eye(N);
FE = tmp(:,1:N-1);
bfe = unif_rnd(N-1,-1,1);

tmp = [xmat*beta_gamma + kron(ones(T,1),FE)*bfe];


y_slx = tmp + eterm;
ted = 0;
[y_slxs,xslx,No,To,Wf]=demeanF(y_slx,xo,N,T,ted,W);

info.iflag = 1;
% info.lflag = 0;
tic;
result1 = lmarginal_static_panel(y_slxs,xslx,Wf,No,To,info); 
toc;

[y_slxs,xslx,No,To,Wf]=demeanF(y_slx,xo,N,T,ted,Wc);

info.iflag = 1;
% info.lflag = 0;
result1b = lmarginal_static_panel(y_slxs,xslx,Wf,No,To,info); 

fprintf(1,'true model is SLX, Wmiles \n');
in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx Wm','sdm Wm','sdem Wm','slx Wc','sdm Wc','sdem Wc');
in.width = 10000;
in.fmt = '%10.4f';
out = [result1.lmarginal result1.probs 
       result1b.lmarginal result1b.probs];

mprint(out,in);

lmarginals = out(:,1);

% find best W-matrix
nmodels = length(lmarginals);
adj = max(lmarginals(:,1));
madj = matsub(lmarginals,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
in2.fmt = strvcat('%16.4f');
in2.rnames = strvcat('model','slx Wm','sdm Wm','sdem Wm','slx Wc','sdm Wc','sdem Wc');
mprint(probs,in2);



y_sdm = F\(tmp + eterm);
ted = 0;
[y_sdms,xsdms,No,To,Wf]=demeanF(y_sdm,xo,N,T,ted,W);

result2 = lmarginal_static_panel(y_sdms,xsdms,Wf,No,To,info); 

[y_sdms,xsdms,No,To,Wf]=demeanF(y_sdm,xo,N,T,ted,Wc);

result2b = lmarginal_static_panel(y_sdms,xsdms,Wf,No,To,info); 

fprintf(1,'true model is SDM \n');

in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx Wm','sdm Wm','sdem Wm','slx Wc','sdm Wc','sdem Wc');
in.width = 10000;
in.fmt = '%10.4f';
out = [result2.lmarginal result2.probs
       result2b.lmarginal result2b.probs];

mprint(out,in);

% find best W-matrix
lmarginals = out(:,1);

nmodels = length(lmarginals);
adj = max(lmarginals(:,1));
madj = matsub(lmarginals,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
in2.fmt = strvcat('%16.4f');
in2.rnames = strvcat('model','slx Wm','sdm Wm','sdem Wm','slx Wc','sdm Wc','sdem Wc');
mprint(probs,in2);


y_sdem = tmp + G\eterm;
[y_sdems,xsdems,No,To,Wf]=demeanF(y_sdem,xo,N,T,ted,W);


result3 = lmarginal_static_panel(y_sdems,xsdems,Wf,No,To,info); 

[y_sdems,xsdems,No,To,Wf]=demeanF(y_sdem,xo,N,T,ted,Wc);


result3b = lmarginal_static_panel(y_sdems,xsdems,Wf,No,To,info); 

fprintf(1,'true model is SDEM \n');

in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx Wm','sdm Wm','sdem Wm','slx Wc','sdm Wc','sdem Wc');
in.width = 10000;
in.fmt = '%10.4f';
out = [result3.lmarginal result3.probs
       result3b.lmarginal result3b.probs];

mprint(out,in);

% find best W-matrix
lmarginals = out(:,1);

nmodels = length(lmarginals);
adj = max(lmarginals(:,1));
madj = matsub(lmarginals,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
in2.fmt = strvcat('%16.4f');
in2.rnames = strvcat('model','slx Wm','sdm Wm','sdem Wm','slx Wc','sdm Wc','sdem Wc');
mprint(probs,in2);

