% simulate model using W and Wc, then calculate log-marginals
% as a test of the 3-dimensional numerical integration

% yt = int + mu + rho*W*yt + phi*ylag + theta*W*ylag + xt*beta +W*xt*gamma + mu + e

% DGP: 
% yt = inv[I_n*t - rho*(I_t \otimes W) - phi*(L \otimes I_n) - theta*(L \otimes W)]
%      x (X*beta + W*X*gamma + I_n*t*alpha + (I_t \otimes iota_n*mu + varepsilon)

% X = [tax, smokers, price, income/pop]
% everything is logged
% W is based on miles of common borders between states
% or Wc is based on 1st-order contiguity

clear all;

load queenW.data;

W = normw(queenW);

n = 49;

T=21; % number of time periods
N=49; % number of regions

Wbig = kron(speye(T-1),W);

rng(86573); % set the random number generator

sigx = 1;
x = randn(N,T)*sqrt(sigx);
taxm = x;
taxt = x(:,2:end);
tax_lagt = x(:,1:end-1);

x = randn(N,T)*sqrt(sigx);
pricem = x;
pricet = x(:,2:end);
price_lagt = x(:,1:end-1);

x = randn(N,T)*sqrt(sigx);
smokersm = x;
smokerst = x(:,2:end);
smokerst_lagt = x(:,1:end-1);

% simulate yt

x1 = (vec(taxt));

x2 = (vec(pricet));

x3 = (vec(smokerst));

k = 3;

x = [x1 x2 x3];

tmp = [-0.5 -1  1];
beta = tmp';

% tmp = [-0.75 -0.75 0.75]; % make sure these don't equal -rho*beta
% gamm = tmp';
% gamm = zeros(3,1);

sige = 1;
sigm = 1;
mu = randn(n,1)*sqrt(sigm); 
int = 1;

iotat = ones(T-1,1);

rho = 0.4;
phi = 0.8;
theta = -0.6; % note theta ~= -rho*phi
 
% we use T-1 to feed the lag for y(t-1)
 e = ones(T-1,1);
 tmp = spdiags([e e zeros(T-1,1)], -1:1, T-1, T-1);
 M = speye(T-1);
 L = tmp - M;
 L(1,1) = (1 - sqrt(1-phi*phi))/phi;
 
IN = speye(N);
IT = speye(T);


eterm = randn(N*(T-1),1)*sqrt(sige);

xmat = [x kron(L,IN)*x kron(M,W)*x kron(L,W)*x  ones(N*(T-1),1)];

%      xmat = [x Wbig*x];
     beta_gamma = [beta
                   0.5*beta
                   0.25*beta
                   0.1*beta
                   int];

      
 F = kron(speye(N),speye(T-1)) - rho*kron(speye(T-1),W) - phi*kron(L,speye(N)) - theta*kron(L,W);
 

 y_sdm = F\(xmat*beta_gamma + ones(N*(T-1),1)*int + eterm);
 
 y_sdem = xmat*beta_gamma + ones(N*(T-1),1)*int + F\eterm;

info.lflag = 0; % use defaults (no approximation)
result1 = lmarginal_dynamic_panel(y_sdm,x,W,N,T,info); 

fprintf(1,'true model is sdmu \n');
in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','dlm','sdmr','sdmu');
in.width = 10000;
in.fmt = '%10.4f';
out = [result1.lmarginal result1.probs];

mprint(out,in);

% true model is sdmu 
% model log-marginal  model probs 
% slx     -3005.5901       0.0000 
% sdm     -3000.7816       0.0000 
% dlm     -2203.5586       0.0000 
% sdmr    -2103.3427       0.0000 
% sdmu    -2030.7219       1.0000 

info.lflag = 1; % use log determinant approximation
result2 = lmarginal_dynamic_panel(y_sdm,x,W,N,T,info); % use defaults (no approximation)

fprintf(1,'true model is sdmu \n');
in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','dlm','sdmr','sdmu');
in.width = 10000;
in.fmt = '%10.4f';
out = [result2.lmarginal result2.probs];

mprint(out,in);
% 
% true model is sdmu 
% model log-marginal  model probs 
% slx     -3005.5901       0.0000 
% sdm     -3001.1822       0.0000 
% dlm     -2203.5586       0.0000 
% sdmr    -2108.6695       0.0000 
% sdmu    -2033.1033       1.0000 

