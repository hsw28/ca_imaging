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

info.lflag = 1;
info.eig = 1;


result1 = lmarginal_cross_section(y_slx,xo,W,info); 
% result1 = lmarginal_cross_section(y_slx,xo,W); 

fprintf(1,'true model is SLX \n');
fprintf(1,'time taken is: %16.4f seconds \n',result1.time);

in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out = [result1.lmarginal result1.probs];

mprint(out,in);


y_sdm = F\(tmp + eterm);

result2 = lmarginal_cross_section(y_sdm,xo,W,info); 
% result2 = lmarginal_cross_section(y_sdm,xo,W); 

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

result3 = lmarginal_cross_section(y_sdem,xo,W,info); 
% result3 = lmarginal_cross_section(y_sdem,xo,W); 

fprintf(1,'true model is SDEM \n');
fprintf(1,'time taken is: %16.4f seconds \n',result3.time);

fprintf('lambda = %10.4f  \n',0.4);

in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out = [result3.lmarginal result3.probs];

mprint(out,in);

% results for info.lflag = 1, info.eig = 1
% true model is SLX 
% time taken is:           0.3590 seconds 
% model log-marginal  model probs 
% slx     -1088.6072       0.3700 
% sdm     -1088.8236       0.2980 
% sdem    -1088.7156       0.3320 
% 
% true model is SDM 
% time taken is:           0.2820 seconds 
% rho =     0.5000  
% model log-marginal  model probs 
% slx     -1177.4566       0.0000 
% sdm     -1104.0876       0.9999 
% sdem    -1112.9241       0.0001 
% 
% true model is SDEM 
% time taken is:           0.2970 seconds 
% lambda =     0.4000  
% model log-marginal  model probs 
% slx     -1124.6364       0.0000 
% sdm     -1102.7602       0.0057 
% sdem    -1097.6032       0.9943 

% ===========================================================
% results for info.lflag = 1, info.eig = 0
% cross_section_demo2
% true model is SLX 
% time taken is:           0.1560 seconds 
% model log-marginal  model probs 
% slx     -1088.6072       0.4447 
% sdm     -1089.1336       0.2627 
% sdem    -1089.0255       0.2927 
% 
% true model is SDM 
% time taken is:           0.1250 seconds 
% rho =     0.5000  
% model log-marginal  model probs 
% slx     -1177.4566       0.0000 
% sdm     -1104.4638       0.9998 
% sdem    -1113.2265       0.0002 
% 
% true model is SDEM 
% time taken is:           0.1410 seconds 
% lambda =     0.4000  
% model log-marginal  model probs 
% slx     -1124.6364       0.0000 
% sdm     -1103.0695       0.0057 
% sdem    -1097.9087       0.9943 
% ===========================================================
% results for info.lflag = 0, info.eig = 1
% cross_section_demo2
% true model is SLX 
% time taken is:           8.2840 seconds 
% model log-marginal  model probs 
% slx     -1088.6072       0.3686 
% sdm     -1088.8197       0.2980 
% sdem    -1088.7077       0.3334 
% 
% true model is SDM 
% time taken is:           8.1870 seconds 
% rho =     0.5000  
% model log-marginal  model probs 
% slx     -1177.4566       0.0000 
% sdm     -1101.6389       0.9996 
% sdem    -1109.5476       0.0004 
% 
% true model is SDEM 
% time taken is:           8.1720 seconds 
% lambda =     0.4000  
% model log-marginal  model probs 
% slx     -1124.6364       0.0000 
% sdm     -1101.6033       0.0034 
% sdem    -1095.9084       0.9966 
% 

% ===========================================================
% results for defaults (no info input arguments)
% uses -1, 1 interval and lndetmc

% true model is SLX 
% time taken is:           0.1720 seconds 
% model log-marginal  model probs 
% slx     -1088.6072       0.4447 
% sdm     -1089.1336       0.2627 
% sdem    -1089.0255       0.2927 
% 
% true model is SDM 
% time taken is:           0.1720 seconds 
% rho =     0.5000  
% model log-marginal  model probs 
% slx     -1177.4566       0.0000 
% sdm     -1104.4638       0.9998 
% sdem    -1113.2265       0.0002 
% 
% true model is SDEM 
% time taken is:           0.1410 seconds 
% lambda =     0.4000  
% model log-marginal  model probs 
% slx     -1124.6364       0.0000 
% sdm     -1103.0695       0.0057 
% sdem    -1097.9087       0.9943 