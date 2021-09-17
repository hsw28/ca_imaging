% PURPOSE: An example of using sarp_g() on a large data set   
%          Gibbs sampling spatial autoregressive probit model                         
%---------------------------------------------------
% USAGE: sarp_gd2 (see sarp_gd for a small data set)
%---------------------------------------------------

clear all;
% Monte Carlo example
n = 1000;             
latt = rand(n,1);
long = rand(n,1);
k = 4;
x = randn(n,k);
W = make_neighborsw(latt,long,18);
% vnames = strvcat('voters','const','educ','homeowners','income');

rho = 0.4;
beta = zeros(k,1);
beta(1:2,1) = -1;
beta(3:4,1) = 1;

y = (speye(n) - rho*W)\(x*beta) + (speye(n) - rho*W)\randn(n,1);
ysave = y;

ndraw = 1200; 
nomit = 200;

% prior.novi = 1;
% result = sar_g(ysave,x,W,ndraw,nomit,prior);
% prt(result,vnames);

prior2.nstep = 1;
y = (y > 0)*1.0; % convert to 0,1 y-values
result2 = sarp_g(y,x,W,ndraw,nomit,prior2);
prt(result2);


% plot densities for comparison
% [h1,f1,y1] = pltdens(result.bdraw(:,1));
[h2,f2,y2] = pltdens(result2.bdraw(:,1));
% [h3,f3,y3] = pltdens(result.bdraw(:,2));
[h4,f4,y4] = pltdens(result2.bdraw(:,2));
% [h5,f5,y5] = pltdens(result.bdraw(:,3));
[h6,f6,y6] = pltdens(result2.bdraw(:,3));

plot(y2,f2,'.g');
xlabel(['true b =' num2str(beta(1,1))]);
pause;
plot(y4,f4,'.g');
xlabel(['true b =' num2str(beta(2,1))]);
pause;
plot(y6,f6,'.g');
xlabel(['true b =' num2str(beta(3,1))]);
pause;


% [h5,f5,y5] = pltdens(result.pdraw);
[h6,f6,y6] = pltdens(result2.pdraw);

plot(y6,f6,'.g');
xlabel(['true rho =' num2str(rho)]);
pause;
