% PURPOSE: An example of using sarp_g() Gibbs sampling
%          spatial autoregressive probit model
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sarp_gd (see also sarp_gd2 for a large data set)
%---------------------------------------------------


clear all;
% NOTE a large data set with 3107 observations
% from Pace and Barry, 
load elect.dat;             % load data on votes
latt = elect(1:400,5);
long = elect(1:400,6);
n = length(latt);


W = make_neighborsw(latt,long,12);


% generated data
defaultStream = RandStream.getDefaultStream;

IN = speye(n); 
rho = 0.8;  % true value of rho
sige = 1;
k = 3;

x = [ones(n,1) randn(n,k)];


beta(1,1) = 0.0;
beta(2,1) = 1.0;
beta(3,1) = -1.0;
beta(4,1) = 1.0;
theta = 0.5*beta(2:end,1);

y = (IN-rho*W)\(x*beta + W*x(:,2:end)*theta) + (IN-rho*W)\randn(n,1);

disp('maximum likelihood estimates based on continuous y');
result = sdm(y,x,W);
prt(result);

z = (y > 0);
z = ones(n,1).*z; % eliminate a logical vector

disp('# of zeros and ones');
[n-sum(z) sum(z)]

% Gibbs sampling function homoscedastic prior
% to maximum likelihood estimates
ndraw = 600;
nomit = 10;

prior2.nsteps = 1;
result2 = sdmp_g(z,x,W,ndraw,nomit,prior2);
prt(result2);


total_obs1 = (result2.total_obs(:,1));

total1 = mean(total_obs1)

tt=1:n;
plot(tt,total_obs1,tt,ones(1,n)*total1,'*');
legend('observation level total effect','mean total effect');
title('variable 1');
pause;

total_obs2 = (result2.total_obs(:,2));

total2 = mean(total_obs2)

tt=1:n;
plot(tt,total_obs2,tt,ones(1,n)*total2,'*');
legend('observation level total effect','mean total effect');
title('variable 2');
pause;

total_obs3 = mean((result2.total_obs(:,3)));

total3 = mean(total_obs3)

tt=1:n;
plot(tt,total_obs3,tt,ones(1,n)*total3,'*');
legend('observation level total effect','mean total effect');
title('variable 3');


