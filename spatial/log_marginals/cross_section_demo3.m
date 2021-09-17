% PURPOSE: An example that compares:
% OLS, SLX, SDEM and SDM model estimates
% using Per Capita Income growth over the period 1987-1993
% for 3,111 US counties and different weight matrices
%---------------------------------------------------

clear all;

load uscounties.data;
% a matrix now exist named uscounties

% the matrix contains 11 columns of county-level data
% col 1  FIPS     
% col 2  LATITUDE 
% col 3  LONGITUDE
% col 4  POP1990  
% col 5  1987_PCI (per capita income)
% col 6  1988_PCI 
% col 7  1989_PCI 
% col 8  1990_PCI 
% col 9 1991_PCI 
% col 10 1992_PCI 
% col 11 1993_PCI 
    
[n,k] = size(uscounties); % find the size of the matrix

pci1987 = uscounties(:,5);  % extract the 5th column from the data matrix 
pci1993 = uscounties(:,11); % creates an n x 1 column vector 
pop1990 = uscounties(:,4);
% calculate growth rate of per capita income over the 1987 to 1993 period
pci_growth = log(pci1993) - log(pci1987);

% make these annualized growth rates
pci_growth = pci_growth;

% do a growth regression
% which involves regressing the growth rate on the (logged) initial level
xmatrix = [log(pci1987) log(pop1990)];

% run SDM model
latt = uscounties(:,2); % extract latt-long coordinates
long = uscounties(:,3);

W5 = make_neighborsw(latt,long,5);
W6 = make_neighborsw(latt,long,6);
W7 = make_neighborsw(latt,long,7);
W8 = make_neighborsw(latt,long,8);
W9 = make_neighborsw(latt,long,9);
W10 = make_neighborsw(latt,long,10);
W11 = make_neighborsw(latt,long,11);
W12 = make_neighborsw(latt,long,12);
W13 = make_neighborsw(latt,long,13);
W14 = make_neighborsw(latt,long,14);
W15 = make_neighborsw(latt,long,15);
W16 = make_neighborsw(latt,long,16);

result5 = lmarginal_cross_section(pci_growth,xmatrix,W5);
result6 = lmarginal_cross_section(pci_growth,xmatrix,W6);
result7 = lmarginal_cross_section(pci_growth,xmatrix,W7);
result8 = lmarginal_cross_section(pci_growth,xmatrix,W8);
result9 = lmarginal_cross_section(pci_growth,xmatrix,W9);
result10 = lmarginal_cross_section(pci_growth,xmatrix,W10);
result11 = lmarginal_cross_section(pci_growth,xmatrix,W11);
result12 = lmarginal_cross_section(pci_growth,xmatrix,W12);
result13 = lmarginal_cross_section(pci_growth,xmatrix,W13);
result14 = lmarginal_cross_section(pci_growth,xmatrix,W14);
result15 = lmarginal_cross_section(pci_growth,xmatrix,W15);
result16 = lmarginal_cross_section(pci_growth,xmatrix,W16);

% extract log-marginals and put them into an 8 by 3 matrix
% where each column is a different model (column 1 = slx, column 2 = sdm column 3 = sdem)

lmarginals = [result5.lmarginal' 
              result6.lmarginal'
              result7.lmarginal'
              result8.lmarginal'
              result9.lmarginal'
              result10.lmarginal'
              result11.lmarginal'
              result12.lmarginal'
              result13.lmarginal'
              result14.lmarginal'
              result15.lmarginal'
              result16.lmarginal'];
             
% calculate model probabilities for each column (model)
nmodels = length(lmarginals);

adj = max(lmarginals);
madj = matsub(lmarginals,adj);

xx = exp(madj);

% compute posterior probabilities
psum = sum(xx);
probs = matdiv(xx,psum);

in.fmt = '%16.4f';

in.cnames = strvcat('slx','sdm','sdem');

nstring = [];
for i=5:16;
 nstring = strvcat(nstring,num2str(i));
end;

in.rnames = strvcat('#neighbors',nstring);

mprint(probs,in);

% #neighbors              slx              sdm             sdem 
% 5                    0.0000           0.0000           0.0000 
% 6                    0.0000           0.0000           0.0000 
% 7                    0.0126           0.0000           0.0000 
% 8                    0.0168           0.0000           0.0000 
% 9                    0.0109           0.0000           0.0000 
% 10                   0.0050           0.0000           0.0000 
% 11                   0.0289           0.0000           0.0000 
% 12                   0.0501           0.0075           0.0000 
% 13                   0.0814           0.0373           0.0024 
% 14                   0.7646           0.6755           0.9719 
% 15                   0.0210           0.0096           0.0005 
% 16                   0.0086           0.2700           0.0252 


% find best model, slx, sdm, sdem
adj2 = max(max(lmarginals));
madj2 = lmarginals-adj2;

% calculate probabilities across all models and weight matrices
xx2 = exp(madj2);

psum2 = sum(vec(xx2));
probs2 = xx2./psum2;

mprint(probs2,in);

% #neighbors              slx              sdm             sdem 
% 5                    0.0000           0.0000           0.0000 
% 6                    0.0000           0.0000           0.0000 
% 7                    0.0000           0.0000           0.0000 
% 8                    0.0000           0.0000           0.0000 
% 9                    0.0000           0.0000           0.0000 
% 10                   0.0000           0.0000           0.0000 
% 11                   0.0000           0.0000           0.0000 
% 12                   0.0000           0.0000           0.0000 
% 13                   0.0000           0.0000           0.0024 
% 14                   0.0000           0.0002           0.9716 
% 15                   0.0000           0.0000           0.0005 
% 16                   0.0000           0.0001           0.0252 



% run growth regressions based on best models
[best_prob,best_index] = max(probs);

W = make_neighborsw(latt,long,best_index(1,1)+4);

result1 = ols(pci_growth,[ones(n,1) xmatrix W*xmatrix]);
vnames = strvcat('income growth','constant','log(pci0)','log(pop0)','W*log(pci0)','W*log(pop0)');
prt(result1,vnames);

% Ordinary Least-squares Estimates 
% Dependent Variable =    income growth 
% R-squared      =    0.2752 
% Rbar-squared   =    0.2743 
% sigma^2        =    0.0003 
% Durbin-Watson  =    1.2492 
% Nobs, Nvars    =   3111,     5 
% ***************************************************************
% Variable           Coefficient      t-statistic    t-probability 
% constant              0.501878        38.514759         0.000000 
% log(pci0)            -0.043279       -17.824268         0.000000 
% log(pop0)             0.001488         4.255509         0.000021 
% W*log(pci0)          -0.005495        -1.800380         0.071898 
% W*log(pop0)           0.005653        11.393887         0.000000 


result2 = sdm(pci_growth,[ones(n,1) xmatrix],W);
vnames2 = strvcat('income growth','constant','log(pci0)','log(pop0)');
prt(result2,vnames2);

% sdm: negative variances from numerical hessian 
% sdm: t-statistics may be inaccurate 
% 
% Spatial Durbin model
% Dependent Variable =    income growth 
% R-squared          =    0.2658   
% Rbar-squared       =    0.2649   
% sigma^2            =    0.0002   
% log-likelihood     =        9948.8791  
% Nobs, Nvars        =   3111,     5 
% # iterations       =     13     
% min and max rho    =   -1.0000,   1.0000 
% total time in secs =    0.7020 
% time for lndet     =    0.1170 
% time for t-stats   =    0.0270 
% time for x-impacts =    0.3700 
% # draws used       =       1000  
% Pace and Barry, 1999 MC lndet approximation used 
% order for MC appr  =     50  
% iter  for MC appr  =     30  
% ***************************************************************
% Variable             Coefficient  Asymptot t-stat    z-probability 
% constant                0.107798        73.615012         0.000000 
% log(pci0)              -0.044240        -4.194162         0.000027 
% log(pop0)               0.001417         0.717889         0.472826 
% W-log(pci0)             0.033603       118.424779         0.000000 
% W-log(pop0)             0.000726         0.292757         0.769708 
% rho                     0.758968      1942.580965         0.000000 
% 
% Direct               Coefficient           t-stat           t-prob         lower 01         upper 99 
% log(pci0)              -0.044248       -23.297646         0.000000        -0.049022        -0.039137 
% log(pop0)               0.001606         5.791287         0.000000         0.000819         0.002281 
% 
% Indirect             Coefficient           t-stat           t-prob         lower 01         upper 99 
% log(pci0)              -0.000026        -0.004208         0.996643        -0.014774         0.016114 
% log(pop0)               0.007372         5.934357         0.000000         0.004196         0.010937 
% 
% Total                Coefficient           t-stat           t-prob         lower 01         upper 99 
% log(pci0)              -0.044273        -7.272205         0.000000        -0.059351        -0.027873 
% log(pop0)               0.008978         7.225427         0.000000         0.006040         0.012553 


result3 = sem(pci_growth,[ones(n,1) xmatrix W*xmatrix],W);
prt(result3,vnames);

% Spatial error Model Estimates 
% Dependent Variable =    income growth 
% R-squared       =    0.5249   
% Rbar-squared    =    0.5243   
% sigma^2         =    0.0002   
% log-likelihood  =        9953.2734  
% Nobs, Nvars     =   3111,     5 
% # iterations    =      0     
% min and max rho =   -0.9900,   0.9900 
% total time in secs =    0.2580 
% time for optimiz   =    0.1040 
% time for lndet     =    0.1130 
% time for t-stats   =    0.0140 
% Pace and Barry, 1999 MC lndet approximation used 
% order for MC appr  =     50  
% iter  for MC appr  =     30  
% ***************************************************************
% Variable           Coefficient  Asymptot t-stat    z-probability 
% constant              0.408188        34.785193         0.000000 
% log(pci0)            -0.044828       -23.847985         0.000000 
% log(pop0)             0.001689         6.146759         0.000000 
% W*log(pci0)           0.009005         8.916597         0.000000 
% W*log(pop0)           0.004685         5.623976         0.000000 
% lambda                0.771000        82.151934         0.000000 

