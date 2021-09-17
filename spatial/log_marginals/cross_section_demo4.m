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

info.lflag = 0;
result5 = lmarginal_cross_section(pci_growth,xmatrix,W5,info);
result6 = lmarginal_cross_section(pci_growth,xmatrix,W6,info);
result7 = lmarginal_cross_section(pci_growth,xmatrix,W7,info);
result8 = lmarginal_cross_section(pci_growth,xmatrix,W8,info);
result9 = lmarginal_cross_section(pci_growth,xmatrix,W9,info);
result10 = lmarginal_cross_section(pci_growth,xmatrix,W10,info);
result11 = lmarginal_cross_section(pci_growth,xmatrix,W11,info);
result12 = lmarginal_cross_section(pci_growth,xmatrix,W12,info);
result13 = lmarginal_cross_section(pci_growth,xmatrix,W13,info);
result14 = lmarginal_cross_section(pci_growth,xmatrix,W14,info);
result15 = lmarginal_cross_section(pci_growth,xmatrix,W15,info);
result16 = lmarginal_cross_section(pci_growth,xmatrix,W16,info);

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

% default settings results (uses MC log det approximation)
% #neighbors              slx              sdm             sdem 
% 5                    0.0000           0.0000           0.0000 
% 6                    0.0000           0.0000           0.0000 
% 7                    0.0126           0.0000           0.0000 
% 8                    0.0168           0.0000           0.0000 
% 9                    0.0109           0.0000           0.0000 
% 10                   0.0050           0.0000           0.0000 
% 11                   0.0289           0.0000           0.0000 
% 12                   0.0501           0.0045           0.0000 
% 13                   0.0814           0.0198           0.0008 
% 14                   0.7646           0.8699           0.9342 
% 15                   0.0210           0.0055           0.0003 
% 16                   0.0086           0.1003           0.0646 

% info.lflag = 0 results (uses actual log det, NO approximation)
% #neighbors              slx              sdm             sdem 
% 5                    0.0000           0.0000           0.0000 
% 6                    0.0000           0.0000           0.0000 
% 7                    0.0126           0.0000           0.0000 
% 8                    0.0168           0.0000           0.0000 
% 9                    0.0109           0.0000           0.0000 
% 10                   0.0050           0.0000           0.0000 
% 11                   0.0289           0.0000           0.0000 
% 12                   0.0501           0.0099           0.0001 
% 13                   0.0814           0.0314           0.0037 
% 14                   0.7646           0.8731           0.9819 
% 15                   0.0210           0.0055           0.0008 
% 16                   0.0086           0.0800           0.0135 

% find best model, slx, sdm, sdem
adj2 = max(max(lmarginals));
madj2 = lmarginals-adj2;

% calculate probabilities across all models and weight matrices
xx2 = exp(madj2);

psum2 = sum(vec(xx2));
probs2 = xx2./psum2;

mprint(probs2,in);

% default settings results (uses MC log det approximation)
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

% info.lflag = 0 results (uses actual log det, NO approximation)
% #neighbors              slx              sdm             sdem 
% 5                    0.0000           0.0000           0.0000 
% 6                    0.0000           0.0000           0.0000 
% 7                    0.0000           0.0000           0.0000 
% 8                    0.0000           0.0000           0.0000 
% 9                    0.0000           0.0000           0.0000 
% 10                   0.0000           0.0000           0.0000 
% 11                   0.0000           0.0000           0.0000 
% 12                   0.0000           0.0000           0.0000 
% 13                   0.0000           0.0000           0.0008 
% 14                   0.0000           0.0004           0.9338 
% 15                   0.0000           0.0000           0.0003 
% 16                   0.0000           0.0000           0.0646 


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
% sigma^2        =    0.0137 
% Durbin-Watson  =    1.2492 
% Nobs, Nvars    =   3111,     5 
% ***************************************************************
% Variable           Coefficient      t-statistic    t-probability 
% constant              3.513147        38.514759         0.000000 
% log(pci0)            -0.302950       -17.824268         0.000000 
% log(pop0)             0.010417         4.255509         0.000021 
% W*log(pci0)          -0.038462        -1.800380         0.071898 
% W*log(pop0)           0.039568        11.393887         0.000000 


result2 = sdm(pci_growth,[ones(n,1) xmatrix],W);
vnames2 = strvcat('income growth','constant','log(pci0)','log(pop0)');
prt(result2,vnames2);

% sdm: negative variances from numerical hessian 
% sdm: t-statistics may be inaccurate 
% 
% Spatial Durbin model
% Dependent Variable =    income growth 
% R-squared          =    0.2667   
% Rbar-squared       =    0.2657   
% sigma^2            =    0.0090   
% log-likelihood     =        3893.4596  
% Nobs, Nvars        =   3111,     5 
% # iterations       =     19     
% min and max rho    =   -1.0000,   1.0000 
% total time in secs =    0.6620 
% time for lndet     =    0.1180 
% time for t-stats   =    0.0270 
% time for x-impacts =    0.3740 
% # draws used       =       1000  
% Pace and Barry, 1999 MC lndet approximation used 
% order for MC appr  =     50  
% iter  for MC appr  =     30  
% ***************************************************************
% Variable             Coefficient  Asymptot t-stat    z-probability 
% constant                0.783569       287.235654         0.000000 
% log(pci0)              -0.309608        -4.219294         0.000025 
% log(pop0)               0.009923         0.723524         0.469358 
% W-log(pci0)             0.232349       119.372040         0.000000 
% W-log(pop0)             0.005445         0.316987         0.751254 
% rho                     0.750994       821.179356         0.000000 
% 
% Direct               Coefficient           t-stat           t-prob         lower 01         upper 99 
% log(pci0)              -0.309916       -22.382194         0.000000        -0.344936        -0.271798 
% log(pop0)               0.011322         5.948718         0.000000         0.006660         0.015853 
% 
% Indirect             Coefficient           t-stat           t-prob         lower 01         upper 99 
% log(pci0)               0.000667         0.015464         0.987663        -0.106179         0.107929 
% log(pop0)               0.050523         6.192081         0.000000         0.029381         0.071623 
% 
% Total                Coefficient           t-stat           t-prob         lower 01         upper 99 
% log(pci0)              -0.309249        -7.509307         0.000000        -0.405279        -0.205298 
% log(pop0)               0.061845         7.579211         0.000000         0.040584         0.084341 
% 


result3 = sem(pci_growth,[ones(n,1) xmatrix W*xmatrix],W);
prt(result3,vnames);

% Spatial error Model Estimates 
% Dependent Variable =    income growth 
% R-squared       =    0.5247   
% Rbar-squared    =    0.5240   
% sigma^2         =    0.0090   
% log-likelihood  =        3898.4853  
% Nobs, Nvars     =   3111,     5 
% # iterations    =      0     
% min and max rho =   -0.9900,   0.9900 
% total time in secs =    0.2570 
% time for optimiz   =    0.1050 
% time for lndet     =    0.1130 
% time for t-stats   =    0.0140 
% Pace and Barry, 1999 MC lndet approximation used 
% order for MC appr  =     50  
% iter  for MC appr  =     30  
% ***************************************************************
% Variable           Coefficient  Asymptot t-stat    z-probability 
% constant              2.864025        33.902642         0.000000 
% log(pci0)            -0.313849       -23.909862         0.000000 
% log(pop0)             0.011825         6.168344         0.000000 
% W*log(pci0)           0.062145         7.553586         0.000000 
% W*log(pop0)           0.032862         5.611711         0.000000 
% lambda                0.769000       256.487945         0.000000 

