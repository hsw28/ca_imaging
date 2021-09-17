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
long = uscounties(:,3); % from the dataset
lmarginals = []; % empty matrix for storing results
neighbors = [];
info.lflag = 0;
info.eig = 1;
tic;
for ii=5:16; % loop over varying nearest neighbors W-matrices
W = make_neighborsw(latt,long,ii);
neighbors = [neighbors
             ii];
result = lmarginal_cross_section(pci_growth,xmatrix,W,info);
% extract log-marginals and put them into a matrix
lmarginals = [lmarginals
              result.lmarginal']; % we transpose the 3 x 1 vector
              % each column is a different model (column 1 = slx, column 2 = sdm column 3 = sdem)
end;
% calculate model probabilities for each column (model)
nmodels = length(lmarginals);
adj = max(lmarginals);
madj = matsub(lmarginals,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [neighbors matdiv(xx,psum)];
in.fmt = strvcat('%10d','%16.4f','%16.4f','%16.4f');
in.cnames = strvcat('#neigbhors','slx','sdm','sdem');
mprint(probs,in);
toc;


% find best model, slx, sdm, sdem
adj2 = max(max(lmarginals));
madj2 = lmarginals-adj2;

% calculate probabilities across all models and weight matrices
xx2 = exp(madj2);

psum2 = sum(vec(xx2));
probs2 = [neighbors xx2./psum2];

mprint(probs2,in);

