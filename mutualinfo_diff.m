function f = mutualinfo_diff(MIcells, MIshuffle)
%finds place cells based on MI> shuffled MI
% SOMEWHAT DEP, USE MI_findgood.mat

diff = MIcells-MIshuffle(:,1:2);
[x,y] = find(diff>0);
f = ([x,y]);

%percent all
l = length(MIcells);
un = unique(f(:,1));
l2 = length(un);
l2/l;
l;

%perceent firing rate >0.01

l = length(MIcells)-length(find(isnan(MIcells)==1));
un = unique(f(:,1));
l2 = length(un);
l2/l;
l;


length(unique(f(:,1)))
