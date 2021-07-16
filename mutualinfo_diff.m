function f = mutualinfo_diff(MIcells, MIshuffle)
%finds place cells based on MI> shuffled MI

diff = MIcells-MIshuffle;
[x,y] = find(diff>0);
f = ([x,y]);

l = length(MIcells);
l2 = length(unique(f(:,1)))-1;
l2/l;

n1 = find(isnan(MIcells(:,1))==1);
n2 = find(isnan(MIcells(:,2))==1);
l = length(MIcells)-length(unique([n1;n2]));
l2 = length(unique(f(:,1)))-1;
l2/l
