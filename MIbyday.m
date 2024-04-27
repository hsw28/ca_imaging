function f = MIbyday(alignmentdata1, alignmentdata2, MI1, MI2, goodcells1, goodcells2)

%function f = fielddistancebyday(alignmentdata1, alignmentdata2, center1, center2, goodcells1, goodcells2)
% takes alignment data and finds the different in fields by individual day


if nargin < 5
    goodcells1 = [1:length(MI1)]; % Example: considering all cells as good
    goodcells2 = [1:length(MI2)];
end


MI1(~goodcells1) = NaN;
MI2(~goodcells2) = NaN;


both = find(alignmentdata1>0 & alignmentdata2>0);
want1 = (alignmentdata1(both));
want2 = (alignmentdata2(both));
MI1 = MI1(want1);
MI2 = MI2(want2);

distance = NaN(length(MI1),1);
for k=1:length(MI1)
  distance(k,1) = abs(MI1(k)-MI2(k));
end

f = distance;
