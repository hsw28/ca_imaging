function [distance centers] = fielddistancebyday(alignmentdata1, alignmentdata2, center1, center2, goodcells1, goodcells2)

%function f = fielddistancebyday(alignmentdata1, alignmentdata2, center1, center2, goodcells1, goodcells2)
% takes alignment data and finds the different in fields by individual day


if nargin < 5
    goodcells1 = [1:length(center1)]; % Example: considering all cells as good
    goodcells2 = [1:length(center2)];
end

center1(~goodcells1) = NaN;
center1(~goodcells2) = NaN;


both = find(alignmentdata1>0 & alignmentdata2>0);
want1 = (alignmentdata1(both));
want2 = (alignmentdata2(both));
center1 = center1(want1,:);
center2 = center2(want2,:);

centers = [];
distance = NaN(length(center1),1);
for k=1:length(center1)

  points = [center1(k,:); center2(k,:)];
  points2 = [center1(k,:), center2(k,:)];
  centers = [centers;points2];
  d = pdist(points, 'euclidean');
  distance(k,1) = d;
end
