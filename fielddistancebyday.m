function f = fielddistancebyday(alignmentdata1, alignmentdata2, center1, center2)
% takes alignment data and finds the different in fields by individual day

both = find(alignmentdata1>0 & alignmentdata2>0);
want1 = (alignmentdata1(both));
want2 = (alignmentdata2(both));
center1 = center1(want1);
center2 = center2(want2);

distance = NaN(length(center1),1);
for k=1:length(center1)
  points = [center1(k,:); center2(k,:)];
  d = pdist(points, 'euclidean');
  distance(k,1) = d;
end

f = distance;
