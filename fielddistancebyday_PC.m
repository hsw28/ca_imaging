function f = fielddistancebyday_PC(alignmentdata1, alignmentdata2, center1, center2, MI1, MI2)

%evaluates field distance by day if the cell in either day is a place cell

if size(MI1,2)>2
  MI1 = MI1(:,3);
  MI2 = MI2(:,3);
end


inboth = find(alignmentdata1>0 & alignmentdata2>0);

want1 = (alignmentdata1(inboth));
want2 = (alignmentdata2(inboth));

PC = [];
for i = 1:length(inboth)
    if MI1(want1(i))>=0.75 || MI2(want2(i))>=0.75
      PC(end+1) = i;
    end
end

want1 = want1(PC);
want2 = want2(PC);

center1 = center1(want1,:);
center2 = center2(want2,:);

distance = NaN(length(center1),1);
for k=1:length(center1)
  points = [center1(k,:); center2(k,:)];
  d = pdist(points, 'euclidean');
  distance(k,1) = d;
end

f = distance;
