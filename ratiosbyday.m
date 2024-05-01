function [ratios1 ratios2] = ratiosbyday(alignmentdata1, alignmentdata2, ratios1, ratios2)

%function f = fielddistancebyday(alignmentdata1, alignmentdata2, ratios1, ratios2)
% takes alignment data and finds the different in fields by individual day


both = find(alignmentdata1>0 & alignmentdata2>0);
want1 = (alignmentdata1(both));
want2 = (alignmentdata2(both));
ratios1 = ratios1(want1,:);
ratios2 = ratios2(want2,:);
