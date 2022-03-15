function f = fielddistancebyday(alignmentdata)
% takes alignment data and finds the different in fields by individual day



dif = NaN;
for k=1:size(alignmentdata,2)-1
  for j=k+1
    difs = abs(alignmentdata(:, j:end)-alignmentdata(:,k:end-1));
    dif = [dif; difs(:)];
  end
end

f = dif*4;
nanmean(dif)*4
nanstd(dif)*4
nanmedian(dif)*4
