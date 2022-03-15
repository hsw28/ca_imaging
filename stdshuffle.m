function f = stdshuffle(alignmentdata, repeatnum)
%takes alignment field data from [repeatnums fwd_field bwd_field] = cellsoverdays(alignmentdata, fieldcentersSTRUCTURE)
%shuffles field locations repeatnum of times and outputs bottom 95 and 99%

alignlinear = alignmentdata(:);
numsindex = find((isnan(alignlinear))==0);
numsamount = length(numsindex);
nums = alignlinear(numsindex);

allstd = NaN(1, repeatnum);
for k = 1:repeatnum
  shuff = randsample(numsamount, numsamount);
  shuff = nums(shuff);
  alignlinear(numsindex) = shuff;
  newmat = reshape(alignlinear, size(alignmentdata));
  newstd = nanstd(newmat');
  newav = nanmean(newstd);
  allstd(:,k) = newav;
end


allstd = sort(allstd');
allstd = allstd';
s95 = ceil(repeatnum*.05);
s99 = ceil(repeatnum*.01);
s995 = ceil(repeatnum*.005);
s999 = ceil(repeatnum*.001);

s95 = allstd(s95);
s99 = allstd(s99);
s995 = allstd(s995);
s999 = allstd(s999);

f = [s95, s99, s995, s999];
