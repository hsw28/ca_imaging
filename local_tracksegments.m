function f = local_tracksegments(localI)
%compare reward end in direction against non reward end in same direction
%compare two reward ends from correct directions against two non reward in both Directions
%compare reward ends with centers

fi = localI(:,1);
bi = localI(:,3);
fpos = localI(:,2).*4;
bpos = localI(:,4).*4;
posall = [fpos;bpos];
Iall = [fi;bi];

[Ni, EDGESi,  BIN] = histcounts(posall);
figure
goodbin = find(BIN>0);
scatter(Iall(goodbin), Ni((BIN(goodbin))));
fitlm(Iall(goodbin), Ni((BIN(goodbin))));


avs = NaN(1,length(Ni));
for k = 1:length(Ni)
  want = find(BIN==k);
  if length(want)>50
  avs(k) = nanmean(Iall(want));
end
end

%figure
%scatter(avs, Ni)
%fitlm(avs, Ni)






div = 5;
sz = (max(fpos)-min(fpos))./div;

frew = find(fpos>(max(fpos)-sz));
fno = find(fpos<(min(fpos)+sz));
fmida = find(fpos>(min(fpos)+(sz*floor(.5*div))));
fmidb = find(fpos<(min(fpos)+(sz*ceil(.5*div))));
fmid = intersect(fmida, fmidb);

bno = find(bpos>(max(bpos)-sz));
brew = find(bpos<(min(bpos)+sz));
bmida = find(bpos>(min(bpos)+(sz*floor(.5*div))));
bmidb = find(bpos<(min(bpos)+(sz*ceil(.5*div))));
bmid = intersect(bmida, bmidb);


[a b] = ttest2(fi(frew), fi(fno)) %
[a b] = ttest2(bi(brew), bi(bno))
[a b] = ttest2([fi(frew); bi(brew)], [fi(fno); bi(bno)]) %
[a b] = ttest2([fi(frew); fi(fno)], fi(fmid)) %
[a b] = ttest2([fi(frew)], fi(fmid))
[a b] = ttest2([fi(fno)], fi(fmid)) %
[a b] = ttest2([bi(brew); bi(bno)], bi(bmid))
[a b] = ttest2([bi(brew)], bi(bmid))
[a b] = ttest2([bi(bno)], bi(bmid))
[a b] = ttest2([fi(frew); bi(brew)], [fi(fmid); bi(bmid)])
[a b] = ttest2([fi(fmid); bi(bmid)], [fi(fno); bi(bno)]) %

%ANOVA REW vs MID vs NO
nanmean([fi(frew); bi(brew)])
nanmean([fi(fmid); bi(bmid)])
nanmean([fi(fno); bi(bno)])
an = ([fi(frew); bi(brew); fi(fmid); bi(bmid); fi(fno); bi(bno)]);
an2 = an;
an2(1:length([fi(frew); bi(brew)])) = 1;
an2(length([fi(frew); bi(brew)])+1:length([fi(frew); bi(brew)])+length([bi(bno); fi(fmid)])) = 2;
an2(length(an)-length([fi(fno); bi(bno)]):end) = 3;
[a b c] = anovan(an, an2, 'display','on')




%frew_ = nanmean(fi(frew))
%fno_ = nanmean(fi(fno))
%brew_ = nanmean(bi(brew))
%bno_ = nanmean(bi(bno))
%rew_ = nanmean(([fi(frew); bi(brew)]))
%no_ = nanmean([fi(fno); bi(bno)])
