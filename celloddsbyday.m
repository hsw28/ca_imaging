function f = celloddsbyday(alignmentdata, daysrecorded)
  %gives days apart and odds you will see a cell
  %days recorded should be a vector like [1, 2, 10, 14, 15, 16]
  %outputs [cell repeaitng, total cells, days apart]

  alignlinear = alignmentdata(:);
  numsindex = find(alignlinear==0);
  alignlinear(numsindex) = NaN; %changing to be NaN if no cell
  alignmentdata = reshape(alignlinear, size(alignmentdata));


output = [];
  for k=1:length(daysrecorded)-1
    for z=k+1:length(daysrecorded)
      inboth = 0;
      allcount = 0;
      for q=1:length(alignmentdata)
        if (alignmentdata(q,k)>0 & alignmentdata(q,z)>0) %cell in both
          inboth = inboth+1;
          allcount = allcount+1;
        elseif (alignmentdata(q,k)>0 & isnan(alignmentdata(q,z))==1) %cell in first day but not second
          allcount = allcount+1;
        end
      end
      daydist = daysrecorded(z) - daysrecorded(k);
      new = [inboth; allcount; daydist];
      output = [output, new];
    end
  end

  f = output';
