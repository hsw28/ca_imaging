function [cellsAll locationAll] = combineCAtraces(trace1, trace2, trace3, trace4, in1, in2, in3, in4)
  %inA = find(alignment_medium(:,14)>0 | alignment_medium(:,15)>0 );
  %inB = find(alignment_medium(:,17)>0 | alignment_medium(:,18)>0 );
  %inAB = intersect(inA, inB);

  %inA21 = alignment_medium(inAB,14);
  %inA22 = alignment_medium(inAB,15);
  %inB24 = alignment_medium(inAB,17);
  %inB25 = alignment_medium(inAB,18);

cells1 = NaN(length(in1), size(trace1,2));
cells2 = NaN(length(in1), size(trace2,2));
cells3 = NaN(length(in1), size(trace3,2));
cells4 = NaN(length(in1), size(trace4,2));

for i=1:length(in1)

  if (in1(i))>0

    cells1(i,:) =  trace1(in1(i),:);
  end
  if (in2(i))>0

    cells2(i,:) = trace2(in2(i),:);
  end
  if (in3(i))>0

    cells3(i,:) = trace3(in3(i),:);
  end
  if (in4(i))>0

    cells4(i,:) = trace4(in4(i),:);
  end

  tracemean = nanmean([cells1(i,:),cells2(i,:),cells3(i,:),cells4(i,:)]);
  if (in1(i))==0

    cells1(i,:) =  ones(1,size(trace1,2)).*tracemean;
  end
  if (in2(i))==0

    cells2(i,:) = ones(1,size(trace2,2)).*tracemean;
  end
  if (in3(i))==0

    cells3(i,:) = ones(1,size(trace3,2)).*tracemean;
  end
  if (in4(i))==0
    
    cells4(i,:) = ones(1,size(trace4,2)).*tracemean;
  end
end



location1 = ones(length(in1), size(trace1,2));
location2 = ones(length(in1), size(trace2,2));
location3 = ones(length(in1), size(trace3,2)).*2;
location4 = ones(length(in1), size(trace4,2)).*2;

%locationAll = [location2, location3];
%cellsAll = [cells2, cells3];
locationAll = [location1, location2, location3, location4];
cellsAll = [cells1, cells2, cells3, cells4];
