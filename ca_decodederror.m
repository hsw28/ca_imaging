function f = ca_decodederror(decoded, pos, decodedinterval)
%returns an error in cm for each decoded time
%decodedinterval is the length of decoding, for ex .5 for half a second

%pos = fixpos(pos);
%pos(:,3) = 0;


if size(decoded,1)>size(decoded,2)
  decoded = decoded';
end



pointstime = decoded(4,:);
X = decoded(1,:);
Y = decoded(2,:);




vel = ca_velocity(pos);
vel = vel(1,:);

alldiffmean = [];
movediffmean = [];
stilldiffmean = [];
velvalue = 0;
velvector = [];
alldiff = [];
movediff = [];
stilldiff = [];
numpoints = [];



realX = [];
realY = [];
realT = [];
predX = [];
predY = [];
predT = [];
realV = [];

decodedinterval = round(15*decodedinterval);

for i=1:length(decoded)

  curtime = pointstime(i)+(decodedinterval);

  %postimes1 = find(pos(:,1)>=(curtime-decodedinterval));
  %postimes2 = find(pos(:,1)<curtime);
  %postimes = intersect(postimes1,postimes2)



    postimes1 = (i-1)*(decodedinterval)+1;
    postimes2 = postimes1+(decodedinterval)-1;
    postimes = [postimes1:1:postimes2];




%  [c index] = min(abs(curtime-pos(:,1)));
%  decinB = round(min(decin, index-1));
%  decinC = round(min(decin, length(vel)-index-1));

  %if mean(vel(index-decinB:index+decinC))>velabove
  if isnan(X(i))==0 & isnan(Y(i))==0

    %pairs = [X(i),Y(i);(pos(index,2)),(pos(index,3))];

    pairs = [X(i),Y(i);nanmean(pos(postimes,2)),nanmean(pos(postimes,3))];
    diff = pdist(pairs,'euclidean');
  else
    diff = NaN;
  end


    alldiff(end+1) = diff;
    %numpoints(end+1) = c;
    realX(end+1) = nanmean(pos(postimes,2));
    realY(end+1) = nanmean(pos(postimes,3));
    realT(end+1) = nanmean(pos(postimes,1));
    %realT(end+1) = pointstime(i);
    predX(end+1) = X(i);
    predY(end+1) = Y(i);
    predT(end+1) = pointstime(i);




end

nanmean(alldiff)./1.000
nanmedian(alldiff)./1.000

size(predX);
size(realX);
size(realT);
f = [realT; predX; predY; realX; realY; alldiff]';
%f = [alldiff/2.5; realT];
