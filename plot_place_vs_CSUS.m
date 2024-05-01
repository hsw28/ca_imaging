function f = plot_place_vs_CSUS(peaks_time, pos, CSUS_id, velthreshold, dim)


  pos = smoothpos(pos);

  goodCSUS = find(CSUS_id(1,:)>0);

  good_CSUStime = pos(goodCSUS,1);
  good_CSUSpos = pos(goodCSUS,:);

  vel = ca_velocity(pos);
  goodvel = find(vel(1,:)>=velthreshold);
  goodtime = pos(goodvel, 1);
  goodpos = pos(goodvel,:);
  goodvel = setdiff(goodvel, goodCSUS);


  mintime = vel(2,1);
  maxtime = vel(2,end);
  tm = vel(2,:);

  highspeedspikes = [];
  CSUSspikes =[];


  for ii=1:length(peaks_time) %finding if in good vel
    [minValue_CSUS,closestIndex] = min(abs(peaks_time(ii)-good_CSUStime));
    [minValue_vel,closestIndex] = min(abs(peaks_time(ii)-goodtime));
    if minValue_CSUS <= 1/15 & isnan(peaks_time(ii))==0 %being CSUS takes precedence
      CSUSspikes(end+1)= peaks_time(ii);
    elseif minValue_vel <= 1/15 & isnan(peaks_time(ii))==0
      highspeedspikes(end+1) = peaks_time(ii);
    end
  end



rate = CA_normalizePosData(highspeedspikes, goodpos, dim, 1.000);
x = isnan(rate);
rate(x) = 0;
rate = imgaussfilt(rate,.75);

subplot(1,2,1)
[nr,nc] = size(rate);
imagesc(rate);
colormap('parula');
%lower and higher three percent of firing sets bounds
numrate = rate(~isnan(rate));
numrate = sort(numrate(:),'descend');
maxratefive = min(numrate(1:ceil(length(numrate)*0.01)));
numrate = sort(numrate(:),'ascend');
minratefive = max(numrate(1:ceil(length(numrate)*0.01)));

imagesc(rate, [maxratefive*.2, maxratefive*1.1]);
colorbar



rate = CA_normalizePosData(CSUSspikes, good_CSUSpos, dim, 1.000);
x = isnan(rate);
[a b] = find(rate>0);
%rate(a,b)=1;
rate(x) = 0;
%rate = imgaussfilt(rate,1);

subplot(1,2,2)
[nr,nc] = size(rate);
colormap('parula');
%lower and higher three percent of firing sets bounds
numrate = rate(~isnan(rate));
numrate = sort(numrate(:),'descend');
maxratefive = min(numrate(1:ceil(length(numrate)*0.01)));
numrate = sort(numrate(:),'ascend');
minratefive = max(numrate(1:ceil(length(numrate)*0.01)));
imagesc(rate, [minratefive*1.5, maxratefive*.7]);
colorbar
