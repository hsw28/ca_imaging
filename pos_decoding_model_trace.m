function model = pos_decoding_model_trace(pos, spike_clusters, tdecode, dim, velthreshold)
%decodes position and outputs decoded x, y, confidence(in percents), and time.
%dim is bin size in cm
%tdecode is decoding in seconds


if isa(spike_clusters,'double')==1
  for k = 1:size(spike_clusters,1)
    name = num2words(k);
    for z = 1:length(name)
      if strfind(name(z),'-')
        name(z) = ' ';
      end
    end
    name = name(find(~isspace(name)));
    curclus = spike_clusters(k,:);
    curclus = curclus(find(~isnan(curclus)));
    cluststruct.(name) = curclus;
  end
  clusters = cluststruct;
end


posData = pos;


timee = pos(:,1);


[cc indexmin] = min(abs(posData(1,1)-timee));
[cc indexmax] = min(abs(posData(end,1)-timee));
decodetimevector = timee(indexmin:indexmax);
if length(decodetimevector)<10
  timevector = timee;
else
  timevector = decodetimevector;
end

if length(peaks_time)>length(pos)
  peaks_time = peaks_time(1:length(pos));
elseif length(peaks_time)<length(pos)
  pos = pos(1:length(peaks_time),:);
end

pos = smoothpos(posData);

tdecodesec = tdecode;

%BIN
psize = 1.000 * dim;

xvals = posData(:,2);
yvals = posData(:,3);
xmin = min(posData(:,2));
ymin = min(posData(:,3));
xmax = max(posData(:,2));
ymax = max(posData(:,3));


xbins = ceil((xmax-xmin)/psize); %number of x
ybins = ceil((ymax-ymin)/psize); %number of y
if ybins ==0
  ybins = 1;
end


xinc = xmin +(0:xbins)*psize; %makes a vectors of all the x values at each increment
yinc = ymin +(0:ybins)*psize; %makes a vector of all the y values at each increment
pos_samp_per_sec = length(posData(:,1))./(posData(end,1)-posData(1,1));

------

peaks_time = spike_clusters;
vel = ca_velocity(pos);
times = vel(2,:);
velocities = vel(1,:);
timeThreshold = 1/15; % second
highVelIndices = find(velocities >= velthreshold);
lowVelIndices = find(velocities < velthreshold);

% Filter out high velocity indices that are too close to low velocities
validHighVelIndices = [];
for ii = 1:length(highVelIndices)
    highVelTime = times(highVelIndices(ii));
    % Find the closest low velocity time
    [~, closestLowVelIndex] = min(abs(highVelTime - times(lowVelIndices)));
    closestLowVelTime = times(lowVelIndices(closestLowVelIndex));

    % Check if the high velocity time is more than 1 second away from the closest low velocity time
    if abs(highVelTime - closestLowVelTime) > timeThreshold
        validHighVelIndices = [validHighVelIndices, highVelIndices(ii)];
    end
end
goodpos = pos(validHighVelIndices,:);
all_highspeedspikes = peaks_time(:,validHighVelIndices);


numunits = size(peaks_time,1);
for k=1:length(names)
  curname = char(names(k));
  now = fxmatrix.(curname);
  highspeedspikes = all_highspeedspikes(k,:);

      if length(highspeedspikes)>0
        [trace_mean occprob] = CA_normalizePosData_trace(highspeedspikes, goodpos, dim, 1.000);

          if (size(trace_mean,1)) < (size(trace_mean,2))
            trace_mean = trace_mean';
          end

          trace_mean(isnan(trace_mean)) = eps;
          fxmatrix.(curname) = trace_mean;

      else
        fxmatrix.(curname) = NaN;
        fprintf('not enough high speed spikes')
      end
end




model = struct('fxmatrix', fxmatrix, 'tdecode', tdecode, 'dim', dim, 'velthreshold', velthreshold, 'xmin', xmin, 'xmax', xmax, 'ymin', ymin, 'ymax', ymax, 'xbins', xbins, 'ybins', ybins)
