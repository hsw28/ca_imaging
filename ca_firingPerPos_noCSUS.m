function f = ca_firingPerPos_noCSUS(pos, clusters, dim, velthreshold, CSUS_id)
%returns firing per position. dim is number of centimeters for binning



spikenames = (fieldnames(clusters));
spikenum = length(spikenames);

psize = 1.000 * dim; %some REAL ratio of pixels to cm




% Find all time points where CSUS_id <=0
goodCSUS = find(CSUS_id(1,:) <= 0);

% Now find the full range of indices to keep
%get vel
vel = ca_velocity(pos);
vel_time = vel(2,:)';
vel_mag  = vel(1,:)';

% Interpolate CSUS labels to velocity timestamps
interp_CSUS = interp1(CSUS_id(2,:), CSUS_id(1,:), vel_time, 'nearest', 0);

[~, uniqueIdx] = unique(pos(:,1), 'stable');
posData = pos(uniqueIdx, :);

% Interpolate X and Y position to velocity timestamps too
x = interp1(posData(:,1), posData(:,2), vel_time, 'linear', NaN);
y = interp1(posData(:,1), posData(:,3), vel_time, 'linear', NaN);


% Now build the posDat used downstream
posData = [vel_time, x, y];
pos_samp_per_sec = length(vel_time)./(posData(end,1)-posData(1,1))


validIdx = (vel_mag >= velthreshold) & (interp_CSUS == 0) & ...
           ~isnan(x) & ~isnan(y);
posDataFast = [vel_time(validIdx), x(validIdx), y(validIdx)];

xmin = min(posDataFast(:,2));
ymin = min(posDataFast(:,3));
xmax = max(posDataFast(:,2));
ymax = max(posDataFast(:,3));

mintime = min(posDataFast(:,1));
maxtime = max(posDataFast(:,1));




xbins = ceil((xmax-xmin)/psize); %number of x
ybins = ceil((ymax-ymin)/psize); %number of y

  timecells = zeros(xbins, ybins); %used to be time
  events = zeros(xbins,ybins);
  xstep = xmax/xbins;
  ystep = ymax/ybins;
  tstep = 1/pos_samp_per_sec;


  xinc = xmin +(0:xbins)*psize; %makes a vectors of all the x values at each increment
  yinc = ymin +(0:ybins)*psize; %makes a vector of all the y values at each increment

%only uses data that is >15cm/s -- first smooths for length of bin


%defiding position
%if length(varargin)>1
%  timecells = cell2mat(varargin)
if 0==1
else
  for x = (1:xbins) %WANT TO PERMUTE THROUGH EACH SQUARE OF SPACE SKIPPING NON OCCUPIED SQUARES. SO EACH BIN SHOULD HAVE TWO COORDINATES
    for y = (1:ybins)
        if x<xbins & y<ybins
          inX = find(posDataFast(:,2)>=xinc(x) & posDataFast(:,2)<xinc(x+1));
          inY = find(posDataFast(:,3)>=yinc(y) & posDataFast(:,3)<yinc(y+1));
        elseif x==xbins & y<ybins
          inX = find(posDataFast(:,2)>=xinc(x));
          inY = find(posDataFast(:,3)>=yinc(y) & posDataFast(:,3)<yinc(y+1));
        elseif x<xbins & y==ybins
          inX = find(posDataFast(:,2)>=xinc(x) & posDataFast(:,2)<xinc(x+1));
          inY = find(posDataFast(:,3)>=yinc(y));
        elseif x==xbins & y==ybins
          inX = find(posDataFast(:,2)>=xinc(x));
          inY = find(posDataFast(:,3)>=yinc(y));
        end

        inboth = intersect(inX, inY);
        timecells(x, y) = length(inboth);
      end
    end
end



for k = 1:spikenum
    events(:)=NaN;
    spikename = char(spikenames(k));
    unit = clusters.(spikename);
    [m firstspike] = min(abs(unit-mintime));
    [m lastspike] = min(abs(unit-maxtime));
    unit = unit(firstspike:lastspike);

    if velthreshold>0
      assvel = assignvelOLD(unit, vel);
      fastspikeindex = find(assvel > velthreshold);
      fastspike = unit(fastspikeindex);
    else
      fastspikeindex = length(unit);
      fastspike = unit;
    end

    ls = placeevent(fastspike, posData); %outputs [event; xposvector; yposvector];
    ls = ls';

    if length(ls)==0
      %myStruct.(spikename) = zeros(xbins, ybins);
      myStruct.(spikename) = 0;
      continue
    end

    minX = min(ls(:,2));
    maxX = max(ls(:,2));
    minY = min(ls(:,3));
    maxY = max(ls(:,3));
    [minValue, minX] = min(abs(minX-xinc));
    minX = minX-1;
    [minValue, maxX] = min(abs(maxX-xinc));
    maxX = maxX+1;
    [minValue, minY] = min(abs(minY-yinc));
    minY = minY-1;
    [minValue, maxY] = min(abs(maxY-yinc));
    maxY = maxY+1;

    if minX<1
      minX = 1;
    end
    if maxX>(xbins)
      maxX = (xbins);
    end
    if minY<1
      minY = 1;
    end
    if maxY>(ybins)
      maxY = (ybins);
    end

    if length(fastspikeindex)>0
    %WILL NEED TO DO THIS FOR ALL CELLS
    for x = minX:maxX
        for y = minY:maxY

          if x<xbins & y<ybins
            inX = find(ls(:,2)>=xinc(x) & ls(:,2)<xinc(x+1));
            inY = find(ls(:,3)>=yinc(y) & ls(:,3)<yinc(y+1));
          elseif x==xbins & y<ybins
            inX = find(ls(:,2)>=xinc(x));
            inY = find(ls(:,3)>=yinc(y) & ls(:,3)<yinc(y+1));
          elseif x<xbins & y==ybins
            inX = find(ls(:,2)>=xinc(x) & ls(:,2)<xinc(x+1));
            inY = find(ls(:,3)>=yinc(y));
          elseif x==xbins & y==ybins
            inX = find(ls(:,2)>=xinc(x));
            inY = find(ls(:,3)>=yinc(y));
          end

          inboth = intersect(inX, inY);
          events(x,y) = length(inboth);
          if events(x,y)>0 & timecells(x,y)==0
            timecells(x,y) = 1;
          end
        end
    end


    rate = events./(timecells*tstep)+eps; %time*tstep is occupancy %want this for all cells
    myStruct.(spikename) = rate;
    else
    rate = zeros(xbins, ybins);
    warning('the cell doesnt have enough points')
    spikename
    end


end

fprintf('firing per complete')
f = myStruct;

%[row,col] = find(rate==Inf)
%for k=1:length(row)
%  [events(row(k),col(k)), timecells(row(k),col(k))];
%end
