function [rate totspikes totstime colorbar spikeprob occprob] = CA_normalizePosData(eventData,posData,dim, varargin)
%TO DO: add figure saving

%This function bins event data based on a user input bin size and
%normalizes based on total time spent in bin
%color maps with range based on highest and lowest three percent of firing
%Args:
%   eventData: A timeseries of cell firings (e.g. the output of abovetheta)
%   posData: The matrix of overall position data with columns [time,x,y]
%   dim: Bin size in cm (only square bins are supported)
%
%
%   ex: map = normalizePosData(lsevents, pos, 4)
%
%Output:
%   rate: A discritized matrix of cell events per second
%   heatmap: A heatmap of the rate matrix



if size(eventData, 1)>size(eventData,2)
	eventData = eventData';
end

%eventData = cutclosest(posData(1,1), posData(end,1), eventData, eventData);
ls = ca_placeevent(eventData,posData);
if size(ls,2)~=3
	ls = ls';
end

if size(varargin)==0
psize = 1 * dim; %some REAL ratio of pixels to cm -- 3.5 for wilson, 2.5 for disterhoft linear, 6.85 for eyeblink
else
psize = cell2mat(varargin) *dim;
end


%only find occupancy map if one hasn't been provided

% Calculate bin edges for position data
xmax = max(posData(:, 2));
xmin = min(posData(:, 2));
ymax = max(posData(:, 3));
ymin = min(posData(:, 3));
xbins = ceil((xmax - xmin) / psize);
ybins = ceil((ymax - ymin) / psize);
xEdges = linspace(xmin, xmax, xbins+1);
yEdges = linspace(ymin, ymax, ybins+1);

% Histogram counts for position data to create occupancy map
[occCounts, ~, ~] = histcounts2(posData(:, 2), posData(:, 3), xEdges, yEdges);
timePerBin = (posData(end, 1) - posData(1, 1)) / length(posData);

% Histogram counts for event data to create events map
[eventCounts, ~, ~] = histcounts2(ls(:, 2), ls(:, 3), xEdges, yEdges);  % Using mapped events data


% Calculate rate, spike probability, and occupancy probability
occupancy = occCounts * timePerBin;

filtWidth = 3;
filtSigma = .5;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);


occupancy(occupancy == 0) = NaN;  % Avoid division by zero to handle empty bins
occupancy2 = nanconv(occupancy,imageFilter, 'edge', 'nanout');

events = nanconv(eventCounts,imageFilter, 'edge', 'nanout');

rate = events./occupancy2;

occprob = occupancy2./nansum(occupancy2);
spikeprob = events./nansum(events);




%[x,y] = find(isinf(rate)==1);
rate(isinf(rate)) = NaN;



totspikes = sum(eventCounts(:));
totstime = sum(occupancy(:));
%rate = rate(:, 15:end-20);



%heat map stuff
%figure


rate = rate./dim;
[nr,nc] = size(rate);
colormap('parula');
%lower and higher three percent of firing sets bounds
numrate = rate(~isnan(rate));
numrate = sort(numrate(:),'descend');
maxratefive = min(numrate(1:ceil(length(numrate)*0.03)));
numrate = sort(numrate(:),'ascend');
minratefive = max(numrate(1:ceil(length(numrate)*0.03)));

%{
pcolor([rate nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
if minratefive ~= maxratefive
		set(gca, 'clim', [minratefive*1.5, maxratefive*.85]);
end
%}

colorbar = [minratefive*1.5, maxratefive*.85];

%axis([20 75, -5 50])
%set(gca,'visible','off')
