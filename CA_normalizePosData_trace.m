function [trace_mean occprob] = CA_normalizePosData_trace(trace_data,posData,dim, varargin)
%TO DO: add figure saving

%%%%same as normalizePosData but finds trace average for location instead of spike rate
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

eventData = trace_data;

if size(eventData, 1)>size(eventData,2)
	eventData = eventData';
end

if size(varargin)==0
psize = 1 * dim; %some REAL ratio of pixels to cm -- 3.5 for wilson, 2.5 for disterhoft linear, 6.85 for eyeblink
else
psize = cell2mat(varargin) *dim;
end



%only find occupancy map if one hasn't been provided

	xmax = max(posData(:,2));
	xmin = min(posData(:,2));
	ymax = max(posData(:,3));
	ymin = min(posData(:,3));
	xbins = ceil((xmax)/psize);
	ybins = ceil((ymax)/psize);

	time = zeros(ybins,xbins);
	events = zeros(ybins,xbins);
	xstep = xmax/xbins;
	ystep = ymax/ybins;
	tstep = (posData(end,1)-posData(1,1))./length(posData);
	%tstep = 1/15;


%occupancy
	for i = 1:xbins
    	for j = 1:ybins
        	A1 = posData(:,2)>((i-1)*xstep) & posData(:,2)<=(i*xstep); %finds all rows that are in the current x axis bin
        	A2 = posData(:,3)>((j-1)*ystep) & posData(:,3)<=(j*ystep); %finds all rows that are in the current y axis bin
        	A = [A1 A2]; %merge results
        	B = sum(A,2); %find the rows that satisfy both previous conditions
					C = B > 1;
        	time(ybins+1-j,i) = sum(C); %set the matrix cell for that bin to the number of rows that satisfy both
					events(ybins+1-j,i) = nanmean(trace_data(C));

			end
		end



%occupancy = imgaussfilt(time*tstep);
%occprob = occupancy./nansum(occupancy);

filtWidth = 3;
filtSigma = .5;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

trace_mean = (events);
trace_mean = nanconv(trace_mean,imageFilter, 'edge', 'nanout');

occupancy = (time*tstep);
occupancy2 = occupancy;
occupancy2(occupancy == 0) = NaN;
occupancy2 = nanconv(occupancy2,imageFilter, 'edge', 'nanout');
occprob = occupancy2./nansum(occupancy2);


%occprob = occupancy./nansum(occupancy);








colorbar = NaN;

%heat map stuff
%figure

%{
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

colorbar = [minratefive*1.5, maxratefive*.85];

%axis([20 75, -5 50])
%set(gca,'visible','off')
%}
