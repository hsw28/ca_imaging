function f = photobleach(ca_shapes, varargin)
  %nargin can be spike_times if you want to nromalize by number of spikes

firstmean = mean(ca_shapes(:,2250:4500)');
secondmean = mean(ca_shapes(:,end-4500:end-2250)');



if length(cell2mat(varargin))>1
  spike_times = cell2mat(varargin);
startend = cellspertime(spike_times);
firstmean = firstmean ./ startend(1,:);
secondmean = secondmean ./ startend(2,:);
end

f = [firstmean; secondmean]';

[h,p,ci,stats] = ttest(firstmean,secondmean, 'Tail','left');




%right means first is bigger

%waves0127_0607
%waves0128_0518

%left second is bigger
%spiketimes0127_0604
%spiketimes0128_0524
%spiketimes0330_0604
