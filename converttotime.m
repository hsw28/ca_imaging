function [spike_val times_sec] = converttotime(train_peak_matrix, timestamps)
%converts output of getspikepeaks to matrix of times.
%can import CSV from timestamps using: timestamps = readtable('file.csv');

if isa(timestamps,'table')
  timestamps = table2array(timestamps);
  timestamps = timestamps(:,2);
end

if size(timestamps,2)==3
  timestamps = timestamps(:,2);
end

if timestamps(5)>2
times_sec = timestamps./1000;
end

[x, y] = find(train_peak_matrix> 0);

spike_val = NaN(size(train_peak_matrix));


for k=1:size(train_peak_matrix,1)
  currclus = train_peak_matrix(k,:);
  index = find(currclus>0);
  wantedstamps = (index)*2;
  wantedstamps = wantedstamps(find(wantedstamps<=length(timestamps))); %this is just checking to make sure your frames arent > timestamps
  times = (timestamps(wantedstamps))./1000;
  spike_val(k, 1:length(times)) = times;
end



[x,y] = (find(isnan(spike_val)==0));

spike_val = spike_val(:,1:max(y));
