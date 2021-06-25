function f = converttotime(train_peak_matrix)
%converts output of getspikepeaks to matrix of times.

length_in_seconds = length(train_peak_matrix)./7.5;


spike_val = NaN(size(train_peak_matrix));
[x, y] = find(train_peak_matrix>0);
for k=1:max(x)
  temp = find(x==k);
  spike_val(k, 1:length(temp)) = y(temp);
end

spike_val = spike_val(:,~all(isnan(spike_val)));
f = spike_val;
frames = size(train_peak_matrix,2);
f = spike_val.*(length_in_seconds./frames); %converting to seconds
