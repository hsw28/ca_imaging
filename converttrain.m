function f = converttrain(train_peak_matrix, length_in_seconds)
%converts spike train to matrix of times. times are in frame count

spike_val = NaN(size(train_peak_matrix));
[x, y] = find(train_peak_matrix>0);
for k=1:max(x)
  temp = find(x==k);
  spike_val(k, 1:length(temp)) = y(temp);
end

spike_val = spike_val(:,~all(isnan(spike_val)));

frames = size(train_peak_matrix,2);
f = spike_val.*(length_in_seconds./frames); %converting to seconds
