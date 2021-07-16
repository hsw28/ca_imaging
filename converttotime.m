function f = converttotime(train_peak_matrix, timestamps)
%converts output of getspikepeaks to matrix of times.



train_peak_matrix = train_peak_matrix.*7.5;

length(timestamps)

spike_val = NaN(size(train_peak_matrix));
size(spike_val)
[x, y] = find(train_peak_matrix>0);

for k=1:max(x)

  index = round(train_peak_matrix(k,:));
  index = index(~isnan(index));

  spike_val(k, 1:length(index)) = timestamps(index)./100;
end

f = spike_val;

%{
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

%}
