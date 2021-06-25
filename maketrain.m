function f = maketrain(sparsepeaks)

spike_val = zeros(size(sparsepeaks));
[x, y] = find(sparsepeaks>0);
mat = [x, y];
spike_val(mat) = 1;

f = spike_val;
