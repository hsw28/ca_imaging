function f = shuffle2columns(data)

sw = 0;
if size(data,2)>size(data,1)
data = data';
sw = 1;
end

sz = length(data);
shuf = randi([1 sz], sz);
f = data(shuf,:);

if sw ==1;
f = f';
end
