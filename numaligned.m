function f = numaligned(alignment)

num = [];
for k=1:length(alignment)
  num(end+1) = length(find(alignment(k,:))>0);
end

mean(num);
std(num);
median(num);

f = num;


tot = 0;
more = 0;
num = [];
for k=1:length(alignment)
  if (alignment(k,1)>0)
      tot = tot+1;
    if length(find(alignment(k,2:6))>0)>0;
      more = more+1;
      num(end+1) = length(find(alignment(k,2:6))>0);
end
end
end

more/tot
mean(num)
std(num)
