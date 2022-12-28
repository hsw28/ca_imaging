function f = moran(cellcenter, fieldcenter, varargin);
  %varargin is weight matrix if you have it

  %  USAGE: result = moran(y,x,W)
  %  where: x = dependent variable vector
  %         y = independent variables matrix
  %         W = contiguity matrix (standardized or unstandardized)

if size(cellcenter,2)~=2
  cellcenter = cellcenter';
end

%{
bicenter = NaN(size(cellcenter,2),1);
for k=1:size(cellcenter,2)
  f1f = fieldcenter(k, 1); %field center forward
  f1b = fieldcenter(k, 2); %field center backward
  if isnan(f1f)==1 && isnan(f1b)==1
    bicenter(k) = NaN;
  elseif isnan(f1f)==1
    bicenter(k) = f1b;
  elseif isnan(f1b)==1
    bicenter(k) = f1f;
  elseif abs(f1f-f1b) <= 15 %same field
    bicenter(k) = mean(f1f, f1b);
  elseif abs(f1f-f1b) > 15
    bicenter(k) = NaN;
  end
end

nans = isnan(bicenter);
bicenter = bicenter(~nans);
cellcenter1 = cellcenter(1,~nans);
cellcenter2 = cellcenter(2,~nans);
cellcenter = [cellcenter1;cellcenter2];

y = cellcenter;
size(cellcenter);
x = bicenter;

[W1 W W3] = xy2cont(cellcenter(1,:), cellcenter(2,:));
W = normw(W);

disweights = dist(cellcenter);
notzero = find(disweights(:)>0);
disweights(notzero) = 1./disweights(notzero);
W = disweights;
%}

y = cellcenter;
x = fieldcenter;


if length(cell2mat(varargin))>1
  W = cell2mat(varargin);

  %notzero = find(W(:)>0);
  %disweights(notzero) = disweights(notzero).^3;
  %W(notzero) = 1./W(notzero);
  %disweights = normw(disweights);

  W = normw(W);
end



  i = y(1,:);
  j = y(2,:);

  n = length(x);
  S = sum(W(:));


  xbar = mean(x);

 %multiplier
   m = n./S;

  %denom
  denom = sum((x-mean(x)).^2);

  %num
  row = x(:)-mean(x);
  col = x(:)-mean(x);
  col = col';
  num = W.*(row*col);
  num = sum(num(:));





I = m*(num./denom);
  f = I;
