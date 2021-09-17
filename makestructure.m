function f = makestructure(varargin)



for k=1:length(varargin)
  name = char(inputname(k));
  myStruct.(name) = cell2mat(varargin(k));

end

f = myStruct;
