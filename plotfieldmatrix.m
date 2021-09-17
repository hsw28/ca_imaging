function f = plotfieldmatrix(cellcenters, characteristic)
%plots numerical location in the location of the cell into a matrix

fieldmatrix = NaN(200,200);
%fieldmatrix = fieldmatrix+eps;

for k=1:length(cellcenters(1,:))
  if isnan(characteristic(k))==0
 fieldmatrix(cellcenters(1,k), cellcenters(2,k)) = characteristic(k);
end
end

f = fieldmatrix;
