function f = makespatialmatrix(spiketimes, spikecenters)
%makes 3d matrix of spike locations where 1 means spiking 0 means none
%matrix is organized as [x coord of field; y coord of field; frame]
%frames aree 7.5hz

%bin edgees
myedges = [min(spiketimes(:)):1/7.5:max(spiketimes(:))];

space = zeros(200,200,length(myedges));
for k = 1:size(spiketimes,1)
  currentcenterX = spikecenters(1,k); %this is for brightest, change 1 and 2 to 3 and 4 for center
  currentcenterY = spikecenters(2,k);
  [N,EDGES,BIN] = histcounts(spiketimes(k,:),myedges);
  BIN = BIN(find(BIN>0));
  space(currentcenterX,currentcenterY,BIN) =  1;
end

f = (space);
%f = ndSparse(space);
