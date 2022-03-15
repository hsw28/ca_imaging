function f = controom(cellcenter, fieldcenter)
  %finnds a contiguity matrix in rook config


nans = isnan(fieldcenter(:,1));
fieldcenterf = fieldcenter(~nans,1);
cellcenter1 = cellcenter(1,~nans);
cellcenter2 = cellcenter(2,~nans);
cellcenterf = [cellcenter1;cellcenter2];

fieldmat = plotfieldmatrix(cellcenterf, ones(length(cellcenterf),1));

[ic,icd] = ixneighbors(fieldmat)

cont = zeros(length(cellcenterf).*length(cellcenterf));


for k=length(ic)
  cont(ic(k),icd(k))=1;
  size(cont)
end


f = sparse(cont);
