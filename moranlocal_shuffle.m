function f = moranlocal_shuffle(cellcenter, fieldcenter)
%shuffles cell locations (NOT field locations) and recomputes local moran's I

%forward
shuffF = randperm(length(fieldcenter)); %returns a row vector containing a random permutation of the integers from 1 to n without repeating elements.
cellcenter = cellcenter(:,shuffF);


f = moranlocal(cellcenter, fieldcenter);
%f = moranlocal_bydist(cellcenter, fieldcenter);
