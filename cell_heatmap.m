function anchorIdx = cell_heatmap(alignment_matrix)
%makes a heat map of shared cells across days, as well as anchor indices for each pair (a matrix of all cells shared in a day)

P = alignment_matrix;

% P is [nCells × nDays] logical
nDays  = size(P,2);
overlapMat = zeros(nDays);          % pair‑wise shared‑cell count

for d1 = 1:nDays
    for d2 = d1+1:nDays
        overlapMat(d1,d2) = nnz(P(:,d1) & P(:,d2));
    end
end

overlapMat = overlapMat + overlapMat';  % symmetric

imagesc(overlapMat); colorbar
title('Shared neurons between day pairs');
xlabel('Day'); ylabel('Day');


nDays = size(alignment_matrix,2)

anchorIdx = cell(nDays, nDays);   % each entry = vector of cell indices
for d1 = 1:nDays
    for d2 = d1+1:nDays
        anchorIdx{d1,d2} = find(P(:,d1) & P(:,d2));
    end
end
