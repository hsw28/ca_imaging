function [matrix1_output, matrix2_output] = findoverlap(matrix1, matrix2, alignmentdata1, alignmentdata2)
%make sure for matrixes, the number of cells is first component

both = find(alignmentdata1>0 & alignmentdata2>0);
want1 = (alignmentdata1(both));
want2 = (alignmentdata2(both));
matrix1_output = matrix1(want1,:);
matrix2_output = matrix2(want2,:);
