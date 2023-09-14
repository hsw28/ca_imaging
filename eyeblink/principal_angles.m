function [degrees PA PB] = principal_angles(coeffA, coeffB, princ_to_analyze)
  %compute the principal angles between two m-dimensional manifolds A and B
  %finds pricipal components between coefficients
  %can specify which components to analyze, ie = principal_angles(coeff_trim22, coeff_trim24, [1:3])


%{
To compute the principal angles between two m-dimensional manifolds A and
B embedded in an n-dimensional neural space, we follow the method by Björck and
Golub48 : consider the corresponding m-dimensional bases Wa and Wa provided
by the m leading PC neural modes, construct their m by m inner product matrix,
and perform a singular value decomposition to obtain

Wa'Wb = PaCPb'

Here Wi, i = a, b are the n by m PC matrices that span the task-specific manifolds
A and B; the corresponding PC neural modes are their column vectors. The
matrices PA and PB, both of dimension m by m, define the new manifold directions
that successively minimize the principal angles. Note that these new projections are
specific to the pair of tasks being compared. The matrix C is a diagonal matrix
whose elements are the ranked cosines of the principal angles
%}

%PC matrices WA and WB for manifolds A and B
Wa = coeffA(:,princ_to_analyze);
Wb = coeffB(:,princ_to_analyze);

% Wa and Wb are n by m matrices, where m is the number of principal components

% Compute the inner product matrix
InnerProductMatrix = Wa' * Wb;

% Perform singular value decomposition
[U, S, V] = svd(InnerProductMatrix);

% Compute the principal angles
principalAngles = acos(diag(S));

% Matrices PA and PB that define the new manifold directions
PA = U;
PBp = V;
PB = V';

% Matrix C with ranked cosines of the principal angles
C = S;
principalAngles = acos(diag(S));
radians = principalAngles;

degrees = rad2deg(radians);
