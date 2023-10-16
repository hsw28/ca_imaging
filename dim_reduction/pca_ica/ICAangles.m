function angles_deg = ICAangles(W1, T1, W2, T2)

% Given outputs from ICA for two datasets
%W1, T1, Zica1, mu1; % from the first dataset
%W2, T2, Zica2, mu2; % from the second dataset

% 1. Obtain ICA basis matrices
B1 = W1 * T1;
B2 = W2 * T2;

% 2. Normalize the basis of the subspaces
B1_normalized = orth(B1);
B2_normalized = orth(B2);

% 3. Compute the SVD of the product of the transposed normalized basis
[U, S, V] = svd(B1_normalized' * B2_normalized);

% Clip the singular values to [-1,1] to handle numerical inaccuracies
S_clipped = min(max(diag(S), -1), 1);

% 4. Obtain the principal angles
angles_rad = acos(S_clipped);
angles_deg = rad2deg(angles_rad);
