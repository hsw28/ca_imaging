function eigenvalues = isomap_eigen(Distance_matrix, k)
  %where D is the distance matrix
  % K =  k-nearest neighbors.
  % For example, for k-nearest neighbors with k = 5:

Distance_matrix;

% Step 1: Construct the k-nearest neighbor graph
   n = size(Distance_matrix, 1);
   [sorted_distances, indices] = sort(Distance_matrix, 2); % Sort distances
   k_neighbors = indices(:, 2:(k+1)); % Get the indices of k-nearest neighbors

   % Step 2: Construct the adjacency matrix W
   W = zeros(n, n);
   for i = 1:n
       W(i, k_neighbors(i, :)) = 1; % Assign 1 to the k-nearest neighbors
   end

   % Step 3: Compute the Laplacian matrix L
   D = diag(sum(W, 2)); % Degree matrix
   L = D - W;

   % Step 4: Compute the eigenvalues of the Laplacian matrix L
   eigenvalues = eig(L);
