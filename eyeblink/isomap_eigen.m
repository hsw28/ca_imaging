  function eigenvalues = isomap_eigen(edge_matrix, k)
    % Input:
    %   edge_matrix: Adjacency matrix representing the neighborhood graph
    %   k: Number of nearest neighbors
    %
    % Output:
    %   eigenvalues: Eigenvalues of the Laplacian matrix

  % Step 1: Compute the Laplacian matrix L from the adjacency matrix
  n = size(edge_matrix, 1);
  D = diag(sum(edge_matrix, 2)); % Degree matrix
  L = D - edge_matrix; % Laplacian matrix

  % Step 2: Compute the eigenvalues of the Laplacian matrix L
  eigenvalues = eig(L);
