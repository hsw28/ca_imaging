function [eigenvalues W D_degree] = isomap_eigen(Distance_matrix, k)
  %where D is the distance matrix
  % K =  k-nearest neighbors.
  % For example, for k-nearest neighbors with k = 5:

D = Distance_matrix;

% Step 1: Construct the Graph Laplacian matrix
n = size(D, 1); % Number of data points
W = zeros(n, n); % Initialize the adjacency matrix


for i = 1:n
    % Find the k-nearest neighbors for each data point
    [~, indices] = sort(D(i, :));
    W(i, indices(1:k)) = 1; % Assign 1 to the k-nearest neighbors
end

% Compute the degree matrix D_degree
D_degree = diag(sum(W, 2));

% Compute the unnormalized Laplacian matrix L
L = D_degree - W;

% Step 2: Compute the eigenvalues of the Laplacian matrix
eigenvalues = eig(L);
