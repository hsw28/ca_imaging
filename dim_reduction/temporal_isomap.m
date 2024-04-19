function f = temporal_isomap(data_matrix, k_size, dimensions)
%DEPRICATED


% Assuming 'time_series_data' is your preprocessed data matrix
% Make sure each row represents a time series

% Compute pairwise Euclidean distances
distance_matrix = pdist(data_matrix);

% Create the neighborhood graph by connecting each point to its k nearest neighbors
%Calculate pairwise distances
distances = pdist2(data_matrix, data_matrix);

% Initialize an adjacency matrix for the KNN graph
num_points = size(data_matrix, 1);
adjacency_matrix = zeros(num_points);

% Find the k-nearest neighbors for each data point
for i = 1:num_points
   [~, sorted_indices] = sort(distances(i, :));
   k_nearest_indices = sorted_indices(2:k_size+1);  % Exclude self, hence starting from 2
   adjacency_matrix(i, k_nearest_indices) = 1;
end


% Now, 'adjacency_matrix' represents the KNN graph

% Calculate the shortest path distances (geodesic distances) using Dijkstra's algorithm
geodesic_distances = graphallshortestpaths(sparse(adjacency_matrix));

f= geodesic_distances;
size(geodesic_distances)
options.dims = dimensions;
%[Y, R, E] = IsomapII(geodesic_distances, 'k', k_size, options);
