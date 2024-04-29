function f = manteltest(data_matrix)
  %can import data directly from outOfFieldFIring

% Data matrix where each row represents a cell
% Columns 1-2: TEBC max firing locations (x, y)
% Columns 3-4: Non-TEBC max firing locations (x, y)

data = data_matrix;

% Separating TEBC and Non-TEBC coordinates
tebc_coords = data(:, 1:2);
non_tebc_coords = data(:, 3:4);

n1 = ~isnan(tebc_coords(:,1));

tebc_coords = tebc_coords(n1,:);
non_tebc_coords = non_tebc_coords(n1,:);

n1 = ~isnan(non_tebc_coords(:,1));

tebc_coords = tebc_coords(n1,:);
non_tebc_coords = non_tebc_coords(n1,:);


% Computing distance matrices
dist_matrix_tebc = squareform(pdist(tebc_coords));
dist_matrix_nontebc = squareform(pdist(non_tebc_coords));

% Mantel test using Pearson correlation of distance matrices
n = numel(dist_matrix_tebc);
mantel_statistic = sum(sum(triu(dist_matrix_tebc, 1) .* triu(dist_matrix_nontebc, 1))) / n;
fprintf('Mantel statistic (correlation between distance matrices): %f\n', mantel_statistic);

% Permutation test to assess significance
%num_permutations = 10000;
num_permutations = 5000;

permutation_stats = zeros(num_permutations, 1);
for i = 1:num_permutations
    permuted_indices = randperm(size(non_tebc_coords, 1));
    permuted_dist_matrix = squareform(pdist(non_tebc_coords(permuted_indices, :)));
    permutation_stats(i) = sum(sum(triu(dist_matrix_tebc, 1) .* triu(permuted_dist_matrix, 1))) / n;
end

% Calculating p-value
p_value = mean(permutation_stats >= mantel_statistic);
fprintf('P-value from permutation test: %f\n', p_value);

% Assuming dist_matrix_tebc and dist_matrix_nontebc are already calculated
% Assuming mantel_statistic and permutation_stats are also calculated

% 1. Distance Matrix Heatmaps
figure;
subplot(1, 3, 1);
imagesc(dist_matrix_tebc);
title('TEBC Distance Matrix');
colorbar;
caxis([0 100]);
subplot(1, 3, 2);
imagesc(dist_matrix_nontebc);
title('Non-TEBC Distance Matrix');
colorbar;
caxis([0 100]);
% Calculate and Plot the Absolute Difference Heatmap
difference_matrix = abs(dist_matrix_tebc - dist_matrix_nontebc);
% Plot Composite Difference Heatmap
subplot(1, 3, 3);
imagesc(difference_matrix);
title('Absolute Difference Heatmap');
caxis([0 75]);
colorbar;

% 2. Scatter Plot of Distances
% Extracting the upper triangular parts of the matrices as vectors
% 2. Scatter Plot of Distances
% Extracting the upper triangular parts of the matrices as vectors
[row, col] = find(triu(ones(size(dist_matrix_tebc)), 1)); % Find indices of upper triangle, excluding diagonal
tebc_distances = dist_matrix_tebc(sub2ind(size(dist_matrix_tebc), row, col));
non_tebc_distances = dist_matrix_nontebc(sub2ind(size(dist_matrix_nontebc), row, col));

figure;
scatter(tebc_distances, non_tebc_distances, 'filled');
xlabel('TEBC Distances');
ylabel('Non-TEBC Distances');
title('Scatter Plot of TEBC vs. Non-TEBC Distances');
grid on;


% 3. Histogram of Permutation Test Results
figure;
histogram(permutation_stats, 'FaceColor', 'b');
hold on;
xline(mantel_statistic, 'LineWidth', 2, 'Color', 'r');
legend('Permutation Statistics', 'Observed Mantel Statistic');
title('Histogram of Permutation Test Results');
xlabel('Mantel Statistic');
ylabel('Frequency');
