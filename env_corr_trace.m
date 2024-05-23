function [overallCorrelation, meanChange, medianChange] = env_corr_trace(traceA, traceB, alignmentdataA, alignmentdataB)
%gets cross correlation of spikes across environments and compares them
%after calculating the cross-correlation matrices for both environments, it computes the absolute differences between these matrices and then calculates an overall Pearson correlation between the flattened matrices to evaluate how similar the co-firing patterns are across the two environments.

%Highpass Filtering: This is used to remove slow drifts that are common in calcium imaging data.
%Normalization: zscore is applied to normalize the activity of each cell to have zero mean and unit variance, which helps to standardize the data and reduce the effect of amplitude variations between cells.
%Cross-Correlation: This part remains the same, but operates on the normalized traces. The xcorr function calculates the cross-correlation and finds the peak value indicating the highest degree of correlation at any lag.

both = find(alignmentdataA>0 & alignmentdataB>0);
want1 = (alignmentdataA(both));
want2 = (alignmentdataB(both));
activityMatrixA = traceA(want1,:);
activityMatrixB = traceB(want2,:);


% Define sampling rate
Fs = 7.5; % Sampling frequency in Hz

% Preprocess the calcium traces (e.g., filtering, detrending)
filteredA = highpass(activityMatrixA, 0.1, Fs);  % High-pass filter to remove slow drifts
filteredB = highpass(activityMatrixB, 0.1, Fs);

% Normalize the traces (optional, depending on data)
normalizedA = zscore(filteredA, [], 2);  % Normalize each cell's activity
normalizedB = zscore(filteredB, [], 2);

% Initialize matrices to store the peak cross-correlation coefficients
n = size(activityMatrixA, 1);  % Number of cells
peakCorrA = NaN(n, n);
peakCorrB = NaN(n, n);

% Calculate peak cross-correlation for each pair of cells in Environment A
for i = 1:n
    for j = i+1:n  % To avoid redundant calculations and self-correlation
        [corrA, lags] = xcorr(normalizedA(i,:), normalizedA(j,:), 'coeff');
        [~, idx] = max(abs(corrA));
        peakCorrA(i, j) = corrA(idx);
        peakCorrA(j, i) = peakCorrA(i, j);  % Symmetric matrix
    end
end

% Calculate peak cross-correlation for each pair of cells in Environment B
for i = 1:n
    for j = i+1:n
        [corrB, lags] = xcorr(normalizedB(i,:), normalizedB(j,:), 'coeff');
        [~, idx] = max(abs(corrB));
        peakCorrB(i, j) = corrB(idx);
        peakCorrB(j, i) = peakCorrB(i, j);  % Symmetric matrix
    end
end

% Compare correlation matrices from Environment A and B
correlationChanges = abs(peakCorrA - peakCorrB);

% Calculate the overall correlation between the matrices to assess stability
overallCorrelation = corr(peakCorrA(:), peakCorrB(:), 'type', 'Pearson', 'Rows', 'complete');

% Display results
fprintf('Overall correlation between the environments: %f\n', overallCorrelation);

% Visualize the peak cross-correlation matrices
figure;
subplot(1,3,1);
imagesc(peakCorrA);
title('Environment A');
xlabel('Neuron Index');
ylabel('Neuron Index');
colorbar;
colormap('jet'); % Color map can be adjusted to preference
axis square;

subplot(1,3,2);
imagesc(peakCorrB);
title('Environment B');
xlabel('Neuron Index');
ylabel('Neuron Index');
colorbar;
colormap('jet');
axis square;

subplot(1,3,3);
imagesc(correlationChanges);
title('Absolute Changes in Correlation');
xlabel('Neuron Index');
ylabel('Neuron Index');
colorbar;
colormap('jet');
axis square;

% Calculate and display mean and median changes
meanChange = mean(correlationChanges(:), 'omitnan');
medianChange = median(correlationChanges(:), 'omitnan');
fprintf('Mean Change in Correlation: %f\n', meanChange)
fprintf('Median Change in Correlation: %f\n', medianChange)

% Add this info to one of the figures
text(0.5, 0.5, sprintf('Mean Change: %f\nMedian Change: %f\nOverall Correlation: %f', meanChange, medianChange, overallCorrelation), ...
    'HorizontalAlignment', 'center');
