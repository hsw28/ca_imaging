function [overallCorrelation, meanChange, medianChange] = env_corr_spikes(spikesA, spikesB, alignmentdataA, alignmentdataB, binSize_in_seconds)
% Gets cross-correlation of spikes across environments and compares them
% After calculating the cross-correlation matrices for both environments, it computes the absolute differences between these matrices and then calculates an overall Pearson correlation between the flattened matrices to evaluate how similar the co-firing patterns are across the two environments.

binSize = binSize_in_seconds;

% Filter neurons present in both datasets
both = find(alignmentdataA > 0 & alignmentdataB > 0);
want1 = alignmentdataA(both);
want2 = alignmentdataB(both);
spikeMatrixA = spikesA(want1,:);
spikeMatrixB = spikesB(want2,:);

% Calculate firing rates based on binned spike counts
firingRateA = calculateFiringRates(spikeMatrixA, binSize);
firingRateB = calculateFiringRates(spikeMatrixB, binSize);

% Normalize firing rates if variability is significant
firingRateA = zscore(firingRateA, [], 2);
firingRateB = zscore(firingRateB, [], 2);

n = size(firingRateA, 1);  % Number of cells

% Initialize matrices to store the peak cross-correlation coefficients
peakCorrA = NaN(n, n);
peakCorrB = NaN(n, n);

% Calculate peak cross-correlation for each pair of cells in Environment A
for i = 1:n
    for j = i+1:n  % To avoid redundant calculations and self-correlation
        [corrA, ~] = xcorr(firingRateA(i,:), firingRateA(j,:), 'coeff');
        peakCorrA(i, j) = max(abs(corrA));
        peakCorrA(j, i) = peakCorrA(i, j);  % Symmetric matrix
    end
end

% Calculate peak cross-correlation for each pair of cells in Environment B
for i = 1:n
    for j = i+1:n
        [corrB, ~] = xcorr(firingRateB(i,:), firingRateB(j,:), 'coeff');
        peakCorrB(i, j) = max(abs(corrB));
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
fprintf('Mean Change in Correlation: %f\n', meanChange);
fprintf('Median Change in Correlation: %f\n', medianChange);
end

function rates = calculateFiringRates(spikeMatrix, binSize)
    % Calculate firing rates based on spike times and a given bin size
    duration = max(max(spikeMatrix, [], 'omitnan'), [], 'omitnan'); % Find latest spike time
    numBins = ceil(duration / binSize);
    rates = zeros(size(spikeMatrix, 1), numBins);
    for i = 1:size(spikeMatrix, 1)
        for j = 1:size(spikeMatrix, 2)
            if ~isnan(spikeMatrix(i, j))
                binIndex = min(numBins, floor(spikeMatrix(i, j) / binSize) + 1);
                rates(i, binIndex) = rates(i, binIndex) + 1;
            end
        end
    end
    rates = rates / binSize; % Convert counts to rates (spikes per second)
end
