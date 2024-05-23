function [meanChange, tauDifference, tauChange] = kendalls_tau(spikesA, spikesB, alignmentdataA, alignmentdataB, binSize_in_seconds)

binSize = binSize_in_seconds;

    both = find(alignmentdataA > 0 & alignmentdataB > 0);
    want1 = alignmentdataA(both);
    want2 = alignmentdataB(both);
    spikesA = spikesA(want1, :);
    spikesB = spikesB(want2, :);

    % Convert spike times to binned spike counts
    binSize = 1; % Define bin size in seconds
    binnedA = spikeTimeToBins(spikesA, binSize);
    binnedB = spikeTimeToBins(spikesB, binSize);

    % Preallocate a matrix to store the Kendall's Tau values
    numCells = size(binnedA, 1);
    tauMatrixA = NaN(numCells, numCells);  % Each element will store the Tau value between two cells for day A
    tauMatrixB = NaN(numCells, numCells);  % Similarly for day B

    % Calculate Kendall's Tau for each pair of cells
    for i = 1:numCells
        for j = i+1:numCells  % Only compute for upper triangular part since the matrix is symmetric
            tauMatrixA(i, j) = corr(binnedA(i,:)', binnedA(j,:)', 'type', 'Kendall', 'Rows', 'pairwise');
            tauMatrixA(j, i) = tauMatrixA(i, j);  % Mirror the value to the lower triangular part
            tauMatrixB(i, j) = corr(binnedB(i,:)', binnedB(j,:)', 'type', 'Kendall', 'Rows', 'pairwise');
            tauMatrixB(j, i) = tauMatrixB(i, j);
        end
    end

    % Calculate the difference and absolute differences between the two correlation matrices
    tauDifference = tauMatrixA - tauMatrixB;
    tauChange = abs(tauDifference);
    meanChange = mean(tauChange, 'all', 'omitnan');

    % Display mean change
    disp(['Mean change in Kendall''s Tau across environments: ', num2str(meanChange)]);

    % Visualize the results
    figure;
    subplot(1,3,1);
    imagesc(tauMatrixA);
    title('Environment A');
    colorbar;

    subplot(1,3,2);
    imagesc(tauMatrixB);
    title('Environment B');
    colorbar;

    subplot(1,3,3);
    imagesc(tauChange);
    title('Absolute Change in Tau');
    colorbar;
end

function binnedData = spikeTimeToBins(spikeMatrix, binSize)
    duration = max(spikeMatrix, [], 'all', 'omitnan'); % Find latest spike time
    numBins = ceil(duration / binSize);
    binnedData = zeros(size(spikeMatrix, 1), numBins);
    for i = 1:size(spikeMatrix, 1)
        for spike = spikeMatrix(i, :)
            if ~isnan(spike)
                bin = min(numBins, floor(spike / binSize) + 1);
                binnedData(i, bin) = binnedData(i, bin) + 1;
            end
        end
    end
end
