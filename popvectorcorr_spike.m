function f = popvectorcorr_spike(spikesA, spikesB, alignmentdataA, alignmentdataB, binSize_in_seconds)
% Population vector correlation using spiking
% Can be normalized or not

binSize = binSize_in_seconds;

% Filter neurons present in both datasets
both = find(alignmentdataA > 0 & alignmentdataB > 0);
want1 = alignmentdataA(both);
want2 = alignmentdataB(both);
spikeMatrixA = spikesA(want1, :);
spikeMatrixB = spikesB(want2, :);

% Convert spike times to firing rates
firingRatesA = convertToRates(spikeMatrixA, binSize);
firingRatesB = convertToRates(spikeMatrixB, binSize);

% Normalize firing rates if desired
firingRatesA = zscore(firingRatesA, 0, 2);
firingRatesB = zscore(firingRatesB, 0, 2);

% Calculate the average firing rate of each neuron
meanRateA = nanmean(firingRatesA, 2);
meanRateB = nanmean(firingRatesB, 2);

% Calculate the population vector correlation between the mean rates
populationVectorCorrelation = corr(meanRateA, meanRateB, 'Rows', 'complete');

% Display the population vector correlation
disp(['Population Vector Correlation between Environment A and B: ', num2str(populationVectorCorrelation)]);

% Visualization of Average Responses
figure;
subplot(1, 2, 1);
bar(meanRateA);
title('Average Firing Rate in Environment A');
xlabel('Neuron Index');
ylabel('Firing Rate (Hz)');

subplot(1, 2, 2);
bar(meanRateB);
title('Average Firing Rate in Environment B');
xlabel('Neuron Index');
ylabel('Firing Rate (Hz)');



% Additional visualization to plot linear regression
validIndices = ~isnan(meanRateA) & ~isnan(meanRateB);
cleanResponseA = meanRateA(validIndices);
cleanResponseB = meanRateB(validIndices);
cleanResponseA = cleanResponseA*100000000;
cleanResponseB = cleanResponseB*100000000;

% Fit a linear model to the cleaned data
f = [cleanResponseA, cleanResponseB];
linearModel = fitlm(cleanResponseA, cleanResponseB);

% Get R-squared and p-value
R2 = linearModel.Rsquared.Ordinary;
pValue = linearModel.Coefficients.pValue(2);  % p-value for the slope

figure;

scatter(cleanResponseA, cleanResponseB, 'filled');
xlabel('Average Firing Rate in Environment A');
ylabel('Average Firing Rate in Environment B');
title('Scatter Plot of Average Responses');
axis square; % Makes the plot's aspect ratio 1:1
grid on; % Adds a grid for easier visualization
hold on;

% Plot the regression line
regressionLine = linearModel.predict([min(cleanResponseA), max(cleanResponseA)]');
plot([min(cleanResponseA), max(cleanResponseA)], regressionLine, 'k--');

% Set plot limits
xlim([min(cleanResponseA), max(cleanResponseA)]);
ylim([min(cleanResponseB), max(cleanResponseB)]);

% Display R^2 and p-value on the plot
text(double(min(cleanResponseA)), double(max(cleanResponseB)), ...
    sprintf('R^2 = %.2f, p = %.3g', R2, pValue), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

hold off;

if var(cleanResponseA) == 0 || var(cleanResponseB) == 0
    disp('Error: One or both of the response variables have zero variance.');
end
end


function rates = convertToRates(spikeMatrix, binSize)
    % Convert spike times to firing rates
    duration = max(max(spikeMatrix)); % Find the latest spike time to define experiment duration
    numBins = ceil(duration / binSize);
    numNeurons = size(spikeMatrix, 1);
    rates = zeros(numNeurons, numBins);

    for i = 1:numNeurons
        % Count spikes in each bin for each neuron
        for j = 1:size(spikeMatrix, 2)
            if ~isnan(spikeMatrix(i, j))
                binIndex = ceil(spikeMatrix(i, j) / binSize);
                rates(i, binIndex) = rates(i, binIndex) + 1;
            end
        end
        % Convert counts to rates
        rates(i, :) = rates(i, :) / binSize;
    end
end
