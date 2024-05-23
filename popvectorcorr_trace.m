function f = popvectorcorr_trace(traceA, traceB, alignmentdataA, alignmentdataB)

%population vector correlation using calcium traces
%Filtering: The data are filtered using a high-pass filter to remove slow drifts, which can obscure true neuronal signals in calcium imaging data.
%Normalization: Each cell’s activity is normalized using z-score normalization, which standardizes the data, making it easier to compare across neurons.
%Mean Response Calculation: The average activity for each neuron across all time points is calculated. This simplifies the data from a time series for each neuron to a single average value per neuron per environment.
%Population Vector Correlation: The correlation between these average activities across the two environments is calculated to determine how similar the overall neural activity patterns are between the two conditions.
%Visualization: Bar graphs show the average activity for each neuron in each environment, and a scatter plot compares these directly, providing visual insights into the correlation and potential linear relationships.
%ex
% rat0222.popcorr_traceA1An = popvectorcorr_trace(rat0222.Ca_traces.CA_traces_2023_05_08, rat0222.Ca_traces.CA_traces_2023_05_09, rat0222.alignment(:,2), rat0222.alignment(:,3));


% Assuming activityMatrixA and activityMatrixB contain your raw calcium trace data
% for Environment A and B respectively
both = find(alignmentdataA>0 & alignmentdataB>0);
want1 = (alignmentdataA(both));
want2 = (alignmentdataB(both));
activityMatrixA = traceA(want1,:);
activityMatrixB = traceB(want2,:);

% Define sampling rate
Fs = 7.5; % Sampling frequency in Hz

% Preprocess the calcium traces (e.g., filtering, detrending)
%activityMatrixA = highpass(activityMatrixA, 0.1, Fs);  % High-pass filter to remove slow drifts
%activityMatrixB = highpass(activityMatrixB, 0.1, Fs);

% Normalize the traces (optional, depending on data)
%activityMatrixA = zscore(activityMatrixA, [], 2);  % Normalize each cell's activity
%activityMatrixB = zscore(activityMatrixB, [], 2);

% Calculate the average response of each neuron over time in each environment
meanResponseA = mean(activityMatrixA, 2);  % Average across columns (time points)
meanResponseB = mean(activityMatrixB, 2);

% Calculate the population vector correlation between the mean responses
populationVectorCorrelation = corr(meanResponseA, meanResponseB);

% Display the population vector correlation
disp(['Population Vector Correlation between Environment A and B: ', num2str(populationVectorCorrelation)]);

% Visualization of Average Responses
figure;
subplot(1, 2, 1);
bar(meanResponseA);
title('Average Neuronal Response in Environment A');
xlabel('Neuron Index');
ylabel('Average Response');

subplot(1, 2, 2);
bar(meanResponseB);
title('Average Neuronal Response in Environment B');
xlabel('Neuron Index');
ylabel('Average Response');

% Visualization of the correlation matrix if desired


% Remove NaN values from both meanResponseA and meanResponseB
validIndices = ~isnan(meanResponseA) & ~isnan(meanResponseB);
cleanResponseA = meanResponseA(validIndices);
cleanResponseB = meanResponseB(validIndices);
%check for constant values
if var(cleanResponseA) == 0 || var(cleanResponseB) == 0
    error('One or both of the variables have zero variance.');
end

% Fit a linear model to the cleaned data
linearModel = fitlm(cleanResponseA, cleanResponseB);

f = [cleanResponseA, cleanResponseB];

% Get R-squared and p-value
R2 = linearModel.Rsquared.Ordinary;
pValue = linearModel.Coefficients.pValue(2);  % p-value for the slope

% Visualization of the correlation matrix if desired
figure;
scatter(cleanResponseA, cleanResponseB, 'filled');
xlabel('Average Response in Environment A');
ylabel('Average Response in Environment B');
title('Scatter Plot of Average Responses');
axis square; % Makes the plot's aspect ratio 1:1
grid on; % Adds a grid for easier visualization
hold on;

% Plot the regression line
regressionLine = linearModel.predict([min(cleanResponseA) max(cleanResponseA)]');
plot([min(cleanResponseA), max(cleanResponseA)], regressionLine, 'k--');

% Set plot limits
xlim([min(cleanResponseA), max(cleanResponseA)]);
ylim([min(cleanResponseB), max(cleanResponseB)]);

% Display R^2 and p-value on the plot
text(double(min(cleanResponseA)), double(max(cleanResponseB)), ...
    sprintf('R^2 = %.2f, p = %.3g', R2, pValue), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Ensure everything is visible
hold off;
