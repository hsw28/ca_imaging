function UMAPtrial_outcome(animalName, win, plot3D, balanceClasses)
% one map for all days, trials compressed
    if nargin < 3, plot3D = false; end
    if nargin < 4, balanceClasses = false; end
    Fs = 7.5;
    animal = evalin('base', animalName);
    dateList = autoDateList(animal);
    G = animal.alignmentALL;
    alignment_data = G;
    nDays = numel(dateList);
    nAligned = size(G, 2);
    if nAligned < nDays
        fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
        dateList = dateList(1:nAligned);
        nDays = nAligned;
    end

    minDayFrac = 0.5;
    neuronDayCounts = sum(G > 0, 2);
    sharedNeurons = find(neuronDayCounts >= round(minDayFrac * nAligned));
    maxNeurons = numel(sharedNeurons);
    fprintf('Using %d shared neurons for embedding.\n', maxNeurons);

    nBins = round(diff(win) * Fs);
    labels = [];
    trialDays = {};
    trialVecs = {};

    for d = 1:nDays
        dateStr = dateList{d};
        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);

        if isempty(X)
            fprintf('Skipping %s: no usable data.\n', dateStr);
            continue;
        end

        fprintf('Day %s: %d trials (correct: %d, incorrect: %d)\n', dateStr, numel(y), sum(y==1), sum(y==0));

        sharedIdx = G(sharedNeurons, d);

        for t = 1:size(X,3)
            trialMat = nan(numel(sharedNeurons), nBins);
            for n = 1:numel(sharedNeurons)
                localIdx = sharedIdx(n);
                if localIdx > 0 && localIdx <= size(X,1)
                    trialMat(n,:) = X(localIdx,:,t);
                end
            end
            trialMat(isnan(trialMat)) = 0;
            vec = double(trialMat(:)');
            trialVecs{end+1,1} = vec;
            labels(end+1,1) = y(t);
            trialDays{end+1,1} = dateStr;
        end
    end

    if isempty(trialVecs)
        error('No consistent trials found.');
    end

    % === Step 1: Determine the length of each trial vector ===
    % trialVecs is a cell array where each cell holds a row vector (a flattened trial matrix)
    vecLens = cellfun(@length, trialVecs);  % returns an array with the length of each trial vector

    % === Step 2: Find the maximum vector length ===
    % This is the longest trial representation — all others will be padded to match this
    maxLen = max(vecLens);

    % === Step 3: Pad all trial vectors to the same length ===
    % This ensures cell2mat works, and all trials have equal-size feature vectors
    for i = 1:numel(trialVecs)
        if length(trialVecs{i}) < maxLen
            % Pad the current trial vector with zeros at the end so it matches maxLen
            trialVecs{i}(end+1:maxLen) = 0;
        end
    end

    % === Step 4: Convert the cell array to a 2D numeric matrix ===
    % Now that all rows are the same length, you can concatenate them vertically
    Xmat = cell2mat(trialVecs);  % size: [#trials × #features]

    % === Step 5: Z-score normalize each feature across all trials ===
    % This centers each column (feature) to mean 0, std 1
    Xmat = zscore(Xmat);

    % === Step 6: Assign labels and day labels to fixed variables ===
    Y = labels;          % trial labels (correct/incorrect, etc.)
    trialDays = trialDays;  % day each trial belongs to

    % === Step 7: Add tiny random noise to all values ===
    % Helps break ties or degenerate structure in downstream analysis (e.g., UMAP/SVM)
    Xmat = Xmat + 1e-6 * randn(size(Xmat));  % jitter each value slightly


    [embedding, umapParams] = run_umap(Xmat, 'n_components', 5, 'n_neighbors', 15, 'min_dist', 0.3, 'metric', 'cosine', 'randomize', true, 'verbose', true);

    if isequal(embedding, -1)
        error('UMAP failed to compute embedding');
    end

    fprintf('Label distribution: %d correct, %d incorrect\n', sum(Y==1), sum(Y==0));
    fprintf('Mean of Xmat: %.3f\n', mean(Xmat(:)));
    fprintf('Std of Xmat: %.3f\n', std(Xmat(:)));

    idx = ~isnan(Y);
    embedding = embedding(idx, :);
    Y = Y(idx);
    Xmat = Xmat(idx, :);
    trialDays = trialDays(idx);

    rng(42);
    part = cvpartition(Y, 'KFold', 5);
    accsXmat = zeros(part.NumTestSets, 1);
    accsEmbed = zeros(part.NumTestSets, 1);

    allTrue = [];
    allPredRaw = [];
    allPredEmbed = [];

    for i = 1:part.NumTestSets
        trainIdx = training(part, i);
        testIdx = test(part, i);

        opts = {'KernelFunction','rbf','Standardize',true};
        if balanceClasses
            opts = [opts, {'ClassNames',[0;1]}];
        end

        mdlRaw = fitcsvm(Xmat(trainIdx,:), Y(trainIdx), opts{:});
        predRaw = predict(mdlRaw, Xmat(testIdx,:));
        accsXmat(i) = sum(predRaw == Y(testIdx)) / numel(testIdx);

        mdlEmbed = fitcsvm(embedding(trainIdx,:), Y(trainIdx), opts{:});
        predEmbed = predict(mdlEmbed, embedding(testIdx,:));
        accsEmbed(i) = sum(predEmbed == Y(testIdx)) / numel(testIdx);

        allTrue = [allTrue; Y(testIdx)];
        allPredRaw = [allPredRaw; predRaw];
        allPredEmbed = [allPredEmbed; predEmbed];
    end

    confRaw = confusionmat(allTrue, allPredRaw);
    confEmbed = confusionmat(allTrue, allPredEmbed);

    fprintf('Linear SVM decoding accuracy in original Xmat space (mean over folds): %.2f%%\n', 100 * mean(accsXmat));
    fprintf('Linear SVM decoding accuracy in UMAP 3D space (mean over folds): %.2f%%\n', 100 * mean(accsEmbed));
    fprintf('UMAP matrix size: [%d trials x %d dims], Shared neurons: %d\n', size(Xmat,1), size(embedding,2), maxNeurons);

    disp('Confusion matrix for Raw model:');
    disp(array2table(confRaw, 'VariableNames', {'Pred_0','Pred_1'}, 'RowNames', {'True_0','True_1'}));

    disp('Confusion matrix for Embed model:');
    disp(array2table(confEmbed, 'VariableNames', {'Pred_0','Pred_1'}, 'RowNames', {'True_0','True_1'}));

    figure;
    hold on;
    scatter(embedding(Y==0,1), embedding(Y==0,2), 30, 'r', 'filled');
    scatter(embedding(Y==1,1), embedding(Y==1,2), 30, 'b', 'filled');
    legend('Incorrect', 'Correct');
    title(sprintf('Trial-Flat UMAP (2D) - %d trials, %d days, %d neurons', size(Xmat,1), numel(unique(trialDays)), maxNeurons));
    xlabel('UMAP 1'); ylabel('UMAP 2'); grid on;

    if plot3D & size(embedding, 2) >= 3
        figure;
        scatter3(embedding(Y==0,1), embedding(Y==0,2), embedding(Y==0,3), 30, 'r', 'filled');
        hold on;
        scatter3(embedding(Y==1,1), embedding(Y==1,2), embedding(Y==1,3), 30, 'b', 'filled');
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3');
        legend('Incorrect', 'Correct');
        title(sprintf('Trial-Flat UMAP (3D) - %d trials, %d days, %d neurons', size(Xmat,1), numel(unique(trialDays)), maxNeurons));
        grid on;
    end

    figure;
    dayLabels = categorical(trialDays);
    uniqueDays = categories(dayLabels);
    cmap = lines(numel(uniqueDays));
    hold on;
    for i = 1:numel(uniqueDays)
        idx = find(dayLabels == uniqueDays{i});
        if all(idx <= size(embedding,1))
            scatter(embedding(idx,1), embedding(idx,2), 30, cmap(i,:), 'filled');
        end
    end
    legend(uniqueDays, 'Location', 'bestoutside');
    title(sprintf('Trial-Flat UMAP Colored by Day (%d trials, %d days)', size(Xmat,1), numel(uniqueDays)));
    xlabel('UMAP 1'); ylabel('UMAP 2'); grid on;
end
