function [accMatrixXmat, accMatrixUMAP, pMatrixXmat, pMatrixUMAP] = UMAPtrial_crossDayDecoder(animalName, win, balanceClasses, doPermutation, numPerms, neighbors, min_dist, metric)
    % trial-flattened Cross-day decoder with aligned cells using alignmentALL
    % UMAP runs on alignedTrain + alignedTest
    % Diagonal uses 5-fold CV
    % Outputs: accMatrixXmat, accMatrixUMAP, p-matrices

    if nargin < 3, balanceClasses = false; end
    if nargin < 4, doPermutation = true; end
    if nargin < 5, numPerms = 100; end
    if nargin < 6, neighbors = 15; end
    if nargin < 7, min_dist = .7; end
    if nargin < 8, metric = 'cosine'; end

    Fs = 7.5;
    animal = evalin('base', animalName);
    dateList = autoDateList(animal);
    G = animal.alignmentALL;

    nDays = numel(dateList);
    nAligned = size(G, 2);
    if nAligned < nDays
        fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
        dateList = dateList(1:nAligned);
        nDays = nAligned;
    end

    nBins = round(diff(win) * Fs);
    minNeurons = 10;

    dataCell = cell(nDays,1);
    labelCell = cell(nDays,1);

    fprintf('Preprocessing all days...\n');
    for d = 1:nDays
        dateStr = dateList{d};
        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);

        if isempty(X) || numel(unique(y(~isnan(y)))) < 2
            fprintf('Skipping %s: insufficient data or only one class.\n', dateStr);
            continue;
        end

        nNeurons = size(X,1);
        if nNeurons < minNeurons
            fprintf('Skipping %s: too few neurons (%d).\n', dateStr, nNeurons);
            continue;
        end

        nTrials = size(X,3);
        trialVecs = zeros(nTrials, nNeurons * nBins);
        for t = 1:nTrials
            mat = X(:,:,t);
            mat(isnan(mat)) = 0;
            trialVecs(t,:) = mat(:)';
        end
        trialVecs = zscore(trialVecs);
        y = y(:);

        dataCell{d} = trialVecs;
        labelCell{d} = y;
    end

    accMatrixXmat = nan(nDays, nDays);
    accMatrixUMAP = nan(nDays, nDays);
    pMatrixXmat   = nan(nDays, nDays);
    pMatrixUMAP   = nan(nDays, nDays);

    opts = {'KernelFunction','linear'};
    if balanceClasses
        opts = [opts, {'ClassNames',[0;1]}];
    end

    fprintf('Starting cross-day decoding with aligned cells and aligned UMAP...\n');
    for trainDay = 1:nDays
        Xtrain_all = dataCell{trainDay};
        ytrain = labelCell{trainDay};
        if isempty(Xtrain_all), continue; end

        for testDay = 1:nDays
            Xtest_all = dataCell{testDay};
            ytest = labelCell{testDay};
            if isempty(Xtest_all), continue; end

            validCells = find(G(:,trainDay) > 0 & G(:,testDay) > 0);
            if isempty(validCells)
                fprintf('Skipping pair Train %d - Test %d: no aligned cells.\n', trainDay, testDay);
                continue;
            end

            trainIdx = G(validCells, trainDay);
            testIdx  = G(validCells, testDay);

            nTrialsTrain = size(Xtrain_all,1);
            nCells = numel(validCells);
            alignedTrain = zeros(nTrialsTrain, nCells * nBins);
            for c = 1:nCells
                idxTrain = trainIdx(c);
                cols = (idxTrain-1)*nBins + (1:nBins);
                alignedTrain(:, (c-1)*nBins + (1:nBins)) = Xtrain_all(:, cols);
            end

            nTrialsTest = size(Xtest_all,1);
            alignedTest = zeros(nTrialsTest, nCells * nBins);
            for c = 1:nCells
                idxTest = testIdx(c);
                cols = (idxTest-1)*nBins + (1:nBins);
                alignedTest(:, (c-1)*nBins + (1:nBins)) = Xtest_all(:, cols);
            end

            % === UMAP on alignedTrain + alignedTest ===
            combined = [alignedTrain; alignedTest];
            [embedding, ~] = run_umap(combined, 'n_components', 3, ...
                'n_neighbors', neighbors, ...
                'min_dist', min_dist, ...
                'metric', metric, ...
                'randomize', true, 'verbose', false);

            Utrain = embedding(1:size(alignedTrain,1), :);
            Utest  = embedding(size(alignedTrain,1)+1:end, :);

            % === Diagonal: 5-fold CV ===
            if trainDay == testDay
                % Xmat CV
                part = cvpartition(ytrain, 'KFold', 5);
                accsX = zeros(part.NumTestSets,1);
                accsU = zeros(part.NumTestSets,1);
                for i = 1:part.NumTestSets
                    trIdx = training(part,i);
                    teIdx = test(part,i);

                    mdl = fitcsvm(alignedTrain(trIdx,:), ytrain(trIdx), opts{:});
                    ypred = predict(mdl, alignedTrain(teIdx,:));
                    accsX(i) = mean(ypred == ytrain(teIdx));

                    mdl = fitcsvm(Utrain(trIdx,:), ytrain(trIdx), opts{:});
                    ypred = predict(mdl, Utrain(teIdx,:));
                    accsU(i) = mean(ypred == ytrain(teIdx));
                end
                accX = mean(accsX);
                accU = mean(accsU);
            else
                % Off-diagonal: train all, test all
                mdlX = fitcsvm(alignedTrain, ytrain, opts{:});
                ypredX = predict(mdlX, alignedTest);
                accX = mean(ypredX == ytest);

                mdlU = fitcsvm(Utrain, ytrain, opts{:});
                ypredU = predict(mdlU, Utest);
                accU = mean(ypredU == ytest);
            end

            accMatrixXmat(trainDay,testDay) = accX;
            accMatrixUMAP(trainDay,testDay) = accU;

            % === Permutation test ===
            if doPermutation
                permAccX = nan(numPerms,1);
                permAccU = nan(numPerms,1);
                for p = 1:numPerms
                    ytrain_shuff = ytrain(randperm(length(ytrain)));

                    if trainDay == testDay
                        % Diagonal: CV
                        part = cvpartition(ytrain_shuff, 'KFold', 5);
                        accsX = zeros(part.NumTestSets,1);
                        accsU = zeros(part.NumTestSets,1);
                        for i = 1:part.NumTestSets
                            trIdx = training(part,i);
                            teIdx = test(part,i);

                            mdl = fitcsvm(alignedTrain(trIdx,:), ytrain_shuff(trIdx), opts{:});
                            ypred = predict(mdl, alignedTrain(teIdx,:));
                            accsX(i) = mean(ypred == ytrain(teIdx));

                            mdl = fitcsvm(Utrain(trIdx,:), ytrain_shuff(trIdx), opts{:});
                            ypred = predict(mdl, Utrain(teIdx,:));
                            accsU(i) = mean(ypred == ytrain(teIdx));
                        end
                        permAccX(p) = mean(accsX);
                        permAccU(p) = mean(accsU);
                    else
                        % Off-diagonal: train all, test all
                        mdlX = fitcsvm(alignedTrain, ytrain_shuff, opts{:});
                        ypredX = predict(mdlX, alignedTest);
                        permAccX(p) = mean(ypredX == ytest);

                        mdlU = fitcsvm(Utrain, ytrain_shuff, opts{:});
                        ypredU = predict(mdlU, Utest);
                        permAccU(p) = mean(ypredU == ytest);
                    end
                end
                % Compute p-values:
                pMatrixXmat(trainDay,testDay) = (sum(permAccX >= accX) + 1) / (numPerms + 1);
                pMatrixUMAP(trainDay,testDay) = (sum(permAccU >= accU) + 1) / (numPerms + 1);
            end
        end
    end
    %    fprintf('Plotting results...\n');
         plotCrossDayMatrix(accMatrixXmat, dateList, 'Xmat Aligned Cross-Day Decoding', false);
         plotCrossDayMatrix(accMatrixUMAP, dateList, 'UMAP Cross-Day Decoding', false);

        if doPermutation
            plotCrossDayMatrix(pMatrixXmat, dateList, 'Xmat P-value Matrix', true);
            plotCrossDayMatrix(pMatrixUMAP, dateList, 'UMAP P-value Matrix', true);
        end

    fprintf('Cross-day decoding done!\n');
end
