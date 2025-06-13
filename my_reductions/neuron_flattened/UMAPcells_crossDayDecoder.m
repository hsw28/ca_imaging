function [accMatrixXmat, accMatrixUMAP, pMatrixXmat, pMatrixUMAP] = UMAPcells_crossDayDecoder(animalName, win, latentDim, n_neighbors, min_dist, metric, numPerms)

    if nargin < 7, numPerms = 50; end
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

    dataCell_corr = cell(nDays,1);
    dataCell_incorr = cell(nDays,1);
    validDays = false(nDays,1);

    %% Preprocess: extract neuron-averaged data
    for d = 1:nDays

        dateStr = dateList{d};
        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
        if isempty(X), continue; end
        if numel(unique(y(~isnan(y)))) < 2, continue; end

        correct_trials = y == 1;
        incorrect_trials = y == 0;
        if sum(correct_trials) < 2 || sum(incorrect_trials) < 2, continue; end

        X_corr = squeeze(nanmean(X(:,:,correct_trials), 3));
        X_incorr = squeeze(nanmean(X(:,:,incorrect_trials), 3));

        X_corr = zscore(X_corr, 0, 2);
        X_incorr = zscore(X_incorr, 0, 2);

        dataCell_corr{d} = X_corr;
        dataCell_incorr{d} = X_incorr;
        validDays(d) = true;
    end

    accMatrixXmat = nan(nDays, nDays);
    accMatrixUMAP = nan(nDays, nDays);
    pMatrixXmat   = nan(nDays, nDays);
    pMatrixUMAP   = nan(nDays, nDays);

    %% Decode across days
    for trainDay = 1:nDays
      trainDay
        if ~validDays(trainDay), continue; end
        for testDay = 1:nDays

            if ~validDays(testDay), continue; end

            shared = find(G(:,trainDay) > 0 & G(:,testDay) > 0);
            if numel(shared) < minNeurons, continue; end

            idxTrain = G(shared, trainDay);
            idxTest  = G(shared, testDay);

            Xtrain = [dataCell_corr{trainDay}(idxTrain,:); dataCell_incorr{trainDay}(idxTrain,:)];
            ytrain = [ones(numel(shared),1); zeros(numel(shared),1)];

            Xtest  = [dataCell_corr{testDay}(idxTest,:); dataCell_incorr{testDay}(idxTest,:)];
            ytest  = [ones(numel(shared),1); zeros(numel(shared),1)];

            % Raw decoding
            if trainDay == testDay
                part = cvpartition(ytrain, 'KFold', 5);
                accsX = zeros(part.NumTestSets,1);
                for i = 1:part.NumTestSets
                    mdl = fitcsvm(Xtrain(part.training(i),:), ytrain(part.training(i)), 'KernelFunction','linear');
                    ypred = predict(mdl, Xtrain(part.test(i),:));
                    accsX(i) = mean(ypred == ytrain(part.test(i)));
                end
                accX = mean(accsX);
            else
                mdl = fitcsvm(Xtrain, ytrain, 'KernelFunction','linear');
                ypred = predict(mdl, Xtest);
                accX = mean(ypred == ytest);
            end
            accMatrixXmat(trainDay, testDay) = accX;

            % UMAP decoding
            XcatTrain = double(Xtrain);
            [Utrain, ~] = run_umap(XcatTrain, 'n_components', latentDim, ...
                'n_neighbors', n_neighbors, 'min_dist', min_dist, 'metric', metric, ...
                'verbose', false, 'cluster_output', 'none');

            if trainDay == testDay
                part = cvpartition(ytrain, 'KFold', 5);
                accsU = zeros(part.NumTestSets,1);
                for i = 1:part.NumTestSets
                    mdl = fitcsvm(Utrain(part.training(i),:), ytrain(part.training(i)), 'KernelFunction','linear');
                    ypred = predict(mdl, Utrain(part.test(i),:));
                    accsU(i) = mean(ypred == ytrain(part.test(i)));
                end
                accU = mean(accsU);
            else
                [Utest, ~] = run_umap(double(Xtest), 'n_components', latentDim, ...
                    'n_neighbors', n_neighbors, 'min_dist', min_dist, 'metric', metric, ...
                    'verbose', false, 'cluster_output', 'none');
                mdl = fitcsvm(Utrain, ytrain, 'KernelFunction','linear');
                ypred = predict(mdl, Utest);
                accU = mean(ypred == ytest);
            end
            accMatrixUMAP(trainDay, testDay) = accU;

            % Permutation test
            accPermX = nan(numPerms,1);
            accPermU = nan(numPerms,1);
            parfor p = 1:numPerms
                yshuff = ytrain(randperm(length(ytrain)));
                if trainDay == testDay
                    accXtmp = zeros(part.NumTestSets,1);
                    accUtmp = zeros(part.NumTestSets,1);
                    for i = 1:part.NumTestSets
                        mdl = fitcsvm(Xtrain(part.training(i),:), yshuff(part.training(i)), 'KernelFunction','linear');
                        accXtmp(i) = mean(predict(mdl, Xtrain(part.test(i),:)) == ytrain(part.test(i)));
                        mdl = fitcsvm(Utrain(part.training(i),:), yshuff(part.training(i)), 'KernelFunction','linear');
                        accUtmp(i) = mean(predict(mdl, Utrain(part.test(i),:)) == ytrain(part.test(i)));
                    end
                    accPermX(p) = mean(accXtmp);
                    accPermU(p) = mean(accUtmp);
                else
                    [Utest_perm, ~] = run_umap(double(Xtest), 'n_components', latentDim, ...
                        'n_neighbors', n_neighbors, 'min_dist', min_dist, 'metric', metric, ...
                        'verbose', false, 'cluster_output', 'none');
                    mdl = fitcsvm(Xtrain, yshuff, 'KernelFunction','linear');
                    accPermX(p) = mean(predict(mdl, Xtest) == ytest);
                    mdl = fitcsvm(Utrain, yshuff, 'KernelFunction','linear');
                    accPermU(p) = mean(predict(mdl, Utest_perm) == ytest);
                end
            end
            pMatrixXmat(trainDay,testDay) = (sum(accPermX >= accX)+1)/(numPerms+1);
            pMatrixUMAP(trainDay,testDay) = (sum(accPermU >= accU)+1)/(numPerms+1);
        end
    end

    % Plot p-value heatmaps
    plotCrossDayMatrix(pMatrixXmat, dateList, 'Neuron Xmat P-value Matrix', true);
    plotCrossDayMatrix(pMatrixUMAP, dateList, 'Neuron UMAP P-value Matrix', true);
end
