function [accMatrixXmat, pMatrixXmat] = XMATcells_crossDayDecoder(animalName, win, Fs, kernel, prep, clf, balanceClasses, numPerms)

  if nargin < 1, error('Must provide animalName'); end
  if nargin < 2, win = [0, 0.75]; end
  if nargin < 3, Fs_vals = [7.5]; end
  if nargin < 4, kernels = {'linear'}; end
  if nargin < 5, preprocessingModes = {'zscore'}; end
  if nargin < 6, classifiers = {'svm'}; end
  if nargin < 7, balanceClasses = false; end
  if nargin < 8, numPerms = 100; end


    animal = evalin('base', animalName);
    dateList = autoDateList(animal);
    G = animal.alignmentALL;
    nDays = numel(dateList);
    nAligned = size(G, 2);
    if nAligned < nDays
      %  fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
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

        switch prep
            case 'zscore'
                X_corr = zscore(X_corr, 0, 2);
                X_incorr = zscore(X_incorr, 0, 2);
            case 'normalize'
                X_corr = normalize(X_corr, 2);
                X_incorr = normalize(X_incorr, 2);
            case 'none'
                % do nothing
            otherwise
                error('Unknown preprocessing method');
        end

        dataCell_corr{d} = X_corr;
        dataCell_incorr{d} = X_incorr;
        validDays(d) = true;
    end

    accMatrixXmat = nan(nDays, nDays);
    pMatrixXmat   = nan(nDays, nDays);

    %% Decode across days
    for trainDay = 1:nDays
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

            if balanceClasses
                minCount = min(sum(ytrain==1), sum(ytrain==0));
                idx1 = find(ytrain==1); idx0 = find(ytrain==0);
                sel = [randsample(idx1, minCount); randsample(idx0, minCount)];
                Xtrain = Xtrain(sel,:);
                ytrain = ytrain(sel);

                idx1 = find(ytest==1); idx0 = find(ytest==0);
                sel = [randsample(idx1, minCount); randsample(idx0, minCount)];
                Xtest = Xtest(sel,:);
                ytest = ytest(sel);
            end

            if strcmp(clf, 'svm')
                decodeFn = @(X, y) fitcsvm(X, y, 'KernelFunction', kernel);
            elseif strcmp(clf, 'nb')
                decodeFn = @(X, y) fitcnb(X, y);
            elseif strcmp(clf, 'tree')
                decodeFn = @(X, y) fitctree(X, y);
            else
                error('Unsupported classifier');
            end

            if trainDay == testDay
                part = cvpartition(ytrain, 'KFold', 5);
                accsX = zeros(part.NumTestSets,1);
                for i = 1:part.NumTestSets
                    mdl = decodeFn(Xtrain(part.training(i),:), ytrain(part.training(i)));
                    ypred = predict(mdl, Xtrain(part.test(i),:));
                    accsX(i) = mean(ypred == ytrain(part.test(i)));
                end
                accX = mean(accsX);
            else
                mdl = decodeFn(Xtrain, ytrain);
                ypred = predict(mdl, Xtest);
                accX = mean(ypred == ytest);
            end
            accMatrixXmat(trainDay, testDay) = accX;

            %% Permutation test
            accPermX = nan(numPerms,1);
            parfor p = 1:numPerms
                yshuff = ytrain(randperm(length(ytrain)));
                if trainDay == testDay
                    accXtmp = zeros(part.NumTestSets,1);
                    for i = 1:part.NumTestSets
                        mdl = decodeFn(Xtrain(part.training(i),:), yshuff(part.training(i)));
                        ypred = predict(mdl, Xtrain(part.test(i),:));
                        accXtmp(i) = mean(ypred == ytrain(part.test(i)));
                    end
                    accPermX(p) = mean(accXtmp);
                else
                    mdl = decodeFn(Xtrain, yshuff);
                    accPermX(p) = mean(predict(mdl, Xtest) == ytest);
                end
            end
            pMatrixXmat(trainDay,testDay) = (sum(accPermX >= accX)+1)/(numPerms+1);
        end
    end

    % Plot p-value heatmaps
    plotCrossDayMatrix(pMatrixXmat, dateList, 'Neuron Xmat P-value Matrix', true);
    end
