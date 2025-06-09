function decodeResults = bayesDecoder_REFgrouped(animalName, trainDays, nPerms, minTrainDays, minTestDays)
% Bayesian decoding with multiple training days, pruning neurons shared across days,
% cross-validation, and permutation testing.
%
% Inputs:
%   animalName  : string, name of data struct in base workspace (e.g. 'rat0314')
%   trainDays   : vector of day indices to train on (e.g. [7 11])
%   nPerms      : number of permutations for p-value estimation (default 1)
%   minTrainDays: minimum days a neuron must appear in training (default 1)
%   minTestDays : minimum days a neuron must appear in testing (default 1)

if nargin < 3, nPerms = 1; end
if nargin < 4, minTrainDays = 1; end
if nargin < 5, minTestDays = 1; end

win = [0, .75];
Fs = 7.5;
nBins = round((win(end) - win(1)) * Fs);

animal = evalin('base', animalName);
G = animal.alignmentALL;
dateList = autoDateList(animal);
nDays = min(numel(dateList), size(G, 2));

decodeResults = struct('testDate', {}, 'nSharedNeurons', {}, ...
    'acc_trial', {}, 'f1_trial', {}, 'pval_acc', {}, 'pval_f1', {}, ...
    'acc_time', {}, 'f1_time', {}, 'pval_time_acc', {}, 'pval_time_f1', {});

% Find neurons present in at least minTrainDays of trainDays
sharedTrain = sum(G(:, trainDays) > 0, 2) >= minTrainDays;
trainNeuronGlobalIDs = find(sharedTrain);
nNeuronsTrain = numel(trainNeuronGlobalIDs);

Xtrain_all = [];
ytrain_all = [];
expectedFeatures = [];

fprintf('Collecting training data...\n');
for i = 1:numel(trainDays)
    d = trainDays(i);
    dateStr = dateList{d};
    [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
    if isempty(X) || numel(unique(y)) < 2
        fprintf('Skipping train day %s: insufficient data or only one class.\n', dateStr);
        continue;
    end

    trainIDs = G(sharedTrain, d);
    validIdx = trainIDs > 0;
    X = X(trainIDs(validIdx), :, :);

    trialMask = squeeze(all(all(~isnan(X),1),2));
    X = X(:,:,trialMask);
    y = y(trialMask);
    valid = ~isnan(y);
    X = X(:,:,valid);
    y = y(valid);

    % Align neurons by global ID
    X_aligned = nan(nNeuronsTrain, size(X, 2), size(X, 3));
    localIDs = G(trainNeuronGlobalIDs, d);
    validAlign = localIDs > 0 & localIDs <= size(X, 1);
    X_aligned(validAlign, :, :) = X(localIDs(validAlign), :, :);

    Xflat = reshape(X_aligned, [], size(X_aligned, 3))';

    if isempty(expectedFeatures)
        expectedFeatures = size(Xflat, 2);
    elseif size(Xflat, 2) ~= expectedFeatures
        fprintf('⚠️ Skipping %s trial block: has %d features, expected %d\n', dateStr, size(Xflat, 2), expectedFeatures);
        continue;
    end

    Xtrain_all = [Xtrain_all; Xflat];
    ytrain_all = [ytrain_all; y];
end

if isempty(Xtrain_all)
    error('No valid training data found after alignment.');
end
fprintf('Initial training data size: %d trials x %d features\n', size(Xtrain_all,1), size(Xtrain_all,2));

for testDay = 1:nDays
    if ismember(testDay, trainDays)
        continue;
    end

    testStr = dateList{testDay};
    [Xtest, ytest] = getDayMatrixFromStruct(animal, testStr, win, nBins, Fs);
    if isempty(Xtest) || numel(unique(ytest)) < 2
        fprintf('Skipping %s: insufficient data or only one class.\n', testStr);
        continue;
    end

    % Find neurons shared between train and test days
    sharedTest = sharedTrain & (G(:, testDay) > 0);
    if sum(sharedTest) < 5
        fprintf('Skipping %s: fewer than 5 neurons shared with training days.\n', testStr);
        continue;
    end

    prunedNeurons = find(sharedTest);
    nNeuronsPruned = numel(prunedNeurons);
    localIDsTest = G(prunedNeurons, testDay);
    XtestLocalIDs = G(G(:, testDay) > 0, testDay);
    XtestRowMap = containers.Map(XtestLocalIDs, 1:numel(XtestLocalIDs));

    % Align test data neurons to global indices
    X_aligned_test = nan(nNeuronsPruned, size(Xtest, 2), size(Xtest, 3));
    for i = 1:nNeuronsPruned
        lid = localIDsTest(i);
        if lid > 0 && isKey(XtestRowMap, lid)
            X_aligned_test(i, :, :) = Xtest(XtestRowMap(lid), :, :);
        end
    end
    Xtest = X_aligned_test;

    % Remove neurons with all NaNs in test data and prune training neurons accordingly
    nanRows = all(all(isnan(Xtest), 3), 2);
    if any(nanRows)
        Xtest(nanRows, :, :) = [];
        keepNeurons = ~nanRows;

        [~, prunedNeuronIdxInTrain] = ismember(prunedNeurons, trainNeuronGlobalIDs);
        fullKeepNeurons = false(length(sharedTrain), 1);
        fullKeepNeurons(prunedNeuronIdxInTrain) = keepNeurons;

        keptNeuronsTrain = sharedTrain;
        neuronsToKeep = keptNeuronsTrain & fullKeepNeurons;

        % Prune training neurons before zero-var removal
        trainNeuronsMask = false(1, nNeuronsTrain);
        trainNeuronsMask(neuronsToKeep) = true;

        % Prune training data to only neurons also present and valid in test
        trainNeuronsIdx = find(trainNeuronsMask);
        featuresPerNeuron = nBins;
        keepTrainFeaturesMask = false(1, length(keepCols));
        for idx = trainNeuronsIdx'
            featureRange = (idx-1)*featuresPerNeuron + (1:featuresPerNeuron);
            keepTrainFeaturesMask(featureRange) = true;
        end

        Xtrain_all_pruned = Xtrain_all(:, keepTrainFeaturesMask);
        Xtest_flat = reshape(Xtest, [], size(Xtest, 3))';

        % Prune test features accordingly (just keep all pruned neurons * bins)
        keepTestFeaturesMask = true(1, size(Xtest_flat, 2));

        % Joint zero variance removal on train and test
        allData = [Xtrain_all_pruned; Xtest_flat(:, keepTestFeaturesMask)];
        jointKeepCols = var(allData, 0, 1) > 0;

        Xtrain_all_clean_pruned = Xtrain_all_pruned(:, jointKeepCols);
        Xtest_flat_clean = Xtest_flat(:, keepTestFeaturesMask);
        Xtest_flat_clean = Xtest_flat_clean(:, jointKeepCols);

    else
        % No neurons removed, just prune zero-variance features on full training data and test data

        Xtest_flat = reshape(Xtest, [], size(Xtest, 3))';

        % Map shared neurons in training data
        [~, testNeuronIdxInTrain] = ismember(prunedNeurons, trainNeuronGlobalIDs);

        % Initialize mask for training features (neurons * bins)
        trainFeatureMask = false(1, size(Xtrain_all, 2));
        keepColsMat = reshape(true(1, size(Xtrain_all, 2)), [], nBins); % initially keep all

        % Build mask for neurons shared in test and train
        neuronsKeptInTest = false(length(prunedNeurons),1);
        for i = 1:length(prunedNeurons)
            trainNeuronIdx = testNeuronIdxInTrain(i);
            if trainNeuronIdx > 0
                neuronsKeptInTest(i) = true;
                featureRange = (trainNeuronIdx-1)*nBins + (1:nBins);
                trainFeatureMask(featureRange) = true;
            end
        end

        % Prune training data
        Xtrain_all_pruned = Xtrain_all(:, trainFeatureMask);

        % Prune test data to only these neurons and bins
        Xtest_flat_pruned = [];
        for i = 1:length(prunedNeurons)
            if neuronsKeptInTest(i)
                featureRange = ((i-1)*nBins + 1) : (i*nBins);
                Xtest_flat_pruned = [Xtest_flat_pruned, Xtest_flat(:, featureRange)];
            end
        end

        % Now concatenate pruned data for joint zero-variance removal
        allData = [Xtrain_all_pruned; Xtest_flat_pruned];
        jointKeepCols = var(allData, 0, 1) > 0;

        Xtrain_all_clean_pruned = Xtrain_all(:, jointKeepCols);
        Xtest_flat_clean = Xtest_flat(:, jointKeepCols);
    end


    if size(Xtest_flat_clean, 2) ~= size(Xtrain_all_clean_pruned, 2)
        error('Mismatch in number of features after pruning test data.');
    end

    accT = crossvalBayesAccuracy(Xtrain_all_clean_pruned, ytrain_all, Xtest_flat_clean, ytest);
    f1T = crossvalBayesF1(Xtrain_all_clean_pruned, ytrain_all, Xtest_flat_clean, ytest);

    perm_acc = nan(nPerms, 1);
    perm_f1 = nan(nPerms, 1);

    if nPerms > 1 && isempty(gcp('nocreate'))
        parpool('local');
    end

    if nPerms > 1
        parfor p = 1:nPerms
            yshuf = ytrain_all(randperm(numel(ytrain_all)));
            perm_acc(p) = crossvalBayesAccuracy(Xtrain_all_clean_pruned, yshuf, Xtest_flat_clean, ytest);
            perm_f1(p) = crossvalBayesF1(Xtrain_all_clean_pruned, yshuf, Xtest_flat_clean, ytest);
        end
    end

    pval_acc = (sum(perm_acc >= accT) + 1) / (nPerms + 1);
    pval_f1 = (sum(perm_f1 >= f1T) + 1) / (nPerms + 1);

    fprintf('Test %s — Trial-flat: Acc = %.2f, F1 = %.3f, pAcc = %.3f, pF1 = %.3f\n', ...
        testStr, accT, f1T, pval_acc, pval_f1);

    decodeResults(end+1) = struct('testDate', testStr, 'nSharedNeurons', nNeuronsPruned, ...
        'acc_trial', accT, 'f1_trial', f1T, 'pval_acc', pval_acc, 'pval_f1', pval_f1, ...
        'acc_time', NaN, 'f1_time', NaN, 'pval_time_acc', NaN, 'pval_time_f1', NaN);
end
end


function [Xclean, keepCols] = removeZeroVarPerClass(X_in, y_in)
% Removes features that have zero variance in any class to prevent Naive Bayes errors
    classes = unique(y_in);
    keepCols = true(1, size(X_in, 2));
    for c = classes(:)'
        idx = y_in == c;
        zeroVarCols = var(X_in(idx,:), 0, 1) == 0;
        keepCols(zeroVarCols) = false;
    end
    Xclean = X_in(:, keepCols);
end

function acc = crossvalBayesAccuracy(Xtrain, ytrain, Xtest, ytest)
    mdlStruct = fitcnbWithVarCheck(Xtrain, ytrain);
    Xtest_clean = Xtest(:, mdlStruct.FeatureSelection);
    yhat = predict(mdlStruct.Model, Xtest_clean);
    acc = mean(yhat == ytest);

  %  if all(yhat == yhat(1))
  %      fprintf('[Warning] All predicted acc labels are identical (%d)!\n', yhat(1));
  %  end
end

function f1 = crossvalBayesF1(Xtrain, ytrain, Xtest, ytest)
    mdlStruct = fitcnbWithVarCheck(Xtrain, ytrain);
    Xtest_clean = Xtest(:, mdlStruct.FeatureSelection);
    yhat = predict(mdlStruct.Model, Xtest_clean);
    f1 = f1score(ytest, yhat);

  %  if all(yhat == yhat(1))
  %      fprintf('[Warning] All predicted f1 labels are identical (%d)!\n', yhat(1));
  %  end
end



function mdlStruct = fitcnbWithVarCheck(X, y)
    % Remove zero variance features per class
    classes = unique(y);
    keepCols = true(1, size(X, 2));
    for c = classes(:)'
        idx = y == c;
        zeroVarCols = var(X(idx,:), 0, 1) == 0;
        keepCols(zeroVarCols) = false;
    end
    Xclean = X(:, keepCols);
%    mdl = fitcnb(Xclean, y);
    mdl = fitcnb(Xclean, y, 'DistributionNames', 'kernel', 'Prior', 'empirical');

    mdlStruct.Model = mdl;
    mdlStruct.FeatureSelection = keepCols;
end

function f1 = f1score(ytrue, ypred)
% Computes F1 score for binary classification
    ytrue = double(ytrue(:));
    ypred = double(ypred(:));
    tp = sum((ytrue == 1) & (ypred == 1));
    fp = sum((ytrue == 0) & (ypred == 1));
    fn = sum((ytrue == 1) & (ypred == 0));
    prec = tp / (tp + fp + eps);
    rec = tp / (tp + fn + eps);
    f1 = 2 * (prec * rec) / (prec + rec + eps);
end
