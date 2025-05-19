function decodeResults = kernelDecoder_REFgrouped(animalName, trainDays, nPerms, minTrainDays, minTestDays)
%% trial flattened decoder with lambda search and class weights

if nargin < 3, nPerms = 1; end
if nargin < 4, minTrainDays = 1; end
if nargin < 5, minTestDays = 1; end

win = [0, 1.3];
Fs = 7.5;
nBins = round((win(end) - win(1)) * Fs);

animal = evalin('base', animalName);
G = animal.alignmentALL;
dateList = autoDateList(animal);
nDays = min(numel(dateList), size(G, 2));

decodeResults = struct('testDate', [], 'nSharedNeurons', [], ...
    'acc_trial', [], 'f1_trial', [], 'pval_acc', [], 'pval_f1', [], ...
    'acc_time', [], 'f1_time', [], 'pval_time_acc', [], 'pval_time_f1', []);

sharedTrain = sum(G(:,trainDays) > 0, 2) >= minTrainDays;

Xtrain_all = [];
ytrain_all = [];
trainNeuronGlobalIDs = [];

expectedFeatureCount = [];

for i = 1:numel(trainDays)
    d = trainDays(i);
    dateStr = dateList{d};
    [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
    if isempty(X) || numel(unique(y)) < 2, continue; end

    trainIDs = G(sharedTrain, d);
    validIdx = trainIDs > 0;
    globalNeuronIDs = find(sharedTrain & G(:, d) > 0);
    trainNeuronGlobalIDs = [trainNeuronGlobalIDs; globalNeuronIDs];

    X = X(trainIDs(validIdx),:,:);


    trialMask = squeeze(all(all(~isnan(X),1),2));
    X = X(:,:,trialMask);
    y = y(trialMask);
    valid = ~isnan(y);
    X = X(:,:,valid);
    y = y(valid);

    fprintf('  Day %s: using %d neurons for training\n', dateStr, sum(validIdx));

    trainNeuronGlobalIDs = unique(trainNeuronGlobalIDs, 'stable');
    nNeuronsTrain = numel(trainNeuronGlobalIDs);

    X_aligned = nan(nNeuronsTrain, size(X,2), size(X,3));
    localIDs = G(trainNeuronGlobalIDs, d);
    validAlign = localIDs > 0 & localIDs <= size(X,1);
    X_aligned(validAlign,:,:) = X(localIDs(validAlign),:,:);

    Xflat = reshape(X_aligned, [], size(X_aligned,3))';

    % Check and enforce consistent feature count
    if isempty(expectedFeatureCount)
        expectedFeatureCount = size(Xflat, 2);
    elseif size(Xflat, 2) ~= expectedFeatureCount
        fprintf('⚠️ Skipping %s trial block: has %d features, expected %d\n', ...
            dateStr, size(Xflat, 2), expectedFeatureCount);
        % Remove added neuron IDs for this skipped day
        trainNeuronGlobalIDs = setdiff(trainNeuronGlobalIDs, globalNeuronIDs, 'stable');
        continue;
    end

    Xtrain_all = [Xtrain_all; Xflat];
    ytrain_all = [ytrain_all; y];
end


trainNeuronGlobalIDs = unique(trainNeuronGlobalIDs, 'stable');
nNeuronsTrain = length(trainNeuronGlobalIDs);

lambdas = logspace(-6, -1, 10);

for testDay = 1:nDays
    if ismember(testDay, trainDays), continue; end

    testStr = dateList{testDay};
    [Xtest, ytest] = getDayMatrixFromStruct(animal, testStr, win, nBins, Fs);
    if isempty(Xtest) || numel(unique(ytest)) < 2
        fprintf('Skipping %s: insufficient data or only one class.\n', testStr);
        continue;
    end

    shared = sum(G(:,trainDays) > 0, 2) >= minTrainDays & G(:,testDay) > 0;
    if sum(shared) < 5
        fprintf('Skipping %s: only %d neurons shared with training days.\n', testStr, sum(shared));
        continue;
    end

    idxTest = G(shared, testDay);
    testNeuronIDs = G(:,testDay);
    [~, testIdxOrder] = ismember(trainNeuronGlobalIDs, testNeuronIDs);
    valid = testIdxOrder > 0;
    testIdxOrder = testIdxOrder(valid);
    alignedTrainNeurons = trainNeuronGlobalIDs(valid);

    localIDs = G(trainNeuronGlobalIDs, testDay);
    XtestLocalIDs = G(G(:, testDay) > 0, testDay);
    XtestRowMap = containers.Map(XtestLocalIDs, 1:numel(XtestLocalIDs));

    X_aligned = nan(nNeuronsTrain, size(Xtest,2), size(Xtest,3));
    for i = 1:numel(trainNeuronGlobalIDs)
        lid = localIDs(i);
        if lid > 0 && isKey(XtestRowMap, lid)
            X_aligned(i,:,:) = Xtest(XtestRowMap(lid),:,:);
        end
    end
    Xtest = X_aligned;

    nanRows = all(all(isnan(Xtest),3),2);
    if any(nanRows)
        Xtest(nanRows,:,:) = [];
        keepNeurons = ~nanRows;


        neuronMask = repmat(keepNeurons(:), nBins, 1);
        neuronMask = neuronMask(:)';

        if length(neuronMask) ~= size(Xtrain_all, 2)
            error('Neuron mask size (%d) does not match training data features (%d).', length(neuronMask), size(Xtrain_all, 2));
        end

        Xtrain_all_test = Xtrain_all(:, neuronMask);
        nNeuronsKeep = sum(keepNeurons);
        Xtrain_reshaped_test = reshape(Xtrain_all_test', nNeuronsKeep, nBins, []);
    else
        Xtrain_all_test = Xtrain_all;
    end

    trialMask = squeeze(all(all(~isnan(Xtest),1),2));
    Xtest = Xtest(:,:,trialMask);
    ytest = ytest(trialMask);

    if isempty(ytest) || size(Xtest,3) == 0
        fprintf('⚠️ Skipping %s: all test trials removed due to NaNs.\n', testStr);
        continue;
    end

    Xtest_flat = reshape(Xtest, [], size(Xtest,3))';

    colMeans = nanmean(Xtrain_all_test, 1);
    nanMask = isnan(Xtrain_all_test);
    Xfill = repmat(colMeans, size(Xtrain_all_test,1), 1);
    Xtrain_all_test(nanMask) = Xfill(nanMask);

    classes = unique(ytrain_all);
    freq = histcounts(ytrain_all, [classes; max(classes)+1]);
    weights = 1 ./ freq;
    classWeights = containers.Map(classes, weights);

    grid_results = table();
    for l = lambdas
        try
            sampleWeights = arrayfun(@(c) classWeights(c), ytrain_all);
            mdl = fitcsvm(Xtrain_all_test, ytrain_all, 'KernelFunction', 'linear', ...
                'BoxConstraint', l, 'Standardize', true, 'ClassNames', classes, 'Weights', sampleWeights);
            yhat = predict(mdl, Xtest_flat);
            acc = mean(yhat == ytest);
            f1 = f1score(ytest, yhat);
            grid_results = [grid_results; table(l, acc, f1)];
        catch
            % skip errors
        end
    end

    if isempty(grid_results)
        warning('No valid models trained on test day %s', testStr);
        continue;
    end

    [~, idxBest] = max(grid_results.acc);
    best_lambda = grid_results.l(idxBest);
    best_acc = grid_results.acc(idxBest);
    best_f1 = grid_results.f1(idxBest);

    sampleWeights = arrayfun(@(c) classWeights(c), ytrain_all);
    mdl = fitcsvm(Xtrain_all_test, ytrain_all, 'KernelFunction', 'linear', ...
        'BoxConstraint', best_lambda, 'Standardize', true, 'ClassNames', classes, 'Weights', sampleWeights);
    yhat = predict(mdl, Xtest_flat);
    acc = mean(yhat == ytest);
    f1 = f1score(ytest, yhat);

    perm_acc = nan(nPerms, 1);
    perm_f1 = nan(nPerms, 1);
    for p = 1:nPerms
        yshuf = ytrain_all(randperm(numel(ytrain_all)));
        mdl_shuf = fitcsvm(Xtrain_all_test, yshuf, 'KernelFunction', 'linear', ...
            'BoxConstraint', best_lambda, 'Standardize', true, 'ClassNames', classes, 'Weights', sampleWeights);
        yhat_shuf = predict(mdl_shuf, Xtest_flat);
        perm_acc(p) = mean(yhat_shuf == ytest);
        perm_f1(p) = f1score(ytest, yhat_shuf);
    end

    pval_acc = (sum(perm_acc >= acc) + 1) / (nPerms + 1);
    pval_f1 = (sum(perm_f1 >= f1) + 1) / (nPerms + 1);

    fprintf('Test %s — Trial-flat: Acc = %.2f, F1 = %.3f, Best Lambda = %.1e, pAcc = %.3f, pF1 = %.3f\n', ...
        testStr, acc, f1, best_lambda, pval_acc, pval_f1);

    decodeResults(end+1) = struct('testDate', testStr, 'nSharedNeurons', nNeuronsTrain, ...
        'acc_trial', acc, 'f1_trial', f1, 'pval_acc', pval_acc, 'pval_f1', pval_f1, ...
        'acc_time', NaN, 'f1_time', NaN, 'pval_time_acc', NaN, 'pval_time_f1', NaN);
end
end

function f1 = f1score(ytrue, ypred)
    ytrue = double(ytrue(:));
    ypred = double(ypred(:));
    tp = sum((ytrue == 1) & (ypred == 1));
    fp = sum((ytrue == 0) & (ypred == 1));
    fn = sum((ytrue == 1) & (ypred == 0));
    prec = tp / (tp + fp + eps);
    rec  = tp / (tp + fn + eps);
    f1 = 2 * (prec * rec) / (prec + rec + eps);
end
