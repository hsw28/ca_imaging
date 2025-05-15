function decodeResults = kernelDecoder_REFgrouped(animalName, trainDays, nPerms, minTrainDays, minTestDays)
% KERNELDECODER_GROUPED
% ----------------------
% Performs kernel-based (RBF SVM) decoding using trial- and timepoint-flattened representations.
% Trains on trials from multiple aligned days and tests on all others.
% Only neurons shared in >= minTrainDays train days and minTestDays test days are used.
%
% INPUTS:
%   animalName    : string, name of struct in base workspace (e.g. 'rat0314')
%   trainDays     : vector of day indices to train on (e.g. [1 2 3])
%   nPerms        : number of permutations for p-value estimation (e.g. 100)
%   minTrainDays  : min number of trainDays a neuron must appear in (default: all)
%   minTestDays   : min number of testDays a neuron must appear in (default: 1)
%
% OUTPUT:
%   decodeResults : struct array with decoding metrics for each test day
  % Kernel-based decoder (RBF SVM) trained on multiple days
  % Trains on combined trials from trainDays and tests on remaining days
  % Uses only shared neurons present in minTrainDays and testDay

  if nargin < 3
    nPerms = 1;
  end
  if nargin < 4
    minTrainDays = numel(trainDays);  % default to strict intersection
  end
  if nargin < 5
    minTestDays = 1;  % lenient default: neuron must appear in test day
  end

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


  sharedTrain = sum(G(:,trainDays) > 0, 2) >= minTrainDays; %get index of cells present in at least minTrainDays

  Xtrain_all = [];
  ytrain_all = [];
  trainNeuronGlobalIDs = [];

  for i = 1:numel(trainDays)
    d = trainDays(i);
    dateStr = dateList{d};
    [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
    if isempty(X) || numel(unique(y)) < 2, continue; end

    % --- Use G to directly get valid neuron indices for this training day ---
    trainIDs = G(sharedTrain, d); %finds shared cells present in said day
    validIdx = trainIDs > 0; %targets the cells that actually exist

    globalNeuronIDs = find(sharedTrain & G(:, d) > 0);  % global neuron indices for this day
    trainNeuronGlobalIDs = [trainNeuronGlobalIDs; globalNeuronIDs];

    X = X(trainIDs(validIdx),:,:);  % Correctly index using global indices, identifies and sorts these cells

    fprintf('  Day %s: using %d neurons for training', dateStr, sum(validIdx));

    %removing NaNs
    trialMask = squeeze(all(all(~isnan(X),1),2));
    X = X(:,:,trialMask);
    y = y(trialMask);
    valid = ~isnan(y);
    X = X(:,:,valid);
    y = y(valid);

% === Gather all global neuron IDs across training days ===
trainNeuronGlobalIDs = [];
for i = 1:numel(trainDays)
    e = trainDays(i);
    globalNeuronIDs = find(sharedTrain & G(:, e) > 0);
    trainNeuronGlobalIDs = [trainNeuronGlobalIDs; globalNeuronIDs];
end
trainNeuronGlobalIDs = unique(trainNeuronGlobalIDs, 'stable');
nNeuronsTrain = numel(trainNeuronGlobalIDs);


% Align training X to full shared neuron set (with padding for missing)
X_aligned = nan(nNeuronsTrain, size(X,2), size(X,3));
localIDs = G(trainNeuronGlobalIDs, d);
validAlign = localIDs > 0 & localIDs <= size(X,1);
X_aligned(validAlign,:,:) = X(localIDs(validAlign),:,:);

% Use the aligned data now
Xflat = reshape(X_aligned, [], size(X_aligned,3))';


% Check feature length AFTER alignment
expectedCols = nNeuronsTrain * nBins;
if size(Xflat, 2) ~= expectedCols
    fprintf('⚠️ Skipping %s trial block: has %d features, expected %d\n', ...
        dateStr, size(Xflat, 2), expectedCols);
    continue;
end


    Xtrain_all = [Xtrain_all; Xflat];
    ytrain_all = [ytrain_all; y];
  end

  trainNeuronGlobalIDs = unique(trainNeuronGlobalIDs, 'stable');


% Enforce neuron order for Xtrain_all before timepoint reshaping
  nNeuronsTrain = length(trainNeuronGlobalIDs);
  Xtrain_all';  % shape: [nFeatures × nTrials] -- want to reshape to [nNeuronsTrain × nBins × nTrials]
  Xtrain_reshaped = reshape(Xtrain_all', nNeuronsTrain, [], size(Xtrain_all,1));
  Xtrain_time = reshape(permute(Xtrain_reshaped, [1 3 2]), nNeuronsTrain, [])';
  ytrain_time = repmat(ytrain_all(:), nBins, 1);

  % Preserve full training timepoint data
    Xtrain_time_full = Xtrain_time;
    ytrain_time_full = ytrain_time;




  for testDay = 1:nDays
    if ismember(testDay, trainDays), continue; end
      % Keep full training matrix for timepoint SVM (do NOT prune by test neuron set)
      Xtrain_time = Xtrain_time_full;
      ytrain_time = ytrain_time_full;
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
    idxTest = testNeuronIDs(testIdxOrder);
    idxTest = idxTest(:);  % ensure column vector
    if isempty(idxTest) || any(idxTest <= 0 | isnan(idxTest))
      error('Invalid or empty indices in idxTest. Some neurons may not be aligned correctly.');
    end
    if any(idxTest <= 0 | isnan(idxTest))
      error('Invalid indices in idxTest. Some neurons may not be aligned correctly.');
    end

    % Align test data to match training neuron order, padding with NaNs if necessary
    localIDs = G(trainNeuronGlobalIDs, testDay);  % local IDs of training neurons on this test day

    % Build a reverse map from localID (as in G) to row index in Xtest
    XtestLocalIDs = G(G(:, testDay) > 0, testDay);  % list of local neuron IDs in Xtest
    XtestRowMap = containers.Map(XtestLocalIDs, 1:numel(XtestLocalIDs));  % localID → row in Xtest

    % Initialize aligned Xtest matrix
    X_aligned = nan(nNeuronsTrain, size(Xtest,2), size(Xtest,3));

    % Fill in aligned matrix using reverse lookup
    for i = 1:numel(trainNeuronGlobalIDs)
        lid = localIDs(i);
        if lid > 0 && isKey(XtestRowMap, lid)
            X_aligned(i,:,:) = Xtest(XtestRowMap(lid),:,:);
        end
    end

    Xtest = X_aligned;

    nanRows = all(all(isnan(Xtest),3),2);
    if any(nanRows)
        fprintf('⚠️ Removing %d test neurons with all-NaN data\n', sum(nanRows));
        Xtest(nanRows,:,:) = [];


        keepNeurons = ~nanRows;  % logical index of neurons to keep

    % Rebuild Xtrain_all to match only these neurons
    neuronMask = reshape(repmat(keepNeurons(:), 1, nBins)', [], 1);  % [neurons*bins x 1]
    Xtrain_all_test = Xtrain_all(:, neuronMask);

    % Also update Xtrain_time accordingly
    nNeuronsKeep = sum(keepNeurons);
    Xtrain_reshaped_test = reshape(Xtrain_all_test', nNeuronsKeep, nBins, []);



    end


    trialMask = squeeze(all(all(~isnan(Xtest),1),2));
    Xtest = Xtest(:,:,trialMask);
    ytest = ytest(trialMask);

    if size(Xtest, 3) == 0
        fprintf('⚠️ Skipping %s: all test trials removed due to NaNs.\n', testStr);
        continue;
    end

    Xtest_flat = reshape(Xtest, [], size(Xtest,3))';


    fprintf('Training size: %d trials × %d features\n', size(Xtrain_all_test,1), size(Xtrain_all_test,2));

    % Simple mean imputation across trials (column-wise)
    colMeans = nanmean(Xtrain_all_test, 1);
    nanMask = isnan(Xtrain_all_test);
    Xfill = repmat(colMeans, size(Xtrain_all_test,1), 1);
    Xtrain_all_test(nanMask) = Xfill(nanMask);



    mdl = fitcsvm(Xtrain_all_test, ytrain_all, 'KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
    yhat = predict(mdl, Xtest_flat);
    acc = mean(yhat == ytest);
    f1  = f1score(ytest, yhat);

  length(find(yhat==1))/numel(yhat) %average
    % Timepoint-flattened
      nNeurons = size(Xtest,1);


    Xtest_time = reshape(permute(Xtest, [1 3 2]), nNeurons, [])';
    ytest_time = repmat(ytest(:), nBins, 1);


    Xtrain_time_test = Xtrain_time_full;
    ytrain_time_test = ytrain_time_full;

    % Impute NaNs in training data using column-wise mean
trainMeans = nanmean(Xtrain_time_test, 1);
nanMaskTrain = isnan(Xtrain_time_test);
XfillTrain = repmat(trainMeans, size(Xtrain_time_test, 1), 1);
Xtrain_time_test(nanMaskTrain) = XfillTrain(nanMaskTrain);
% Drop only test rows with NaNs (e.g., full missing neuron)
validTest  = ~any(isnan(Xtest_time), 2);
Xtest_time = Xtest_time(validTest,:);
ytest_time = ytest_time(validTest);

nTrainCols = size(Xtrain_time_test, 2);
nTestCols = size(Xtest_time, 2);
if nTestCols < nTrainCols
    padSize = nTrainCols - nTestCols;
    Xtest_time = [Xtest_time, nan(size(Xtest_time,1), padSize)];
elseif nTestCols > nTrainCols
    Xtest_time = Xtest_time(:, 1:nTrainCols); % truncate to match
end



    fprintf('  Pre-prune timepoint size: Xtrain_time [%d×%d], Xtest_time [%d×%d]\n', ...
        size(Xtrain_time,1), size(Xtrain_time,2), size(Xtest_time,1), size(Xtest_time,2));
        fprintf('  Post-prune rows: %d train, %d test\n', size(Xtrain_time_test,1), sum(validTest));

    if isempty(Xtrain_time) || isempty(Xtest_time)
        fprintf('⚠️ Skipping %s: no valid timepoint-aligned data\n', testStr);
        continue;
    end


    mdl_time = fitcsvm(Xtrain_time_test, ytrain_time_test, 'KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
    yhat_time = predict(mdl_time, Xtest_time);
    acc_time = mean(yhat_time == ytest_time);
    f1_time  = f1score(ytest_time, yhat_time);

    yhat_time'

    length(find(yhat_time==1))/numel(yhat_time) %average

    % Permutation
    perm_acc = nan(nPerms,1);
    perm_f1  = nan(nPerms,1);
    perm_acc_time = nan(nPerms,1);
    perm_f1_time  = nan(nPerms,1);


    % Mean imputation for training trial-level data
    Xtrain_trial = Xtrain_all_test;  % This should match the test neuron subset
    colMeans_trial = nanmean(Xtrain_trial, 1);
    nanMask_trial = isnan(Xtrain_trial);
    Xfill_trial = repmat(colMeans_trial, size(Xtrain_trial, 1), 1);
    Xtrain_trial(nanMask_trial) = Xfill_trial(nanMask_trial);


    % Predefine fixed training and test data for permutations
    Xtrain_trial_perm   = Xtrain_all_test;
    ytrain_all_perm     = ytrain_all;

    Xtrain_time_perm    = Xtrain_time_test;
    ytrain_time_perm    = ytrain_time_test;

    Xtest_flat_perm     = Xtest_flat;
    Xtest_time_perm     = Xtest_time;
    ytest_perm          = ytest;
    ytest_time_perm     = ytest_time;


    if nPerms > 1 && isempty(gcp('nocreate'))
        parpool('local');
    end

    if nPerms > 1
      parfor p = 1:nPerms
        try
          % Shuffle labels ONLY (not rows of X)
          yshuf = ytrain_all(randperm(numel(ytrain_all)));
          mdl_shuf = fitcsvm(Xtrain_all_test, yshuf, 'KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
          yhat_shuf = predict(mdl_shuf, Xtest_flat);
          perm_acc(p) = mean(yhat_shuf == ytest);
          perm_f1(p)  = f1score(ytest, yhat_shuf);

          % Timepoint-level (shuffle only labels)
          yshuf_time = ytrain_time_test(randperm(numel(ytrain_time_test)));
          mdl_shuf_time = fitcsvm(Xtrain_time_test, yshuf_time, 'KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
          yhat_shuf_time = predict(mdl_shuf_time, Xtest_time);
          perm_acc_time(p) = mean(yhat_shuf_time == ytest_time);
          perm_f1_time(p)  = f1score(ytest_time, yhat_shuf_time);

        catch err
          warning('Permutation %d failed: %s', p, err.message);
          perm_acc(p) = NaN;
          perm_f1(p) = NaN;
          perm_acc_time(p) = NaN;
          perm_f1_time(p) = NaN;
        end
      end
    end


    fprintf('  Observed accuracy: %.2f\n', acc);
    fprintf('  Permutation acc mean: %.2f, max: %.2f\n', mean(perm_acc), max(perm_acc));

    fprintf('  Valid permutations: %d of %d\n', sum(~isnan(perm_acc)), nPerms);

    validPerms = ~isnan(perm_acc);
    % Include the observed value in the null to avoid p = 0
    pval_acc = (sum(perm_acc(validPerms) >= acc) + 1) / (sum(validPerms) + 1);
    pval_f1  = (sum(perm_f1(validPerms) >= f1) + 1) / (sum(validPerms) + 1);

    validPermsTime = ~isnan(perm_acc_time);
    pval_time_acc = (sum(perm_acc_time(validPermsTime) >= acc_time) + 1) / (sum(validPermsTime) + 1);
    pval_time_f1  = (sum(perm_f1_time(validPermsTime) >= f1_time) + 1) / (sum(validPermsTime) + 1);




    decodeResults(end+1) = struct('testDate', testStr, 'nSharedNeurons', sum(shared), ...
      'acc_trial', acc, 'f1_trial', f1, 'pval_acc', pval_acc, 'pval_f1', pval_f1, ...
      'acc_time', acc_time, 'f1_time', f1_time, 'pval_time_acc', pval_time_acc, 'pval_time_f1', pval_time_f1);

    fprintf('[Train: %s → Test: %s] Neurons: %d\n', strjoin(dateList(trainDays), ','), testStr, sum(shared));
    fprintf('  Label balance: %d correct, %d incorrect\n', sum(ytest == 1), sum(ytest == 0));
    fprintf('  Trial-flat:     Acc = %.2f%% (p=%.3f), F1 = %.3f (p=%.3f)\n', acc*100, mean(perm_acc >= acc, 'omitnan'), f1, mean(perm_f1 >= f1, 'omitnan'));
    fprintf('  Timepoint-flat: Acc = %.2f%% (p=%.3f), F1 = %.3f (p=%.3f)\n', acc_time*100, mean(perm_acc_time >= acc_time, 'omitnan'), f1_time, mean(perm_f1_time >= f1_time, 'omitnan'));
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
