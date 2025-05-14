function decodeResults = kernelDecoder_REF(animalName, trainDay, nPerms)
    % Kernel-based decoder using RBF kernel (nonlinear SVM)
    % Trains on one day and tests on all others using aligned neurons
    % Performs both trial-flattened and timepoint-flattened decoding

    if nargin < 3
      nPerms = 1;
    end

    win = [0, 1.3];
    Fs = 7.5;
    nBins = round((win(end) - win(1)) * Fs);

    animal = evalin('base', animalName);
    G = animal.alignmentALL;
    dateList = autoDateList(animal);
    nDays = min(numel(dateList), size(G, 2));

    if ischar(trainDay) && (strcmpi(trainDay, 'end') || strcmpi(trainDay, 'last'))
      trainDay = nDays;
    end

    decodeResults = struct('testDate', [], 'nSharedNeurons', [], ...
        'acc_trial', [], 'f1_trial', [], 'pval_acc', [], 'pval_f1', [], ...
        'acc_time', [], 'f1_time', [], 'pval_time_acc', [], 'pval_time_f1', []);

    trainStr = dateList{trainDay};

    for testDay = 1:nDays
      if testDay == trainDay
          continue;
      end

      [Xtrain, ytrain] = getDayMatrixFromStruct(animal, trainStr, win, nBins, Fs);
      testStr = dateList{testDay};
      [Xtest, ytest] = getDayMatrixFromStruct(animal, testStr, win, nBins, Fs);
      if isempty(Xtrain) || numel(unique(ytrain)) < 2 || isempty(Xtest) || numel(unique(ytest)) < 2
          fprintf('Skipping %s: insufficient data or only one class.\n', testStr);
          continue;
      end

      shared = G(:,trainDay) > 0 & G(:,testDay) > 0;
      idxTrain = G(shared, trainDay);
      idxTest  = G(shared, testDay);

      if numel(idxTrain) < 5
          fprintf('Skipping %s: only %d shared neurons.\n', testStr, numel(idxTrain));
          continue;
      end

      % === TRIAL-FLATTENED ===
      Xtrain_flat = reshape(Xtrain(idxTrain,:,:), [], size(Xtrain,3))';
      trialMaskTrain = all(~isnan(Xtrain_flat), 2);
      Xtrain_flat = Xtrain_flat(trialMaskTrain, :);
      ytrain = ytrain(trialMaskTrain);

      Xtest_flat = reshape(Xtest(idxTest,:,:), [], size(Xtest,3))';
      trialMaskTest = all(~isnan(Xtest_flat), 2);
      Xtest_flat = Xtest_flat(trialMaskTest, :);
      ytest_trial = ytest(trialMaskTest);

      mdl = fitcsvm(Xtrain_flat, ytrain, 'KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
      yhat = predict(mdl, Xtest_flat);
      accT = mean(yhat == ytest_trial);
      f1T = f1score(ytest_trial, yhat);

      % === TIMEPOINT-FLATTENED ===
      [nC, nBins, ~] = size(Xtrain);
      Xtrain_sub = Xtrain(idxTrain,:,:);
      [nC_train, ~, ~] = size(Xtrain_sub);
      Xtrain_time = reshape(permute(Xtrain_sub, [1 3 2]), nC_train, [])';
      ytrain_time = repmat(ytrain(:), nBins, 1);
      validTrain = ~any(isnan(Xtrain_time), 2);
      Xtrain_time = Xtrain_time(validTrain,:);
      ytrain_time = ytrain_time(validTrain);

      Xtest_sub = Xtest(idxTest,:,:);
      [nC_test, ~, ~] = size(Xtest_sub);
      Xtest_time = reshape(permute(Xtest_sub, [1 3 2]), nC_test, [])';
      ytest_time = repmat(ytest(:), nBins, 1);
      validTest = ~any(isnan(Xtest_time), 2);
      Xtest_time = Xtest_time(validTest,:);
      ytest_time = ytest_time(validTest);

      mdl_time = fitcsvm(Xtrain_time, ytrain_time, 'KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
      yhat_time = predict(mdl_time, Xtest_time);
      accTime = mean(yhat_time == ytest_time);
      f1Time = f1score(ytest_time, yhat_time);

      % === Permutation testing for trial-level ===
      perm_acc = nan(nPerms,1);
      perm_f1  = nan(nPerms,1);
      perm_acc_time = nan(nPerms,1);
      perm_f1_time  = nan(nPerms,1);

      if nPerms > 1 && isempty(gcp('nocreate'))
          parpool('local');
      end
      if nPerms > 1
          parfor p = 1:nPerms
              try
                % Trial-level shuffle
                yshuf = ytrain(randperm(numel(ytrain)));
                mdl_shuf = fitcsvm(Xtrain_flat, yshuf, 'KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
                yhat_shuf = predict(mdl_shuf, Xtest_flat);
                perm_acc(p) = mean(yhat_shuf == ytest_trial);
                perm_f1(p)  = f1score(ytest_trial, yhat_shuf);

                % Timepoint-level shuffle
                yshuf_time = ytrain_time(randperm(numel(ytrain_time)));
                mdl_shuf_time = fitcsvm(Xtrain_time, yshuf_time, 'KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
                yhat_shuf_time = predict(mdl_shuf_time, Xtest_time);
                perm_acc_time(p) = mean(yhat_shuf_time == ytest_time);
                perm_f1_time(p)  = f1score(ytest_time, yhat_shuf_time);
              catch
                perm_acc(p) = NaN;
                perm_f1(p)  = NaN;
                perm_acc_time(p) = NaN;
                perm_f1_time(p)  = NaN;
              end
          end
      end

      pval_acc = mean(perm_acc >= accT, 'omitnan');
      pval_f1  = mean(perm_f1  >= f1T, 'omitnan');
      pval_time_acc = mean(perm_acc_time >= accTime, 'omitnan');
      pval_time_f1  = mean(perm_f1_time  >= f1Time, 'omitnan');

      decodeResults(end+1) = struct('testDate', testStr, 'nSharedNeurons', numel(idxTrain), ...
          'acc_trial', accT, 'f1_trial', f1T, 'pval_acc', pval_acc, 'pval_f1', pval_f1, ...
          'acc_time', accTime, 'f1_time', f1Time, 'pval_time_acc', pval_time_acc, 'pval_time_f1', pval_time_f1);

          fprintf('[Train: %s → Test: %s]  Neurons: %d\n', trainStr, testStr, numel(idxTrain));
          fprintf('  Label balance: %d correct, %d incorrect\n', sum(ytest == 1), sum(ytest == 0));

      fprintf('  Trial-flat:     Acc = %.2f%% (p=%.3f), F1 = %.3f (p=%.3f)\n', accT*100, pval_acc, f1T, pval_f1);
      fprintf('  Timepoint-flat: Acc = %.2f%% (p=%.3f), F1 = %.3f (p=%.3f)\n', accTime*100, pval_time_acc, f1Time, pval_time_f1);
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
