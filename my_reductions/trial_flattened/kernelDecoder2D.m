function decodeResults = kernelDecoder2D(animalName, nPerms)
    %~~~~decodes its own day
  % 2D kernel-based decoder using RBF kernel (nonlinear SVM)
  % Input: animalName (string), nPerms (int) - number of permutations
  % Usage:
  %   decodeResults = kernelDecoder('rat0314', 500);

%classifying in a 2D projected space
%Low-dimensional exploratory analysis
%Visualizing decision boundaries
%Seeing whether separability is present in the first two PCs



  if nargin < 2
    nPerms = 1;
  end

  win = [0, 1.3];
  Fs = 7.5;
  nBins = round((win(end) - win(1)) * Fs);

  animal = evalin('base', animalName);
  G = animal.alignmentALL;
  dateList = autoDateList(animal);
  nDays = min(numel(dateList), size(G, 2));

  decodeResults = struct('date', [], ...
      'acc_trial', [], 'f1_trial', [], ...
      'acc_time', [], 'f1_time', [], ...
      'perm_acc_trial', [], 'perm_f1_trial', [], ...
      'perm_acc_time', [], 'perm_f1_time', []);

  for d = 1:nDays
      dateStr = dateList{d};
      [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);

      if isempty(X) || numel(unique(y)) < 2
          fprintf('Skipping %s: insufficient data or only one class.\n', dateStr);
          continue;
      end

      trialMask = squeeze(all(all(~isnan(X),1),2));
      X = X(:,:,trialMask);
      y = y(trialMask);
      valid = ~isnan(y);
      X = X(:,:,valid);
      y = y(valid);

      Xtrial_raw = reshape(X, [], size(X,3))';
      [coeff_trial, Xtrial] = pca(Xtrial_raw, 'NumComponents', 2);
      y = y(~any(isnan(Xtrial),2));
      Xtrial = Xtrial(~any(isnan(Xtrial),2), :);

      [nC, nTime, nTrials] = size(X);
      Xtime_raw = reshape(permute(X,[1 3 2]), nC, [])';
      ytime = repmat(y(:), nTime, 1);
      [coeff_time, Xtime] = pca(Xtime_raw, 'NumComponents', 2);
      ytime = ytime(~any(isnan(Xtime),2));
      Xtime = Xtime(~any(isnan(Xtime),2), :);

      accT = crossvalKernelAccuracy(Xtrial, y);
      f1T  = crossvalKernelF1(Xtrial, y);

      accTime = crossvalKernelAccuracy(Xtime, ytime);
      f1Time  = crossvalKernelF1(Xtime, ytime);

      perm_acc_trial = nan(nPerms,1);
      perm_f1_trial  = nan(nPerms,1);
      perm_acc_time  = nan(nPerms,1);
      perm_f1_time   = nan(nPerms,1);

      if nPerms > 1 && isempty(gcp('nocreate'))
          parpool('local');
      end

      if nPerms > 1
          parfor p = 1:nPerms
              yshuf_trial = y(randperm(numel(y)));
              perm_acc_trial(p) = crossvalKernelAccuracy(Xtrial, yshuf_trial);
              perm_f1_trial(p)  = crossvalKernelF1(Xtrial, yshuf_trial);

              yshuf_time = ytime(randperm(numel(ytime)));
              perm_acc_time(p) = crossvalKernelAccuracy(Xtime, yshuf_time);
              perm_f1_time(p)  = crossvalKernelF1(Xtime, yshuf_time);
          end
      end

      decodeResults(end+1) = struct('date', dateStr, ...
          'acc_trial', accT, 'f1_trial', f1T, ...
          'acc_time', accTime, 'f1_time', f1Time, ...
          'perm_acc_trial', perm_acc_trial, 'perm_f1_trial', perm_f1_trial, ...
          'perm_acc_time', perm_acc_time, 'perm_f1_time', perm_f1_time);

      pval_acc_trial = mean(perm_acc_trial >= accT);
      pval_f1_trial  = mean(perm_f1_trial  >= f1T);
      pval_acc_time  = mean(perm_acc_time  >= accTime);
      pval_f1_time   = mean(perm_f1_time   >= f1Time);

      fprintf('[%s]\n', dateStr);
      fprintf('  Label balance: %d correct, %d incorrect\n', sum(y==1), sum(y==0));
      fprintf('  Trial-flat:     accuracy = %.2f%% (shuff p=%.3f), f1 = %.3f (f1 p=%.3f)\n', ...
          accT*100, pval_acc_trial, f1T, pval_f1_trial);
      fprintf('  Timepoint-flat: accuracy = %.2f%% (shuff p=%.3f), f1 = %.3f (f1 p=%.3f)\n', ...
          accTime*100, pval_acc_time, f1Time, pval_f1_time);
  end
end

function acc = crossvalKernelAccuracy(X, y)
    n = size(X, 1); yhat = nan(n,1);
    for i = 1:n
        train = setdiff(1:n, i);
        mdl = fitcsvm(X(train,:), y(train), 'KernelFunction','rbf', ...
            'KernelScale','auto', 'Standardize',true);
        yhat(i) = predict(mdl, X(i,:));
    end
    acc = mean(yhat == y);
end

function f1 = crossvalKernelF1(X, y)
    n = size(X, 1); yhat = nan(n,1);
    for i = 1:n
        train = setdiff(1:n, i);
        mdl = fitcsvm(X(train,:), y(train), 'KernelFunction','rbf', ...
            'KernelScale','auto', 'Standardize',true);
        yhat(i) = predict(mdl, X(i,:));
    end
    f1 = f1score(y, yhat);
end
