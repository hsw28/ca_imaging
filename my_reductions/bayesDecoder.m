function decodeResults = bayesDecoder(animalName, nPerms)
    %~~~~decodes its own day
  % Bayesian classifier using Naive Bayes on full-dimensional data
  % Input: animalName (string), nPerms (int) - number of permutations
  % Usage:
  %   decodeResults = bayesDecoder('rat0314', 500);

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

      Xtrial = reshape(X, [], size(X,3))';
      valid_trial = all(~isnan(Xtrial), 2);
      Xtrial = Xtrial(valid_trial, :);
      y = y(valid_trial);

      [nC, nTime, ~] = size(X);
      Xtime = reshape(permute(X,[1 3 2]), nC, [])';
      ytime = repmat(y(:), nTime, 1);
      timeMask = ~any(isnan(Xtime), 2);
      Xtime = Xtime(timeMask, :);
      ytime = ytime(timeMask);

      accT = crossvalBayesAccuracy(Xtrial, y);
      f1T  = crossvalBayesF1(Xtrial, y);

      accTime = crossvalBayesAccuracy(Xtime, ytime);
      f1Time  = crossvalBayesF1(Xtime, ytime);

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
              perm_acc_trial(p) = crossvalBayesAccuracy(Xtrial, yshuf_trial);
              perm_f1_trial(p)  = crossvalBayesF1(Xtrial, yshuf_trial);

              yshuf_time = ytime(randperm(numel(ytime)));
              perm_acc_time(p) = crossvalBayesAccuracy(Xtime, yshuf_time);
              perm_f1_time(p)  = crossvalBayesF1(Xtime, yshuf_time);
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

function acc = crossvalBayesAccuracy(X, y)
    n = size(X, 1); yhat = nan(n,1);
    for i = 1:n
        train = setdiff(1:n, i);
        mdl = fitcnb(X(train,:), y(train));
        yhat(i) = predict(mdl, X(i,:));
    end
    acc = mean(yhat == y);
end

function f1 = crossvalBayesF1(X, y)
    n = size(X, 1); yhat = nan(n,1);
    for i = 1:n
        train = setdiff(1:n, i);
        mdl = fitcnb(X(train,:), y(train));
        yhat(i) = predict(mdl, X(i,:));
    end
    f1 = f1score(y, yhat);
end
