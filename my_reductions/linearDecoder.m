function decodeResults = linearDecoder(animalName, nPerms)
  %input animal name and number of permutations to run shuffling.
  %if you dont put in number of perms shuffling will be skipped
  % ie  linearDecoder('rat0314')
  %     linearDecoder('rat0314', 500)    <-- this is slow
  % outputs accuracy, p score compared to shuffled, f1 score and f1 p value:
        % Trial-flat:     accuracy = 65.52% (shuff p=1.000), f1 = 0.000 (f1 p=1.000)
        % Timepoint-flat: accuracy = 65.52% (shuff p=1.000), f1 = 0.000 (f1 p=1.000)


  if nargin < 2
    nPerms = 1;
  end


  %win = [0, .93];
  %Fs = 7.5;

  win = [0, 1.3];
  Fs = 7.5;

  nBins = round((win(end)-win(1)) * 7.5);

  % Access the structured data from the base workspace
  animal = evalin('base', animalName);
  G = animal.alignmentALL;    % cell-to-day global index matrix

  % Automatically get list of date strings
  dateList = autoDateList(animal);
  nDays = numel(dateList);

  % Only analyze aligned sessions

  nAligned = size(G, 2);
  if nAligned < nDays
      fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
      dateList = dateList(1:nAligned);
      nDays = nAligned;
  end


    dateList = autoDateList(animal);
    decodeResults = struct('date', [], ...
        'acc_trial', [], 'f1_trial', [], ...
        'acc_time', [], 'f1_time', [], ...
        'perm_acc_trial', [], 'perm_f1_trial', [], ...
        'perm_acc_time', [], 'perm_f1_time', []);

  for d = 1:nDays
        dateStr = dateList{d};
        nBins = round((win(2) - win(1)) * Fs);
        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);

        if isempty(X) || numel(unique(y)) < 2
            fprintf('Skipping %s: insufficient data or only one class.\n', dateStr);
            continue;
        end

        % Clean valid trials
        size(y)
        trialMask = squeeze(all(all((X),1),2));

        X = X(:,:,trialMask);
        y = y(trialMask);
        
        valid = ~isnan(y);
        X = X(:,:,valid);
        y = y(valid);

        % === TRIAL-FLATTENED ===
        Xtrial = reshape(X, [], size(X,3))';  % trials × features
        accT = crossvalAccuracy(Xtrial, y);
        f1T  = crossvalF1(Xtrial, y);

        % === TIMEPOINT-FLATTENED ===
        [nC, nTime, nTrials] = size(X);
        Xtime = reshape(permute(X,[1 3 2]), nC, [])';  % (trials*time) × neurons
        ytime = repmat(y(:), nTime, 1);

        timeMask = ~any(isnan(Xtime), 2);
        Xtime = Xtime(timeMask, :);
        ytime = ytime(timeMask);

        accTime = crossvalAccuracy(Xtime, ytime);
        f1Time  = crossvalF1(Xtime, ytime);

        % === PERMUTATIONS ===
        perm_acc_trial = nan(nPerms,1);
        perm_f1_trial  = nan(nPerms,1);
        perm_acc_time  = nan(nPerms,1);
        perm_f1_time   = nan(nPerms,1);

        if nPerms > 1 && isempty(gcp('nocreate'))
            parpool('local');
        end

        if nPerms>1
          parfor p = 1:nPerms
              % Trial-level
              yshuf_trial = y(randperm(numel(y)));
              perm_acc_trial(p) = crossvalAccuracy(Xtrial, yshuf_trial);
              perm_f1_trial(p)  = crossvalF1(Xtrial, yshuf_trial);

              % Timepoint-level
              yshuf_time = ytime(randperm(numel(ytime)));
              perm_acc_time(p) = crossvalAccuracy(Xtime, yshuf_time);
              perm_f1_time(p)  = crossvalF1(Xtime, yshuf_time);
          end
        end

        % === Store Results ===
        decodeResults(end+1) = struct('date', dateStr, ...
            'acc_trial', accT, 'f1_trial', f1T, ...
            'acc_time', accTime, 'f1_time', f1Time, ...
            'perm_acc_trial', perm_acc_trial, 'perm_f1_trial', perm_f1_trial, ...
            'perm_acc_time', perm_acc_time, 'perm_f1_time', perm_f1_time);

            % === Empirical p-values ===
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

function acc = crossvalAccuracy(X, y)
    n = size(X, 1); yhat = nan(n,1);
    for i = 1:n
        train = setdiff(1:n,i);
        mdl = fitclinear(X(train,:), y(train), 'Learner','logistic', 'Regularization','lasso');
        yhat(i) = predict(mdl, X(i,:));
    end
    acc = mean(yhat == y);
end


function f1 = crossvalF1(X, y)
    n = size(X, 1); yhat = nan(n,1);
    for i = 1:n
        train = setdiff(1:n,i);
        mdl = fitclinear(X(train,:), y(train), 'Learner','logistic', 'Regularization','lasso');
        yhat(i) = predict(mdl, X(i,:));
    end
    f1 = f1score(y, yhat);
end


function f = f1score(ytrue, ypred)
    ytrue = double(ytrue(:));
    ypred = double(ypred(:));
    tp = sum((ytrue == 1) & (ypred == 1));
    fp = sum((ytrue == 0) & (ypred == 1));
    fn = sum((ytrue == 1) & (ypred == 0));
    prec = tp / (tp + fp + eps);
    rec  = tp / (tp + fn + eps);
    f = 2 * (prec * rec) / (prec + rec + eps);
end
