function decodeResults = linearDecoder(animalName)


  %win = [0, .93];
  %Fs = 7.5;
  %nBins = 7;

  win = [0, 1.3];
  Fs = 7.5;
  nBins = 10;

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


    decodeResults = struct('date', [], 'acc_trial_flat', [], 'acc_time_flat', []);



    for d = 1:numel(dateList)
        dateStr = dateList{d};
        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
        if isempty(X) || numel(unique(y)) < 2
            fprintf('Skipping %s: Not enough usable trials or only one class.\n', dateStr);
            continue;
        end

        % Clean data: remove trials with NaNs
        trialMask = squeeze(all(all(~isnan(X), 1), 2)); % trials without NaNs
        X = X(:,:,trialMask);
        y = y(trialMask);

        % ---------- Trial-Flattened ----------
        X_trial_flat = reshape(X, [], size(X,3))'; % trials × (neurons*time)
        y_trial = y(:);

        acc_trial = crossvalAccuracy(X_trial_flat, y_trial);

        % ---------- Timepoint-Flattened ----------
        [nC, nTime, nTrials] = size(X);
        X_time_flat = reshape(permute(X, [1 3 2]), nC, [])';  % (trials*time) × neurons
        y_time = repmat(y(:), nTime, 1);

        % Remove timepoints with NaNs
        nanMask = any(isnan(X_time_flat), 2);
        X_time_flat(nanMask,:) = [];
        y_time(nanMask) = [];

        acc_time = crossvalAccuracy(X_time_flat, y_time);

        % Store
        decodeResults(end+1) = struct('date', dateStr, ...
                                      'acc_trial_flat', acc_trial, ...
                                      'acc_time_flat', acc_time);

        fprintf('[%s] Trial-flattened: %.2f%% | Timepoint-flattened: %.2f%%\n', ...
            dateStr, acc_trial*100, acc_time*100);
    end
end

function acc = crossvalAccuracy(X, y)
    n = size(X, 1);
    yhat = nan(n, 1);
    for i = 1:n
        Xtrain = X([1:i-1, i+1:end], :);
        ytrain = y([1:i-1, i+1:end]);
        Xtest  = X(i,:);

        model = fitclinear(Xtrain, ytrain, 'Learner', 'logistic', 'Regularization', 'lasso');
        yhat(i) = predict(model, Xtest);
    end
    acc = mean(yhat == y);
end
