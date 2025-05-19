function baselineCheck(animalName, trainDays)
    win = [0, 1.3];
    Fs = 7.5;
    nBins = round((win(end) - win(1)) * Fs);

    animal = evalin('base', animalName);
    G = animal.alignmentALL;
    dateList = autoDateList(animal);

    unionNeurons = any(G(:, trainDays) > 0, 2);
    neuronIDs = find(unionNeurons);
    nNeurons = numel(neuronIDs);

    allX = [];
    ally = [];

    for d = trainDays
        [X, y] = getDayMatrixFromStruct(animal, dateList{d}, win, nBins, Fs);
        if isempty(X) || numel(unique(y)) < 2
            fprintf('Skipping day %s due to insufficient data or one class.\n', dateList{d});
            continue;
        end

        localIDs = G(neuronIDs, d);
        X_aligned = nan(nNeurons, size(X,2), size(X,3));
        validNeurons = localIDs > 0 & localIDs <= size(X,1);
        X_aligned(validNeurons,:,:) = X(localIDs(validNeurons),:,:);

        Xflat = reshape(X_aligned, nNeurons*nBins, size(X_aligned,3))';

        % Impute NaNs with neuron-wise mean across trials (column-wise)
        nanMask = isnan(Xflat);
        colMean = nanmean(Xflat, 1);
        colMean(isnan(colMean)) = 0; % replace NaN means with zero if any

        % Replace NaNs column-wise
        for c = 1:size(Xflat, 2)
            Xflat(nanMask(:, c), c) = colMean(c);
        end


        allX = [allX; Xflat];
        ally = [ally; y];
    end

    validLabels = ~isnan(ally);
    allX = allX(validLabels, :);
    ally = ally(validLabels);

    if isempty(allX)
        error('No data available after alignment and cleaning.');
    end

    mdl = fitclinear(allX, ally, 'Learner', 'logistic', 'Regularization', 'ridge');
    loss = manualCrossVal(allX, ally, 5);
    fprintf('Manual 5-fold CV misclassification rate: %.4f\n', loss);

end

function loss = manualCrossVal(X, y, k)
    if nargin < 3, k = 5; end
    n = size(X, 1);
    idx = crossvalind('Kfold', n, k);  % Generate fold indices
    losses = zeros(k,1);
    for fold = 1:k
        testIdx = (idx == fold);
        trainIdx = ~testIdx;
        Xtrain = X(trainIdx, :);
        ytrain = y(trainIdx);
        Xtest = X(testIdx, :);
        ytest = y(testIdx);

        mdl = fitclinear(Xtrain, ytrain, 'Learner', 'logistic', 'Regularization', 'ridge');
        yhat = predict(mdl, Xtest);
        losses(fold) = mean(yhat ~= ytest);  % misclassification rate
    end
    loss = mean(losses);
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
