function baselineResults =  linearDecoder_REF(animalName, trainDays, testDays)
    win = [0, .75];
    Fs = 7.5;
    nBins = round((win(end) - win(1)) * Fs);

    animal = evalin('base', animalName);
    G = animal.alignmentALL;
    dateList = autoDateList(animal);

    % Prepare training data (aggregate trainDays)
    Xtrain_all = [];
    ytrain_all = [];
    for i = 1:numel(trainDays)
        d = trainDays(i);
        dateStr = dateList{d};
        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
        if isempty(X) || numel(unique(y)) < 2, continue; end

        % Flatten to summary as in your decoder
        Xflat = reshape(X, size(X,1)*size(X,2), size(X,3))';
        Xtrain_all = [Xtrain_all; Xflat];
        ytrain_all = [ytrain_all; y];
    end

    % Prepare test data (aggregate testDays)
    Xtest_all = [];
    ytest_all = [];
    for i = 1:numel(testDays)
        d = testDays(i);
        dateStr = dateList{d};
        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
        if isempty(X) || numel(unique(y)) < 2, continue; end

        Xflat = reshape(X, size(X,1)*size(X,2), size(X,3))';
        Xtest_all = [Xtest_all; Xflat];
        ytest_all = [ytest_all; y];
    end

    % Train logistic regression
    mdl = fitclinear(Xtrain_all, ytrain_all, 'Learner', 'logistic', 'Regularization', 'ridge');

    % Predict test labels
    yhat = predict(mdl, Xtest_all);

    % Compute accuracy and F1 score
    acc = mean(yhat == ytest_all);
    f1 = f1score(ytest_all, yhat);

    fprintf('Logistic Regression Baseline — Acc: %.3f, F1: %.3f\n', acc, f1);

    baselineResults = struct('Accuracy', acc, 'F1Score', f1, 'Model', mdl);
end

function f1 = f1score(ytrue, ypred)
    ytrue = double(ytrue(:));
    ypred = double(ypred(:));
    tp = sum((ytrue == 1) & (ypred == 1));
    fp = sum((ytrue == 0) & (ypred == 1));
    fn = sum((ytrue == 1) & (ypred == 0));
    prec = tp / (tp + fp + eps);
    rec = tp / (tp + fn + eps);
    f1 = 2 * (prec * rec) / (prec + rec + eps);
end
