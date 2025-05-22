function decodeResults = bayesDecoder(animalName, nPerms)
  %within day decoding
if nargin < 2
    nPerms = 1;
end

win = [0, .75];
Fs = 7.5;
nBins = round((win(end) - win(1)) * Fs);

animal = evalin('base', animalName);
G = animal.alignmentALL;
dateList = autoDateList(animal);
nDays = min(numel(dateList), size(G, 2));

decodeResults = struct('date', {}, ...
    'acc_trial', {}, 'f1_trial', {}, ...
    'acc_time', {}, 'f1_time', {}, ...
    'perm_acc_trial', {}, 'perm_f1_trial', {}, ...
    'perm_acc_time', {}, 'perm_f1_time', {});

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
    ytrial = y(valid_trial);

    % Remove zero variance features per class to avoid fitcnb errors
    Xtrial = removeZeroVarFeaturesPerClass(Xtrial, ytrial);

    [nC, nTime, ~] = size(X);
    Xtime = reshape(permute(X,[1 3 2]), nC, [])';
    ytime = repmat(y(:), nTime, 1);
    timeMask = ~any(isnan(Xtime), 2);
    Xtime = Xtime(timeMask, :);
    ytime = ytime(timeMask);

    Xtime = removeZeroVarFeaturesPerClass(Xtime, ytime);

    accT = crossvalBayesAccuracy(Xtrial, ytrial);
    f1T  = crossvalBayesF1(Xtrial, ytrial);

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
            yshuf_trial = ytrial(randperm(numel(ytrial)));
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

function Xclean = removeZeroVarFeaturesPerClass(X, y)
classes = unique(y);
colsToKeep = true(1, size(X,2));
for c = classes(:)'
    idx = y==c;
    classData = X(idx,:);
    zeroVarCols = var(classData, 0, 1) == 0;
    colsToKeep(zeroVarCols) = false;
end
Xclean = X(:, colsToKeep);
end


function colsToKeep = removeZeroVarFeatures(X, y)
    % Remove features with zero variance in any class to avoid Naive Bayes errors
    classes = unique(y);
    keep = true(1, size(X,2));
    for c = classes'
        idx = y == c;
        classVar = var(X(idx, :), 0, 1);
        keep = keep & (classVar > 0);
    end
    colsToKeep = find(keep);
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


function acc = crossvalBayesAccuracy(X, y, k)
    if nargin < 3, k = 5; end
    cvp = cvpartition(length(y), 'KFold', k);
    yhat = nan(size(y));
    for i = 1:k
        trainIdx = training(cvp, i);
        testIdx = test(cvp, i);

        Xtrain = X(trainIdx, :);
        ytrain = y(trainIdx);
        Xtest = X(testIdx, :);

        % Binarize or discretize features as counts for multinomial
        Xtrain_disc = discretizeFeaturesForMultinomial(Xtrain);
        Xtest_disc = discretizeFeaturesForMultinomial(Xtest);

        if numel(unique(ytrain)) < 2
            continue; % skip fold if only one class present
        end

        Xtrain_disc = discretizeFeaturesForMultinomial(Xtrain);
        Xtest_disc = discretizeFeaturesForMultinomial(Xtest);
        mdl = fitcnb(Xtrain_disc, ytrain, 'DistributionNames', 'mn');
        yhat(testIdx) = predict(mdl, Xtest_disc);

    end
    valid = ~isnan(yhat);
    acc = mean(yhat(valid) == y(valid));
end

function Xdisc = discretizeFeaturesForMultinomial(X)
    nBins = 10;
    Xdisc = zeros(size(X));
    for f = 1:size(X,2)
        colData = X(:, f);
        minVal = min(colData);
        maxVal = max(colData);
        if minVal == maxVal
            Xdisc(:, f) = 1;  % single bin if no variance
        else
            edges = linspace(minVal, maxVal, nBins+1);
            Xdisc(:, f) = discretize(colData, edges);
            Xdisc(isnan(Xdisc(:, f)), f) = nBins;  % assign max bin to NaNs
        end
    end
end



function f1 = crossvalBayesF1(X, y, k)
    if nargin < 3, k = 5; end
    cvp = cvpartition(length(y), 'KFold', k);
    yhat = nan(size(y));
    for i = 1:k
        trainIdx = training(cvp, i);
        testIdx = test(cvp, i);

        Xtrain = X(trainIdx, :);
        ytrain = y(trainIdx);
        Xtest = X(testIdx, :);

        varTrain = var(Xtrain, 0, 1);
        keepCols = varTrain > 1e-12;
        if all(~keepCols) || numel(unique(ytrain)) < 2
            continue;
        end

        Xtrain = Xtrain(:, keepCols);
        Xtest = Xtest(:, keepCols);


        Xtrain_disc = discretizeFeaturesForMultinomial(Xtrain);
        Xtest_disc = discretizeFeaturesForMultinomial(Xtest);
        mdl = fitcnb(Xtrain_disc, ytrain, 'DistributionNames', 'mn');
        yhat(testIdx) = predict(mdl, Xtest_disc);

    end
    valid = ~isnan(yhat);
    f1 = f1score(y(valid), yhat(valid));
end
