
function decodeResults = kernelDecoder_REFgrouped_binned(animalName, trainDays, nPerms, minTrainDays, minTestDays)
%time point flattened, multiple training days, bins time points into 3 bins

if nargin < 3, nPerms = 1; end
if nargin < 4, minTrainDays = numel(trainDays); end
if nargin < 5, minTestDays = 1; end

win = [0, 1.3];
Fs = 7.5;
nBins = round((win(end) - win(1)) * Fs);

animal = evalin('base', animalName);
G = animal.alignmentALL;
dateList = autoDateList(animal);
nDays = min(numel(dateList), size(G, 2));

decodeResults = struct('testDate', {}, 'nSharedNeurons', {}, ...
  'acc_trial', {}, 'f1_trial', {}, 'pval_acc', {}, 'pval_f1', {}, ...
  'acc_time', {}, 'f1_time', {}, 'pval_time_acc', {}, 'pval_time_f1', {});

sharedTrain = sum(G(:,trainDays) > 0, 2) >= minTrainDays;
trainNeuronGlobalIDs = find(sharedTrain);
nNeuronsTrain = numel(trainNeuronGlobalIDs);

Xtrain_all = [];
ytrain_all = [];
actualTrainDays = [];

for i = 1:numel(trainDays)
  d = trainDays(i);
  dateStr = dateList{d};
  [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
  if isempty(X) || numel(unique(y)) < 2, continue; end

  trainIDs = G(sharedTrain, d);
  validIdx = trainIDs > 0;
  localIDs = trainIDs(validIdx);
  X = X(localIDs,:,:) ;

  trialMask = squeeze(all(all(~isnan(X),1),2));
  X = X(:,:,trialMask);
  y = y(trialMask);
  valid = ~isnan(y);
  X = X(:,:,valid);
  y = y(valid);

  fprintf('  Day %s: using %d neurons for training\n', dateStr, sum(validIdx));

  X_aligned = nan(nNeuronsTrain, size(X,2), size(X,3));
  gid_map = G(trainNeuronGlobalIDs, d);
  good = gid_map > 0 & gid_map <= size(X,1);
  X_aligned(good,:,:) = X(gid_map(good),:,:);

  Xflat = reshape(X_aligned, [], size(X_aligned,3))';
  if isempty(Xtrain_all)
    Xtrain_all = Xflat;
  else
    if size(Xflat,2) ~= size(Xtrain_all,2)
      fprintf('\n⚠️ Inconsistent feature count on %s. Skipping day.\n', dateStr);
      continue;
    end
    Xtrain_all = [Xtrain_all; Xflat];
  end
  ytrain_all = [ytrain_all; y];
  actualTrainDays = [actualTrainDays, d];
end

if isempty(Xtrain_all)
  warning('No valid training data found. Exiting.');
  return;
end

nTrainTrials = size(Xtrain_all, 1);
Xtrain_time_tensor = permute(reshape(Xtrain_all', [], nBins, nTrainTrials), [1 3 2]);
nNeuronsTrain = size(Xtrain_time_tensor, 1);

win1 = 1:round(nBins/3);
win2 = round(nBins/3)+1:round(2*nBins/3);
win3 = round(2*nBins/3)+1:nBins;

Xtrain_summary = zeros(nTrainTrials, 3 * nNeuronsTrain);
for t = 1:nTrainTrials
  Z1 = mean(Xtrain_time_tensor(:, t, win1), 3);
  Z2 = mean(Xtrain_time_tensor(:, t, win2), 3);
  Z3 = mean(Xtrain_time_tensor(:, t, win3), 3);
  Xtrain_summary(t, :) = [Z1; Z2; Z3]';
end

nanMask_train = isnan(Xtrain_summary);
colMeans_train = nanmean(Xtrain_summary, 1);
Xfill_train = repmat(colMeans_train, size(Xtrain_summary, 1), 1);
Xtrain_summary(nanMask_train) = Xfill_train(nanMask_train);

ytrain_all_use = ytrain_all;

for testDay = 1:nDays
  testStr = dateList{testDay};
  [Xtest, ytest] = getDayMatrixFromStruct(animal, testStr, win, nBins, Fs);

  if isempty(Xtest), continue; end
  if numel(unique(ytest)) < 2, continue; end

  trialMask = squeeze(all(all(~isnan(Xtest),1),2));
  Xtest = Xtest(:,:,trialMask);
  ytest = ytest(trialMask);
  if isempty(ytest), continue; end

  Xtest_time_tensor = permute(Xtest, [1 3 2]);
  nNeuronsTest = size(Xtest_time_tensor, 1);
  nTestTrials = size(Xtest_time_tensor, 2);
  Xtest_summary = nan(nTestTrials, 3 * nNeuronsTrain);

  for t = 1:nTestTrials
    Z1 = mean(Xtest_time_tensor(:, t, win1), 3);
    Z2 = mean(Xtest_time_tensor(:, t, win2), 3);
    Z3 = mean(Xtest_time_tensor(:, t, win3), 3);
    xsum = [Z1; Z2; Z3]';
    xsum = padOrTrimToLength(xsum, size(Xtest_summary,2));
    Xtest_summary(t, :) = xsum;
  end

  nanMask_test = isnan(Xtest_summary);
  colMeans_test = nanmean(Xtest_summary, 1);
  colMeans_test(isnan(colMeans_test)) = 0;
  Xfill_test = repmat(colMeans_test, size(Xtest_summary, 1), 1);
  Xtest_summary(nanMask_test) = Xfill_test(nanMask_test);

  ytest_use = ytest;

  mdl = fitcsvm(Xtrain_summary, ytrain_all_use, 'KernelFunction','linear', 'Standardize', true);
  yhat = predict(mdl, Xtest_summary);
  acc = mean(yhat == ytest_use);
  f1 = f1score(ytest_use, yhat);
  av = length(find(yhat==1))/numel(yhat); %average
  if av==1
    fprintf('model guessed all 1s')
end

  perm_acc = nan(nPerms, 1);
  perm_f1  = nan(nPerms, 1);

  for p = 1:nPerms
    yshuf = ytrain_all_use(randperm(numel(ytrain_all_use)));
    mdl_shuf = fitcsvm(Xtrain_summary, yshuf, 'KernelFunction','linear', 'Standardize', true);
    yhat_shuf = predict(mdl_shuf, Xtest_summary);
    perm_acc(p) = mean(yhat_shuf == ytest_use);
    perm_f1(p)  = f1score(ytest_use, yhat_shuf);
  end

  pval_time_acc = (sum(perm_acc >= acc) + 1) / (nPerms + 1);
  pval_time_f1  = (sum(perm_f1 >= f1) + 1) / (nPerms + 1);

  fprintf('Test %s — Timepoint-flat: Acc = %.2f, F1 = %.2f, pAcc = %.3f, pF1 = %.3f\n', testStr, acc, f1, pval_time_acc, pval_time_f1);

  decodeResults(end+1) = struct('testDate', testStr, 'nSharedNeurons', nNeuronsTest, ...
    'acc_trial', NaN, 'f1_trial', NaN, 'pval_acc', NaN, 'pval_f1', NaN, ...
    'acc_time', acc, 'f1_time', f1, 'pval_time_acc', pval_time_acc, 'pval_time_f1', pval_time_f1);
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

function padded = padOrTrimToLength(vec, targetLen)
  if numel(vec) < targetLen
    padded = nan(1, targetLen);
    padded(1:numel(vec)) = vec;
  else
    padded = vec(1:targetLen);
  end
end
