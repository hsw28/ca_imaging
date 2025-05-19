function resultsTable = testTrainingDayCombos(animalName, days, nPerms, minTrainDays, minTestDays)
resultsTable = table('Size',[0 4], 'VariableTypes',{'string','double','double','double'}, ...
    'VariableNames', {'TrainDays', 'Loss', 'Accuracy', 'F1'});

combos = {};
for k = 1:numel(days)
    c = nchoosek(days,k);
    for i = 1:size(c,1)
        combos{end+1} = c(i,:);
    end
end

for i = 1:length(combos)
    trainSet = combos{i};
    fprintf('Testing training days: %s\n', mat2str(trainSet));
    try
        [loss, accuracy, f1] = baselineCheck(animalName, trainSet);
        resultsTable = [resultsTable; {mat2str(trainSet), loss, accuracy, f1}];
    catch ME
        warning('Failed on days %s: %s', mat2str(trainSet), ME.message);
        resultsTable = [resultsTable; {mat2str(trainSet), NaN, NaN, NaN}];
    end
end
end





function [allX, ally] = getDataForDays(animalName, trainDays, minTrainDays)
% Helper to load and concatenate data for given days, aligning neurons across days

animal = evalin('base', animalName);
G = animal.alignmentALL;
dateList = autoDateList(animal);
win = [0, 1.3];
Fs = 7.5;
nBins = round((win(end) - win(1)) * Fs);

sharedTrain = sum(G(:,trainDays) > 0, 2) >= minTrainDays;
trainNeuronGlobalIDs = find(sharedTrain);

allX = [];
ally = [];

for d = trainDays
    dateStr = dateList{d};
    [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);
    if isempty(X) || numel(unique(y)) < 2
        continue;
    end

    trainIDs = G(sharedTrain, d);
    validIdx = trainIDs > 0;

    X = X(trainIDs(validIdx),:,:);
    trialMask = squeeze(all(all(~isnan(X),1),2));
    X = X(:,:,trialMask);
    y = y(trialMask);
    valid = ~isnan(y);
    X = X(:,:,valid);
    y = y(valid);

    % Align to full neuron set
    nNeuronsTrain = numel(trainNeuronGlobalIDs);
    X_aligned = nan(nNeuronsTrain, size(X,2), size(X,3));
    localIDs = G(trainNeuronGlobalIDs, d);
    validAlign = localIDs > 0 & localIDs <= size(X,1);
    X_aligned(validAlign,:,:) = X(localIDs(validAlign),:,:);

    Xflat = reshape(X_aligned, [], size(X_aligned,3))';
    allX = [allX; Xflat];
    ally = [ally; y];
end

end
