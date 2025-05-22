function UMAPtrial_dayByDay(animalName, win, plot3D, balanceClasses)
  %trial compressed
    if nargin < 3, plot3D = false; end
    if nargin < 4, balanceClasses = false; end
    Fs = 7.5;
    animal = evalin('base', animalName);
    dateList = autoDateList(animal);
    G = animal.alignmentALL;
    nDays = numel(dateList);
    nAligned = size(G, 2);
    if nAligned < nDays
        fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
        dateList = dateList(1:nAligned);
        nDays = nAligned;
    end

    nBins = round(diff(win) * Fs);
    minNeurons = 10;

    nCols = 5;
    nRows = ceil(nDays / nCols);
    figure;
    plotIdx = 1;

    for d = 1:nDays
        dateStr = dateList{d};
        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);

        if isempty(X) || numel(unique(y)) < 2
            fprintf('Skipping %s: insufficient data or only one class.\n', dateStr);
            continue;
        end

        nNeurons = size(X,1);
        if nNeurons < minNeurons
            fprintf('Skipping %s: too few neurons (%d).\n', dateStr, nNeurons);
            continue;
        end

        nTrials = size(X,3);
        trialVecs = zeros(nTrials, nNeurons * nBins);
        for t = 1:nTrials
            mat = X(:,:,t);
            mat(isnan(mat)) = 0;
            trialVecs(t,:) = mat(:)';
        end

        trialVecs = zscore(trialVecs);
        y = y(:);

        try
            subplot(nRows, nCols, plotIdx);
            [embedding, ~] = run_umap(trialVecs, 'n_components', 2, 'n_neighbors', 15, 'min_dist', 0.3, 'metric', 'cosine', 'randomize', true, 'verbose', false);
            scatter(embedding(y==0,1), embedding(y==0,2), 30, 'r', 'filled');
            hold on;
            scatter(embedding(y==1,1), embedding(y==1,2), 30, 'b', 'filled');
            title(sprintf('%s (%d)', dateStr, nTrials));
            xlabel('UMAP 1'); ylabel('UMAP 2'); grid on;
            axis tight;
            plotIdx = plotIdx + 1;
        catch
            fprintf('UMAP failed on %s\n', dateStr);
            continue;
        end

        % SVM decoding in original and UMAP space
        rng(42);
        part = cvpartition(y, 'KFold', 5);
        accRaw = zeros(part.NumTestSets, 1);
        accUMAP = zeros(part.NumTestSets, 1);
        for i = 1:part.NumTestSets
            trainIdx = training(part,i);
            testIdx = test(part,i);
            opts = {'KernelFunction','linear'};
            if balanceClasses
                opts = [opts, {'ClassNames',[0;1]}];
            end
            mdlX = fitcsvm(trialVecs(trainIdx,:), y(trainIdx), opts{:});
            mdlU = fitcsvm(embedding(trainIdx,:), y(trainIdx), opts{:});
            accRaw(i) = mean(predict(mdlX, trialVecs(testIdx,:)) == y(testIdx));
            accUMAP(i) = mean(predict(mdlU, embedding(testIdx,:)) == y(testIdx));
        end

        fprintf('Day %s: Xmat Acc = %.2f%% | UMAP Acc = %.2f%%\n', dateStr, 100*mean(accRaw), 100*mean(accUMAP));
    end
end
