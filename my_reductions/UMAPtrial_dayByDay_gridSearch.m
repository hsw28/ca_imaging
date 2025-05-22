function [accMatrix, paramList] = UMAPtrial_dayByDay_gridSearch(animalName, win, plot3D, balanceClasses, componentsList, neighborsList, minDistList, metricList)
    if nargin < 3, plot3D = false; end
    if nargin < 4, balanceClasses = false; end
    if nargin < 5, componentsList = [2, 3, 5]; end
    if nargin < 6, neighborsList = [10, 15, 30]; end
    if nargin < 7, minDistList = [0.1, 0.3, 0.5]; end
    if nargin < 8, metricList = {'cosine', 'euclidean'}; end

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

    accMatrix = [];
    paramList = [];

    for nc = componentsList
        for nn = neighborsList
            for md = minDistList
                for mt = 1:length(metricList)
                    metric = metricList{mt};
                    fprintf('\n=== Params: n_comp=%d, n_neighbors=%d, min_dist=%.2f, metric=%s ===\n', nc, nn, md, metric);
                    accUMAP_allDays = nan(1, nDays); % Collect UMAP accuracy per day


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
                            [embedding, ~] = run_umap(trialVecs, 'n_components', nc, 'n_neighbors', nn, ...
                                'min_dist', md, 'metric', metric, 'randomize', true, 'verbose', false);
                        catch
                            fprintf('UMAP failed on %s\n', dateStr);
                            continue;
                        end
                        % SVM decoding in UMAP space
                        rng(42);
                        part = cvpartition(y, 'KFold', 5);
                        accUMAP = zeros(part.NumTestSets, 1);
                        for i = 1:part.NumTestSets
                            trainIdx = training(part,i);
                            testIdx = test(part,i);
                            opts = {'KernelFunction','linear'};
                            if balanceClasses
                                opts = [opts, {'ClassNames',[0;1]}];
                            end
                            mdlU = fitcsvm(embedding(trainIdx,:), y(trainIdx), opts{:});
                            accUMAP(i) = mean(predict(mdlU, embedding(testIdx,:)) == y(testIdx));
                        end
                        accUMAP_allDays(d) = mean(accUMAP);
                        fprintf('Day %s: UMAP Acc = %.2f%%\n', dateStr, 100*mean(accUMAP));
                    end
                    accMatrix = [accMatrix; accUMAP_allDays];
                    paramList = [paramList; nc, nn, md, mt];  % mt = index into metricList
                end
            end
        end
    end
end
