function results = UMAPtrial_dayByDay(animalName, win, plot3D, balanceClasses, showPermPlots)
    %trial-flattened
%    👉 Split this day's trials into 5 folds:
%In each fold:
%Train on 80% of trials from this day
%Test on 20% of trials from this day
%Repeat for all 5 folds → average the accuracy.

    if nargin < 3, plot3D = false; end
    if nargin < 4, balanceClasses = false; end
    if nargin < 5, showPermPlots = true; end  % NEW: option to disable permutation plots

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

    % Initialize main figure for UMAP scatter subplots:
    figure;
    plotIdx = 1;

    % Initialize results matrix:
    % Columns: [Xmat Acc | Xmat p | UMAP Acc | UMAP p]
    results = nan(nDays, 4);

    for d = 1:nDays
        dateStr = dateList{d};
        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs);

        if isempty(X) || numel(unique(y(~isnan(y)))) < 2
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
            % UMAP embedding:
            [embedding, ~] = run_umap(trialVecs, 'n_components', 3, 'n_neighbors', 50, 'min_dist', 0.3, 'metric', 'euclidean', 'randomize', true, 'verbose', false);

            % UMAP scatter plot (in main figure):
            subplot(nRows, nCols, plotIdx);
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
        opts = {'KernelFunction','linear'};
        if balanceClasses
            opts = [opts, {'ClassNames',[0;1]}];
        end
        for i = 1:part.NumTestSets
            trainIdx = training(part,i);
            testIdx = test(part,i);
            mdlX = fitcsvm(trialVecs(trainIdx,:), y(trainIdx), opts{:});
            mdlU = fitcsvm(embedding(trainIdx,:), y(trainIdx), opts{:});
            accRaw(i) = mean(predict(mdlX, trialVecs(testIdx,:)) == y(testIdx));
            accUMAP(i) = mean(predict(mdlU, embedding(testIdx,:)) == y(testIdx));
        end

        % Final real accuracies:
        real_acc_Xmat = mean(accRaw);
        real_acc_UMAP = mean(accUMAP);

        fprintf('Day %s: Xmat Acc = %.2f%% | UMAP Acc = %.2f%%\n', dateStr, 100*real_acc_Xmat, 100*real_acc_UMAP);

        %% PERMUTATION TESTS:
        num_perms = 100;
        fprintf('Running %d permutations for %s...\n', num_perms, dateStr);

        % Perm test for Xmat:
        perm_acc_Xmat = zeros(num_perms,1);
        for p = 1:num_perms
            y_shuff = y(randperm(length(y)));
            accPerm = zeros(part.NumTestSets, 1);
            for i = 1:part.NumTestSets
                trainIdx = training(part,i);
                testIdx = test(part,i);
                mdl = fitcsvm(trialVecs(trainIdx,:), y_shuff(trainIdx), opts{:});
                accPerm(i) = mean(predict(mdl, trialVecs(testIdx,:)) == y_shuff(testIdx));
            end
            perm_acc_Xmat(p) = mean(accPerm);
        end
        pval_Xmat = (sum(perm_acc_Xmat >= real_acc_Xmat) + 1) / (num_perms + 1);

        % Perm test for UMAP:
        perm_acc_UMAP = zeros(num_perms,1);
        for p = 1:num_perms
            y_shuff = y(randperm(length(y)));
            accPerm = zeros(part.NumTestSets, 1);
            for i = 1:part.NumTestSets
                trainIdx = training(part,i);
                testIdx = test(part,i);
                mdl = fitcsvm(embedding(trainIdx,:), y_shuff(trainIdx), opts{:});
                accPerm(i) = mean(predict(mdl, embedding(testIdx,:)) == y_shuff(testIdx));
            end
            perm_acc_UMAP(p) = mean(accPerm);
        end
        pval_UMAP = (sum(perm_acc_UMAP >= real_acc_UMAP) + 1) / (num_perms + 1);

        % Report:
        fprintf('Permutation p-values: Xmat = %.5f | UMAP = %.5f\n', pval_Xmat, pval_UMAP);

        %% Save results:
        results(d,:) = [real_acc_Xmat, pval_Xmat, real_acc_UMAP, pval_UMAP];

        %% Optionally plot UMAP permutation histogram:
        %if showPermPlots
        %    figure;
        %    histogram(perm_acc_UMAP, 'Normalization','probability');
        %    hold on;
        %    yL = ylim;
        %    plot([real_acc_UMAP real_acc_UMAP], yL, 'r-', 'LineWidth',2);
        %    xlabel('Accuracy');
        %    ylabel('Probability');
        %    title(sprintf('%s UMAP Perm p=%.5f', dateStr, pval_UMAP));
        %end
    end

    % Final message:
    fprintf('==== All done! Returning results matrix: ====\n');
    disp(results);
end
