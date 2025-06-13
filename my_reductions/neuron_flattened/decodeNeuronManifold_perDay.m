function [accRaw, pRaw, accUMAP, pUMAP] = decodeNeuronManifold_perDay(animalName, win, latentDim, n_neighbors, min_dist, metric, numPerms)

    if nargin < 7, numPerms = 100; end
    Fs = 7.5;

    animal = evalin('base', animalName);  % extract once, not inside loop
    dateList = autoDateList(animal);
    nDays = numel(dateList);

    G = animal.alignmentALL;
    alignment_data = G;
    nAligned = size(G, 2);
    if nAligned < nDays
        fprintf('Only %d days aligned. Will stop analysis at that point.\n', nAligned);
        dateList = dateList(1:nAligned);
        nDays = nAligned;
    end

    accRaw = nan(1,nDays);
    pRaw   = nan(1,nDays);
    accUMAP = nan(1,nDays);
    pUMAP   = nan(1,nDays);

    for d = 1:nDays
        dateStr = dateList{d};
        if isBadDay(animal, dateStr), continue; end

        [X, y] = getDayMatrixFromStruct(animal, dateStr, win, round(diff(win)*Fs), Fs);
        if isempty(X), continue; end

        correct_trials = y == 1;
        incorrect_trials = y == 0;
        if sum(correct_trials) < 2 || sum(incorrect_trials) < 2, continue; end

        X_corr = squeeze(nanmean(X(:,:,correct_trials), 3));
        X_incorr = squeeze(nanmean(X(:,:,incorrect_trials), 3));
        X_corr = zscore(X_corr, 0, 2);
        X_incorr = zscore(X_incorr, 0, 2);

        valid_corr = all(isfinite(X_corr), 2);
        valid_incorr = all(isfinite(X_incorr), 2);
        X_corr_valid = X_corr(valid_corr, :);
        X_incorr_valid = X_incorr(valid_incorr, :);

        Xcat = double([X_corr_valid; X_incorr_valid]);
        labels = [ones(size(X_corr_valid,1),1); zeros(size(X_incorr_valid,1),1)];
        part = cvpartition(labels, 'KFold', 5);

        % Raw decoding
        accs = nan(part.NumTestSets,1);
        for i = 1:part.NumTestSets
            tr = training(part,i); te = test(part,i);
            mdl = fitcsvm(Xcat(tr,:), labels(tr), 'KernelFunction','linear');
            accs(i) = mean(predict(mdl, Xcat(te,:)) == labels(te));
        end
        accRaw(d) = mean(accs);

        % Raw permutation
        accPerm = nan(numPerms,1);
        parfor p = 1:numPerms
            shuffLabels = labels(randperm(length(labels)));
            accs_p = nan(part.NumTestSets,1);
            for i = 1:part.NumTestSets
                tr = training(part,i); te = test(part,i);
                mdl = fitcsvm(Xcat(tr,:), shuffLabels(tr), 'KernelFunction','linear');
                accs_p(i) = mean(predict(mdl, Xcat(te,:)) == shuffLabels(te));
            end
            accPerm(p) = mean(accs_p);
        end
        pRaw(d) = (sum(accPerm >= accRaw(d)) + 1) / (numPerms + 1);

        % UMAP decoding
        [embedding, ~] = run_umap(Xcat, 'n_components', latentDim, ...
            'n_neighbors', n_neighbors, 'min_dist', min_dist, ...
            'metric', metric, 'verbose', false, 'cluster_output', 'none');

        accs = nan(part.NumTestSets,1);
        for i = 1:part.NumTestSets
            tr = training(part,i); te = test(part,i);
            mdl = fitcsvm(embedding(tr,:), labels(tr), 'KernelFunction','linear');
            accs(i) = mean(predict(mdl, embedding(te,:)) == labels(te));
        end
        accUMAP(d) = mean(accs);

        % UMAP permutation
        accPerm = nan(numPerms,1);
        for p = 1:numPerms
            shuffLabels = labels(randperm(length(labels)));
            accs_p = nan(part.NumTestSets,1);
            for i = 1:part.NumTestSets
                tr = training(part,i); te = test(part,i);
                mdl = fitcsvm(embedding(tr,:), shuffLabels(tr), 'KernelFunction','linear');
                accs_p(i) = mean(predict(mdl, embedding(te,:)) == shuffLabels(te));
            end
            accPerm(p) = mean(accs_p);
        end
        pUMAP(d) = (sum(accPerm >= accUMAP(d)) + 1) / (numPerms + 1);
    end
end
