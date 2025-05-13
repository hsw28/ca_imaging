
function results = quantifyUMAPManifolds(Z_corr_aligned, Z_incorr_aligned, validDays)
% Quantify UMAP manifolds across days and trial outcomes
% For each day: measure separation between trial types, within-manifold structure, and alignment across days

    % Inputs:
    %   Z_corr_aligned: cell array of [latentDim x nCells] UMAP coords for correct trials
    %   Z_incorr_aligned: same for incorrect trials
    %   validDays: list of day indices that were successfully aligned

    nDays = numel(validDays);
    latentDim = size(Z_corr_aligned{validDays(1)},1);

    results = struct();
    results.centroidDist = nan(1, nDays);
    results.withinVar = struct('correct', nan(1,nDays), 'incorrect', nan(1,nDays));
    results.procrustesErr = nan(1, nDays);

    refDay = validDays(end);
    Zref = Z_corr_aligned{refDay}';

    for i = 1:nDays
        d = validDays(i);
        Zc = Z_corr_aligned{d}';
        Zi = Z_incorr_aligned{d}';

        % 1. Distance between correct and incorrect centroids
        mu_c = mean(Zc, 1);
        mu_i = mean(Zi, 1);
        results.centroidDist(i) = norm(mu_c - mu_i);

        % 2. Within-manifold variance (mean squared distance to centroid)
        results.withinVar.correct(i) = mean(vecnorm(Zc - mu_c, 2, 2).^2);
        results.withinVar.incorrect(i) = mean(vecnorm(Zi - mu_i, 2, 2).^2);

        % 3. Procrustes error vs reference day
        if d == refDay
            results.procrustesErr(i) = 0;
        else
            shared = min(size(Zref,1), size(Zc,1));
            err = procrustes(Zref(1:shared,:), Zc(1:shared,:), 'Scaling', false, 'Reflection', false);
            results.procrustesErr(i) = err;           
        end
    end

    % Plot summary metrics
    figure;
    subplot(3,1,1); plot(validDays, results.centroidDist, '-o'); ylabel('Centroid Distance'); title('Correct vs Incorrect Separation');
    subplot(3,1,2); plot(validDays, results.withinVar.correct, '-ob'); hold on;
    plot(validDays, results.withinVar.incorrect, '-or'); ylabel('Within-Manifold Variance'); legend('Correct','Incorrect');
    subplot(3,1,3); plot(validDays, results.procrustesErr, '-ok'); ylabel('Procrustes Err'); xlabel('Day'); title('Drift from Reference');
end
