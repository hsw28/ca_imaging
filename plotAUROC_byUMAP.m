function plotAUROC_byUMAP(ratName, refDaySpec, pcList)
% plotAUROC_byPC(ratName, refDaySpec, pcList)
% Runs PCA alignment + decoding for a given rat and plots AUROC curves
% across days using different numbers of principal components.
% Also visualizes trial-level UMAPs colored by CR vs no-CR.
%
% Inputs:
%   ratName     - string, e.g. 'rat0314'
%   refDaySpec  - optional: index (e.g. 2) or 'last' (default)
%   pcList      - vector of integers, e.g. [2 4 6 8 10]

% Runs UMAP alignment + decoding for a given rat and plots AUROC curves
% across days using different numbers of components (UMAP dims).
% Also visualizes trial-level UMAPs colored by CR vs no-CR when 2D.
%
% Inputs:
%   ratName     - string, e.g. 'rat0314'
%   refDaySpec  - optional: index (e.g. 2) or 'last' (default)
%   pcList      - vector of integers, e.g. [2 4 6 8 10]

if nargin < 3
    refDaySpec = 'last';
end

if ~ischar(ratName) && ~isstring(ratName)
    error('First input must be rat name as string, e.g. "rat0314"');
end

% Get rat struct from base workspace
if ~evalin('base', sprintf('exist("%s", "var")', ratName))
    error('Variable "%s" not found in base workspace.', ratName);
end
animal = evalin('base', ratName);

% Preallocate AUROC matrix [nPCs × nDays]
nPCs = numel(pcList);
AUROCs = [];
legLabels = cell(1,nPCs);
Z_trials_all = cell(1,nPCs);
labels_all = cell(1,nPCs);

% Run for each UMAP dimensionality
for i = 1:nPCs
    latentDim = pcList(i);
    [AUROC, Z_trials, labels] = runUMAP_fromStruct(ratName, refDaySpec, latentDim);
    AUROCs(i,:) = AUROC;
    legLabels{i} = sprintf('%d D', latentDim);

    % Save trials for visualization
    if latentDim == 2
        Z_trials_all{i} = Z_trials;
        labels_all{i} = labels;
    end
end

% Plot AUROC vs Day
figure; hold on;
plot(AUROCs', '-o', 'LineWidth', 1.5);
legend(legLabels, 'Location','southeast');
xlabel('Day index'); ylabel('AUROC');
title(sprintf('Cross-day decoder performance (%s, UMAP)', ratName), 'Interpreter','none');
ylim([0 1]); grid on;
hline(0.5);

% Optional: Plot UMAP projections if 2D
if any(pcList == 2)
    idx = find(pcList == 2, 1);
    Z = Z_trials_all{idx};
    lbl = labels_all{idx};
    figure('Name', sprintf('%s: UMAP for CR vs no-CR', ratName));
    nDays = numel(Z);
    for d = 1:nDays
        if isempty(Z{d}) || isempty(lbl{d}), continue; end
        subplot(ceil(nDays/4), 4, d);
        scatter(Z{d}(1,:), Z{d}(2,:), 30, lbl{d}, 'filled');
        title(sprintf('Day %d', d));
        axis tight; axis square;
        colormap([0.6 0.6 0.6; 0 0.5 1]);
        set(gca,'XTick',[],'YTick',[]);
    end
    sgtitle(sprintf('%s — UMAP Projections by Trial Type (2D)', ratName), 'Interpreter','none');
end
end
