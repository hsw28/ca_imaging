function plotAUROC_byPC(ratName, refDaySpec, pcList)
% plotAUROC_byPC(ratName, pcList, refDaySpec)
% Runs PCA alignment + decoding for a given rat and plots AUROC curves
% across days using different numbers of principal components.
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

% Get rat struct from workspace
if ~evalin('base', sprintf('exist("%s", "var")', ratName))
    error('Variable "%s" not found in base workspace.', ratName);
end
animal = evalin('base', ratName);

% Preallocate AUROC matrix [nPCs × nDays]
nPCs = numel(pcList);
AUROCs = [];

% Store for legend
legLabels = cell(1,nPCs);


% Run for each PC count
for i = 1:nPCs
    latentDim = pcList(i);
    AUROC = runPCA_fromStruct(ratName, refDaySpec, latentDim);
    AUROCs(i,:) = AUROC;
    legLabels{i} = sprintf('%d PCs', latentDim);
end

% Plot
figure; hold on;
plot(AUROCs', '-o', 'LineWidth', 1.5);
legend(legLabels, 'Location','southeast');
xlabel('Day index'); ylabel('AUROC');
title(sprintf('Cross-day decoder performance (%s)', ratName), 'Interpreter','none');
ylim([0 1]); grid on;

if exist('dateList','var') && ~isempty(dateList)
    xticks(1:numel(dateList));
    xticklabels(dateList);
    xtickangle(45);
end
end
