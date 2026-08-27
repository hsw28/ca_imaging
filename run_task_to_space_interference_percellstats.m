function OUT = run_task_to_space_interference_percellstats(inputArg, varargin)
% run_task_to_space_interference_percellstats
%
% Quantifies per-cell Fig. 4g spatial-map similarity values without changing
% the Fig. 4g calculation code. The Pearson r values are extracted from:
%   R(i).C.SC.byDay.withTask.z{d}
% and converted back to r with tanh().
%
% Categories:
%   Preserved        r >= 0.5
%   Weakly preserved 0 <= r < 0.5
%   Reorganized      r < 0
%
% Usage:
%   OUT = run_task_to_space_interference_percellstats();
%   OUT = run_task_to_space_interference_percellstats({'rat0222','rat0307','rat0313','rat0314','rat0816'});
%   OUT = run_task_to_space_interference_percellstats(R);
%   OUT = run_task_to_space_interference_percellstats(ratNames, 'DoPlots', false);

defaultRats = {'rat0222','rat0307','rat0313','rat0314','rat0816'};

if nargin < 1 || isempty(inputArg)
    R = run_task_to_space_interference(defaultRats, varargin{:});
elseif isstruct(inputArg)
    R = inputArg;
elseif iscell(inputArg) || isstring(inputArg) || ischar(inputArg)
    R = run_task_to_space_interference(cellstr(inputArg), varargin{:});
else
    error('Input must be empty, a rat-name list, or an existing R struct from run_task_to_space_interference.');
end

[cellTable, perRat] = extractFig4gPerCellR(R);

allR = cellTable.r;
cats = {'Preserved', 'Weakly preserved', 'Reorganized'};
countsPooled = [sum(allR >= 0.5), sum(allR >= 0 & allR < 0.5), sum(allR < 0)];
if isempty(allR)
    pctPooled = [NaN NaN NaN];
else
    pctPooled = 100 .* countsPooled ./ numel(allR);
end

perRatCounts = vertcat(perRat.counts);
perRatPct = vertcat(perRat.percentages);
grandMeanPct = mean(perRatPct, 1, 'omitnan');
medianR = median(allR, 'omitnan');
iqrR = prctile(allR, [25 75]);

ratLabels = string({perRat.rat}).';
countTable = array2table(perRatCounts, ...
    'VariableNames', {'Preserved_n', 'WeaklyPreserved_n', 'Reorganized_n'});
countTable.Rat = ratLabels;
countTable.N = sum(perRatCounts, 2);
countTable = movevars(countTable, {'Rat','N'}, 'Before', 1);

pctTable = array2table(perRatPct, ...
    'VariableNames', {'Preserved_pct', 'WeaklyPreserved_pct', 'Reorganized_pct'});
pctTable.Rat = ratLabels;
pctTable = movevars(pctTable, 'Rat', 'Before', 1);

pooledTable = table( ...
    ["Pooled"; "Grand mean across rats"], ...
    [numel(allR); NaN], ...
    [countsPooled(1); NaN], ...
    [pctPooled(1); grandMeanPct(1)], ...
    [countsPooled(2); NaN], ...
    [pctPooled(2); grandMeanPct(2)], ...
    [countsPooled(3); NaN], ...
    [pctPooled(3); grandMeanPct(3)], ...
    'VariableNames', {'Group','N','Preserved_n','Preserved_pct', ...
                      'WeaklyPreserved_n','WeaklyPreserved_pct', ...
                      'Reorganized_n','Reorganized_pct'});

F = figure('Color', 'w', 'Name', 'Fig 4g per-cell spatial similarity categories');
tlo = tiledlayout(F, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axHist = nexttile(tlo, 1);
histogram(allR, 'BinLimits', [-1 1], 'BinWidth', 0.05, ...
    'Normalization', 'percentage', ...
    'FaceColor', [0.35 0.55 0.80], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
xline(0, 'k-', 'LineWidth', 1);
xline(0.5, 'k--', 'LineWidth', 1);
xlim([-1 1]);
xlabel('Fig. 4g per-cell spatial similarity (Pearson r)');
ylabel('Neurons (%)');
title('Per-cell tEBC-to-baseline spatial similarity');
set(gca, 'Box', 'off', 'LineWidth', 1);
drawGrandMeanInset(axHist, grandMeanPct);

nexttile(tlo, 2);
barData = [perRatPct; grandMeanPct];
b = bar(barData, 'stacked', 'EdgeColor', 'none');
b(1).FaceColor = [0.20 0.60 0.35];
b(2).FaceColor = [0.85 0.65 0.25];
b(3).FaceColor = [0.75 0.30 0.25];
ylim([0 100]);
ylabel('Neurons (%)');
xticks(1:(numel(perRat) + 1));
xticklabels([ratLabels; "Grand mean"]);
xtickangle(35);
legend(cats, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Spatial similarity categories by rat');
set(gca, 'Box', 'off', 'LineWidth', 1);

OUT = struct();
OUT.R = R;
OUT.figure = F;
OUT.cellTable = cellTable;
OUT.perRat = perRat;
OUT.countsByRat = countTable;
OUT.percentagesByRat = pctTable;
OUT.pooledAndGrandMean = pooledTable;
OUT.categoryThresholds = struct('Preserved', 'r >= 0.5', ...
                                'WeaklyPreserved', '0 <= r < 0.5', ...
                                'Reorganized', 'r < 0');
OUT.summary = struct('totalNeurons', numel(allR), ...
                     'medianPearsonR', medianR, ...
                     'iqrPearsonR', iqrR, ...
                     'percentReorganized_r_lt_0', pctPooled(3), ...
                     'percentPreserved_r_ge_0p5', pctPooled(1));

printPlainTextStats(perRat, grandMeanPct, numel(allR), medianR, iqrR, pctPooled);

end


function [cellTable, perRat] = extractFig4gPerCellR(R)
rat_all = strings(0,1);
day_all = strings(0,1);
dayIndex_all = [];
cellIndex_all = [];
z_all = [];
r_all = [];

perRat = repmat(struct('rat', '', 'r', [], 'counts', [], 'percentages', []), numel(R), 1);

for i = 1:numel(R)
    ratName = getRatName(R, i);
    zByDay = getWithTaskZByDay(R(i));
    dayNames = getDayNames(R, i, numel(zByDay));

    thisR = [];
    for d = 1:numel(zByDay)
        z = zByDay{d};
        if isempty(z), continue; end
        z = z(:);
        r = tanh(z);
        keep = isfinite(r);
        if ~any(keep), continue; end

        z = z(keep);
        r = r(keep);
        cellIdx = find(keep);

        n = numel(r);
        rat_all = [rat_all; repmat(string(ratName), n, 1)]; %#ok<AGROW>
        day_all = [day_all; repmat(string(dayNames{d}), n, 1)]; %#ok<AGROW>
        dayIndex_all = [dayIndex_all; d*ones(n,1)]; %#ok<AGROW>
        cellIndex_all = [cellIndex_all; cellIdx(:)]; %#ok<AGROW>
        z_all = [z_all; z(:)]; %#ok<AGROW>
        r_all = [r_all; r(:)]; %#ok<AGROW>
        thisR = [thisR; r(:)]; %#ok<AGROW>
    end

    counts = [sum(thisR >= 0.5), sum(thisR >= 0 & thisR < 0.5), sum(thisR < 0)];
    if isempty(thisR)
        pct = [NaN NaN NaN];
    else
        pct = 100 .* counts ./ numel(thisR);
    end

    perRat(i).rat = ratName;
    perRat(i).r = thisR;
    perRat(i).counts = counts;
    perRat(i).percentages = pct;
end

category = strings(numel(r_all), 1);
category(r_all >= 0.5) = "Preserved";
category(r_all >= 0 & r_all < 0.5) = "Weakly preserved";
category(r_all < 0) = "Reorganized";

cellTable = table(rat_all, day_all, dayIndex_all, cellIndex_all, z_all, r_all, category, ...
    'VariableNames', {'Rat','Day','DayIndex','CellIndex','z','r','Category'});
end


function drawGrandMeanInset(axHist, grandMeanPct)
colors = [0.20 0.60 0.35;
          0.85 0.65 0.25;
          0.75 0.30 0.25];
labels = {'Pres.', 'Weak', 'Reorg.'};

drawnow;
oldUnits = axHist.Units;
axHist.Units = 'normalized';
p = axHist.Position;
axHist.Units = oldUnits;

insetPos = [p(1) + 0.54*p(3), p(2) + 0.68*p(4), 0.38*p(3), 0.18*p(4)];
axInset = axes('Position', insetPos);
hold(axInset, 'on');
x0 = 0;
for k = 1:3
    w = grandMeanPct(k);
    if ~isfinite(w), w = 0; end
    rectangle(axInset, 'Position', [x0 0 w 1], ...
        'FaceColor', colors(k,:), 'EdgeColor', 'none');
    if w >= 12
        text(axInset, x0 + w/2, 0.5, sprintf('%.1f%%', grandMeanPct(k)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 8, 'Color', 'w', 'FontWeight', 'bold');
    end
    x0 = x0 + w;
end
xlim(axInset, [0 100]);
ylim(axInset, [0 1]);
set(axInset, 'XTick', [], 'YTick', [], 'Box', 'on', 'LineWidth', 0.75);
title(axInset, 'Grand mean', 'FontSize', 8, 'FontWeight', 'normal');

for k = 1:3
    text(axInset, 0, -0.32 - 0.28*(k-1), labels{k}, ...
        'Units', 'normalized', 'Color', colors(k,:), 'FontSize', 7, ...
        'FontWeight', 'bold');
end
hold(axInset, 'off');
axes(axHist);
end


function printPlainTextStats(perRat, grandMeanPct, totalN, medianR, iqrR, pctPooled)
fprintf('\nFig. 4g per-cell spatial similarity categories\n\n');

for i = 1:numel(perRat)
    ratLabel = regexprep(perRat(i).rat, '^rat', 'Rat ');
    c = perRat(i).counts;
    p = perRat(i).percentages;

    fprintf('%s:\n', ratLabel);
    fprintf('Preserved: %d neurons (%.1f%%)\n', c(1), p(1));
    fprintf('Weakly preserved: %d neurons (%.1f%%)\n', c(2), p(2));
    fprintf('Reorganized: %d neurons (%.1f%%)\n\n', c(3), p(3));
end

fprintf('Grand mean across rats (equal weighting across rats):\n\n');
fprintf('Preserved: %.1f%%\n', grandMeanPct(1));
fprintf('Weakly preserved: %.1f%%\n', grandMeanPct(2));
fprintf('Reorganized: %.1f%%\n\n', grandMeanPct(3));

fprintf('Total neurons analyzed: %d\n', totalN);
fprintf('Median Pearson r: %.4f\n', medianR);
fprintf('IQR: %.4f to %.4f\n', iqrR(1), iqrR(2));
fprintf('Percentage of neurons with r < 0: %.1f%%\n', pctPooled(3));
fprintf('Percentage of neurons with r >= 0.5: %.1f%%\n\n', pctPooled(1));
end


function zByDay = getWithTaskZByDay(Ri)
if isfield(Ri, 'C') && isfield(Ri.C, 'SC') && isfield(Ri.C.SC, 'byDay') && ...
        isfield(Ri.C.SC.byDay, 'withTask') && isfield(Ri.C.SC.byDay.withTask, 'z')
    zByDay = Ri.C.SC.byDay.withTask.z;
    return;
end

error(['Could not find Fig. 4g per-cell values at ' ...
       'R(i).C.SC.byDay.withTask.z. Re-run run_task_to_space_interference ' ...
       'with DoSingleCell=true.']);
end


function ratName = getRatName(R, i)
if isfield(R(i), 'animal') && ~isempty(R(i).animal)
    ratName = char(R(i).animal);
else
    ratName = sprintf('rat%d', i);
end
end


function dayNames = getDayNames(R, i, nDays)
if isfield(R(i), 'meta') && isfield(R(i).meta, 'days') && numel(R(i).meta.days) >= nDays
    dayNames = R(i).meta.days;
else
    dayNames = arrayfun(@(d) sprintf('day%d', d), 1:nDays, 'UniformOutput', false);
end
end
