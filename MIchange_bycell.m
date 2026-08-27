function OUT = MIchange_bycell(RatList, tol)
% MIchange_bycell  Neuron-level change in spatial MI after excluding tEBC.
%
% Uses the same already-computed Fig. 4a MI arrays and inclusion mask:
%   WITH tEBC    -> rat.MI_wCSUS_all
%   WITHOUT tEBC -> rat.MI_noCSUS15_all
%   SPEED MATCHED -> rat.MI_noCSUS15_controlSpeed_all(:,2)
%   include      -> rat.ratemask_all == 1
%
% If the pooled *_all fields are not present, this function builds the same
% vectors from the last 3 days ending at rat.An, matching the helper logic
% used elsewhere for Fig. 4a summaries. It does not recompute MI.
%
% Main deltas:
%   delta_MI             = MI_without_tEBC - MI_with_tEBC
%   delta_MI_speedmatch  = MI_without_tEBC - MI_without_tEBC_speed_matched
%
% Usage:
%   OUT = MIchange_bycell();
%   OUT = MIchange_bycell({'rat0222','rat0307','rat0313','rat0314','rat0816'});
%   OUT = MIchange_bycell([], 1e-12);

if nargin < 1 || isempty(RatList)
    RatList = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end
if nargin < 2 || isempty(tol)
    tol = 1e-12;
end

categories = {'Increased spatial information', 'No change', 'Decreased spatial information'};
shortCats = {'Increased', 'No change', 'Decreased'};
nRats = numel(RatList);

perRat = repmat(struct( ...
    'rat', '', ...
    'days', {{}}, ...
    'withMI', [], ...
    'withoutMI', [], ...
    'speedMatchedMI', [], ...
    'deltaMI', [], ...
    'deltaMI_speedmatch', [], ...
    'counts', [], ...
    'percentages', [], ...
    'counts_speedmatch', [], ...
    'percentages_speedmatch', []), nRats, 1);

allDelta = [];
allRat = {};
allDeltaSpeed = [];
allRatSpeed = {};

for r = 1:nRats
    ratName = RatList{r};
    if ~evalin('base', sprintf('exist(''%s'', ''var'')', ratName))
        error('Variable %s not found in base workspace.', ratName);
    end

    rat = evalin('base', ratName);
    [withMI, withoutMI, speedMatchedMI, keepMask, daysUsed] = getFig4aVectors(rat);

    withMI = withMI(:);
    withoutMI = withoutMI(:);
    speedMatchedMI = speedMatchedMI(:);
    keepMask = keepMask(:);

    if numel(withMI) ~= numel(withoutMI) || numel(withMI) ~= numel(speedMatchedMI) || ...
            numel(withMI) ~= numel(keepMask)
        error(['%s has mismatched vector lengths: MI_wCSUS=%d, ' ...
               'MI_noCSUS15=%d, MI_noCSUS15_controlSpeed=%d, ratemask=%d.'], ...
              ratName, numel(withMI), numel(withoutMI), numel(speedMatchedMI), numel(keepMask));
    end

    valid = keepMask == 1 & isfinite(withMI) & isfinite(withoutMI);
    withUse = withMI(valid);
    withoutUse = withoutMI(valid);
    delta = withoutUse - withUse;

    validSpeed = keepMask == 1 & isfinite(withoutMI) & isfinite(speedMatchedMI);
    withoutUseSpeed = withoutMI(validSpeed);
    speedMatchedUse = speedMatchedMI(validSpeed);
    deltaSpeed = withoutUseSpeed - speedMatchedUse;

    counts = [sum(delta > tol), sum(abs(delta) <= tol), sum(delta < -tol)];
    if isempty(delta)
        pct = [NaN NaN NaN];
    else
        pct = 100 .* counts ./ numel(delta);
    end

    countsSpeed = [sum(deltaSpeed > tol), sum(abs(deltaSpeed) <= tol), sum(deltaSpeed < -tol)];
    if isempty(deltaSpeed)
        pctSpeed = [NaN NaN NaN];
    else
        pctSpeed = 100 .* countsSpeed ./ numel(deltaSpeed);
    end

    perRat(r).rat = ratName;
    perRat(r).days = daysUsed;
    perRat(r).withMI = withUse;
    perRat(r).withoutMI = withoutUse;
    perRat(r).speedMatchedMI = speedMatchedUse;
    perRat(r).deltaMI = delta;
    perRat(r).deltaMI_speedmatch = deltaSpeed;
    perRat(r).counts = counts;
    perRat(r).percentages = pct;
    perRat(r).counts_speedmatch = countsSpeed;
    perRat(r).percentages_speedmatch = pctSpeed;

    allDelta = [allDelta; delta(:)]; %#ok<AGROW>
    allRat = [allRat; repmat({ratName}, numel(delta), 1)]; %#ok<AGROW>
    allDeltaSpeed = [allDeltaSpeed; deltaSpeed(:)]; %#ok<AGROW>
    allRatSpeed = [allRatSpeed; repmat({ratName}, numel(deltaSpeed), 1)]; %#ok<AGROW>
end

pooledCounts = [sum(allDelta > tol), sum(abs(allDelta) <= tol), sum(allDelta < -tol)];
if isempty(allDelta)
    pooledPct = [NaN NaN NaN];
else
    pooledPct = 100 .* pooledCounts ./ numel(allDelta);
end

perRatCounts = vertcat(perRat.counts);
perRatPct = vertcat(perRat.percentages);
grandMeanPct = mean(perRatPct, 1, 'omitnan');

pooledCountsSpeed = [sum(allDeltaSpeed > tol), sum(abs(allDeltaSpeed) <= tol), sum(allDeltaSpeed < -tol)];
if isempty(allDeltaSpeed)
    pooledPctSpeed = [NaN NaN NaN];
else
    pooledPctSpeed = 100 .* pooledCountsSpeed ./ numel(allDeltaSpeed);
end

perRatCountsSpeed = vertcat(perRat.counts_speedmatch);
perRatPctSpeed = vertcat(perRat.percentages_speedmatch);
grandMeanPctSpeed = mean(perRatPctSpeed, 1, 'omitnan');

ratLabels = strrep(RatList(:), 'rat', 'Rat ');
countTable = array2table(perRatCounts, ...
    'VariableNames', {'Increased_n', 'NoChange_n', 'Decreased_n'});
countTable.Rat = ratLabels;
countTable.N = sum(perRatCounts, 2);
countTable = movevars(countTable, {'Rat','N'}, 'Before', 1);

pctTable = array2table(perRatPct, ...
    'VariableNames', {'Increased_pct', 'NoChange_pct', 'Decreased_pct'});
pctTable.Rat = ratLabels;
pctTable = movevars(pctTable, 'Rat', 'Before', 1);

countTableSpeed = array2table(perRatCountsSpeed, ...
    'VariableNames', {'Increased_n', 'NoChange_n', 'Decreased_n'});
countTableSpeed.Rat = ratLabels;
countTableSpeed.N = sum(perRatCountsSpeed, 2);
countTableSpeed = movevars(countTableSpeed, {'Rat','N'}, 'Before', 1);

pctTableSpeed = array2table(perRatPctSpeed, ...
    'VariableNames', {'Increased_pct', 'NoChange_pct', 'Decreased_pct'});
pctTableSpeed.Rat = ratLabels;
pctTableSpeed = movevars(pctTableSpeed, 'Rat', 'Before', 1);

pooledTable = table( ...
    ["Pooled"; "Grand mean across rats"], ...
    [numel(allDelta); NaN], ...
    [pooledCounts(1); NaN], ...
    [pooledPct(1); grandMeanPct(1)], ...
    [pooledCounts(2); NaN], ...
    [pooledPct(2); grandMeanPct(2)], ...
    [pooledCounts(3); NaN], ...
    [pooledPct(3); grandMeanPct(3)], ...
    'VariableNames', {'Group','N','Increased_n','Increased_pct', ...
                      'NoChange_n','NoChange_pct','Decreased_n','Decreased_pct'});

pooledTableSpeed = table( ...
    ["Pooled"; "Grand mean across rats"], ...
    [numel(allDeltaSpeed); NaN], ...
    [pooledCountsSpeed(1); NaN], ...
    [pooledPctSpeed(1); grandMeanPctSpeed(1)], ...
    [pooledCountsSpeed(2); NaN], ...
    [pooledPctSpeed(2); grandMeanPctSpeed(2)], ...
    [pooledCountsSpeed(3); NaN], ...
    [pooledPctSpeed(3); grandMeanPctSpeed(3)], ...
    'VariableNames', {'Group','N','Increased_n','Increased_pct', ...
                      'NoChange_n','NoChange_pct','Decreased_n','Decreased_pct'});

cellTable = table( ...
    string(allRat), allDelta, ...
    classifyDelta(allDelta, tol), ...
    'VariableNames', {'Rat','delta_MI','Category'});

cellTableSpeed = table( ...
    string(allRatSpeed), allDeltaSpeed, ...
    classifyDelta(allDeltaSpeed, tol), ...
    'VariableNames', {'Rat','delta_MI_without_minus_speedmatch','Category'});

F = figure('Color', 'w', 'Name', 'Neuron-level MI change after excluding tEBC');
tlo = tiledlayout(F, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(tlo, 1);
histogram(allDelta, 'Normalization', 'count', ...
    'FaceColor', [0.35 0.55 0.80], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
xline(0, 'k-', 'LineWidth', 1);
xlabel('\DeltaMI = MI without tEBC - MI with tEBC');
ylabel('Number of neurons');
title('No tEBC - with tEBC');
set(gca, 'Box', 'off', 'LineWidth', 1);

nexttile(tlo, 2);
barData = [perRatPct; grandMeanPct];
b = bar(barData, 'stacked', 'EdgeColor', 'none');
b(1).FaceColor = [0.20 0.60 0.35];
b(2).FaceColor = [0.70 0.70 0.70];
b(3).FaceColor = [0.75 0.30 0.25];
ylim([0 100]);
ylabel('Neurons (%)');
xticks(1:(nRats + 1));
xticklabels([ratLabels; {'Grand mean'}]);
xtickangle(35);
legend(shortCats, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('No tEBC - with tEBC');
set(gca, 'Box', 'off', 'LineWidth', 1);

nexttile(tlo, 3);
histogram(allDeltaSpeed, 'Normalization', 'count', ...
    'FaceColor', [0.35 0.55 0.80], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
xline(0, 'k-', 'LineWidth', 1);
xlabel('\DeltaMI = MI without tEBC - MI speed-matched no tEBC');
ylabel('Number of neurons');
title('No tEBC - speed-matched no tEBC');
set(gca, 'Box', 'off', 'LineWidth', 1);

nexttile(tlo, 4);
barDataSpeed = [perRatPctSpeed; grandMeanPctSpeed];
b = bar(barDataSpeed, 'stacked', 'EdgeColor', 'none');
b(1).FaceColor = [0.20 0.60 0.35];
b(2).FaceColor = [0.70 0.70 0.70];
b(3).FaceColor = [0.75 0.30 0.25];
ylim([0 100]);
ylabel('Neurons (%)');
xticks(1:(nRats + 1));
xticklabels([ratLabels; {'Grand mean'}]);
xtickangle(35);
legend(shortCats, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('No tEBC - speed-matched no tEBC');
set(gca, 'Box', 'off', 'LineWidth', 1);

OUT = struct();
OUT.figure = F;
OUT.tolerance = tol;
OUT.categories = categories;
OUT.perRat = perRat;
OUT.countsByRat = countTable;
OUT.percentagesByRat = pctTable;
OUT.pooledAndGrandMean = pooledTable;
OUT.cellTable = cellTable;
OUT.deltaMI_all = allDelta;
OUT.speedmatch.countsByRat = countTableSpeed;
OUT.speedmatch.percentagesByRat = pctTableSpeed;
OUT.speedmatch.pooledAndGrandMean = pooledTableSpeed;
OUT.speedmatch.cellTable = cellTableSpeed;
OUT.speedmatch.deltaMI_all = allDeltaSpeed;

fprintf('\n=== Neuron-level MI change: MI_without_tEBC - MI_with_tEBC ===\n');
fprintf('Tolerance for no change: abs(delta_MI) <= %.3g\n\n', tol);
disp(countTable);
disp(pctTable);
fprintf('\nPooled neurons:\n');
fprintf('  Increased: %d / %d (%.6f%%)\n', pooledCounts(1), numel(allDelta), pooledPct(1));
fprintf('  No change: %d / %d (%.6f%%)\n', pooledCounts(2), numel(allDelta), pooledPct(2));
fprintf('  Decreased: %d / %d (%.6f%%)\n', pooledCounts(3), numel(allDelta), pooledPct(3));
fprintf('\nGrand mean across rats:\n');
fprintf('  Increased: %.6f%%\n', grandMeanPct(1));
fprintf('  No change: %.6f%%\n', grandMeanPct(2));
fprintf('  Decreased: %.6f%%\n\n', grandMeanPct(3));

fprintf('\n=== Neuron-level MI change: MI_without_tEBC - MI_speed_matched_without_tEBC ===\n');
fprintf('Tolerance for no change: abs(delta_MI) <= %.3g\n\n', tol);
disp(countTableSpeed);
disp(pctTableSpeed);
fprintf('\nPooled neurons:\n');
fprintf('  Increased: %d / %d (%.6f%%)\n', pooledCountsSpeed(1), numel(allDeltaSpeed), pooledPctSpeed(1));
fprintf('  No change: %d / %d (%.6f%%)\n', pooledCountsSpeed(2), numel(allDeltaSpeed), pooledPctSpeed(2));
fprintf('  Decreased: %d / %d (%.6f%%)\n', pooledCountsSpeed(3), numel(allDeltaSpeed), pooledPctSpeed(3));
fprintf('\nGrand mean across rats:\n');
fprintf('  Increased: %.6f%%\n', grandMeanPctSpeed(1));
fprintf('  No change: %.6f%%\n', grandMeanPctSpeed(2));
fprintf('  Decreased: %.6f%%\n\n', grandMeanPctSpeed(3));

end


function [withMI, withoutMI, speedMatchedMI, keepMask, daysUsed] = getFig4aVectors(rat)
if isfield(rat, 'MI_wCSUS_all') && isfield(rat, 'MI_noCSUS15_all') && ...
        isfield(rat, 'MI_noCSUS15_controlSpeed_all') && isfield(rat, 'ratemask_all')
    withMI = rat.MI_wCSUS_all;
    withoutMI = rat.MI_noCSUS15_all;
    speedMatchedMI = pickSpeedMatchedColumn(rat.MI_noCSUS15_controlSpeed_all);
    keepMask = rat.ratemask_all;
    daysUsed = {};
    return;
end

requiredFields = {'MI_wCSUS', 'MI_noCSUS15', 'MI_noCSUS15_controlSpeed', 'ratemask'};
missing = requiredFields(~cellfun(@(f) isfield(rat, f), requiredFields));
if ~isempty(missing)
    error('Missing fields needed to build Fig. 4a vectors: %s', strjoin(missing, ', '));
end

dateList = autoDateList(rat);
idx = find(strcmp(dateList, rat.An), 1);
if isempty(idx) || idx < 3
    daysUsed = dateList(max(1, numel(dateList) - 2):numel(dateList));
else
    daysUsed = dateList(idx-2:idx);
end

withMI = collectDayVector(rat.MI_wCSUS, 'MI_', daysUsed);
withoutMI = collectDayVector(rat.MI_noCSUS15, 'MI_', daysUsed);
speedMatchedMI = pickSpeedMatchedColumn(collectDayMatrix(rat.MI_noCSUS15_controlSpeed, 'MIspeedmatch_', daysUsed));
keepMask = collectDayVector(rat.ratemask, 'ratemask_', daysUsed);
end


function values = collectDayVector(dayStruct, prefix, days)
values = [];
for d = 1:numel(days)
    fieldName = [prefix days{d}];
    if isfield(dayStruct, fieldName)
        v = dayStruct.(fieldName);
        values = [values; v(:)]; %#ok<AGROW>
    else
        warning('Missing field %s; skipping this day.', fieldName);
    end
end
end


function values = collectDayMatrix(dayStruct, prefix, days)
values = [];
for d = 1:numel(days)
    fieldName = [prefix days{d}];
    if isfield(dayStruct, fieldName)
        v = dayStruct.(fieldName);
        values = [values; v]; %#ok<AGROW>
    else
        warning('Missing field %s; skipping this day.', fieldName);
    end
end
end


function values = pickSpeedMatchedColumn(M)
if isvector(M)
    values = M(:);
    return;
end
if size(M, 2) < 2
    error('MI_noCSUS15_controlSpeed must have at least 2 columns to match Fig. 4a column-2 usage.');
end
values = M(:, 2);
end


function cat = classifyDelta(delta, tol)
cat = strings(numel(delta), 1);
cat(delta > tol) = "Increased spatial information";
cat(abs(delta) <= tol) = "No change";
cat(delta < -tol) = "Decreased spatial information";
end
