function results = drift(rats, outputDir, varargin)
%DRIFT Quantify longitudinal stability of place- and task-cell identity.
%
% results = drift(rats, outputDir)
% results = drift(..., 'Dates', explicitDates)
%
% rats may be a cell array of rat structures or base-workspace variable
% names. outputDir may be empty to suppress file export. Automatic date
% mapping follows the manuscript convention used by plotProportionModulated:
% alignment columns 1:3 are the three chronological sessions ending at
% RAT.An (An-2, An-1, An). The mapping is accepted only when all three days
% exist in both MI and modulation results and indices are size-compatible.
%
% Explicit dates may be:
%   * an nRat-by-3 cell array;
%   * a 1-by-nRat cell array whose entries are 1-by-3 cell arrays; or
%   * a struct with fields matching rat names, each containing three dates.
%
% Task-enhanced classification exactly follows plotProportionModulated: a
% stored right-tailed permutation p-value is significant when p < TaskAlpha
% (defaults: TaskTail='right', TaskAlpha=0.05). This matches the manuscript's
% task-enhanced population; it should not be described as all task-modulated
% cells. Struct outputs with p_right/p_left/p_two retain task-down and
% two-sided results for optional secondary analyses.

p = inputParser;
p.addRequired('rats', @(x) iscell(x) && ~isempty(x));
p.addRequired('outputDir', @(x) isempty(x) || ischar(x) || isstring(x));
p.addParameter('Dates', [], @(x) isempty(x) || iscell(x) || isstruct(x));
p.addParameter('RatNames', {}, @(x) iscell(x) || isstring(x));
p.addParameter('NShuffles', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('NBootstrap', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('Seed', 20260728, @(x) isnumeric(x) && isscalar(x));
p.addParameter('PlaceThreshold', 0.95, @(x) isnumeric(x) && isscalar(x));
p.addParameter('TaskAlpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('TaskTail', 'right', @(x) any(strcmpi(string(x), ["right","left","two","two-sided"])));
p.addParameter('MinTracked', 5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('SaveFigures', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(rats, outputDir, varargin{:});
opt = p.Results;
opt.NShuffles = round(opt.NShuffles);
opt.NBootstrap = round(opt.NBootstrap);
outputDir = char(outputDir);

% Resolve inputs without modifying the source structures.
[ratData, ratNames] = resolveRats(rats, opt.RatNames);
nRats = numel(ratData);
rng(opt.Seed, 'twister');

if ~isempty(outputDir) && ~isfolder(outputDir)
    mkdir(outputDir);
end

pairCols = [1 2; 1 3; 2 3];
pairLabels = ["Day 1-Day 2"; "Day 1-Day 3"; "Day 2-Day 3"];
codingTypes = ["Place"; "TaskEnhanced"];
rows = repmat(emptyRow(), 0, 1);
dateMapRows = repmat(struct('Rat',"",'Day1',"",'Day2',"",'Day3',"",'Source',""), 0, 1);

for r = 1:nRats
    rat = ratData{r};
    ratName = string(ratNames{r});

    if ~isfield(rat, 'alignment') || ~isnumeric(rat.alignment) || size(rat.alignment,2) < 3
        warning('drift:MissingAlignment', ...
            '[%s] Missing numeric alignment with at least three columns; emitting no pair rows for this rat.', ratName);
        continue
    end

    explicit = explicitDatesForRat(opt.Dates, r, ratNames{r}, nRats);
    [dates, dateSource] = mapCriterionDates(rat, ratName, explicit);
    dateMapRows(end+1) = struct('Rat',ratName,'Day1',string(dates{1}), ...
        'Day2',string(dates{2}),'Day3',string(dates{3}),'Source',string(dateSource)); %#ok<AGROW>

    % Read classifications once per day. Empty values are retained and
    % produce explicit NaN result rows rather than silently dropping data.
    place = cell(1,3);
    task = cell(1,3);
    for d = 1:3
        place{d} = readPlace(rat, dates{d}, ratName, opt.PlaceThreshold);
        task{d} = readTask(rat, dates{d}, ratName, opt.TaskAlpha, opt.TaskTail);
        if place{d}.n > 0 && task{d}.n > 0 && place{d}.n ~= task{d}.n
            warning('drift:MismatchedCellCounts', ...
                '[%s %s] MI has %d cells but modulation has %d cells; indices are validated separately.', ...
                ratName, dates{d}, place{d}.n, task{d}.n);
        end
    end

    for q = 1:3
        a = pairCols(q,1);
        b = pairCols(q,2);
        for c = 1:2
            if c == 1
                classA = place{a}; classB = place{b};
            else
                classA = task{a}; classB = task{b};
            end
            row = analyzePair(rat.alignment(:,a), rat.alignment(:,b), ...
                classA, classB, ratName, dates{a}, dates{b}, ...
                pairLabels(q), codingTypes(c), opt);
            rows(end+1) = row; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    tidy = struct2table(repmat(emptyRow(),0,1));
else
    tidy = struct2table(rows);
end
dateMapping = struct2table(dateMapRows);

% The rat is the statistical unit: average the three session pairs within
% each rat, then perform a paired t-test across rats. This avoids treating
% the 15 rat-by-pair observations as independent.
[summary, stats] = ratSummary(tidy, ratNames, opt.NBootstrap);
[sessionStats] = sessionPairStatistics(tidy, ratNames, opt.NBootstrap);
[fig1, fig2] = makeFigures(tidy, summary, stats, opt.Visible);
if logical(opt.Verbose)
    printStatistics(summary, stats, sessionStats);
end

if ~isempty(outputDir)
    writetable(tidy, fullfile(outputDir, 'drift_pairwise_results.csv'));
    writetable(summary, fullfile(outputDir, 'drift_rat_summary.csv'));
    writetable(sessionStats, fullfile(outputDir, 'drift_session_pair_statistics.csv'));
    writetable(dateMapping, fullfile(outputDir, 'drift_date_mapping.csv'));
    if logical(opt.SaveFigures)
        exportgraphics(fig1, fullfile(outputDir, 'drift_rat_level_comparison.pdf'), 'ContentType','vector');
        exportgraphics(fig1, fullfile(outputDir, 'drift_rat_level_comparison.png'), 'Resolution',300);
        exportgraphics(fig2, fullfile(outputDir, 'drift_by_session_pair.pdf'), 'ContentType','vector');
        exportgraphics(fig2, fullfile(outputDir, 'drift_by_session_pair.png'), 'Resolution',300);
    end
end

results = struct('pairwiseTable',tidy, 'ratSummaryTable',summary, ...
    'dateMappingTable',dateMapping, 'statistics',stats, ...
    'sessionPairStatistics',sessionStats, ...
    'ratLevelFigure',fig1, 'sessionPairFigure',fig2, 'options',opt);
end

function s = emptyRow()
s = struct('Rat',"",'SessionA',"",'SessionB',"",'SessionPair',"", ...
    'CodingType',"",'NTracked',NaN,'NClassifiedA',NaN,'NClassifiedB',NaN, ...
    'NOverlap',NaN,'Recurrence',NaN,'Jaccard',NaN, ...
    'NullOverlapMean',NaN,'CorrectedOverlap',NaN,'OverlapZ',NaN,'OverlapP',NaN, ...
    'NullRecurrenceMean',NaN,'CorrectedRecurrence',NaN,'RecurrenceZ',NaN,'RecurrenceP',NaN, ...
    'NullJaccardMean',NaN,'CorrectedJaccard',NaN,'JaccardZ',NaN,'JaccardP',NaN, ...
    'NTaskUpA',NaN,'NTaskUpB',NaN,'NTaskUpOverlap',NaN,'TaskUpRecurrence',NaN, ...
    'NTaskDownA',NaN,'NTaskDownB',NaN,'NTaskDownOverlap',NaN,'TaskDownRecurrence',NaN, ...
    'NRecurrentTaskDirectional',NaN,'NDirectionPreserved',NaN,'DirectionPreservedFraction',NaN);
end

function [data, names] = resolveRats(rats, suppliedNames)
n = numel(rats);
data = cell(1,n);
names = cell(1,n);
suppliedNames = cellstr(suppliedNames);
for i = 1:n
    if ischar(rats{i}) || (isstring(rats{i}) && isscalar(rats{i}))
        names{i} = char(rats{i});
        try
            data{i} = evalin('base', names{i});
        catch ME
            error('drift:RatNotFound', 'Could not read base-workspace variable %s: %s', names{i}, ME.message);
        end
    elseif isstruct(rats{i}) && isscalar(rats{i})
        data{i} = rats{i};
        if i <= numel(suppliedNames) && ~isempty(suppliedNames{i})
            names{i} = suppliedNames{i};
        elseif isfield(rats{i}, 'name') && (ischar(rats{i}.name) || isstring(rats{i}.name))
            names{i} = char(rats{i}.name);
        else
            names{i} = sprintf('rat%02d', i);
            warning('drift:GeneratedRatName', ...
                'Structure input %d has no name; using %s. Supply ''RatNames'' for manuscript labels.', i, names{i});
        end
    else
        error('drift:BadRatInput', 'rats{%d} must be a scalar structure or variable name.', i);
    end
end
end

function dates = explicitDatesForRat(spec, idx, name, nRats)
dates = {};
if isempty(spec), return; end
if isstruct(spec)
    if isfield(spec, name), dates = spec.(name); end
elseif iscell(spec)
    if size(spec,1) == nRats && size(spec,2) == 3
        dates = spec(idx,:);
    elseif numel(spec) == nRats && iscell(spec{idx})
        dates = spec{idx};
    elseif nRats == 1 && numel(spec) == 3
        dates = spec;
    end
end
if ~isempty(dates)
    dates = cellstr(string(dates));
    if numel(dates) ~= 3
        error('drift:BadExplicitDates', '[%s] Explicit date list must contain exactly three dates.', name);
    end
end
end

function [dates, source] = mapCriterionDates(rat, ratName, explicit)
if ~isempty(explicit)
    dates = normalizeDates(explicit);
    validateDateSizes(rat, dates, ratName);
    source = 'explicit';
    return
end

if ~isfield(rat,'An') || isempty(rat.An)
    error('drift:DatesRequired', ...
        '[%s] Cannot infer criterion dates safely because RAT.An is missing. Supply ''Dates''.', ratName);
end
try
    allDates = cellstr(autoDateList(rat));
catch ME
    error('drift:DatesRequired', ...
        '[%s] autoDateList failed (%s). Supply ''Dates''.', ratName, ME.message);
end
allDates = normalizeDates(allDates);
an = normalizeDate(rat.An);
iAn = find(strcmp(allDates, an), 1);
if isempty(iAn) || iAn < 3
    error('drift:DatesRequired', ...
        '[%s] Could not identify three dates ending at RAT.An=%s. Supply ''Dates''.', ratName, an);
end
dates = allDates(iAn-2:iAn);
validateDateAvailability(rat, dates, ratName);
validateDateSizes(rat, dates, ratName);
source = 'automatic: chronological An-2 through An';
end

function validateDateAvailability(rat, dates, ratName)
for d = 1:3
    miField = ['MI_' dates{d}];
    modField = ['mod_' dates{d}];
    if ~isfield(rat,'MI_noCSUS15_shuff') || ~isfield(rat.MI_noCSUS15_shuff,miField) || ...
            ~isfield(rat,'mod') || ~isfield(rat.mod,modField)
        warning('drift:MissingDateField', ...
            '[%s] Date %s is missing from MI and/or modulation results.', ratName, dates{d});
        error('drift:DatesRequired', ...
            '[%s] Date %s is not present in both MI and modulation results. Supply a verified explicit ''Dates'' mapping or populate missing results.', ...
            ratName, dates{d});
    end
end
end

function validateDateSizes(rat, dates, ratName)
for d = 1:3
    miField = ['MI_' dates{d}];
    modField = ['mod_' dates{d}];
    if ~isfield(rat,'MI_noCSUS15_shuff') || ~isfield(rat.MI_noCSUS15_shuff,miField) || ...
            ~isfield(rat,'mod') || ~isfield(rat.mod,modField)
        % Explicit mappings are allowed when a result is missing so that
        % readPlace/readTask can warn and emit a visible NaN row.
        continue
    end
    ids = rat.alignment(:,d);
    ids = ids(isfinite(ids) & ids > 0 & ids == fix(ids));
    if isempty(ids), continue; end
    nMI = size(rat.MI_noCSUS15_shuff.(miField),1);
    nMod = modLength(rat.mod.(modField));
    if max(ids) > nMI || (~isnan(nMod) && max(ids) > nMod)
        error('drift:DatesRequired', ...
            ['[%s] Automatic/explicit date %s is size-incompatible with alignment column %d ' ...
             '(max index %d, MI rows %d, modulation rows %g). Verify and supply corrected ''Dates''.'], ...
             ratName, dates{d}, d, max(ids), nMI, nMod);
    end
end
end

function n = modLength(x)
n = NaN;
if isnumeric(x) || islogical(x)
    if isvector(x), n = numel(x); else, n = size(x,1); end
elseif isstruct(x)
    f = {'p_right','p_left','p_two','p_two_sided','p'};
    for i = 1:numel(f)
        if isfield(x,f{i}) && isnumeric(x.(f{i}))
            n = numel(x.(f{i})); return
        end
    end
end
end

function out = readPlace(rat, date, ratName, threshold)
out = emptyClass();
fld = ['MI_' date];
if ~isfield(rat,'MI_noCSUS15_shuff') || ~isfield(rat.MI_noCSUS15_shuff,fld)
    warning('drift:MissingMI','[%s %s] Missing MI_noCSUS15_shuff.%s.',ratName,date,fld); return
end
x = rat.MI_noCSUS15_shuff.(fld);
if ~isnumeric(x) || size(x,2) < 3
    warning('drift:BadMI','[%s %s] MI result must be numeric with at least three columns.',ratName,date); return
end
score = x(:,3);
out.any = isfinite(score) & score >= threshold;
out.valid = isfinite(score);
out.n = numel(score);
end

function out = readTask(rat, date, ratName, alpha, tail)
out = emptyClass();
fld = ['mod_' date];
if ~isfield(rat,'mod') || ~isfield(rat.mod,fld)
    warning('drift:MissingMod','[%s %s] Missing mod.%s.',ratName,date,fld); return
end
x = rat.mod.(fld);
if isnumeric(x) && isvector(x)
    if any(strcmpi(string(tail), ["two","two-sided"]))
        warning('drift:LegacyModCannotBeTwoSided', ...
            ['[%s %s] mod.%s is a legacy selected-tail numeric vector and cannot define all ' ...
             'task-modulated cells. Rerun plotProportionModulated with SaveMod=true to save ' ...
             'p_two/p_right/p_left, or explicitly request TaskTail=''right''/''left''.'], ...
             ratName,date,fld);
        return
    end
    pv = x(:);
    out.any = isfinite(pv) & pv < alpha;
    out.valid = isfinite(pv);
    out.n = numel(pv);
    % Numeric vectors contain only the p-value for the selected Tail. The
    % default plotProportionModulated Tail is right, so direction is known
    % only when the caller confirms a one-sided tail.
    if strcmpi(tail,'right')
        out.up = out.any; out.down = false(size(out.any)); out.hasDirection = true;
    elseif strcmpi(tail,'left')
        out.down = out.any; out.up = false(size(out.any)); out.hasDirection = true;
    end
elseif isstruct(x)
    pr = getVectorField(x, {'p_right'});
    pl = getVectorField(x, {'p_left'});
    pt = getVectorField(x, {'p_two','p_two_sided','p'});
    switch lower(char(tail))
        case 'right', pv = pr;
        case 'left', pv = pl;
        otherwise, pv = pt;
    end
    if isempty(pv)
        warning('drift:BadModFormat', ...
            '[%s %s] Struct mod output lacks the p-value field required for TaskTail=%s.',ratName,date,tail); return
    end
    pv = pv(:);
    out.any = isfinite(pv) & pv < alpha;
    out.valid = isfinite(pv);
    out.n = numel(pv);
    if ~isempty(pr) && ~isempty(pl) && numel(pr)==out.n && numel(pl)==out.n
        out.up = isfinite(pr(:)) & pr(:) < alpha;
        out.down = isfinite(pl(:)) & pl(:) < alpha;
        out.hasDirection = true;
    end
else
    warning('drift:BadModFormat', ...
        '[%s %s] Cannot interpret mod.%s. Expected a numeric p-value vector or p-value struct.',ratName,date,fld);
end
end

function x = getVectorField(s, choices)
x = [];
for i = 1:numel(choices)
    if isfield(s,choices{i}) && isnumeric(s.(choices{i})) && isvector(s.(choices{i}))
        x = s.(choices{i}); return
    end
end
end

function x = emptyClass()
x = struct('any',[],'valid',[],'n',0,'up',[],'down',[],'hasDirection',false);
end

function row = analyzePair(colA, colB, aClass, bClass, ratName, dateA, dateB, pairLabel, codingType, opt)
row = emptyRow();
row.Rat = ratName; row.SessionA = string(dateA); row.SessionB = string(dateB);
row.SessionPair = pairLabel; row.CodingType = codingType;
if isempty(aClass.any) || isempty(bClass.any), return; end

integerA = isfinite(colA) & colA == fix(colA);
integerB = isfinite(colB) & colB == fix(colB);
bad = (colA ~= 0 & ~integerA) | (colB ~= 0 & ~integerB) | colA < 0 | colB < 0;
if any(bad)
    warning('drift:InvalidAlignmentIndex', ...
        '[%s %s] %d rows contain negative, nonfinite, or noninteger alignment indices and were excluded.', ...
        ratName, pairLabel, nnz(bad));
end
detected = integerA & integerB & colA > 0 & colB > 0;
inRange = detected & colA <= aClass.n & colB <= bClass.n;
if nnz(detected & ~inRange) > 0
    warning('drift:IndexOutOfRange', ...
        '[%s %s %s] %d tracked rows exceed classification array bounds and were excluded.', ...
        ratName, pairLabel, codingType, nnz(detected & ~inRange));
end
idxA = colA(inRange); idxB = colB(inRange);
valid = aClass.valid(idxA) & bClass.valid(idxB);
if nnz(~valid) > 0
    warning('drift:UnanalyzedCells', ...
        '[%s %s %s] %d aligned cells have NaN/unavailable classifications and were excluded.', ...
        ratName, pairLabel, codingType, nnz(~valid));
end
idxA = idxA(valid); idxB = idxB(valid);
va = logical(aClass.any(idxA)); vb = logical(bClass.any(idxB));
row.NTracked = numel(va);
if row.NTracked < opt.MinTracked
    warning('drift:TooFewTracked', '[%s %s %s] Only %d usable tracked cells.', ...
        ratName, pairLabel, codingType, row.NTracked);
end
if isempty(va), return; end

[overlap, recurrence, jaccard] = metrics(va,vb);
row.NClassifiedA = nnz(va); row.NClassifiedB = nnz(vb); row.NOverlap = overlap;
row.Recurrence = recurrence; row.Jaccard = jaccard;

null = nan(opt.NShuffles,3);
for k = 1:opt.NShuffles
    [null(k,1),null(k,2),null(k,3)] = metrics(va, vb(randperm(numel(vb))));
end
[row.NullOverlapMean,row.CorrectedOverlap,row.OverlapZ,row.OverlapP] = nullStats(overlap,null(:,1));
[row.NullRecurrenceMean,row.CorrectedRecurrence,row.RecurrenceZ,row.RecurrenceP] = nullStats(recurrence,null(:,2));
[row.NullJaccardMean,row.CorrectedJaccard,row.JaccardZ,row.JaccardP] = nullStats(jaccard,null(:,3));

if codingType == "TaskEnhanced" && aClass.hasDirection && bClass.hasDirection
    upA = logical(aClass.up(idxA)); upB = logical(bClass.up(idxB));
    dnA = logical(aClass.down(idxA)); dnB = logical(bClass.down(idxB));
    row.NTaskUpA=nnz(upA); row.NTaskUpB=nnz(upB); row.NTaskUpOverlap=nnz(upA&upB);
    row.TaskUpRecurrence = safeDivide(row.NTaskUpOverlap,row.NTaskUpA);
    row.NTaskDownA=nnz(dnA); row.NTaskDownB=nnz(dnB); row.NTaskDownOverlap=nnz(dnA&dnB);
    row.TaskDownRecurrence = safeDivide(row.NTaskDownOverlap,row.NTaskDownA);
    recurrent = va & vb;
    directional = recurrent & (upA | dnA) & (upB | dnB);
    preserved = directional & ((upA & upB) | (dnA & dnB));
    row.NRecurrentTaskDirectional = nnz(directional);
    row.NDirectionPreserved = nnz(preserved);
    row.DirectionPreservedFraction = safeDivide(nnz(preserved),nnz(directional));
end
end

function [o,r,j] = metrics(a,b)
o = nnz(a & b);
r = safeDivide(o,nnz(a));
j = safeDivide(o,nnz(a | b));
end

function y = safeDivide(a,b)
if b == 0, y = NaN; else, y = a/b; end
end

function [mu,corrected,z,pv] = nullStats(observed,null)
null = null(isfinite(null));
if isempty(null) || ~isfinite(observed)
    mu=NaN; corrected=NaN; z=NaN; pv=NaN; return
end
mu = mean(null);
corrected = observed-mu;
sd = std(null,0);
if sd > 0, z=corrected/sd; else, z=NaN; end
pv = (nnz(abs(null-mu) >= abs(observed-mu)) + 1) / (numel(null)+1);
end

function [summary,stats] = ratSummary(tidy, ratNames, nBootstrap)
n = numel(ratNames);
Rat = string(ratNames(:));
PlaceRawRecurrence = nan(n,1);
PlaceNullRecurrence = nan(n,1);
PlaceCorrectedRecurrence = nan(n,1);
TaskRawRecurrence = nan(n,1);
TaskNullRecurrence = nan(n,1);
TaskCorrectedRecurrence = nan(n,1);
for i=1:n
    isPlace=tidy.Rat==Rat(i) & tidy.CodingType=="Place";
    isTask=tidy.Rat==Rat(i) & tidy.CodingType=="TaskEnhanced";
    PlaceRawRecurrence(i)=mean(tidy.Recurrence(isPlace),'omitnan');
    PlaceNullRecurrence(i)=mean(tidy.NullRecurrenceMean(isPlace),'omitnan');
    PlaceCorrectedRecurrence(i)=mean(tidy.CorrectedRecurrence(isPlace),'omitnan');
    TaskRawRecurrence(i)=mean(tidy.Recurrence(isTask),'omitnan');
    TaskNullRecurrence(i)=mean(tidy.NullRecurrenceMean(isTask),'omitnan');
    TaskCorrectedRecurrence(i)=mean(tidy.CorrectedRecurrence(isTask),'omitnan');
end
DifferenceTaskMinusPlace = TaskCorrectedRecurrence-PlaceCorrectedRecurrence;
summary = table(Rat,PlaceRawRecurrence,PlaceNullRecurrence,PlaceCorrectedRecurrence, ...
    TaskRawRecurrence,TaskNullRecurrence,TaskCorrectedRecurrence,DifferenceTaskMinusPlace);
keep = isfinite(PlaceCorrectedRecurrence) & isfinite(TaskCorrectedRecurrence);
stats = struct('method','Rat-level means across three pairs; rat is the statistical unit', ...
    'nRats',nnz(keep),'meanDifference',NaN,'bootstrapCI',[NaN NaN], ...
    'nBootstrap',nBootstrap,'permutationP',NaN,'nExactPermutations',NaN, ...
    'ttestH',NaN,'ttestP',NaN,'ttestCI',[NaN NaN],'tstat',NaN,'df',NaN, ...
    'signrankP',NaN,'signrankSignedRank',NaN);
if nnz(keep)>=2
    differences=DifferenceTaskMinusPlace(keep);
    stats.meanDifference=mean(differences);

    % Nonparametric rat-level bootstrap confidence interval for the mean
    % paired task-minus-place difference.
    bootMeans=nan(nBootstrap,1);
    for b=1:nBootstrap
        bootMeans(b)=mean(differences(randi(numel(differences),numel(differences),1)));
    end
    stats.bootstrapCI=prctile(bootMeans,[2.5 97.5]);

    % Exact paired randomization test: enumerate every sign flip of the
    % rat-level differences (2^5 = 32 permutations for the full dataset).
    nPerm=2^numel(differences);
    permMeans=nan(nPerm,1);
    for b=0:nPerm-1
        signs=2*bitget(b,1:numel(differences))-1;
        permMeans(b+1)=mean(differences(:)'.*signs);
    end
    stats.nExactPermutations=nPerm;
    stats.permutationP=mean(abs(permMeans)>=abs(stats.meanDifference));

    [stats.ttestH,stats.ttestP,stats.ttestCI,t] = ...
        ttest(TaskCorrectedRecurrence(keep),PlaceCorrectedRecurrence(keep));
    stats.tstat=t.tstat; stats.df=t.df;

    % Sensitivity analysis; exact behavior depends on MATLAB's signrank
    % implementation and handling of zero paired differences.
    if exist('signrank','file')==2
        [stats.signrankP,~,sr]=signrank(TaskCorrectedRecurrence(keep),PlaceCorrectedRecurrence(keep));
        if isfield(sr,'signedrank'), stats.signrankSignedRank=sr.signedrank; end
    end
else
    warning('drift:TooFewRats','Fewer than two rats have paired rat-level data; inferential tests not run.');
end
end

function out = sessionPairStatistics(tidy,ratNames,nBootstrap)
% Compare place versus task-enhanced recurrence separately for each pair.
pairs=["Day 1-Day 2";"Day 1-Day 3";"Day 2-Day 3"];
nPairs=numel(pairs);
SessionPair=pairs;
NRats=zeros(nPairs,1);
PlaceRawMean=nan(nPairs,1); PlaceNullMean=nan(nPairs,1); PlaceCorrectedMean=nan(nPairs,1);
TaskRawMean=nan(nPairs,1); TaskNullMean=nan(nPairs,1); TaskCorrectedMean=nan(nPairs,1);
MeanDifference=nan(nPairs,1); BootstrapCILow=nan(nPairs,1); BootstrapCIHigh=nan(nPairs,1);
PermutationP=nan(nPairs,1); PermutationPHolm=nan(nPairs,1);
TStat=nan(nPairs,1); DF=nan(nPairs,1); TTestP=nan(nPairs,1);
TTestCILow=nan(nPairs,1); TTestCIHigh=nan(nPairs,1);
SignrankP=nan(nPairs,1); SignedRank=nan(nPairs,1);

for q=1:nPairs
    placeRaw=nan(numel(ratNames),1); placeNull=placeRaw; placeCorr=placeRaw;
    taskRaw=placeRaw; taskNull=placeRaw; taskCorr=placeRaw;
    for r=1:numel(ratNames)
        rat=string(ratNames{r});
        ip=tidy.Rat==rat & tidy.SessionPair==pairs(q) & tidy.CodingType=="Place";
        it=tidy.Rat==rat & tidy.SessionPair==pairs(q) & tidy.CodingType=="TaskEnhanced";
        if nnz(ip)==1
            placeRaw(r)=tidy.Recurrence(ip); placeNull(r)=tidy.NullRecurrenceMean(ip);
            placeCorr(r)=tidy.CorrectedRecurrence(ip);
        end
        if nnz(it)==1
            taskRaw(r)=tidy.Recurrence(it); taskNull(r)=tidy.NullRecurrenceMean(it);
            taskCorr(r)=tidy.CorrectedRecurrence(it);
        end
    end
    keep=isfinite(placeCorr)&isfinite(taskCorr);
    NRats(q)=nnz(keep);
    PlaceRawMean(q)=mean(placeRaw(keep),'omitnan');
    PlaceNullMean(q)=mean(placeNull(keep),'omitnan');
    PlaceCorrectedMean(q)=mean(placeCorr(keep),'omitnan');
    TaskRawMean(q)=mean(taskRaw(keep),'omitnan');
    TaskNullMean(q)=mean(taskNull(keep),'omitnan');
    TaskCorrectedMean(q)=mean(taskCorr(keep),'omitnan');
    if NRats(q)<2, continue; end

    differences=taskCorr(keep)-placeCorr(keep);
    MeanDifference(q)=mean(differences);
    bootMeans=nan(nBootstrap,1);
    for b=1:nBootstrap
        bootMeans(b)=mean(differences(randi(NRats(q),NRats(q),1)));
    end
    ci=prctile(bootMeans,[2.5 97.5]);
    BootstrapCILow(q)=ci(1); BootstrapCIHigh(q)=ci(2);

    nPerm=2^NRats(q);
    permMeans=nan(nPerm,1);
    for b=0:nPerm-1
        signs=2*bitget(b,1:NRats(q))-1;
        permMeans(b+1)=mean(differences(:)'.*signs);
    end
    PermutationP(q)=mean(abs(permMeans)>=abs(MeanDifference(q)));

    [~,TTestP(q),ciT,t]=ttest(taskCorr(keep),placeCorr(keep));
    TStat(q)=t.tstat; DF(q)=t.df; TTestCILow(q)=ciT(1); TTestCIHigh(q)=ciT(2);
    if exist('signrank','file')==2
        [SignrankP(q),~,sr]=signrank(taskCorr(keep),placeCorr(keep));
        if isfield(sr,'signedrank'), SignedRank(q)=sr.signedrank; end
    end
end

% Holm adjustment controls family-wise error across the three pair-specific
% place-versus-task-enhanced tests.
valid=find(isfinite(PermutationP));
if ~isempty(valid)
    [sortedP,order]=sort(PermutationP(valid));
    adjusted=cummax((numel(valid)-(1:numel(valid))'+1).*sortedP);
    adjusted=min(adjusted,1);
    tmp=nan(numel(valid),1); tmp(order)=adjusted;
    PermutationPHolm(valid)=tmp;
end

out=table(SessionPair,NRats,PlaceRawMean,PlaceNullMean,PlaceCorrectedMean, ...
    TaskRawMean,TaskNullMean,TaskCorrectedMean,MeanDifference, ...
    BootstrapCILow,BootstrapCIHigh,PermutationP,PermutationPHolm, ...
    TStat,DF,TTestP,TTestCILow,TTestCIHigh,SignrankP,SignedRank);
end

function [fig1,fig2] = makeFigures(tidy,summary,stats,visibility)
fig1=figure('Color','w','Visible',char(visibility),'Name','Rat-level corrected recurrence');
hold on
valid=isfinite(summary.PlaceCorrectedRecurrence)&isfinite(summary.TaskCorrectedRecurrence);
for i=find(valid)'
    plot([1 2],[summary.PlaceCorrectedRecurrence(i) summary.TaskCorrectedRecurrence(i)],'-o', ...
        'Color',[.65 .65 .65],'MarkerFaceColor',[.85 .85 .85]);
end
vals=[summary.PlaceCorrectedRecurrence summary.TaskCorrectedRecurrence];
mu=mean(vals(valid,:),1,'omitnan'); se=std(vals(valid,:),0,1,'omitnan')/sqrt(max(nnz(valid),1));
errorbar([1 2],mu,se,'k','LineStyle','none','LineWidth',2,'CapSize',10);
scatter([1 2],mu,75,'k','filled');
xlim([.6 2.4]); xticks([1 2]); xticklabels({'Place','Task-enhanced'}); ylabel('Chance-corrected recurrence');
title(sprintf('Mean task-enhanced minus place = %.3f, exact permutation p=%.3g', ...
    stats.meanDifference,stats.permutationP));
box off

fig2=figure('Color','w','Visible',char(visibility),'Name','Recurrence by session pair');
tl=tiledlayout(fig2,2,1,'TileSpacing','compact','Padding','compact');
pairs=["Day 1-Day 2","Day 1-Day 3","Day 2-Day 3"];
colors=[.20 .45 .85; .85 .35 .25];

% Top: the intuitive quantities on their native probability scale.
ax1=nexttile(tl); hold(ax1,'on');
for q=1:3
    for c=1:2
        typ=["Place","TaskEnhanced"]; x=q+(-.12+0.24*(c-1));
        use=tidy.SessionPair==pairs(q)&tidy.CodingType==typ(c);
        raw=tidy.Recurrence(use); null=tidy.NullRecurrenceMean(use);
        raw=raw(isfinite(raw)); null=null(isfinite(null));
        scatter(ax1,x-.025+zeros(size(raw)),raw,28,colors(c,:),'filled','MarkerFaceAlpha',.55);
        scatter(ax1,x+.025+zeros(size(null)),null,30,colors(c,:),'d','MarkerFaceColor','w', ...
            'MarkerEdgeColor',colors(c,:),'LineWidth',1);
        if ~isempty(raw)
            errorbar(ax1,x-.025,mean(raw),std(raw)/sqrt(numel(raw)),'Color',colors(c,:), ...
                'LineWidth',2,'CapSize',7);
        end
    end
end
xlim(ax1,[.5 3.5]); xticks(ax1,1:3); xticklabels(ax1,{'Day 1-Day 2','Day 1-Day 3','Day 2-Day 3'});
ylabel(ax1,'Recurrence probability');
title(ax1,'Raw recurrence (filled) and shuffled expectation (open diamonds)');
box(ax1,'off');

% Bottom: observed minus shuffled expectation.
ax2=nexttile(tl); hold(ax2,'on');
for q=1:3
    for c=1:2
        typ=["Place","TaskEnhanced"]; x=q+(-.12+0.24*(c-1));
        y=tidy.CorrectedRecurrence(tidy.SessionPair==pairs(q)&tidy.CodingType==typ(c));
        y=y(isfinite(y));
        scatter(ax2,x+zeros(size(y)),y,34,colors(c,:),'filled','MarkerFaceAlpha',.65);
        if ~isempty(y)
            errorbar(ax2,x,mean(y),std(y)/sqrt(numel(y)),'Color',colors(c,:), ...
                'LineWidth',2,'CapSize',8);
        end
    end
end
yline(ax2,0,':','Color',[.4 .4 .4]);
xlim(ax2,[.5 3.5]); xticks(ax2,1:3); xticklabels(ax2,{'Day 1-Day 2','Day 1-Day 3','Day 2-Day 3'});
ylabel(ax2,'Observed - shuffled'); title(ax2,'Chance-corrected recurrence');
legend(ax2,{'Place rats','Place mean','Task-enhanced rats','Task-enhanced mean'},'Location','best');
box(ax2,'off');
end

function printStatistics(summary,stats,sessionStats)
% Print manuscript-relevant descriptive and inferential statistics.
fprintf('\n=== Longitudinal classification stability ===\n');
fprintf('Rat-level means across the three session pairs:\n');
fprintf('%-12s  %9s %9s %9s  %9s %9s %9s  %9s\n', ...
    'Rat','PlaceRaw','PlaceNull','PlaceCorr','TaskRaw','TaskNull','TaskCorr','Task-Place');
for i=1:height(summary)
    fprintf('%-12s  %9.3f %9.3f %9.3f  %9.3f %9.3f %9.3f  %9.3f\n', ...
        char(summary.Rat(i)),summary.PlaceRawRecurrence(i),summary.PlaceNullRecurrence(i), ...
        summary.PlaceCorrectedRecurrence(i),summary.TaskRawRecurrence(i), ...
        summary.TaskNullRecurrence(i),summary.TaskCorrectedRecurrence(i), ...
        summary.DifferenceTaskMinusPlace(i));
end

valid=isfinite(summary.PlaceCorrectedRecurrence)&isfinite(summary.TaskCorrectedRecurrence);
fprintf('\nGroup comparison (rat is the statistical unit; n=%d):\n',stats.nRats);
if any(valid)
    fprintf('  Place corrected recurrence:         %.3f +/- %.3f SEM\n', ...
        mean(summary.PlaceCorrectedRecurrence(valid)), ...
        std(summary.PlaceCorrectedRecurrence(valid))/sqrt(nnz(valid)));
    fprintf('  Task-enhanced corrected recurrence: %.3f +/- %.3f SEM\n', ...
        mean(summary.TaskCorrectedRecurrence(valid)), ...
        std(summary.TaskCorrectedRecurrence(valid))/sqrt(nnz(valid)));
end
fprintf('  Mean paired difference (task-enhanced - place): %.3f\n',stats.meanDifference);
fprintf('  Rat-level bootstrap 95%% CI: [%.3f, %.3f] (%d resamples)\n', ...
    stats.bootstrapCI(1),stats.bootstrapCI(2),stats.nBootstrap);
fprintf('  Exact paired sign-permutation: p=%.4g (%d permutations)\n', ...
    stats.permutationP,stats.nExactPermutations);
fprintf('  Paired t-test: t(%g)=%.3f, p=%.4g, 95%% CI=[%.3f, %.3f]\n', ...
    stats.df,stats.tstat,stats.ttestP,stats.ttestCI(1),stats.ttestCI(2));
if isfinite(stats.signrankP)
    fprintf('  Wilcoxon signed-rank sensitivity: p=%.4g, signed rank=%.3f\n', ...
        stats.signrankP,stats.signrankSignedRank);
else
    fprintf('  Wilcoxon signed-rank sensitivity: unavailable\n');
end

fprintf('\nSession-pair-specific comparisons across rats:\n');
for q=1:height(sessionStats)
    s=sessionStats(q,:);
    fprintf('\n  %s (n=%d rats)\n',char(s.SessionPair),s.NRats);
    fprintf('    Place:         raw=%.3f, shuffled=%.3f, corrected=%.3f\n', ...
        s.PlaceRawMean,s.PlaceNullMean,s.PlaceCorrectedMean);
    fprintf('    Task-enhanced: raw=%.3f, shuffled=%.3f, corrected=%.3f\n', ...
        s.TaskRawMean,s.TaskNullMean,s.TaskCorrectedMean);
    fprintf('    Mean paired difference: %.3f, bootstrap 95%% CI=[%.3f, %.3f]\n', ...
        s.MeanDifference,s.BootstrapCILow,s.BootstrapCIHigh);
    fprintf('    Exact sign-permutation: p=%.4g, Holm-adjusted p=%.4g\n', ...
        s.PermutationP,s.PermutationPHolm);
    fprintf('    Paired t-test: t(%g)=%.3f, p=%.4g, 95%% CI=[%.3f, %.3f]\n', ...
        s.DF,s.TStat,s.TTestP,s.TTestCILow,s.TTestCIHigh);
    if isfinite(s.SignrankP)
        fprintf('    Wilcoxon signed-rank: p=%.4g, signed rank=%.3f\n', ...
            s.SignrankP,s.SignedRank);
    else
        fprintf('    Wilcoxon signed-rank: unavailable\n');
    end
end
fprintf('===============================================\n\n');
end

function dates = normalizeDates(dates)
dates = cellfun(@normalizeDate, cellstr(string(dates)), 'UniformOutput', false);
end

function d = normalizeDate(d)
d = char(string(d));
d = regexprep(d, '^(MI_|mod_|CA_peaks_|CA_traces_|CS_)', '');
d = strrep(d,'-','_');
end
