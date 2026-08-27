function OUT = run_deconv_event_threshold_robustness(ratNames, varargin)
%RUN_DECONV_EVENT_THRESHOLD_ROBUSTNESS Re-run key panels with stricter CNMF-E events.
%
%   OUT = run_deconv_event_threshold_robustness()
%
%   Builds conservative copies of the rat structs in the base workspace by
%   thresholding CNMF-E deconvolved event amplitudes per cell/session, then
%   temporarily swaps those copies into the base workspace and runs the figure
%   commands that consume rat.Ca_peaks.
%
%   Default event criterion:
%       amplitude > 50th percentile of nonzero event amplitudes
%   computed separately for each cell and session. This retains the largest
%   50% of inferred events for each cell/session.
%
%   The original rat variables are restored automatically when this function
%   exits, even if a panel errors.
%
%   Example when each rat already contains rat.Ca_deconvolved:
%       OUT = run_deconv_event_threshold_robustness();
%
%   Example rebuilding rat.Ca_deconvolved from original CNMF-E files:
%       OUT = run_deconv_event_threshold_robustness([], ...
%           'CnmfeSearchRoot', ['/Users/Hannah/Library/CloudStorage/' ...
%           'OneDrive-NorthwesternUniversity/Desktop/videos/eyeblink']);
%
%   Example with explicit event-amplitude structs:
%       deconvByRat.rat0314 = rat0314_deconvolved;
%       OUT = run_deconv_event_threshold_robustness([], ...
%           'EventEstimateByRat', deconvByRat);

if nargin < 1 || isempty(ratNames)
    ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end
ratNames = cellstr(ratNames);

p = inputParser;
addParameter(p, 'ThresholdMode', 'percentile', @(s) ischar(s) || isstring(s));
addParameter(p, 'Threshold', 50, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'EventEstimateByRat', struct(), @isstruct);
addParameter(p, 'EventEstimateRowsByRat', struct(), @isstruct);
addParameter(p, 'EventEstimateRATField', 'Ca_deconvolved', @(s) ischar(s) || isstring(s));
addParameter(p, 'CnmfeSearchRoot', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'RatFolderMap', defaultRatFolderMap(), @isstruct);
addParameter(p, 'RunFigures', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'RunSupplement', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'RunShuffleControls', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SaveMatPath', '', @(s) ischar(s) || isstring(s));
parse(p, varargin{:});
opt = p.Results;
opt.RunFigures = logical(opt.RunFigures);
opt.RunSupplement = logical(opt.RunSupplement);
opt.RunShuffleControls = logical(opt.RunShuffleControls);
opt.ThresholdMode = validatestring(opt.ThresholdMode, {'percentile', 'sd'}, ...
    'run_deconv_event_threshold_robustness', 'ThresholdMode');

OUT = struct();
OUT.params = opt;
OUT.ratNames = ratNames;
OUT.thresholdedRats = struct();
OUT.eventStats = table();
OUT.panels = struct('name', {}, 'status', {}, 'message', {});

fprintf('\nBuilding conservative event-thresholded rat structs...\n');
for r = 1:numel(ratNames)
    ratName = ratNames{r};
    if ~evalin('base', sprintf('exist(''%s'', ''var'')', ratName))
        error('run_deconv_event_threshold_robustness:MissingRat', ...
            'Base workspace variable %s was not found.', ratName);
    end

    rat = evalin('base', ratName);
    if ~isfield(opt.EventEstimateByRat, ratName) && ...
            ~isfield(rat, char(opt.EventEstimateRATField)) && ...
            ~isempty(opt.CnmfeSearchRoot)
        fprintf('Rebuilding %s.%s from CNMF-E files...\n', ...
            ratName, char(opt.EventEstimateRATField));
        ratSearchRoot = searchRootForRat(char(opt.CnmfeSearchRoot), ratName, opt.RatFolderMap);
        [rat.(char(opt.EventEstimateRATField)), rowIdxByDay, fileInfo] = ...
            buildCnmfeDeconvolvedFromFiles(rat, ...
            'SearchRoot', ratSearchRoot, ...
            'TraceFields', last3TraceFields(rat), ...
            'Verbose', true);
        OUT.rebuiltDeconvolved.(ratName).rowIdxByDay = rowIdxByDay;
        OUT.rebuiltDeconvolved.(ratName).fileInfo = fileInfo;
        assignin('base', ratName, rat);
    end

    eventArgs = {};
    if isfield(opt.EventEstimateByRat, ratName)
        eventArgs = [eventArgs, {'eventEstimateStruct', opt.EventEstimateByRat.(ratName)}]; %#ok<AGROW>
    else
        eventArgs = [eventArgs, {'eventEstimateRATField', char(opt.EventEstimateRATField)}]; %#ok<AGROW>
    end
    if isfield(opt.EventEstimateRowsByRat, ratName)
        eventArgs = [eventArgs, {'eventEstimateRowsStruct', opt.EventEstimateRowsByRat.(ratName)}]; %#ok<AGROW>
    end

    ratThresholded = rat;
    [ratThresholded, info] = eventsfromtrace(ratThresholded, ...
        'thresholds', opt.Threshold, ...
        'sourceMode', 'deconvolvedAmplitude', ...
        'thresholdMode', opt.ThresholdMode, ...
        'thresholdBasis', 'positiveStd', ...
        'thresholdCenter', 'mean', ...
        'daysMode', 'last3toAn', ...
        'outputFormat', 'times', ...
        'outputPrefix', outputPrefixForThreshold(opt), ...
        eventArgs{:});

    outputStruct = info.outputStructs{1};
    ratThresholded.Ca_peaks = ratThresholded.(outputStruct);
    OUT.thresholdedRats.(ratName) = ratThresholded;
    OUT.eventStats = [OUT.eventStats; summarizeEventRetention(ratName, rat, ratThresholded, info)]; %#ok<AGROW>
end

OUT.overallRetention = summarizeOverallRetention(OUT.eventStats);
printRetentionSummary(OUT.overallRetention, opt);

if opt.RunFigures
    backups = swapBaseRats(ratNames, OUT.thresholdedRats);
    cleanupObj = onCleanup(@() restoreBaseRats(backups)); %#ok<NASGU>
    OUT.panels = runRobustnessPanels(ratNames, opt.RunSupplement, opt.RunShuffleControls);
end

if ~isempty(opt.SaveMatPath)
    save(opt.SaveMatPath, 'OUT', '-v7.3');
    fprintf('Saved robustness output to %s\n', opt.SaveMatPath);
end
end


function stats = summarizeEventRetention(ratName, ratOriginal, ratThresholded, info)
rows = {};
for i = 1:numel(info.traceFields)
    traceField = info.traceFields{i};
    day = regexprep(traceField, '^CA_traces_', '');
    peakField = ['CA_peaks_' day];

    originalEvents = countEventTimes(ratOriginal.Ca_peaks.(peakField));
    retainedEvents = countEventTimes(ratThresholded.Ca_peaks.(peakField));
    removedEvents = originalEvents - retainedEvents;
    if originalEvents > 0
        removedPct = 100 .* removedEvents ./ originalEvents;
        retainedPct = 100 .* retainedEvents ./ originalEvents;
    else
        removedPct = NaN;
        retainedPct = NaN;
    end

    rows(end + 1, :) = {ratName, day, originalEvents, retainedEvents, ...
        removedEvents, removedPct, retainedPct}; %#ok<AGROW>
end

stats = cell2table(rows, 'VariableNames', { ...
    'Rat', 'Day', 'OriginalEvents', 'RetainedEvents', ...
    'RemovedEvents', 'RemovedPct', 'RetainedPct'});
end


function n = countEventTimes(eventMatrix)
if isempty(eventMatrix)
    n = 0;
    return;
end
eventMatrix = double(eventMatrix);
n = sum(isfinite(eventMatrix(:)) & eventMatrix(:) > 0);
end


function overall = summarizeOverallRetention(eventStats)
original = sum(eventStats.OriginalEvents, 'omitnan');
retained = sum(eventStats.RetainedEvents, 'omitnan');
removed = original - retained;
overall = struct();
overall.originalEvents = original;
overall.retainedEvents = retained;
overall.removedEvents = removed;
overall.removedPct = 100 .* removed ./ original;
overall.retainedPct = 100 .* retained ./ original;
end


function printRetentionSummary(overall, opt)
fprintf(['\nConservative CNMF-E amplitude threshold retained %d/%d events ' ...
    '(%.1f%% retained; %.1f%% removed).\n'], ...
    overall.retainedEvents, overall.originalEvents, ...
    overall.retainedPct, overall.removedPct);
switch char(opt.ThresholdMode)
    case 'percentile'
        fprintf(['Suggested wording: We retained only the largest %.0f%% of ' ...
            'CNMF-E inferred events for each cell/session. This removed %.1f%% ' ...
            'of inferred events while preserving the principal results.\n\n'], ...
            100 - opt.Threshold, overall.removedPct);
    case 'sd'
        fprintf(['Suggested wording: This removed %.1f%% of inferred events while ' ...
            'preserving the principal results.\n\n'], overall.removedPct);
end
end


function outputPrefix = outputPrefixForThreshold(opt)
switch char(opt.ThresholdMode)
    case 'percentile'
        outputPrefix = 'Ca_peaks_ampPct';
    case 'sd'
        outputPrefix = 'Ca_peaks_ampStd';
end
end


function folderMap = defaultRatFolderMap()
folderMap = struct();
folderMap.rat0222 = '022223';
folderMap.rat0307 = '030723';
folderMap.rat0313 = '031323';
folderMap.rat0314 = '031423';
folderMap.rat0816 = '081622';
end


function ratSearchRoot = searchRootForRat(searchRoot, ratName, folderMap)
ratSearchRoot = searchRoot;
if isempty(searchRoot) || ~isfield(folderMap, ratName)
    return;
end

candidate = fullfile(searchRoot, folderMap.(ratName));
if isfolder(candidate)
    ratSearchRoot = candidate;
    fprintf('Using rat-specific CNMF-E search root for %s: %s\n', ...
        ratName, ratSearchRoot);
end
end


function backups = swapBaseRats(ratNames, thresholdedRats)
backups = struct();
for r = 1:numel(ratNames)
    ratName = ratNames{r};
    backups.(ratName) = evalin('base', ratName);
    assignin('base', ratName, thresholdedRats.(ratName));
end
fprintf('Temporarily replaced base-workspace rat structs with thresholded copies.\n');
end


function restoreBaseRats(backups)
ratNames = fieldnames(backups);
for r = 1:numel(ratNames)
    assignin('base', ratNames{r}, backups.(ratNames{r}));
end
fprintf('Restored original base-workspace rat structs.\n');
end


function panels = runRobustnessPanels(ratNames, runSupplement, runShuffleControls)
panels = struct('name', {}, 'status', {}, 'message', {});

panels(end + 1) = runPanel('Fig2b_task_modulated_cells', ...
    @() plotProportionModulated('SaveMod', true)); %#ok<AGROW>

panels(end + 1) = runPanel('Fig2c_epoch_modulation_main', ...
    @() epochModulation_venn(ratNames, 'PlaceOnly', false, ...
        'CollapsePreIntoOne', true)); %#ok<AGROW>

panels(end + 1) = runPanel('Fig2g_h_speed_matched_event_rates', ...
    @() run_speedBinMatched(ratNames)); %#ok<AGROW>

panels(end + 1) = runPanel('Fig3a_b_population_temporal_structure', ...
    @() populationOrthogonality_group(ratNames, 'Mode', 'demean', ...
        'NSplits', 15)); %#ok<AGROW>

panels(end + 1) = runPanel('Fig4a_MI_without_tEBC_events', ...
    @() runFig4MI(ratNames)); %#ok<AGROW>

panels(end + 1) = runPanel('Fig4f_g_task_to_space_interference', ...
    @() run_task_to_space_interference(ratNames)); %#ok<AGROW>

panels(end + 1) = runPanel('Fig4h_task_space_stability', ...
    @() run_task_space_stability_2x2(ratNames)); %#ok<AGROW>

if runSupplement
    panels(end + 1) = runPanel('Supp_Fig2_epoch_modulation_uncollapsed', ...
        @() epochModulation_venn(ratNames, 'PlaceOnly', false, ...
            'CollapsePreIntoOne', false)); %#ok<AGROW>

    panels(end + 1) = runPanel('Supp_Fig2_epoch_modulation_dep', ...
        @() epochModulation_venn_dep(ratNames, 'CollapsePreIntoOne', false)); %#ok<AGROW>

    panels(end + 1) = runPanel('Supp_Fig2_spikes_per_bin', ...
        @() SpikesPerBin(ratNames)); %#ok<AGROW>

    panels(end + 1) = runPanel('Supp_Fig3_speed_summaries', ...
        @() runFig3Supplement(ratNames)); %#ok<AGROW>

    panels(end + 1) = runPanel('Supp_Fig4_MI_controls', ...
        @() runFig4MIControls(ratNames, runShuffleControls)); %#ok<AGROW>

    panels(end + 1) = runPanel('Supp_Fig4_space_task_interference', ...
        @() run_space_to_task_interference(ratNames)); %#ok<AGROW>
end
end


function panel = runPanel(name, fn)
panel = struct('name', name, 'status', 'ok', 'message', '');
fprintf('\n--- Running %s ---\n', name);
try
    fn();
catch ME
    panel.status = 'error';
    panel.message = ME.message;
    warning('run_deconv_event_threshold_robustness:PanelFailed', ...
        '%s failed: %s', name, ME.message);
end
end


function runFig4MI(ratNames)
buildCsus15(ratNames);
for r = 1:numel(ratNames)
    ratName = ratNames{r};
    rat = evalin('base', ratName);
    rat.MI_noCSUS15 = mutualinfo_openfield_noCSUS( ...
        rat.Ca_peaks, rat.pos, 4, 2.5, rat.Ca_ts, rat.csus15);
    rat.MI_wCSUS = mutualinfo_openfield_wCSUS( ...
        rat.Ca_peaks, rat.pos, 4, 2.5, rat.Ca_ts, rat.csus15);
    rat = attachLast3AllFields(rat);
    assignin('base', ratName, rat);
end
plotMIperAnimal(ratNames);
if allRatsHaveField(ratNames, 'MI_noCSUS15_shuff')
    MIstats;
else
    fprintf('Skipping MIstats because MI_noCSUS15_shuff is not present for all rats.\n');
end
end


function runFig4MIControls(ratNames, runShuffleControls)
buildCsus15(ratNames);
if runShuffleControls && exist('run_MI_noCSUS_Shuffle', 'file') == 2
    run_MI_noCSUS_Shuffle(15);
end
run_MI_control_matchSpikes(15);
attachAllFieldsForBaseRats(ratNames);
MI_distribution;
if exist('run_MI_control_matchSpikesSpeed', 'file') == 2
    run_MI_control_matchSpikesSpeed(15);
end
end


function runFig3Supplement(ratNames)
if exist('plotFRvsSpeedSummary', 'file') == 2
    plotFRvsSpeedSummary('fdr');
end
if exist('plotFRvsSpeedWithinTrials', 'file') == 2
    plotFRvsSpeedWithinTrials('correction', 'fdr');
end
if exist('run_accelBinMatched', 'file') == 2
    run_accelBinMatched(ratNames);
end
if exist('plotTraceFRDist_div', 'file') == 2
    plotTraceFRDist_div;
end
end


function attachAllFieldsForBaseRats(ratNames)
for r = 1:numel(ratNames)
    ratName = ratNames{r};
    rat = evalin('base', ratName);
    rat = attachLast3AllFields(rat);
    assignin('base', ratName, rat);
end
end


function rat = attachLast3AllFields(rat)
dateList = autoDateList(rat);
idx = find(strcmp(dateList, rat.An), 1);
if isempty(idx) || idx < 3
    days = dateList(max(1, numel(dateList) - 2):numel(dateList));
else
    days = dateList(idx-2:idx);
end

rat.ratemask_all = collectDayVector(rat.ratemask, 'ratemask_', days);
if isfield(rat, 'MI_noCSUS15')
    rat.MI_noCSUS15_all = collectDayVector(rat.MI_noCSUS15, 'MI_', days);
end
if isfield(rat, 'MI_wCSUS')
    rat.MI_wCSUS_all = collectDayVector(rat.MI_wCSUS, 'MI_', days);
end
if isfield(rat, 'MI_noCSUS15_control')
    rat.MI_noCSUS15_control_all = collectDayMatrix(rat.MI_noCSUS15_control, 'MI_', days);
end
if isfield(rat, 'MI_noCSUS15_controlSpeed')
    rat.MI_noCSUS15_controlSpeed_all = collectDayMatrix(rat.MI_noCSUS15_controlSpeed, 'MIspeedmatch_', days);
end
end


function values = collectDayVector(dayStruct, prefix, days)
values = [];
for d = 1:numel(days)
    fieldName = [prefix days{d}];
    if isfield(dayStruct, fieldName)
        v = dayStruct.(fieldName);
        values = [values; v(:)]; %#ok<AGROW>
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
    end
end
end


function tf = allRatsHaveField(ratNames, fieldName)
tf = true;
for r = 1:numel(ratNames)
    rat = evalin('base', ratNames{r});
    tf = tf && isfield(rat, fieldName);
end
end


function buildCsus15(ratNames)
for r = 1:numel(ratNames)
    ratName = ratNames{r};
    rat = evalin('base', ratName);
    if ~isfield(rat, 'csus15') || ~isstruct(rat.csus15)
        rat.csus15 = BULKconverttoframe15(rat.US_times, rat.Ca_ts);
        assignin('base', ratName, rat);
    end
end
end


function traceFields = last3TraceFields(rat)
dateList = autoDateList(rat);
idx = find(strcmp(dateList, rat.An), 1);
if isempty(idx) || idx < 3
    days = dateList(max(1, numel(dateList) - 2):numel(dateList));
else
    days = dateList(idx-2:idx);
end
traceFields = cellfun(@(d) ['CA_traces_' d], days, 'UniformOutput', false);
traceFields = traceFields(isfield(rat.Ca_traces, traceFields));
end
