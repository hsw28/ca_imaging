function [RAT, eventInfo] = eventsfromtrace(RAT, varargin)
%EVENTSFROMTRACE Recompute or threshold calcium events from saved matrices.
%
%   [RAT, eventInfo] = eventsfromtrace(RAT, 'thresholds', [3 5])
%
%   Reads traces from the last three RAT.Ca_traces.CA_traces_* fields ending
%   at RAT.An, then writes analogous peak matrices to
%   RAT.Ca_peaks_std3.CA_peaks_* and RAT.Ca_peaks_std5.CA_peaks_* by default.
%   By default, the function thresholds CNMF-E deconvolved event estimates
%   (the inferred event/spike-amplitude matrix) rather than re-detecting
%   local maxima from the fluorescence trace.
%
%   Example:
%       rat0222 = eventsfromtrace(rat0222, 'thresholds', [3 5], ...
%           'eventEstimateStruct', rat0222.Ca_deconvolved);
%       peaks3 = rat0222.Ca_peaks_std3.CA_peaks_2023_05_09;
%       peaks5 = rat0222.Ca_peaks_std5.CA_peaks_2023_05_09;
%
%   Options:
%       thresholds      Numeric vector of SD thresholds. Default: [3 5]
%       sourceMode      'deconvolvedAmplitude' or 'tracePeaks'. Default:
%                       'deconvolvedAmplitude'. In deconvolvedAmplitude mode,
%                       thresholds are applied to CNMF-E inferred event
%                       amplitudes. In tracePeaks mode, CIAtah
%                       computeSignalPeaks is applied to Ca_traces.
%       eventEstimateStruct
%                       Struct containing CNMF-E deconvolved event estimates
%                       with fields matching CA_traces_* or CA_peaks_* dates.
%                       Required for deconvolvedAmplitude mode unless
%                       eventEstimateRATField is present in RAT.
%       eventEstimateRATField
%                       RAT field containing eventEstimateStruct. Default:
%                       'Ca_deconvolved'.
%       eventEstimateRows
%                       Optional rows to keep from the deconvolved estimate
%                       before thresholding, e.g. find(validCNMFE == 1) when
%                       the source matrix still contains all CNMF-E candidates.
%       eventEstimateRowsStruct
%                       Optional struct of per-day row indices. Field names may
%                       match CA_traces_DATE, CA_peaks_DATE, or DATE.
%       thresholdBasis  'allStd' or 'positiveStd'. Default: 'allStd'
%       thresholdCenter 'zero' or 'mean'. Default: 'mean'. With 'mean',
%                       events are retained when amplitude exceeds
%                       mean(basis) + threshold*std(basis), computed
%                       separately for each cell/session.
%       thresholdMode   'sd' or 'percentile'. Default: 'sd'. In percentile
%                       mode, thresholds are percentiles of each cell's
%                       nonzero deconvolved event-amplitude distribution,
%                       e.g. threshold 75 retains the largest 25% of
%                       nonzero inferred events for that cell/session.
%       detectMethod    CIAtah computeSignalPeaks method. Default: 'diff'
%       minTimeBtEvents Minimum event spacing in frames. Default: 8
%       daysMode        'last3toAn' or 'all'. Default: 'last3toAn'
%       daysToUse       Optional date strings, e.g. {'2023_05_07','2023_05_09'}.
%                       Overrides daysMode when supplied.
%       outputFormat    'times' or 'binary'. Default: 'times'
%       outputPrefix    Output struct prefix. Default: 'Ca_peaks_std'
%       overwrite       If true and one threshold is supplied, write to
%                       RAT.Ca_peaks. Default: false
%       waitbarOn       Passed to computeSignalPeaks. Default: 0
%       makePlots       Passed to computeSignalPeaks. Default: 0
%       parallel        Passed to computeSignalPeaks. Default: 0
%       makeTraceFigure If true, make raw trace/event overlay figure. Default: false
%       figureThreshold Threshold to show in trace figure. Default: first threshold
%       figureField     Trace field to show. Default: last selected CA_traces_* field
%       figureCells     Cell indices to show. Default: auto-select active cells
%       figureWindow    Frame window [start end] for plot. Default: full trace
%       frameRate       Frames per second for x-axis. Default: 7.5
%       saveFigurePath  Optional path for saving figure. Default: ''

opts.thresholds = [3 5];
opts.sourceMode = 'deconvolvedAmplitude';
opts.eventEstimateStruct = [];
opts.eventEstimateRATField = 'Ca_deconvolved';
opts.eventEstimateRows = [];
opts.eventEstimateRowsStruct = [];
opts.thresholdBasis = 'allStd';
opts.thresholdCenter = 'mean';
opts.thresholdMode = 'sd';
opts.detectMethod = 'diff';
opts.minTimeBtEvents = 8;
opts.daysMode = 'last3toAn';
opts.daysToUse = {};
opts.outputFormat = 'times';
opts.outputPrefix = 'Ca_peaks_std';
opts.overwrite = false;
opts.waitbarOn = 0;
opts.makePlots = 0;
opts.parallel = 0;
opts.makeTraceFigure = false;
opts.figureThreshold = [];
opts.figureField = '';
opts.figureCells = [];
opts.figureWindow = [];
opts.frameRate = 7.5;
opts.saveFigurePath = '';
opts = parseInputs(opts, varargin{:});

if ~isfield(RAT, 'Ca_traces') || ~isstruct(RAT.Ca_traces)
    error('eventsfromtrace:MissingTraces', ...
        'Input RAT must contain a struct field named Ca_traces.');
end

traceFields = fieldnames(RAT.Ca_traces);
traceFields = traceFields(startsWith(traceFields, 'CA_traces_'));

if isempty(traceFields)
    error('eventsfromtrace:NoTraceFields', ...
        'No fields beginning with CA_traces_ were found in RAT.Ca_traces.');
end

allTraceFields = traceFields;
traceFields = selectTraceFieldsForDays(RAT, traceFields, opts);

eventInfo = struct();
eventInfo.thresholds = opts.thresholds;
eventInfo.sourceMode = opts.sourceMode;
eventInfo.eventEstimateRATField = opts.eventEstimateRATField;
eventInfo.thresholdBasis = opts.thresholdBasis;
eventInfo.thresholdCenter = opts.thresholdCenter;
eventInfo.thresholdMode = opts.thresholdMode;
eventInfo.detectMethod = opts.detectMethod;
eventInfo.minTimeBtEvents = opts.minTimeBtEvents;
eventInfo.daysMode = opts.daysMode;
eventInfo.daysToUse = traceFieldsToDates(traceFields);
eventInfo.availableTraceFields = allTraceFields;
eventInfo.traceFields = traceFields;
eventInfo.outputStructs = cell(size(opts.thresholds));
eventInfo.outputFields = cell(numel(opts.thresholds), numel(traceFields));
eventInfo.eventCounts = cell(numel(opts.thresholds), numel(traceFields));

for thresholdNo = 1:numel(opts.thresholds)
    threshold = opts.thresholds(thresholdNo);

    if opts.overwrite && numel(opts.thresholds) == 1
        outputStruct = 'Ca_peaks';
    else
        outputStruct = sprintf('%s%s', opts.outputPrefix, thresholdLabel(threshold));
    end

    if ~isfield(RAT, outputStruct) || ~isstruct(RAT.(outputStruct))
        RAT.(outputStruct) = struct();
    end

    eventInfo.outputStructs{thresholdNo} = outputStruct;

    for fieldNo = 1:numel(traceFields)
        traceField = traceFields{fieldNo};
        peakField = regexprep(traceField, '^CA_traces_', 'CA_peaks_');
        traces = RAT.Ca_traces.(traceField);

        if isempty(traces)
            RAT.(outputStruct).(peakField) = [];
            eventInfo.outputFields{thresholdNo, fieldNo} = peakField;
            eventInfo.eventCounts{thresholdNo, fieldNo} = [];
            continue;
        end

        traces = double(traces);
        if size(traces, 1) > size(traces, 2)
            warning('eventsfromtrace:PossibleTransposedTraces', ...
                ['%s has more rows than columns. Expected ' ...
                'cells x frames; leaving orientation unchanged.'], traceField);
        end

        switch opts.sourceMode
            case 'tracePeaks'
                [peaksBinary, peaksArray] = computeSignalPeaks(traces, ...
                    'makePlots', opts.makePlots, ...
                    'makeSummaryPlots', 0, ...
                    'waitbarOn', opts.waitbarOn, ...
                    'outputInfo', 0, ...
                    'parallel', opts.parallel, ...
                    'numStdsForThresh', threshold, ...
                    'detectMethod', opts.detectMethod, ...
                    'minTimeBtEvents', opts.minTimeBtEvents);

            case 'deconvolvedAmplitude'
                eventEstimates = getEventEstimateMatrix(RAT, traceField, opts);
                eventEstimates = double(eventEstimates);
                rowIdx = getEventEstimateRows(traceField, opts);
                if ~isempty(rowIdx)
                    eventEstimates = eventEstimates(rowIdx, :);
                end
                if ~isequal(size(eventEstimates), size(traces))
                    error('eventsfromtrace:TraceEventSizeMismatch', ...
                        ['%s is %d x %d, but its deconvolved event estimate is %d x %d. ' ...
                        'Expected manually accepted cells x frames in the same order.'], ...
                        traceField, size(traces, 1), size(traces, 2), ...
                        size(eventEstimates, 1), size(eventEstimates, 2));
                end
                [peaksBinary, peaksArray] = thresholdDeconvolvedEvents( ...
                    eventEstimates, threshold, opts);
        end

        peaks = formatPeakOutput(peaksBinary, RAT, traceField, opts);

        RAT.(outputStruct).(peakField) = peaks;
        RAT.(outputStruct).([peakField '_binary']) = peaksBinary;
        RAT.(outputStruct).([peakField '_array']) = peaksArray;
        eventInfo.outputFields{thresholdNo, fieldNo} = peakField;
        eventInfo.eventCounts{thresholdNo, fieldNo} = sum(peaksBinary, 2);

        fprintf('%s.%s: threshold %.3g SD, %d cells, %d total events\n', ...
            outputStruct, peakField, threshold, size(peaksBinary, 1), sum(peaksBinary(:)));
    end
end

if opts.makeTraceFigure
    eventInfo.traceFigure = makeTraceOverlayFigure(RAT, opts, eventInfo);
else
    eventInfo.traceFigure = [];
end
end

function opts = parseInputs(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('eventsfromtrace:BadInputs', 'Options must be name/value pairs.');
end

for i = 1:2:numel(varargin)
    name = varargin{i};
    value = varargin{i + 1};
    if ~ischar(name) && ~isstring(name)
        error('eventsfromtrace:BadOptionName', 'Option names must be strings.');
    end
    name = char(name);
    if ~isfield(opts, name)
        error('eventsfromtrace:UnknownOption', 'Unknown option: %s', name);
    end
    opts.(name) = value;
end

opts.thresholds = opts.thresholds(:)';
if isempty(opts.thresholds) || ~isnumeric(opts.thresholds)
    error('eventsfromtrace:BadThresholds', 'thresholds must be a numeric vector.');
end

opts.daysMode = validatestring(opts.daysMode, {'last3toAn', 'all'}, ...
    'eventsfromtrace', 'daysMode');
opts.outputFormat = validatestring(opts.outputFormat, {'times', 'binary'}, ...
    'eventsfromtrace', 'outputFormat');
opts.sourceMode = validatestring(opts.sourceMode, {'deconvolvedAmplitude', 'tracePeaks'}, ...
    'eventsfromtrace', 'sourceMode');
opts.thresholdBasis = validatestring(opts.thresholdBasis, {'allStd', 'positiveStd'}, ...
    'eventsfromtrace', 'thresholdBasis');
opts.thresholdCenter = validatestring(opts.thresholdCenter, {'zero', 'mean'}, ...
    'eventsfromtrace', 'thresholdCenter');
opts.thresholdMode = validatestring(opts.thresholdMode, {'sd', 'percentile'}, ...
    'eventsfromtrace', 'thresholdMode');
end


function eventEstimates = getEventEstimateMatrix(RAT, traceField, opts)
eventStruct = opts.eventEstimateStruct;
if isempty(eventStruct)
    if isfield(RAT, opts.eventEstimateRATField) && isstruct(RAT.(opts.eventEstimateRATField))
        eventStruct = RAT.(opts.eventEstimateRATField);
    else
        error('eventsfromtrace:MissingEventEstimateStruct', ...
            ['sourceMode ''deconvolvedAmplitude'' requires CNMF-E deconvolved ' ...
            'event estimates with amplitudes. Pass them as eventEstimateStruct, ' ...
            'or store them in RAT.%s. Timestamp-only Ca_peaks fields cannot be ' ...
            'amplitude-thresholded.'], opts.eventEstimateRATField);
    end
end

dateLabel = regexprep(traceField, '^CA_traces_', '');
candidateFields = { ...
    traceField, ...
    ['CA_peaks_' dateLabel], ...
    ['CA_deconvolved_' dateLabel], ...
    ['CA_events_' dateLabel], ...
    ['extractedPeaks_' dateLabel]};

eventField = '';
for i = 1:numel(candidateFields)
    if isfield(eventStruct, candidateFields{i})
        eventField = candidateFields{i};
        break;
    end
end

if isempty(eventField)
    error('eventsfromtrace:MissingEventEstimateField', ...
        ['Could not find a deconvolved event estimate field for %s. Looked for: %s. ' ...
        'Use fields named CA_traces_DATE, CA_peaks_DATE, CA_deconvolved_DATE, ' ...
        'CA_events_DATE, or extractedPeaks_DATE.'], ...
        traceField, strjoin(candidateFields, ', '));
end

eventEstimates = eventStruct.(eventField);
end


function rowIdx = getEventEstimateRows(traceField, opts)
rowIdx = [];
if ~isempty(opts.eventEstimateRowsStruct)
    dateLabel = regexprep(traceField, '^CA_traces_', '');
    candidateFields = { ...
        traceField, ...
        ['CA_peaks_' dateLabel], ...
        ['CA_deconvolved_' dateLabel], ...
        ['date_' dateLabel]};

    for i = 1:numel(candidateFields)
        if isfield(opts.eventEstimateRowsStruct, candidateFields{i})
            rowIdx = opts.eventEstimateRowsStruct.(candidateFields{i});
            rowIdx = rowIdx(:)';
            return;
        end
    end
end

if ~isempty(opts.eventEstimateRows)
    rowIdx = opts.eventEstimateRows(:)';
end
end


function [peaksBinary, peaksArray] = thresholdDeconvolvedEvents(eventEstimates, threshold, opts)
if size(eventEstimates, 1) > size(eventEstimates, 2)
    warning('eventsfromtrace:PossibleTransposedEvents', ...
        ['Deconvolved event matrix has more rows than columns. Expected cells x frames; ' ...
        'leaving orientation unchanged.']);
end

peaksBinary = false(size(eventEstimates));
peaksArray = cell(size(eventEstimates, 1), 1);

for cellNo = 1:size(eventEstimates, 1)
    eventTrace = eventEstimates(cellNo, :);
    cutoff = deconvolvedEventCutoff(eventTrace, threshold, opts);
    if ~isfinite(cutoff) || cutoff <= 0
        peaksArray{cellNo} = [];
        continue;
    end

    aboveThreshold = eventTrace > cutoff;
    candidateFrames = chooseEventMaximaFromThresholdRuns(eventTrace, aboveThreshold);
    candidateFrames = enforceMinimumEventSpacing( ...
        eventTrace, candidateFrames, opts.minTimeBtEvents);
    peaksBinary(cellNo, candidateFrames) = true;
    peaksArray{cellNo} = candidateFrames(:)';
end
end


function cutoff = deconvolvedEventCutoff(eventTrace, threshold, opts)
eventTrace = eventTrace(:)';
eventTrace = eventTrace(isfinite(eventTrace));
if isempty(eventTrace)
    cutoff = NaN;
    return;
end

switch opts.thresholdBasis
    case 'allStd'
        basis = eventTrace;
    case 'positiveStd'
        basis = eventTrace(eventTrace > 0);
        if numel(basis) < 2
            basis = eventTrace;
        end
end

switch opts.thresholdMode
    case 'percentile'
        positiveAmplitudes = eventTrace(eventTrace > 0);
        if isempty(positiveAmplitudes)
            cutoff = NaN;
        else
            cutoff = prctile(positiveAmplitudes, threshold);
        end
        return;
end

sigma = std(basis, 0, 'omitnan');
if ~isfinite(sigma) || sigma == 0
    sigma = std(eventTrace, 0, 'omitnan');
end

switch opts.thresholdCenter
    case 'zero'
        center = 0;
    case 'mean'
        center = mean(basis, 'omitnan');
end
cutoff = center + threshold .* sigma;
end


function eventFrames = chooseEventMaximaFromThresholdRuns(eventTrace, aboveThreshold)
eventFrames = [];
runStarts = find(diff([false aboveThreshold]) == 1);
runStops = find(diff([aboveThreshold false]) == -1);

for runNo = 1:numel(runStarts)
    frames = runStarts(runNo):runStops(runNo);
    [~, localIdx] = max(eventTrace(frames));
    eventFrames(end + 1) = frames(localIdx); %#ok<AGROW>
end
end


function keptFrames = enforceMinimumEventSpacing(eventTrace, eventFrames, minSpacing)
if isempty(eventFrames)
    keptFrames = eventFrames;
    return;
end

eventFrames = sort(eventFrames(:)');
keptFrames = [];
for frameNo = eventFrames
    if isempty(keptFrames) || frameNo - keptFrames(end) >= minSpacing
        keptFrames(end + 1) = frameNo; %#ok<AGROW>
    elseif eventTrace(frameNo) > eventTrace(keptFrames(end))
        keptFrames(end) = frameNo;
    end
end
end


function peaks = formatPeakOutput(peaksBinary, RAT, traceField, opts)
switch opts.outputFormat
    case 'binary'
        peaks = peaksBinary;
    case 'times'
        timestampField = regexprep(traceField, '^CA_traces_', 'CA_time_');
        if isfield(RAT, 'Ca_ts') && isstruct(RAT.Ca_ts) && isfield(RAT.Ca_ts, timestampField)
            peaks = converttotime(peaksBinary, RAT.Ca_ts.(timestampField));
        else
            warning('eventsfromtrace:MissingTimestamps', ...
                ['Could not find RAT.Ca_ts.%s. Writing event times from frameRate ' ...
                'instead of timestamp file.'], timestampField);
            peaks = binaryPeaksToFrameTimes(peaksBinary, opts.frameRate);
        end
end
end


function peakTimes = binaryPeaksToFrameTimes(peaksBinary, frameRate)
peakTimes = NaN(size(peaksBinary));
for cellNo = 1:size(peaksBinary, 1)
    peakFrames = find(peaksBinary(cellNo, :) > 0);
    if isempty(peakFrames)
        continue;
    end
    times = (peakFrames - 1) ./ frameRate;
    peakTimes(cellNo, 1:numel(times)) = times;
end

validCols = find(any(~isnan(peakTimes), 1));
if isempty(validCols)
    peakTimes = NaN(size(peaksBinary, 1), 0);
else
    peakTimes = peakTimes(:, 1:max(validCols));
end
end


function traceFields = selectTraceFieldsForDays(RAT, traceFields, opts)
if ~isempty(opts.daysToUse)
    daysToUse = normalizeDayList(opts.daysToUse);
else
    switch opts.daysMode
        case 'all'
            return;
        case 'last3toAn'
            if ~isfield(RAT, 'An') || isempty(RAT.An)
                error('eventsfromtrace:MissingAn', ...
                    'RAT.An is required when daysMode is ''last3toAn''.');
            end

            traceDates = traceFieldsToDates(traceFields);
            [traceDates, sortIdx] = sort(traceDates);
            traceFields = traceFields(sortIdx);

            anDate = normalizeDayString(RAT.An);
            anIdx = find(strcmp(traceDates, anDate), 1);
            if isempty(anIdx)
                error('eventsfromtrace:AnNotFound', ...
                    'RAT.An (%s) was not found among RAT.Ca_traces fields.', anDate);
            end

            dayIdx = max(1, anIdx - 2):anIdx;
            traceFields = traceFields(dayIdx);
            return;
    end
end

traceDates = traceFieldsToDates(traceFields);
keep = ismember(traceDates, daysToUse);
traceFields = traceFields(keep);
if isempty(traceFields)
    error('eventsfromtrace:NoSelectedTraceFields', ...
        'No CA_traces fields matched the requested daysToUse.');
end

[~, sortIdx] = sort(traceFieldsToDates(traceFields));
traceFields = traceFields(sortIdx);
end


function days = normalizeDayList(days)
if ischar(days) || isstring(days)
    days = cellstr(days);
elseif ~iscell(days)
    error('eventsfromtrace:BadDaysToUse', ...
        'daysToUse must be a string, string array, or cell array of date strings.');
end

days = cellfun(@normalizeDayString, days, 'UniformOutput', false);
end


function day = normalizeDayString(day)
day = char(day);
day = regexprep(day, '^CA_traces_', '');
day = regexprep(day, '^CA_peaks_', '');
day = strrep(day, '-', '_');
day = strtrim(day);
end


function dates = traceFieldsToDates(traceFields)
dates = regexprep(traceFields, '^CA_traces_', '');
end


function figHandle = makeTraceOverlayFigure(RAT, opts, eventInfo)
if isempty(opts.figureThreshold)
    plotThreshold = opts.thresholds(1);
else
    plotThreshold = opts.figureThreshold;
end

if isempty(opts.figureField)
    traceField = eventInfo.traceFields{end};
else
    traceField = opts.figureField;
end

if ~isfield(RAT.Ca_traces, traceField)
    error('eventsfromtrace:BadFigureField', ...
        'figureField %s was not found in RAT.Ca_traces.', traceField);
end

outputStruct = sprintf('%s%s', opts.outputPrefix, thresholdLabel(plotThreshold));
if opts.overwrite && numel(opts.thresholds) == 1
    outputStruct = 'Ca_peaks';
end
peakField = regexprep(traceField, '^CA_traces_', 'CA_peaks_');

if ~isfield(RAT, outputStruct) || ~isfield(RAT.(outputStruct), peakField)
    error('eventsfromtrace:MissingFigurePeaks', ...
        'Could not find %s.%s. Include figureThreshold in thresholds.', outputStruct, peakField);
end

traces = double(RAT.Ca_traces.(traceField));
binaryPeakField = [peakField '_binary'];
if isfield(RAT.(outputStruct), binaryPeakField)
    peaksBinary = RAT.(outputStruct).(binaryPeakField);
else
    peaksBinary = RAT.(outputStruct).(peakField);
end

nCells = size(traces, 1);
nFrames = size(traces, 2);
if isempty(opts.figureCells)
    eventCounts = sum(peaksBinary, 2);
    validCells = find(eventCounts > 0);
    if isempty(validCells)
        validCells = 1:min(nCells, 5);
    else
        [~, sortIdx] = sort(eventCounts(validCells), 'descend');
        validCells = validCells(sortIdx);
        validCells = validCells(1:min(numel(validCells), 5));
    end
    cellsToPlot = validCells(:)';
else
    cellsToPlot = opts.figureCells(:)';
end
cellsToPlot = cellsToPlot(cellsToPlot >= 1 & cellsToPlot <= nCells);
if isempty(cellsToPlot)
    error('eventsfromtrace:NoFigureCells', 'No valid figureCells were supplied.');
end

if isempty(opts.figureWindow)
    frameWindow = [1 nFrames];
else
    frameWindow = opts.figureWindow;
end
frameWindow(1) = max(1, round(frameWindow(1)));
frameWindow(2) = min(nFrames, round(frameWindow(2)));
if frameWindow(2) <= frameWindow(1)
    error('eventsfromtrace:BadFigureWindow', 'figureWindow must be [startFrame endFrame].');
end

frames = frameWindow(1):frameWindow(2);
timeSeconds = (frames - 1) ./ opts.frameRate;
figHandle = figure('Color', 'w', 'Name', sprintf('%s %s %.3g SD', traceField, peakField, plotThreshold));

nPlot = numel(cellsToPlot);
for plotNo = 1:nPlot
    cellNo = cellsToPlot(plotNo);
    subplot(nPlot, 1, plotNo);
    trace = traces(cellNo, frames);
    plot(timeSeconds, trace, 'k', 'LineWidth', 1); hold on;

    ylabel(sprintf('Cell %d', cellNo));
    box off;
    set(gca, 'TickDir', 'out');

    yLimits = ylim;
    yRange = diff(yLimits);
    if yRange == 0
        yRange = max(abs(yLimits(1)), 1);
    end

    peakFrames = find(peaksBinary(cellNo, :) > 0);
    peakFrames = peakFrames(peakFrames >= frameWindow(1) & peakFrames <= frameWindow(2));
    if ~isempty(peakFrames)
        peakTimes = (peakFrames - 1) ./ opts.frameRate;
        peakY = traces(cellNo, peakFrames) + 0.08 .* yRange;
        scatter(peakTimes, peakY, 28, [0.85 0.1 0.1], 'v', 'filled');
        ylim([yLimits(1), max(yLimits(2), max(peakY) + 0.05 .* yRange)]);
    end

    if plotNo == 1
        title(sprintf('%s, %.3g SD threshold, %s', traceField, plotThreshold, strrep(outputStruct, '_', '\_')));
    end
    if plotNo < nPlot
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time (s)');
    end
end

if ~isempty(opts.saveFigurePath)
    saveas(figHandle, opts.saveFigurePath);
    fprintf('Saved trace overlay figure to %s\n', opts.saveFigurePath);
end
end

function label = thresholdLabel(threshold)
label = regexprep(sprintf('%.12g', threshold), '\.', 'p');
end
