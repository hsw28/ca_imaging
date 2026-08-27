function figHandle = plotCnmfeEventExamples(sessionFolder, varargin)
%PLOTCNMFEEVENTEXAMPLES Plot CNMF-E traces with inferred events.
%
%   figHandle = plotCnmfeEventExamples(sessionFolder)
%
%   Makes a three-column panel for selected manually accepted cells:
%       1) CNMF-E extracted fluorescence trace with inferred event times.
%       2) Corresponding CNMF-E deconvolved activity estimate.
%       3) Event-triggered extracted fluorescence snippets.
%
%   Cell indices are accepted-cell indices, matching rows in rat.Ca_traces.
%   The function maps them back to the original CNMF-E row indices using
%   validCNMFE from the sorted file.
%
%   The function expects sessionFolder, or a folder beneath it, to contain:
%       * *_cnmfeAnalysis.mat
%       * *_cnmfeAnalysisSorted*.mat
%       * timeStamps.csv
%
%   Example:
%       plotCnmfeEventExamples(['/Users/Hannah/Library/CloudStorage/' ...
%           'OneDrive-NorthwesternUniversity/Desktop/videos/eyeblink/' ...
%           '031423/2023_05_19/14_45_47'], ...
%           'windowSeconds', [120 210], ...
%           'saveFigurePath', 'rat0314_2023_05_19_cnmfe_events.png');
%
%   Useful options:
%       cells           Accepted-cell indices to plot. Default: auto-select
%                       moderate event-rate cells near the 40th to 70th
%                       percentiles.
%       windowSeconds   Time window [start end] in seconds. Default: auto-select
%                       a 90 s window with a moderate number of events.
%       taskTimes       Optional task onset times in seconds for light shading.
%       usTimes         Optional US/shock times in seconds for tick marks.

opts.cells = [];
opts.window = [];
opts.windowSeconds = [];
opts.defaultWindowSec = 90;
opts.frameRate = 7.5;
opts.maxCells = 5;
opts.autoPercentiles = [40 47.5 55 62.5 70];
opts.minWindowEvents = 4;
opts.preferredWindowEvents = [4 15];
opts.traceSource = 'extractedSignalsEst';
opts.eventSource = 'extractedPeaks';
opts.analysisFile = '';
opts.sortedFile = '';
opts.timestampFile = '';
opts.maxEventMarkers = Inf;
opts.snippetWindowSec = [-2 4];
opts.maxSnippets = 5;
opts.taskTimes = [];
opts.taskWindowSec = [0 2];
opts.usTimes = [];
opts.eventMarkerStyle = 'tick';
opts.eventMarkerColor = [0.9 0 0];
opts.eventMarkerSizeFrac = 0.08;
opts.showAggregate = true;
opts.aggregateCells = 'figure';
opts.aggregateShowSubset = 0;
opts.aggregateSubsetSeed = 1;
opts.aggregateSemColor = [0.75 0.75 0.75];
opts.aggregateSemAlpha = 0.45;
opts.maxTraceRobustZ = 8;
opts.maxDerivativeRobustZ = 12;
opts.saveFigurePath = '';
opts.vectorFormat = 'pdf';
opts.preserveEditableText = true;
opts = parseInputs(opts, varargin{:});

[analysisFile, sortedFile, timestampFile] = findSessionFiles(sessionFolder, opts);
analysisData = load(analysisFile);
sortedData = load(sortedFile);
timestamps = readSessionTimestamps(timestampFile);

if ~isfield(analysisData, 'cnmfeAnalysisOutput')
    error('plotCnmfeEventExamples:MissingCnmfeOutput', ...
        '%s does not contain cnmfeAnalysisOutput.', analysisFile);
end
if ~isfield(sortedData, 'validCNMFE')
    error('plotCnmfeEventExamples:MissingValidCnmfe', ...
        '%s does not contain validCNMFE.', sortedFile);
end

cnmfe = analysisData.cnmfeAnalysisOutput;
if ~isfield(cnmfe, opts.traceSource)
    error('plotCnmfeEventExamples:MissingTraceSource', ...
        'cnmfeAnalysisOutput.%s was not found.', opts.traceSource);
end
if ~isfield(cnmfe, opts.eventSource)
    error('plotCnmfeEventExamples:MissingEventSource', ...
        'cnmfeAnalysisOutput.%s was not found.', opts.eventSource);
end

acceptedOriginalIdx = find(sortedData.validCNMFE == 1);
if isempty(acceptedOriginalIdx)
    error('plotCnmfeEventExamples:NoAcceptedCells', ...
        'No accepted cells were found in validCNMFE.');
end

tracesAll = double(cnmfe.(opts.traceSource));
eventsAll = double(cnmfe.(opts.eventSource));

if size(tracesAll, 1) ~= size(eventsAll, 1)
    error('plotCnmfeEventExamples:TraceEventSizeMismatch', ...
        'Trace and event matrices have different cell counts.');
end

acceptedOriginalIdx = acceptedOriginalIdx(acceptedOriginalIdx <= size(tracesAll, 1));
traces = tracesAll(acceptedOriginalIdx, :);
events = eventsAll(acceptedOriginalIdx, :);

if isempty(opts.cells)
    cellsToPlot = selectModerateEventCells(events, opts);
else
    cellsToPlot = opts.cells(:)';
end
cellsToPlot = cellsToPlot(cellsToPlot >= 1 & cellsToPlot <= size(traces, 1));
if isempty(cellsToPlot)
    error('plotCnmfeEventExamples:NoValidCells', ...
        'No valid accepted-cell indices were supplied.');
end

nFrames = min(size(traces, 2), size(events, 2));
frameWindow = chooseFrameWindow(events, cellsToPlot, timestamps, nFrames, opts);
frameWindow(1) = max(1, round(frameWindow(1)));
frameWindow(2) = min(nFrames, round(frameWindow(2)));
if frameWindow(2) <= frameWindow(1)
    error('plotCnmfeEventExamples:BadWindow', ...
        'window must be [startFrame endFrame].');
end

if isempty(opts.cells)
    cellsToPlot = selectModerateWindowEventCells(traces, events, frameWindow, opts);
    cellsToPlot = cellsToPlot(cellsToPlot >= 1 & cellsToPlot <= size(traces, 1));
    if isempty(cellsToPlot)
        error('plotCnmfeEventExamples:NoWindowMatchedCells', ...
            'No cells had at least %d events in the plotted window.', opts.minWindowEvents);
    end
end

frames = frameWindow(1):frameWindow(2);
timeSeconds = frameToSeconds(frames, timestamps, opts.frameRate);
nExampleRows = numel(cellsToPlot);
nRows = nExampleRows + double(opts.showAggregate);

figHandle = figure('Color', 'w', 'Name', 'CNMF-E traces and inferred events', ...
    'Position', [100 100 1300 230 * nRows]);
if opts.preserveEditableText
    set(figHandle, 'Renderer', 'painters');
end

for plotNo = 1:nExampleRows
    cellNo = cellsToPlot(plotNo);

    trace = traces(cellNo, frames);
    eventTrace = events(cellNo, frames);
    eventFramesAll = find(events(cellNo, :) > 0);
    eventFrames = eventFramesAll(eventFramesAll >= frameWindow(1) & eventFramesAll <= frameWindow(2));
    nEventsInWindow = numel(eventFrames);

    subplot(nRows, 3, (plotNo - 1) * 3 + 1);
    plot(timeSeconds, trace, 'k', 'LineWidth', 1); hold on;
    ylabel('\DeltaF/F');
    addRowLabel(gca, sprintf('Cell %d (n = %d events)', cellNo, nEventsInWindow));
    markTaskTimes(gca, opts);

    yLimits = ylim;
    if ~isempty(eventFrames)
        eventTimes = frameToSeconds(eventFrames, timestamps, opts.frameRate);
        [eventTimes, eventFrames] = thinEvents(eventTimes, eventFrames, opts.maxEventMarkers);
        eventY = traces(cellNo, eventFrames);
        plotTraceEventMarkers(eventTimes, eventY, yLimits, opts);
    end
    ylim(yLimits);

    box off;
    set(gca, 'TickDir', 'out');
    if plotNo == 1
        title('CNMF-E extracted \DeltaF/F');
    end
    if plotNo < nExampleRows
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time (s)');
    end

    subplot(nRows, 3, (plotNo - 1) * 3 + 2);
    plot(timeSeconds, eventTrace, 'Color', [0.85 0.1 0.1], 'LineWidth', 1);
    ylabel('Deconvolved activity');
    markTaskTimes(gca, opts);
    box off;
    set(gca, 'TickDir', 'out');
    if plotNo == 1
        title('CNMF-E deconvolved activity');
    end
    if plotNo < nExampleRows
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time (s)');
    end

    subplot(nRows, 3, (plotNo - 1) * 3 + 3);
    plotEventTriggeredSnippets(traces(cellNo, :), eventFrames, opts, plotNo == nExampleRows);
    if plotNo == 1
        title('Event-triggered \DeltaF/F');
    end
end

if opts.showAggregate
    subplot(nRows, 3, (nExampleRows * 3 + 1):(nExampleRows * 3 + 3));
    plotPopulationEventTriggered(traces, events, cellsToPlot, opts);
end

if ~isempty(opts.saveFigurePath)
    saveFigureOutputs(figHandle, opts.saveFigurePath, opts.vectorFormat);
end
end


function opts = parseInputs(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('plotCnmfeEventExamples:BadInputs', ...
        'Options must be name/value pairs.');
end

for i = 1:2:numel(varargin)
    name = varargin{i};
    value = varargin{i + 1};
    if ~ischar(name) && ~isstring(name)
        error('plotCnmfeEventExamples:BadOptionName', ...
            'Option names must be strings.');
    end
    name = char(name);
    if ~isfield(opts, name)
        error('plotCnmfeEventExamples:UnknownOption', ...
            'Unknown option: %s', name);
    end
    opts.(name) = value;
end

opts.traceSource = validatestring(opts.traceSource, ...
    {'extractedSignalsEst', 'extractedSignals'}, ...
    'plotCnmfeEventExamples', 'traceSource');
opts.eventSource = validatestring(opts.eventSource, ...
    {'extractedPeaks'}, ...
    'plotCnmfeEventExamples', 'eventSource');
opts.aggregateCells = validatestring(opts.aggregateCells, ...
    {'figure', 'allAccepted'}, ...
    'plotCnmfeEventExamples', 'aggregateCells');
opts.vectorFormat = validatestring(opts.vectorFormat, ...
    {'pdf', 'svg'}, ...
    'plotCnmfeEventExamples', 'vectorFormat');
opts.preserveEditableText = logical(opts.preserveEditableText);
end


function cellsToPlot = selectModerateEventCells(events, opts)
eventCounts = sum(events > 0, 2);
activeCells = find(eventCounts > 0);
if isempty(activeCells)
    cellsToPlot = 1:min(size(events, 1), opts.maxCells);
    return;
end

counts = eventCounts(activeCells);
percentiles = opts.autoPercentiles(:)';
percentiles = percentiles(1:min(numel(percentiles), opts.maxCells));
targets = prctile(counts, percentiles);

cellsToPlot = [];
for target = targets
    [~, order] = sort(abs(counts - target), 'ascend');
    for idx = order(:)'
        candidate = activeCells(idx);
        if ~ismember(candidate, cellsToPlot)
            cellsToPlot(end + 1) = candidate; %#ok<AGROW>
            break;
        end
    end
end

if numel(cellsToPlot) < opts.maxCells
    medianTarget = median(counts);
    [~, order] = sort(abs(counts - medianTarget), 'ascend');
    for idx = order(:)'
        candidate = activeCells(idx);
        if ~ismember(candidate, cellsToPlot)
            cellsToPlot(end + 1) = candidate; %#ok<AGROW>
        end
        if numel(cellsToPlot) >= opts.maxCells
            break;
        end
    end
end
end


function cellsToPlot = selectModerateWindowEventCells(traces, events, frameWindow, opts)
windowCounts = sum(events(:, frameWindow(1):frameWindow(2)) > 0, 2);
eligibleCells = find(windowCounts >= opts.minWindowEvents);

if isempty(eligibleCells)
    cellsToPlot = [];
    return;
end

cleanCells = eligibleCells(isCleanTraceWindow(traces, eligibleCells, frameWindow, opts));
if ~isempty(cleanCells)
    eligibleCells = cleanCells;
end

preferred = eligibleCells( ...
    windowCounts(eligibleCells) >= opts.preferredWindowEvents(1) & ...
    windowCounts(eligibleCells) <= opts.preferredWindowEvents(2));
if isempty(preferred)
    preferred = eligibleCells;
end

counts = windowCounts(preferred);
percentiles = opts.autoPercentiles(:)';
percentiles = percentiles(1:min(numel(percentiles), opts.maxCells));
targets = prctile(counts, percentiles);

cellsToPlot = [];
for target = targets
    [~, order] = sort(abs(counts - target), 'ascend');
    for idx = order(:)'
        candidate = preferred(idx);
        if ~ismember(candidate, cellsToPlot)
            cellsToPlot(end + 1) = candidate; %#ok<AGROW>
            break;
        end
    end
end

if numel(cellsToPlot) < opts.maxCells
    target = median(opts.preferredWindowEvents);
    [~, order] = sort(abs(windowCounts(eligibleCells) - target), 'ascend');
    for idx = order(:)'
        candidate = eligibleCells(idx);
        if ~ismember(candidate, cellsToPlot)
            cellsToPlot(end + 1) = candidate; %#ok<AGROW>
        end
        if numel(cellsToPlot) >= opts.maxCells
            break;
        end
    end
end
end


function cleanMask = isCleanTraceWindow(traces, cellIdx, frameWindow, opts)
cleanMask = true(size(cellIdx));
frames = frameWindow(1):frameWindow(2);

for i = 1:numel(cellIdx)
    trace = traces(cellIdx(i), frames);
    trace = trace(:)';
    center = median(trace, 'omitnan');
    scale = 1.4826 .* median(abs(trace - center), 'omitnan');
    if ~isfinite(scale) || scale == 0
        scale = std(trace, 0, 'omitnan');
    end
    if ~isfinite(scale) || scale == 0
        continue;
    end

    robustZ = abs((trace - center) ./ scale);
    dTrace = diff(trace);
    dCenter = median(dTrace, 'omitnan');
    dScale = 1.4826 .* median(abs(dTrace - dCenter), 'omitnan');
    if ~isfinite(dScale) || dScale == 0
        dScale = std(dTrace, 0, 'omitnan');
    end
    if ~isfinite(dScale) || dScale == 0
        dRobustZ = 0;
    else
        dRobustZ = max(abs((dTrace - dCenter) ./ dScale), [], 'omitnan');
    end

    cleanMask(i) = max(robustZ, [], 'omitnan') <= opts.maxTraceRobustZ && ...
        dRobustZ <= opts.maxDerivativeRobustZ;
end
end


function frameWindow = chooseFrameWindow(events, cellsToPlot, timestamps, nFrames, opts)
if ~isempty(opts.window)
    frameWindow = opts.window;
    return;
end

if ~isempty(opts.windowSeconds)
    frameWindow = secondsToFrames(opts.windowSeconds, timestamps, opts.frameRate);
    return;
end

windowFrames = max(2, round(opts.defaultWindowSec .* opts.frameRate));
if windowFrames >= nFrames
    frameWindow = [1 nFrames];
    return;
end

stepFrames = max(1, round(windowFrames ./ 4));
starts = 1:stepFrames:(nFrames - windowFrames + 1);
targetEvents = max(3, 5 .* numel(cellsToPlot));
bestScore = Inf;
bestStart = starts(1);

for startFrame = starts
    stopFrame = startFrame + windowFrames - 1;
    eventChunk = events(cellsToPlot, startFrame:stopFrame) > 0;
    eventCount = sum(eventChunk(:));
    if eventCount == 0
        score = Inf;
    else
        score = abs(eventCount - targetEvents);
    end
    if score < bestScore
        bestScore = score;
        bestStart = startFrame;
    end
end

frameWindow = [bestStart, bestStart + windowFrames - 1];
end


function frames = secondsToFrames(secondsWindow, timestamps, frameRate)
secondsWindow = secondsWindow(:)';
if numel(secondsWindow) ~= 2
    error('plotCnmfeEventExamples:BadWindowSeconds', ...
        'windowSeconds must be [startSeconds endSeconds].');
end

if ~isempty(timestamps)
    frames = zeros(1, 2);
    for i = 1:2
        [~, timestampIdx] = min(abs(timestamps - secondsWindow(i)));
        frames(i) = max(1, round(timestampIdx ./ 2));
    end
else
    frames = round(secondsWindow .* frameRate) + 1;
end
end


function [eventTimes, eventFrames] = thinEvents(eventTimes, eventFrames, maxEvents)
if numel(eventTimes) <= maxEvents
    return;
end

keepIdx = round(linspace(1, numel(eventTimes), maxEvents));
eventTimes = eventTimes(keepIdx);
eventFrames = eventFrames(keepIdx);
end


function plotTraceEventMarkers(eventTimes, eventY, yLimits, opts)
markerStyle = validatestring(opts.eventMarkerStyle, {'tick', 'circle'}, ...
    'plotCnmfeEventExamples', 'eventMarkerStyle');

switch markerStyle
    case 'circle'
        plot(eventTimes, eventY, 'o', 'Color', opts.eventMarkerColor, ...
            'MarkerFaceColor', opts.eventMarkerColor, 'MarkerSize', 4, ...
            'LineStyle', 'none');
    case 'tick'
        yRange = diff(yLimits);
        if yRange == 0
            yRange = max(abs(yLimits(1)), 1);
        end
        halfHeight = opts.eventMarkerSizeFrac .* yRange ./ 2;
        for eventNo = 1:numel(eventTimes)
            plot([eventTimes(eventNo) eventTimes(eventNo)], ...
                [eventY(eventNo) - halfHeight, eventY(eventNo) + halfHeight], ...
                'Color', opts.eventMarkerColor, 'LineWidth', 1.6);
        end
end
end


function addRowLabel(ax, labelText)
text(ax, -0.32, 0.5, labelText, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Rotation', 90, ...
    'FontWeight', 'normal', ...
    'Clipping', 'off');
end


function markTaskTimes(ax, opts)
axes(ax); %#ok<LAXES>
hold on;
yLimits = ylim(ax);

if ~isempty(opts.taskTimes)
    taskTimes = opts.taskTimes(:)';
    for t = taskTimes
        shadePatch = patch(ax, t + [opts.taskWindowSec(1) opts.taskWindowSec(2) ...
            opts.taskWindowSec(2) opts.taskWindowSec(1)], ...
            [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
            [0.55 0.7 1.0], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
        try
            uistack(shadePatch, 'bottom');
        catch
        end
    end
end

if ~isempty(opts.usTimes)
    usTimes = opts.usTimes(:)';
    for t = usTimes
        xline(ax, t, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.75);
    end
end
ylim(ax, yLimits);
end


function plotEventTriggeredSnippets(trace, eventFrames, opts, showXLabel)
if isempty(eventFrames)
    text(0.5, 0.5, 'No events in window', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

eventFrames = eventFrames(:)';
if numel(eventFrames) > opts.maxSnippets
    keepIdx = round(linspace(1, numel(eventFrames), opts.maxSnippets));
    eventFrames = eventFrames(keepIdx);
end

[snippets, snippetTime] = extractEventSnippets(trace, eventFrames, opts);
plotEventTriggeredMatrix(snippets, snippetTime, showXLabel);
end


function [snippets, snippetTime] = extractEventSnippets(trace, eventFrames, opts)
preFrames = round(abs(opts.snippetWindowSec(1)) .* opts.frameRate);
postFrames = round(opts.snippetWindowSec(2) .* opts.frameRate);
snippetFrames = -preFrames:postFrames;
snippetTime = snippetFrames ./ opts.frameRate;

snippets = NaN(numel(eventFrames), numel(snippetFrames));
for i = 1:numel(eventFrames)
    idx = eventFrames(i) + snippetFrames;
    valid = idx >= 1 & idx <= numel(trace);
    snippets(i, valid) = trace(idx(valid));

    baselineIdx = snippetTime < 0;
    baseline = median(snippets(i, baselineIdx), 'omitnan');
    snippets(i, :) = snippets(i, :) - baseline;
end
end


function plotEventTriggeredMatrix(snippets, snippetTime, showXLabel)
plot(snippetTime, snippets', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.75); hold on;
plot(snippetTime, mean(snippets, 1, 'omitnan'), 'k', 'LineWidth', 1.5);
xline(0, '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 1);
box off;
set(gca, 'TickDir', 'out');
ylabel('\DeltaF/F');
if showXLabel
    xlabel('Time from event (s)');
else
    set(gca, 'XTickLabel', []);
end
end


function plotPopulationEventTriggered(traces, events, cellsToPlot, opts)
switch opts.aggregateCells
    case 'figure'
        aggregateCells = cellsToPlot;
    case 'allAccepted'
        aggregateCells = 1:size(traces, 1);
end

allSnippets = [];
for cellNo = aggregateCells(:)'
    eventFrames = find(events(cellNo, :) > 0);
    if isempty(eventFrames)
        continue;
    end
    [snippets, snippetTime] = extractEventSnippets(traces(cellNo, :), eventFrames, opts);
    snippets = zscoreSnippetsToBaseline(snippets, snippetTime);
    allSnippets = [allSnippets; snippets]; %#ok<AGROW>
end

if isempty(allSnippets)
    text(0.5, 0.5, 'No aggregate events', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

semTrace = std(allSnippets, 0, 1, 'omitnan') ./ sqrt(size(allSnippets, 1));
meanTrace = mean(allSnippets, 1, 'omitnan');

if opts.aggregateShowSubset > 0
    rng(opts.aggregateSubsetSeed);
    nSubset = min(opts.aggregateShowSubset, size(allSnippets, 1));
    subsetIdx = randperm(size(allSnippets, 1), nSubset);
    plot(snippetTime, allSnippets(subsetIdx, :)', 'Color', [0.9 0.9 0.9], ...
        'LineWidth', 0.25); hold on;
else
    hold on;
end
fill([snippetTime fliplr(snippetTime)], ...
    [meanTrace - semTrace, fliplr(meanTrace + semTrace)], ...
    opts.aggregateSemColor, 'FaceAlpha', opts.aggregateSemAlpha, 'EdgeColor', 'none');
plot(snippetTime, meanTrace, 'k', 'LineWidth', 1.8);
xline(0, '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.2);
box off;
set(gca, 'TickDir', 'out');
ylabel('\DeltaF/F z-scored to pre-event baseline');
xlabel('Time from event (s)');
title(sprintf('Population event-triggered \\DeltaF/F (n = %d events)', size(allSnippets, 1)));
end


function zSnippets = zscoreSnippetsToBaseline(snippets, snippetTime)
baselineIdx = snippetTime >= -2 & snippetTime < 0;
zSnippets = NaN(size(snippets));

for i = 1:size(snippets, 1)
    baseline = snippets(i, baselineIdx);
    mu = mean(baseline, 'omitnan');
    sigma = std(baseline, 0, 'omitnan');
    if ~isfinite(sigma) || sigma == 0
        sigma = 1;
    end
    zSnippets(i, :) = (snippets(i, :) - mu) ./ sigma;
end
end


function saveFigureOutputs(figHandle, saveFigurePath, vectorFormat)
[outDir, baseName, ext] = fileparts(saveFigurePath);
if isempty(outDir)
    outDir = pwd;
end
if isempty(ext)
    ext = '.png';
end

pngPath = fullfile(outDir, [baseName '.png']);
vectorPath = fullfile(outDir, [baseName '.' char(vectorFormat)]);

try
    exportgraphics(figHandle, pngPath, 'Resolution', 300);
catch
    saveas(figHandle, pngPath);
end

try
    set(figHandle, 'Renderer', 'painters');
    switch char(vectorFormat)
        case 'pdf'
            print(figHandle, vectorPath, '-dpdf', '-painters', '-bestfit');
        case 'svg'
            print(figHandle, vectorPath, '-dsvg', '-painters');
    end
catch
    try
        exportgraphics(figHandle, vectorPath, 'ContentType', 'vector');
    catch
        saveas(figHandle, vectorPath);
    end
end

if ~strcmpi(ext, '.png') && ~strcmpi(saveFigurePath, vectorPath)
    try
        saveas(figHandle, saveFigurePath);
    catch
    end
end

fprintf('Saved CNMF-E event example figure to %s and %s\n', pngPath, vectorPath);
end


function [analysisFile, sortedFile, timestampFile] = findSessionFiles(sessionFolder, opts)
if ~isfolder(sessionFolder)
    error('plotCnmfeEventExamples:MissingSessionFolder', ...
        'Session folder was not found: %s', sessionFolder);
end

if isempty(opts.analysisFile)
    analysisFiles = dir(fullfile(sessionFolder, '**', '*_cnmfeAnalysis.mat'));
    analysisFiles = analysisFiles(~contains({analysisFiles.name}, 'Sorted'));
    analysisFile = chooseLargestFile(sessionFolder, analysisFiles, 'CNMF-E analysis');
else
    analysisFile = opts.analysisFile;
end

if isempty(opts.sortedFile)
    sortedFiles = dir(fullfile(sessionFolder, '**', '*_cnmfeAnalysisSorted*.mat'));
    sortedFile = chooseLargestFile(sessionFolder, sortedFiles, 'CNMF-E sorted');
else
    sortedFile = opts.sortedFile;
end

if isempty(opts.timestampFile)
    timestampFiles = dir(fullfile(sessionFolder, '**', 'timeStamps.csv'));
    timestampFile = chooseLargestFile(sessionFolder, timestampFiles, 'timestamp');
else
    timestampFile = opts.timestampFile;
end
end


function filePath = chooseLargestFile(folderPath, files, label)
if isempty(files)
    error('plotCnmfeEventExamples:MissingFile', ...
        'No %s file was found in %s.', label, folderPath);
end

[~, idx] = max([files.bytes]);
if isfield(files, 'folder') && ~isempty(files(idx).folder)
    filePath = fullfile(files(idx).folder, files(idx).name);
else
    filePath = fullfile(folderPath, files(idx).name);
end
end


function timestamps = readSessionTimestamps(timestampFile)
timestamps = readtable(timestampFile);
timestamps = table2array(timestamps);
if size(timestamps, 2) >= 2
    timestamps = timestamps(:, 2);
end
timestamps = timestamps(:);

if numel(timestamps) >= 5 && timestamps(5) > 2
    timestamps = timestamps ./ 1000;
end
end


function timeSeconds = frameToSeconds(frames, timestamps, frameRate)
frames = frames(:)';
timestampFrames = frames .* 2;
if ~isempty(timestamps) && max(timestampFrames) <= numel(timestamps)
    timeSeconds = timestamps(timestampFrames);
else
    timeSeconds = (frames - 1) ./ frameRate;
end
end
