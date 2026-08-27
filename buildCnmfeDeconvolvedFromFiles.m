function [Ca_deconvolved, rowIdxByDay, fileInfo] = buildCnmfeDeconvolvedFromFiles(RAT, varargin)
%BUILDCNMFECONVOLVEDFROMFILES Rebuild accepted-cell CNMF-E event amplitudes.
%
%   [Ca_deconvolved, rowIdxByDay, fileInfo] =
%       buildCnmfeDeconvolvedFromFiles(RAT, 'SearchRoot', rootFolder)
%
%   For each RAT.Ca_traces.CA_traces_DATE field, this searches the original
%   session folders for the matching DATE, loads *_cnmfeAnalysis.mat and
%   *_cnmfeAnalysisSorted*.mat, and stores accepted-cell rows of
%   cnmfeAnalysisOutput.extractedPeaks.
%
%   Output fields are named CA_peaks_DATE so they can be paired directly with
%   RAT.Ca_peaks.CA_peaks_DATE.

p = inputParser;
addParameter(p, 'SearchRoot', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'TraceFields', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
opts = p.Results;
opts.Verbose = logical(opts.Verbose);

if isempty(opts.SearchRoot)
    error('buildCnmfeDeconvolvedFromFiles:MissingSearchRoot', ...
        'Pass the parent video folder as SearchRoot.');
end
searchRoot = char(opts.SearchRoot);
if ~isfolder(searchRoot)
    error('buildCnmfeDeconvolvedFromFiles:MissingSearchRoot', ...
        'SearchRoot was not found: %s', searchRoot);
end

if ~isfield(RAT, 'Ca_traces') || ~isstruct(RAT.Ca_traces)
    error('buildCnmfeDeconvolvedFromFiles:MissingTraces', ...
        'RAT must contain RAT.Ca_traces.');
end

if isempty(opts.TraceFields)
    traceFields = fieldnames(RAT.Ca_traces);
    traceFields = traceFields(startsWith(traceFields, 'CA_traces_'));
else
    traceFields = cellstr(opts.TraceFields);
end

Ca_deconvolved = struct();
rowIdxByDay = struct();
fileInfo = struct();

for i = 1:numel(traceFields)
    traceField = traceFields{i};
    dateLabel = regexprep(traceField, '^CA_traces_', '');
    peakField = ['CA_peaks_' dateLabel];

    expectedSize = size(RAT.Ca_traces.(traceField));
    if opts.Verbose
        fprintf('Searching for %s under %s...\n', dateLabel, searchRoot);
    end
    [analysisFile, sortedFile] = findMatchingCnmfeFiles( ...
        searchRoot, dateLabel, expectedSize, opts.Verbose);

    analysisData = load(analysisFile);
    sortedData = load(sortedFile);

    if ~isfield(analysisData, 'cnmfeAnalysisOutput') || ...
            ~isfield(analysisData.cnmfeAnalysisOutput, 'extractedPeaks')
        error('buildCnmfeDeconvolvedFromFiles:MissingExtractedPeaks', ...
            '%s does not contain cnmfeAnalysisOutput.extractedPeaks.', analysisFile);
    end
    if ~isfield(sortedData, 'validCNMFE')
        error('buildCnmfeDeconvolvedFromFiles:MissingValidCNMFE', ...
            '%s does not contain validCNMFE.', sortedFile);
    end

    eventAmplitudes = double(analysisData.cnmfeAnalysisOutput.extractedPeaks);
    goodRows = find(sortedData.validCNMFE == 1);
    goodRows = goodRows(goodRows <= size(eventAmplitudes, 1));
    eventAmplitudes = eventAmplitudes(goodRows, :);

    if ~isequal(size(eventAmplitudes), expectedSize)
        error('buildCnmfeDeconvolvedFromFiles:SizeMismatch', ...
            ['%s rebuilt to %d x %d, but %s is %d x %d. ' ...
            'Check that the matched CNMF-E file is the correct session.'], ...
            peakField, size(eventAmplitudes, 1), size(eventAmplitudes, 2), ...
            traceField, expectedSize(1), expectedSize(2));
    end

    Ca_deconvolved.(peakField) = eventAmplitudes;
    rowIdxByDay.(peakField) = goodRows(:)';
    rowIdxByDay.(['date_' dateLabel]) = goodRows(:)';
    fileInfo.(peakField) = struct( ...
        'analysisFile', analysisFile, ...
        'sortedFile', sortedFile, ...
        'nAccepted', numel(goodRows), ...
        'nFrames', size(eventAmplitudes, 2));

    if opts.Verbose
        fprintf('Loaded %s: %d accepted cells x %d frames\n', ...
            peakField, size(eventAmplitudes, 1), size(eventAmplitudes, 2));
    end
end
end


function [analysisFile, sortedFile] = findMatchingCnmfeFiles(searchRoot, dateLabel, expectedSize, verbose)
dateDirs = findCandidateDateDirs(searchRoot, dateLabel, verbose);

if isempty(dateDirs)
    error('buildCnmfeDeconvolvedFromFiles:DateNotFound', ...
        'Could not find a session date folder named %s under %s.', dateLabel, searchRoot);
end

candidateErrors = {};
for i = 1:numel(dateDirs)
    datePath = fullfile(dateDirs(i).folder, dateDirs(i).name);
    if verbose
        fprintf('  [%d/%d] Checking date folder: %s\n', i, numel(dateDirs), datePath);
    end
    analysisFiles = dir(fullfile(datePath, '**', '*_cnmfeAnalysis.mat'));
    analysisFiles = analysisFiles(~contains({analysisFiles.name}, 'Sorted'));
    sortedFiles = dir(fullfile(datePath, '**', '*_cnmfeAnalysisSorted*.mat'));

    if verbose
        fprintf('      Found %d analysis files and %d sorted files\n', ...
            numel(analysisFiles), numel(sortedFiles));
    end

    if isempty(analysisFiles) || isempty(sortedFiles)
        continue;
    end

    [~, analysisOrder] = sort([analysisFiles.bytes], 'descend');
    [~, sortedOrder] = sort([sortedFiles.bytes], 'descend');

    for a = analysisOrder(:)'
        for s = sortedOrder(:)'
            thisAnalysis = fullfile(analysisFiles(a).folder, analysisFiles(a).name);
            thisSorted = fullfile(sortedFiles(s).folder, sortedFiles(s).name);
            if verbose
                fprintf('      Testing %s with %s...\n', ...
                    analysisFiles(a).name, sortedFiles(s).name);
            end
            try
                if cnmfeFilesMatchExpectedSize(thisAnalysis, thisSorted, expectedSize)
                    analysisFile = thisAnalysis;
                    sortedFile = thisSorted;
                    if verbose
                        fprintf('      Matched CNMF-E files for %s\n', dateLabel);
                    end
                    return;
                end
            catch ME
                candidateErrors{end + 1} = sprintf('%s: %s', thisAnalysis, ME.message); %#ok<AGROW>
            end
        end
    end
end

if isempty(candidateErrors)
    error('buildCnmfeDeconvolvedFromFiles:NoMatchingFiles', ...
        'Found %s folders, but no CNMF-E analysis/sorted file pair matched %d x %d.', ...
        dateLabel, expectedSize(1), expectedSize(2));
else
    error('buildCnmfeDeconvolvedFromFiles:NoMatchingFiles', ...
        'No CNMF-E file pair matched %s (%d x %d). First error: %s', ...
        dateLabel, expectedSize(1), expectedSize(2), candidateErrors{1});
end
end


function dateDirs = findCandidateDateDirs(searchRoot, dateLabel, verbose)
dateDirs = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, 'isdir', {}, 'datenum', {});

% Fast paths first. These cover the usual layouts:
%   root/DATE
%   root/RAT_FOLDER/DATE
%   root/RAT_FOLDER/SESSION/DATE
patterns = { ...
    fullfile(searchRoot, dateLabel), ...
    fullfile(searchRoot, '*', dateLabel), ...
    fullfile(searchRoot, '*', '*', dateLabel)};

for p = 1:numel(patterns)
    if verbose
        fprintf('  Fast search pattern %d/%d: %s\n', p, numel(patterns), patterns{p});
    end
    hits = dir(patterns{p});
    hits = hits([hits.isdir]);
    dateDirs = appendUniqueDirs(dateDirs, hits);
    if ~isempty(dateDirs)
        return;
    end
end

% Fallback: manual recursive search with progress output. This is slower,
% but unlike dir(root/**/DATE), it tells you where MATLAB is spending time.
if verbose
    fprintf('  Fast search found no date folder. Starting recursive search...\n');
end
dateDirs = recursiveFindDateDirs(searchRoot, dateLabel, verbose, 0);
end


function out = appendUniqueDirs(out, hits)
for i = 1:numel(hits)
    if ~hits(i).isdir
        continue;
    end
    fullPath = fullfile(hits(i).folder, hits(i).name);
    alreadyHave = false;
    for j = 1:numel(out)
        if strcmp(fullfile(out(j).folder, out(j).name), fullPath)
            alreadyHave = true;
            break;
        end
    end
    if ~alreadyHave
        out(end + 1) = hits(i); %#ok<AGROW>
    end
end
end


function matches = recursiveFindDateDirs(folderPath, dateLabel, verbose, depth)
matches = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, 'isdir', {}, 'datenum', {});
if depth > 5
    return;
end

if verbose
    fprintf('    scanning depth %d: %s\n', depth, folderPath);
end

items = dir(folderPath);
items = items([items.isdir]);
items = items(~ismember({items.name}, {'.', '..'}));

for i = 1:numel(items)
    if strcmp(items(i).name, dateLabel)
        matches = appendUniqueDirs(matches, items(i));
        continue;
    end
    childPath = fullfile(items(i).folder, items(i).name);
    childMatches = recursiveFindDateDirs(childPath, dateLabel, verbose, depth + 1);
    matches = appendUniqueDirs(matches, childMatches);
    if ~isempty(matches)
        return;
    end
end
end


function tf = cnmfeFilesMatchExpectedSize(analysisFile, sortedFile, expectedSize)
analysisData = load(analysisFile, 'cnmfeAnalysisOutput');
sortedData = load(sortedFile, 'validCNMFE');

if ~isfield(analysisData, 'cnmfeAnalysisOutput') || ...
        ~isfield(analysisData.cnmfeAnalysisOutput, 'extractedPeaks') || ...
        ~isfield(sortedData, 'validCNMFE')
    tf = false;
    return;
end

eventSize = size(analysisData.cnmfeAnalysisOutput.extractedPeaks);
goodRows = find(sortedData.validCNMFE == 1);
goodRows = goodRows(goodRows <= eventSize(1));
tf = numel(goodRows) == expectedSize(1) && eventSize(2) == expectedSize(2);
end
