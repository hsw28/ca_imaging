function [QC, BinSummary, BinOccupancy] = speedmatch_diagnostics(R, varargin)
% speedmatch_diagnostics  Summarize QC for speed/acceleration matched output.
%
% Usage:
%   R = run_speedBinMatched({'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%       'SpeedBinWidth',2,'MinDurPerBin',0.1,'MinBins',1,'win',[0 2],'Plot',false);
%   [QC, BinSummary, BinOccupancy] = speedmatch_diagnostics(R);
%
%   Racc = run_accelBinMatched({'rat0222','rat0307','rat0313','rat0314','rat0816'});
%   [QCacc, AccelSummary, AccelOccupancy] = speedmatch_diagnostics(Racc, ...
%       'AnalysisType','acceleration');
%
% Optional name/value:
%   'WriteCSV'      false
%   'QCFile'        'speed_bin_matching_QC.csv'
%   'BinSummaryFile' 'matched_bin_summary.csv'
%   'BinOccupancyFile' 'matched_bin_occupancy.csv'
%   'AnalysisType'  'auto'  % 'auto' | 'speed' | 'acceleration'
%   'Print'         true

ip = inputParser;
ip.addParameter('WriteCSV', false, @islogical);
ip.addParameter('QCFile', 'speed_bin_matching_QC.csv', @(x) ischar(x) || isstring(x));
ip.addParameter('BinSummaryFile', 'matched_bin_summary.csv', @(x) ischar(x) || isstring(x));
ip.addParameter('BinOccupancyFile', 'matched_bin_occupancy.csv', @(x) ischar(x) || isstring(x));
ip.addParameter('SpeedBinsFile', '', @(x) ischar(x) || isstring(x)); % legacy alias
ip.addParameter('AnalysisType', 'auto', @(s) any(validatestring(lower(char(s)), {'auto','speed','acceleration'})));
ip.addParameter('Print', true, @islogical);
ip.parse(varargin{:});
P = ip.Results;
analysisType = lower(char(P.AnalysisType));
if strlength(string(P.SpeedBinsFile)) > 0
    P.BinSummaryFile = P.SpeedBinsFile;
end

if nargin < 1 || isempty(R)
    R = evalin('base', 'R');
end

qcRows = {};
summaryRows = {};
occupancyRows = {};

for r = 1:numel(R)
    if ~isfield(R, 'rat') || isempty(R(r).rat)
        ratName = sprintf('rat_%d', r);
    else
        ratName = char(R(r).rat);
    end

    if ~isfield(R(r), 'perDay') || isempty(R(r).perDay)
        continue;
    end

    for d = 1:numel(R(r).perDay)
        S = R(r).perDay{d};
        if isempty(S) || ~isstruct(S)
            continue;
        end

        thisType = infer_analysis_type(S, analysisType);
        binLabel = bin_label_for_type(thisType);

        D = getDropStruct(S);
        nCellsTotal = getfield_default(D, 'nCells', numel(getfield_default(S, 'pVal', [])));
        nMasked = getfield_default(D, 'nMasked', NaN);
        nEmpty = getfield_default(D, 'nEmpty', NaN);
        nBelowMinBins = getfield_default(D, 'nBelowMinBins', NaN);
        nTested = getfield_default(D, 'nTested', sum(isfinite(getfield_default(S, 'pVal', []))));
        matchedSpeedBinsAvailable = getfield_default(D, 'pairedBinsAvail', NaN);
        minBinsRequired = getfield_default(D, 'MinBins', NaN);
        minDurPerBinSec = getfield_default(D, 'MinDurPerBin', NaN);

        nExcluded = nCellsTotal - nMasked - nTested;
        if ~isfinite(nExcluded)
            nExcluded = NaN;
        end

        nBinsUsedAll = getfield_default(S, 'nBinsUsed', []);
        nBinsUsed = nBinsUsedAll(nBinsUsedAll > 0 & isfinite(nBinsUsedAll));

        medianBinsUsed = nanScalar(@median, nBinsUsed);
        meanBinsUsed = nanScalar(@mean, nBinsUsed);
        minBinsUsed = nanScalar(@min, nBinsUsed);
        maxBinsUsed = nanScalar(@max, nBinsUsed);
        p25BinsUsed = nanPercentile(nBinsUsed, 25);
        p75BinsUsed = nanPercentile(nBinsUsed, 75);

        pVal = getfield_default(S, 'pVal', []);
        sigDiffFDR = getfield_default(S, 'sigDiffFDR', false(size(pVal)));
        tested = isfinite(pVal);
        nSigFDR = sum(sigDiffFDR(tested));
        pctSigFDR = getfield_default(S, 'pctSigDiffFDR', 100 * mean(sigDiffFDR(tested), 'omitnan'));

        DL = getfield_default(S, 'dayLevel', struct());
        dayLevelNCells = getfield_default(DL, 'nCells', NaN);
        dayLevelDeltaHz = getfield_default(DL, 'delta', NaN);
        dayLevelP = getfield_default(DL, 'p', NaN);

        binCenters = collect_bin_centers(S, thisType);
        binMin = nanScalar(@min, binCenters);
        binMax = nanScalar(@max, binCenters);
        binMedian = nanScalar(@median, binCenters);
        binP25 = nanPercentile(binCenters, 25);
        binP75 = nanPercentile(binCenters, 75);
        nCellBinPairs = numel(binCenters);

        M = get_matching_struct(S);
        [matchedCenterMin, matchedCenterMax, matchedCenterMedian, matchedCenterP25, matchedCenterP75] = ...
            summarize_matching_centers(M, binCenters);
        [trialDurTotal, nonTrialDurTotal, durRatioMedian, durRatioP25, durRatioP75, ...
            medianAbsLog2DurRatio, overlap, l1Distance, nOccupancyBins] = summarize_matching_quality(M);
        [nTotalBins, nTrialOccupiedBins, nNonTrialOccupiedBins, nEitherOccupiedBins, ...
            nBothOccupiedBins, nRetainedBins] = summarize_bin_denominators(M, matchedSpeedBinsAvailable);

        qcRows(end+1,:) = { ...
            ratName, d, thisType, binLabel, nCellsTotal, nMasked, nEmpty, nBelowMinBins, ...
            nTested, nExcluded, nSigFDR, pctSigFDR, ...
            matchedSpeedBinsAvailable, minBinsRequired, minDurPerBinSec, ...
            nTotalBins, nTrialOccupiedBins, nNonTrialOccupiedBins, nEitherOccupiedBins, ...
            nBothOccupiedBins, nRetainedBins, ...
            medianBinsUsed, meanBinsUsed, minBinsUsed, maxBinsUsed, ...
            p25BinsUsed, p75BinsUsed, dayLevelNCells, dayLevelDeltaHz, dayLevelP, ...
            binMin, binMax, binMedian, binP25, binP75, nCellBinPairs, ...
            matchedCenterMin, matchedCenterMax, matchedCenterMedian, matchedCenterP25, matchedCenterP75, ...
            trialDurTotal, nonTrialDurTotal, durRatioMedian, durRatioP25, durRatioP75, ...
            medianAbsLog2DurRatio, overlap, l1Distance, nOccupancyBins}; %#ok<AGROW>

        if ~isempty(binCenters)
            summaryRows(end+1,:) = { ...
                ratName, d, thisType, binLabel, binMin, binMax, binMedian, binP25, binP75, ...
                nCellBinPairs, nTotalBins, nTrialOccupiedBins, nNonTrialOccupiedBins, ...
                nEitherOccupiedBins, nBothOccupiedBins, nRetainedBins, ...
                matchedCenterMin, matchedCenterMax, matchedCenterMedian, ...
                matchedCenterP25, matchedCenterP75}; %#ok<AGROW>
        end

        if ~isempty(M) && isfield(M, 'binCenters')
            centers = M.binCenters(:);
            trialDur = getfield_default(M, 'trialDurSec', nan(size(centers)));
            nonDur = getfield_default(M, 'nonTrialDurSec', nan(size(centers)));
            ratio = getfield_default(M, 'nonTrialToTrialDurRatio', nonDur ./ max(trialDur, eps));
            trialProb = getfield_default(M, 'trialProb', nan(size(centers)));
            nonProb = getfield_default(M, 'nonTrialProb', nan(size(centers)));
            for b = 1:numel(centers)
                occupancyRows(end+1,:) = { ...
                    ratName, d, thisType, binLabel, b, centers(b), trialDur(b), ...
                    nonDur(b), ratio(b), log2(ratio(b)), trialProb(b), nonProb(b)}; %#ok<AGROW>
            end
        end
    end
end

qcNames = { ...
    'rat', 'sessionIndex', 'analysisType', 'binLabel', ...
    'nCellsTotal', 'nMasked', 'nEmpty', ...
    'nBelowMinBins', 'nTested', 'nExcluded', 'nSigFDR', 'pctSigFDR', ...
    'matchedBinsAvailable', 'minBinsRequired', 'minDurPerBinSec', ...
    'nTotalBins', 'nTrialOccupiedBins', 'nNonTrialOccupiedBins', ...
    'nEitherOccupiedBins', 'nBothOccupiedBins', 'nRetainedBins', ...
    'medianBinsUsedPerCell', 'meanBinsUsedPerCell', ...
    'minBinsUsedPerCell', 'maxBinsUsedPerCell', ...
    'p25BinsUsedPerCell', 'p75BinsUsedPerCell', ...
    'dayLevelNCells', 'dayLevelDeltaHz', 'dayLevelP', ...
    'cellUsedBinMin', 'cellUsedBinMax', 'cellUsedBinMedian', ...
    'cellUsedBinP25', 'cellUsedBinP75', 'nCellBinPairs', ...
    'matchedBinMin', 'matchedBinMax', 'matchedBinMedian', ...
    'matchedBinP25', 'matchedBinP75', ...
    'trialDurTotalSec', 'nonTrialDurTotalSec', ...
    'medianNonTrialToTrialDurRatio', 'p25NonTrialToTrialDurRatio', ...
    'p75NonTrialToTrialDurRatio', 'medianAbsLog2DurRatio', ...
    'occupancyOverlap', 'occupancyL1Distance', 'nOccupancyBins'};

summaryNames = { ...
    'rat', 'sessionIndex', 'analysisType', 'binLabel', ...
    'cellUsedBinMin', 'cellUsedBinMax', 'cellUsedBinMedian', ...
    'cellUsedBinP25', 'cellUsedBinP75', 'nCellBinPairs', ...
    'nTotalBins', 'nTrialOccupiedBins', 'nNonTrialOccupiedBins', ...
    'nEitherOccupiedBins', 'nBothOccupiedBins', 'nRetainedBins', ...
    'matchedBinMin', 'matchedBinMax', 'matchedBinMedian', ...
    'matchedBinP25', 'matchedBinP75'};

occupancyNames = { ...
    'rat', 'sessionIndex', 'analysisType', 'binLabel', 'matchedBinIndex', ...
    'binCenter', 'trialDurSec', 'nonTrialDurSec', ...
    'nonTrialToTrialDurRatio', 'log2NonTrialToTrialDurRatio', ...
    'trialProb', 'nonTrialProb'};

if isempty(qcRows)
    QC = cell2table(cell(0, numel(qcNames)), 'VariableNames', qcNames);
else
    QC = cell2table(qcRows, 'VariableNames', qcNames);
end

if isempty(summaryRows)
    BinSummary = cell2table(cell(0, numel(summaryNames)), 'VariableNames', summaryNames);
else
    BinSummary = cell2table(summaryRows, 'VariableNames', summaryNames);
end

if isempty(occupancyRows)
    BinOccupancy = cell2table(cell(0, numel(occupancyNames)), 'VariableNames', occupancyNames);
else
    BinOccupancy = cell2table(occupancyRows, 'VariableNames', occupancyNames);
end

if P.Print
    print_summary(QC);
end

if P.WriteCSV
    writetable(QC, P.QCFile);
    writetable(BinSummary, P.BinSummaryFile);
    writetable(BinOccupancy, P.BinOccupancyFile);
    fprintf('Wrote %s, %s, and %s\n', P.QCFile, P.BinSummaryFile, P.BinOccupancyFile);
end
end

function analysisType = infer_analysis_type(S, requested)
if ~strcmp(requested, 'auto')
    analysisType = requested;
    return;
end
if isfield(S, 'accelBinCenters') || isfield(S, 'AccelEdges')
    analysisType = 'acceleration';
else
    analysisType = 'speed';
end
end

function label = bin_label_for_type(analysisType)
switch analysisType
    case 'acceleration'
        label = 'acceleration';
    otherwise
        label = 'speed';
end
end

function D = getDropStruct(S)
if isfield(S, 'drop') && isstruct(S.drop)
    D = S.drop;
else
    D = struct();
end
end

function val = getfield_default(S, fieldName, defaultVal)
if isstruct(S) && isfield(S, fieldName)
    val = S.(fieldName);
else
    val = defaultVal;
end
end

function x = collect_bin_centers(S, analysisType)
x = [];
fieldName = 'speedBinCenters';
if strcmp(analysisType, 'acceleration')
    fieldName = 'accelBinCenters';
end

if ~isfield(S, fieldName) || isempty(S.(fieldName))
    return;
end

binCenters = S.(fieldName);
if iscell(binCenters)
    for i = 1:numel(binCenters)
        xi = binCenters{i};
        x = [x; xi(:)]; %#ok<AGROW>
    end
else
    x = binCenters(:);
end
x = x(isfinite(x));
end

function M = get_matching_struct(S)
if isfield(S, 'matching') && isstruct(S.matching)
    M = S.matching;
else
    M = [];
end
end

function [mn, mx, med, p25, p75] = summarize_matching_centers(M, fallbackCenters)
if ~isempty(M) && isfield(M, 'binCenters')
    centers = M.binCenters(:);
else
    centers = fallbackCenters(:);
end
centers = centers(isfinite(centers));
mn = nanScalar(@min, centers);
mx = nanScalar(@max, centers);
med = nanScalar(@median, centers);
p25 = nanPercentile(centers, 25);
p75 = nanPercentile(centers, 75);
end

function [trialTotal, nonTotal, ratioMed, ratioP25, ratioP75, medAbsLog2, overlap, l1, nBins] = summarize_matching_quality(M)
if isempty(M)
    [trialTotal, nonTotal, ratioMed, ratioP25, ratioP75, medAbsLog2, overlap, l1, nBins] = deal(NaN);
    return;
end

trialDur = getfield_default(M, 'trialDurSec', []);
nonDur = getfield_default(M, 'nonTrialDurSec', []);
ratio = getfield_default(M, 'nonTrialToTrialDurRatio', nonDur ./ max(trialDur, eps));
log2Ratio = getfield_default(M, 'log2NonTrialToTrialDurRatio', log2(ratio));

trialTotal = sum(trialDur, 'omitnan');
nonTotal = sum(nonDur, 'omitnan');
ratioMed = nanScalar(@median, ratio);
ratioP25 = nanPercentile(ratio, 25);
ratioP75 = nanPercentile(ratio, 75);
medAbsLog2 = nanScalar(@median, abs(log2Ratio(isfinite(log2Ratio))));
overlap = getfield_default(M, 'overlap', NaN);
l1 = getfield_default(M, 'l1Distance', NaN);
nBins = numel(getfield_default(M, 'binCenters', []));
end

function [nTotal, nTrial, nNon, nEither, nBoth, nRetained] = summarize_bin_denominators(M, fallbackRetained)
if isempty(M)
    nTotal = NaN;
    nTrial = NaN;
    nNon = NaN;
    nEither = NaN;
    nBoth = NaN;
    nRetained = fallbackRetained;
    return;
end

nTotal = getfield_default(M, 'nTotalBins', NaN);
nTrial = getfield_default(M, 'nTrialOccupiedBins', NaN);
nNon = getfield_default(M, 'nNonTrialOccupiedBins', NaN);
nEither = getfield_default(M, 'nEitherOccupiedBins', NaN);
nBoth = getfield_default(M, 'nBothOccupiedBins', NaN);
nRetained = getfield_default(M, 'nRetainedBins', fallbackRetained);
end

function y = nanScalar(funHandle, x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = funHandle(x);
end
end

function y = nanPercentile(x, p)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = prctile(x, p);
end
end

function print_summary(QC)
if isempty(QC)
    fprintf('No matched-bin diagnostics found.\n');
    return;
end

fprintf('\n=== Matched-bin diagnostics by rat/session ===\n');
for i = 1:height(QC)
    fprintf(['%s session %d: tested=%d/%d, excluded=%d, belowMinBins=%d, ' ...
        'matchedBins=%g/%g total (%g both-occupied), bins/cell median=%.1f [%.1f %.1f], ' ...
        'matched %s=%.2f-%.2f, occupancy overlap=%.3f, L1=%.3f, ' ...
        'dur ratio median=%.2f, sig=%d (%.1f%%)\n'], ...
        QC.rat{i}, QC.sessionIndex(i), QC.nTested(i), QC.nCellsTotal(i), ...
        QC.nExcluded(i), QC.nBelowMinBins(i), QC.matchedBinsAvailable(i), ...
        QC.nTotalBins(i), QC.nBothOccupiedBins(i), ...
        QC.medianBinsUsedPerCell(i), QC.minBinsUsedPerCell(i), ...
        QC.maxBinsUsedPerCell(i), QC.binLabel{i}, QC.matchedBinMin(i), QC.matchedBinMax(i), ...
        QC.occupancyOverlap(i), QC.occupancyL1Distance(i), ...
        QC.medianNonTrialToTrialDurRatio(i), ...
        QC.nSigFDR(i), QC.pctSigFDR(i));
end

totalTested = sum(QC.nTested, 'omitnan');
totalCells = sum(QC.nCellsTotal, 'omitnan');
totalExcluded = sum(QC.nExcluded, 'omitnan');
totalSig = sum(QC.nSigFDR, 'omitnan');

fprintf('\nPooled: tested=%d/%d cells, excluded=%d, sig=%d (%.1f%% of tested)\n', ...
    totalTested, totalCells, totalExcluded, totalSig, ...
    100 * totalSig / max(totalTested, 1));
end
