function R = taskPlaceExamples(ratNames, varargin)
% taskPlaceExamples  Reviewer figure: non-task place map vs trace activity.
%
% Each row is one neuron; columns are the non-tEBC rate map, task/trace
% spatial map with the established non-tEBC field boundary, and PETH.
% The implementation follows the conventions used by CA_normalizePosData,
% plotCellComposite, plotRateMaskSummary, and tracecellactivity:
%   * spatial bins = 2.5 cm; map display smoothing sigma = 0.75 bins
%   * non-task samples are running >= 4 cm/s and outside CS+[0,2] s
%   * comparison maps use CS+[0,2] s by default, or CS+[0.25,0.75] s
%     when TraceOnly=true
%   * the non-task mask is rate > mean(rate) + N*SD(rate), default N = 1
%   * PETH is [-1,2] s, 50-ms bins, shown as absolute event rate (Hz)
% Spatial maps are peak-normalized within each condition to visualize the
% activity location; absolute event rates are quantified separately in Fig. Sx.
%
% Automatic selection searches the three days ending at rat.An and ranks
% spatially significant cells into outside-field, inside-field, and mixed categories.
% Manual selection bypasses automatic selection:
%   S = table(["rat0314";"rat0816"], ["2023_05_25";"2022_11_04"], ...
%             [12;2], 'VariableNames', {'Animal','Session','CellID'});
%   R = taskPlaceExamples([], 'SelectedNeurons', S);
% Fig. S20 starting points can be reproduced/inspected with:
%   S = {'rat0314','2023_05_22',296; 'rat0314','2023_05_22',305; ...
%        'rat0816','2022_11_03',14};
%   R = taskPlaceExamples([], 'SelectedNeurons', S);
%
% Useful outputs:
%   R.summary     selected-cell summary table
%   R.candidates ranked candidate table
%   R.figure      figure handle

% Name/value options:
%   SelectedNeurons  table/struct/cell array with animal, session, cell ID
%   ExcludeNeurons   table/struct/N-by-3 cell array of cells to skip during
%                    automatic selection
%   NumNeurons       total examples, distributed across three spatial
%                    relationships (positive integer; default 3)
%   Days             optional cellstr of dates to use for every animal
%   MapBin           spatial bin size in cm (default 2.5)
%   SmoothSigma      display smoothing in bins (default 0.75)
%   MapNormalizationPercentile robust positive-rate ceiling used for map
%                    display normalization (default 100 = actual maximum)
%   VelThresh        non-task running threshold in cm/s (default 4)
%   MaskNStd         Fig. S20 threshold in SD above map mean (default 1)
%   TraceOnly        false: map the full 0--2 s task period (default);
%                    true: map only the 0.25--0.75 s trace interval
%   PETHBin          event-rate bin width in seconds (default 0.05)
%   MinFieldBins     minimum contiguous bins in a clear place field (default 5)
%   MaxPlaceFields   maximum clear fields for automatic selection (default 2)
%   MinCSUSFraction  minimum fraction of 0--2 s task events occurring
%                    before US onset at 0.75 s (default 0.5)
%   OutsideMaxOverlap maximum primary-field event fraction for the
%                    outside-field archetype (default 0.15)
%   OverlapMinOverlap minimum primary-field event fraction for the
%                    inside-field archetype (default 0.50)
%   OutsideMaxFieldPeak maximum within-field map intensity relative to the
%                    comparison-map peak for an outside example (default 0.20)
%   OverlapMinFieldPeak minimum within-field relative map intensity for an
%                    inside-field example (default 0.70)
%   InsideMinMapMass minimum comparison-map activity mass inside the field
%                    for the inside-field example (default 0.40)
%   MixedMinOverlap  lower inside fraction for mixed examples (default 0.25)
%   MixedMaxOverlap  upper inside fraction for mixed examples (default 0.75)
%   SavePath         optional .pdf, .svg, .png, or .fig output path
%   SummaryPath      optional .csv summary output path

% R = taskPlaceExamples({'rat0222','rat0307','rat0313', ...
%                                'rat0314','rat0816'});

% See also CA_normalizePosData, plotRateMaskSummary, tracecellactivity.


if nargin < 1 || isempty(ratNames)
    ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end

p = inputParser;
p.addParameter('SelectedNeurons', [], @(x) isempty(x) || istable(x) || isstruct(x) || iscell(x));
p.addParameter('ExcludeNeurons', [], @(x) isempty(x) || istable(x) || isstruct(x) || iscell(x));
p.addParameter('NumNeurons', 3, @(x) isnumeric(x) && isscalar(x) && ...
    x>=1 && x==round(x));
p.addParameter('Days', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('MapBin', 2.5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('SmoothSigma', 0.75, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MapNormalizationPercentile', 100, ...
    @(x) isnumeric(x) && isscalar(x) && x>=90 && x<=100);
p.addParameter('VelThresh', 4, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MaskNStd', 1, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('TaskExclusionWin', [0 2], @(x) isnumeric(x) && numel(x)==2 && x(2)>x(1));
p.addParameter('TraceWin', [0.25 0.75], @(x) isnumeric(x) && numel(x)==2 && x(2)>x(1));
p.addParameter('TraceOnly', false, @(x) islogical(x) && isscalar(x));
p.addParameter('PETHWindow', [-1 2], @(x) isnumeric(x) && numel(x)==2 && x(2)>x(1));
p.addParameter('PETHBin', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('PETHSmoothBins', 1, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MIThreshold', 0.95, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinNonTaskEvents', 10, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MinTraceEvents', 3, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MinFieldBins', 5, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('MaxPlaceFields', 2, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('MinCSUSFraction', 0.5, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('OutsideMaxOverlap', 0.15, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('OverlapMinOverlap', 0.50, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('OutsideMaxFieldPeak', 0.20, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('OverlapMinFieldPeak', 0.70, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('InsideMinMapMass', 0.40, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('MixedMinOverlap', 0.25, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('MixedMaxOverlap', 0.75, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('SavePath', '', @(x) ischar(x) || isstring(x));
p.addParameter('SummaryPath', '', @(x) ischar(x) || isstring(x));
p.addParameter('Visible', 'on', @(x) any(strcmpi(x,{'on','off'})));
p.parse(varargin{:});
opt = p.Results;
opt.ExcludeTable = normalizeSelection(opt.ExcludeNeurons);
if opt.TraceOnly
    opt.SpatialComparisonWin = opt.TraceWin;
else
    opt.SpatialComparisonWin = [0 2];
end
ratNames = cellstr(ratNames);

manual = normalizeSelection(opt.SelectedNeurons);
if ~isempty(manual)
    selected = buildManualCandidates(manual, opt);
    selected = labelArchetypesByRow(selected);
    candidates = selected;
else
    candidates = collectCandidates(ratNames, opt);
    if isempty(candidates)
        error('taskPlaceExamples:NoCandidates', ...
            ['No cells passed the spatial/place-field screening gates. Check that the rat ' ...
             'structures (including MI_noCSUS15_shuff and traceneurons) are loaded, ' ...
             'or pass SelectedNeurons to inspect chosen cells directly.']);
    end
    selected = chooseDiverseCandidates(candidates, opt.NumNeurons, opt);
end

[fig, summary] = plotSelectedCandidates(selected, opt);
R.figure = fig;
R.summary = summary;
R.candidates = candidateTable(candidates);
R.selected = selected;
R.options = opt;

disp(summary);
if strlength(string(opt.SummaryPath)) > 0
    writetable(summary, char(opt.SummaryPath));
end
if strlength(string(opt.SavePath)) > 0
    saveReviewerFigure(fig, char(opt.SavePath));
end
end


function C = collectCandidates(ratNames, opt)
C = emptyCandidateStruct();
for r = 1:numel(ratNames)
    ratName = ratNames{r};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratName))
        warning('taskPlaceExamples:MissingAnimal', '%s is not loaded; skipping.', ratName);
        continue
    end
    rat = evalin('base', ratName);
    days = analysisDays(rat, opt.Days);
    for d = 1:numel(days)
        day = days{d};
        peakFld = ['CA_peaks_' day];
        if ~hasDayInputs(rat, day), continue, end
        daySpikes = rat.Ca_peaks.(peakFld);
        nCells = size(daySpikes,1);
        fprintf('[taskPlaceExamples] Screening %s, %s (%d cells)\n', ...
            ratName, day, nCells);
        csQuick = rat.CS_times.(['CS_' day]);
        for cellID = 1:nCells
            if isExcludedCell(opt.ExcludeTable,ratName,day,cellID), continue, end
            if ~fetchRateMask(rat,day,cellID), continue, end
            mi = fetchSpatialMI(rat, day, cellID);
            traceFlag = fetchTraceFlag(rat, day, cellID);
            if ~(isfinite(mi) && mi >= opt.MIThreshold), continue, end
            if isfinite(traceFlag) && traceFlag~=1, continue, end

            % Cheap upper-bound event-count screen before velocity,
            % occupancy maps, field labeling, and trial alignment.
            quickEvents = daySpikes(cellID,:);
            quickEvents = quickEvents(isfinite(quickEvents) & quickEvents>0);
            quickTraceCount = sum(inWindows(quickEvents,csQuick,opt.TraceWin));
            quickFullTaskCount = sum(inWindows(quickEvents,csQuick,[0 2]));
            quickCSUSCount = sum(inWindows(quickEvents,csQuick,[0 opt.TraceWin(2)]));
            minQuickTrace = 1;
            if traceFlag==1, minQuickTrace=opt.MinTraceEvents; end
            if quickTraceCount < minQuickTrace
                continue
            end
            if quickFullTaskCount==0 || ...
                    quickCSUSCount/quickFullTaskCount < opt.MinCSUSFraction
                continue
            end
            if sum(~inWindows(quickEvents,csQuick,opt.TaskExclusionWin)) < opt.MinNonTaskEvents
                continue
            end
            try
                c = computeCandidate(ratName, rat, day, cellID, opt, true);
            catch
                % Sparse/poorly occupied cells are expected during the
                % automatic screen; reject them without flooding output.
                continue
            end
            if c.nonTaskEventCount < opt.MinNonTaskEvents
                continue
            end
            if c.placeFieldCount < 1 || c.placeFieldCount > opt.MaxPlaceFields
                continue
            end
            C(end+1) = c; %#ok<AGROW>
        end
    end
end
if isempty(C), return, end
score = [C.rankScore];
[~, ord] = sort(score, 'descend');
C = C(ord);
end


function C = buildManualCandidates(S, opt)
C = emptyCandidateStruct();
for k = 1:height(S)
    ratName = char(S.Animal(k));
    day = char(S.Session(k));
    cellID = S.CellID(k);
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratName))
        error('taskPlaceExamples:MissingAnimal', '%s is not loaded.', ratName);
    end
    rat = evalin('base', ratName);
    if ~hasDayInputs(rat, day)
        error('taskPlaceExamples:MissingSession', '%s %s lacks position, events, or CS times.', ratName, day);
    end
    C(end+1) = computeCandidate(ratName, rat, day, cellID, opt, false); %#ok<AGROW>
end
end


function c = computeCandidate(ratName, rat, day, cellID, opt, enforceClearField)
peakFld = ['CA_peaks_' day];
posFld = ['pos_' day];
csFld = ['CS_' day];
spikeMat = rat.Ca_peaks.(peakFld);
if cellID < 1 || cellID > size(spikeMat,1) || cellID ~= round(cellID)
    error('Cell ID %g is outside 1:%d.', cellID, size(spikeMat,1));
end
pos = double(rat.pos.(posFld));
csData = rat.CS_times.(csFld);
cs = double(csData(:));
events = double(spikeMat(cellID,:));
events = events(isfinite(events) & events>0);
fullTaskEventCount=sum(inWindows(events,cs,[0 2]));
csUsEventCount=sum(inWindows(events,cs,[0 opt.TraceWin(2)]));
if fullTaskEventCount>0
    csUsEventFraction=csUsEventCount/fullTaskEventCount;
else
    csUsEventFraction=NaN;
end

[nonPos, nonEvents, tracePos, traceEvents] = splitSpatialData(pos, events, cs, opt);
if size(nonPos,1)<5 || size(tracePos,1)<5
    error('Insufficient non-task or trace-period position samples.');
end
if isempty(nonEvents) || isempty(traceEvents)
    error('Insufficient non-task or trace-period events.');
end

nonRaw = CA_normalizePosData(nonEvents(:), nonPos, opt.MapBin, 1);
nonMap = displaySmooth(nonRaw, opt.SmoothSigma);
nonEdges = mapEdges(nonPos, size(nonRaw));

maskRate = nonRaw;
maskRate(maskRate==0) = NaN;
maskCutoff = mean(maskRate(:),'omitnan') + opt.MaskNStd*std(maskRate(:),0,'omitnan');
nonMask = maskRate > maskCutoff;
if ~any(nonMask(:)), error('Non-task spatial mask is empty.'), end
primaryFieldMask = dominantFieldMask(nonMask,3);
[~,nonTaskPeakIdx]=max(nonRaw(:),[],'omitnan');
primaryFieldContainsPeak=~isempty(nonTaskPeakIdx) && ...
    primaryFieldMask(nonTaskPeakIdx);
nonTaskGlobalPeak=max(nonRaw(:),[],'omitnan');
if any(primaryFieldMask(:)) && isfinite(nonTaskGlobalPeak) && nonTaskGlobalPeak>0
    primaryNonTaskPeakRatio=max(nonRaw(primaryFieldMask),[],'omitnan')/nonTaskGlobalPeak;
else
    primaryNonTaskPeakRatio=0;
end
[fieldCount, fieldMass, peakContrast, fieldClarity] = ...
    placeFieldClarity(maskRate, nonMask, opt.MinFieldBins);
if enforceClearField && (fieldCount < 1 || fieldCount > opt.MaxPlaceFields)
    error('Cell does not have one or two clear non-task place fields.');
end

% Only cells surviving the inexpensive metadata/count checks and the
% non-task field check incur trace-map and PETH computation.
traceRaw = CA_normalizePosData(traceEvents(:), tracePos, opt.MapBin, 1);
traceMap = displaySmooth(traceRaw, opt.SmoothSigma);
traceEdges = mapEdges(tracePos, size(traceRaw));
traceMaskRate = traceRaw;
traceMaskRate(traceMaskRate==0) = NaN;
traceCutoff = mean(traceMaskRate(:),'omitnan') + std(traceMaskRate(:),0,'omitnan');
traceHighMask = traceMaskRate > traceCutoff;
[traceFieldCount,~,~,~] = placeFieldClarity(traceMaskRate,traceHighMask,2);

[insideFrac, outsideFrac] = eventFractions(traceEvents, pos, nonMask, nonEdges);
[primaryInsideFrac, primaryOutsideFrac] = ...
    eventFractions(traceEvents, pos, primaryFieldMask, nonEdges);
tracePeakInsidePrimary = mapPeakInsideMask( ...
    traceMap,traceEdges,primaryFieldMask,nonEdges);
[primaryMapMassFraction, primaryMapPeakRatio] = mapActivityInMask( ...
    traceMap,traceEdges,primaryFieldMask,nonEdges);
[pethTime, trialPETH, pethMean, pethSEM] = eventPETH(events, cs, opt);
[traceT, traceP, traceDelta] = pairedTraceModulation(events, cs, opt.TraceWin);
[arenaX, arenaY, arenaSource] = arenaBoundary(rat, day, pos);

validNon = isfinite(nonRaw);
validTrace = isfinite(traceRaw);
validFraction = mean(validNon(:)) * mean(validTrace(:));
baselinePETH = trialPETH(:,pethTime<0);
pethStrength = abs(traceDelta) / (std(baselinePETH(:),0,'omitnan') + eps);
rankScore = 4 + log1p(numel(nonEvents)) + 1.5*log1p(numel(traceEvents)) + ...
    min(abs(traceT),8) + 2*validFraction + min(pethStrength,4) + ...
    20*fieldClarity + 8*csUsEventFraction;

c = struct();
c.animal = ratName;
c.session = day;
c.cellID = cellID;
c.spatialMI = fetchSpatialMI(rat, day, cellID);
c.placeSignificant = isfinite(c.spatialMI) && c.spatialMI >= opt.MIThreshold;
c.traceFlag = fetchTraceFlag(rat, day, cellID);
c.traceStatistic = traceT;
c.traceP = traceP;
c.traceModulated = (c.traceFlag==1) || ...
    (~isfinite(c.traceFlag) && isfinite(traceP) && traceP<0.05);
c.traceRateDelta = traceDelta;
c.csUsEventFraction = csUsEventFraction;
c.nonTaskPeakRate = max(nonMap(:),[],'omitnan');
c.tracePeakRate = max(traceMap(:),[],'omitnan');
c.traceInsideFraction = insideFrac;
c.traceOutsideFraction = outsideFrac;
c.primaryInsideFraction = primaryInsideFrac;
c.primaryOutsideFraction = primaryOutsideFrac;
c.tracePeakInsidePrimary = tracePeakInsidePrimary;
c.primaryMapMassFraction = primaryMapMassFraction;
c.primaryMapPeakRatio = primaryMapPeakRatio;
c.primaryFieldContainsPeak = primaryFieldContainsPeak;
c.primaryNonTaskPeakRatio = primaryNonTaskPeakRatio;
c.nonTaskEventCount = numel(nonEvents);
c.traceEventCount = numel(traceEvents);
c.placeFieldCount = fieldCount;
c.placeFieldMass = fieldMass;
c.peakContrast = peakContrast;
c.fieldClarity = fieldClarity;
c.traceFieldCount = traceFieldCount;
c.archetype = '';
c.rankScore = rankScore;
c.nonMap = nonMap;
c.traceMap = traceMap;
c.nonMask = nonMask;
c.primaryFieldMask = primaryFieldMask;
c.nonEdges = nonEdges;
c.traceEdges = traceEdges;
c.pethTime = pethTime;
c.trialPETH = trialPETH;
c.pethMean = pethMean;
c.pethSEM = pethSEM;
c.arenaX = arenaX;
c.arenaY = arenaY;
c.arenaSource = arenaSource;
c.maskCutoff = maskCutoff;
end


function primaryMask = dominantFieldMask(mask,minArea)
% Display only the neuron's primary established field. The complete mask is
% retained separately for classification and event-fraction calculations.
primaryMask=false(size(mask));
cc=bwconncomp(mask,8);
if cc.NumObjects==0, return, end
componentArea=cellfun(@numel,cc.PixelIdxList);
[largestArea,idx]=max(componentArea);
if largestArea>=minArea
    primaryMask(cc.PixelIdxList{idx})=true;
end
end


function [massFraction,peakRatio] = mapActivityInMask(rateMap,rateEdges,mask,maskEdges)
% Quantify overlap exactly in the display coordinate system. peakRatio asks
% whether any strong comparison-period hotspot remains inside the boundary,
% even when the global peak or most individual events occur elsewhere.
[X,Y]=meshgrid(rateEdges.xCenters,rateEdges.yCenters);
bx=discretize(X,maskEdges.x); by=discretize(Y,maskEdges.y);
valid=isfinite(bx) & isfinite(by) & bx>=1 & by>=1 & ...
    bx<=size(mask,2) & by<=size(mask,1) & isfinite(rateMap);
inside=false(size(rateMap));
if any(valid(:))
    validLinear=find(valid);
    maskLinear=sub2ind(size(mask),by(valid),bx(valid));
    inside(validLinear)=mask(maskLinear);
end
positiveMap=max(rateMap,0);
totalMass=sum(positiveMap(valid),'omitnan');
if totalMass>0
    massFraction=sum(positiveMap(inside & valid),'omitnan')/totalMass;
else
    massFraction=NaN;
end
globalPeak=max(positiveMap(valid),[],'omitnan');
if globalPeak>0 && any(inside(:))
    peakRatio=max(positiveMap(inside),[],'omitnan')/globalPeak;
else
    peakRatio=0;
end
end


function peakInside = mapPeakInsideMask(rateMap,rateEdges,mask,maskEdges)
peakInside=false;
if isempty(rateMap) || ~any(isfinite(rateMap(:))), return, end
[~,linearIdx]=max(rateMap(:),[],'omitnan');
[row,col]=ind2sub(size(rateMap),linearIdx);
x=rateEdges.xCenters(col); y=rateEdges.yCenters(row);
bx=discretize(x,maskEdges.x); by=discretize(y,maskEdges.y);
if isfinite(bx) && isfinite(by) && bx>=1 && by>=1 && ...
        bx<=size(mask,2) && by<=size(mask,1)
    peakInside=mask(by,bx);
end
end


function [nFields, fieldMass, peakContrast, clarity] = placeFieldClarity(rate, mask, minBins)
% Count only substantial 8-connected mask components. This matches the
% existing place-field convention of requiring at least five adjacent bins
% while allowing the Fig. S20 mean+N*SD threshold to define the mask itself.
cc = bwconncomp(mask,8);
componentSize = cellfun(@numel,cc.PixelIdxList);
major = find(componentSize >= minBins);
nFields = numel(major);

validRate = rate(isfinite(rate) & rate>0);
if isempty(validRate)
    fieldMass=0; peakContrast=0; clarity=0; return
end
majorIdx = [];
for k=major(:)'
    majorIdx = [majorIdx; cc.PixelIdxList{k}(:)]; %#ok<AGROW>
end
if isempty(majorIdx)
    fieldMass=0;
else
    fieldMass=sum(rate(majorIdx),'omitnan')/sum(validRate,'omitnan');
end
peakContrast=max(validRate)/max(mean(validRate,'omitnan'),eps);

% Concentrated rate mass and high peak/background contrast both improve
% clarity. One field receives a small preference over two; >2 is gated in
% automatic selection but remains available through manual overrides.
fieldNumberWeight = double(nFields==1) + 0.85*double(nFields==2);
clarity = fieldNumberWeight * fieldMass * min(log2(max(peakContrast,1)),3);
end


function [nonPos, nonEvents, tracePos, traceEvents] = splitSpatialData(pos, events, cs, opt)
t = pos(:,1);
xy = pos(:,2:3);
vel = ca_velocity(pos);
velTime = vel(2,:)';
velMag = vel(1,:)';
vx = interp1(t, xy(:,1), velTime, 'linear', NaN);
vy = interp1(t, xy(:,2), velTime, 'linear', NaN);

inTrialVel = inWindows(velTime, cs, opt.TaskExclusionWin);
keep = velMag >= opt.VelThresh & ~inTrialVel & isfinite(vx) & isfinite(vy);
nonPos = [velTime(keep), vx(keep), vy(keep)];
eventVel = interp1(velTime, velMag, events, 'linear', NaN);
nonEventKeep = ~inWindows(events, cs, opt.TaskExclusionWin) & eventVel >= opt.VelThresh;
nonEvents = events(nonEventKeep);

traceKeep = inWindows(t, cs, opt.SpatialComparisonWin) & all(isfinite(xy),2);
tracePos = pos(traceKeep,1:3);
traceEvents = events(inWindows(events, cs, opt.SpatialComparisonWin));
end


function tf = inWindows(t, anchors, win)
tf = false(size(t));
for k = 1:numel(anchors)
    tf = tf | (t >= anchors(k)+win(1) & t < anchors(k)+win(2));
end
end


function rate = displaySmooth(rate, sigma)
rate(isnan(rate)) = 0;
if sigma > 0
    rate = imgaussfilt(rate, sigma);
end
end


function edges = mapEdges(pos, mapSize)
edges.x = linspace(min(pos(:,2)), max(pos(:,2)), mapSize(2)+1);
edges.y = linspace(min(pos(:,3)), max(pos(:,3)), mapSize(1)+1);
edges.xCenters = (edges.x(1:end-1)+edges.x(2:end))/2;
edges.yCenters = (edges.y(1:end-1)+edges.y(2:end))/2;
end


function [inside, outside] = eventFractions(events, pos, mask, edges)
if isempty(events)
    inside = NaN; outside = NaN; return
end
sx = interp1(pos(:,1), pos(:,2), events, 'linear', NaN);
sy = interp1(pos(:,1), pos(:,3), events, 'linear', NaN);
bx = discretize(sx, edges.x);
by = discretize(sy, edges.y);
valid = isfinite(bx) & isfinite(by) & bx>=1 & by>=1 & bx<=size(mask,2) & by<=size(mask,1);
if ~any(valid)
    inside = NaN; outside = NaN; return
end
in = mask(sub2ind(size(mask), by(valid), bx(valid)));
inside = mean(in);
outside = 1-inside;
end


function [t, trialActivity, m, sem] = eventPETH(events, cs, opt)
edges = opt.PETHWindow(1):opt.PETHBin:opt.PETHWindow(2);
t = edges(1:end-1) + opt.PETHBin/2;
trialRate = zeros(numel(cs),numel(t));
for tr = 1:numel(cs)
    rel = events(events>=cs(tr)+edges(1) & events<cs(tr)+edges(end))-cs(tr);
    trialRate(tr,:) = histcounts(rel,edges)/opt.PETHBin;
end
if opt.PETHSmoothBins>0
    w = ones(1,2*round(opt.PETHSmoothBins)+1);
    w = w/sum(w);
    for tr = 1:size(trialRate,1)
        trialRate(tr,:) = conv(trialRate(tr,:),w,'same');
    end
end
trialActivity = trialRate;
m = mean(trialActivity,1,'omitnan');
n = sum(isfinite(trialActivity),1);
sem = std(trialActivity,0,1,'omitnan')./sqrt(max(n,1));
end


function [tstat, pval, deltaRate] = pairedTraceModulation(events, cs, traceWin)
dur = diff(traceWin);
preWin = [-dur 0];
traceRate = zeros(numel(cs),1);
preRate = zeros(numel(cs),1);
for tr = 1:numel(cs)
    traceRate(tr) = sum(events>=cs(tr)+traceWin(1) & events<cs(tr)+traceWin(2))/dur;
    preRate(tr) = sum(events>=cs(tr)+preWin(1) & events<cs(tr)+preWin(2))/dur;
end
deltaRate = mean(traceRate-preRate,'omitnan');
if numel(cs)>=2 && any(traceRate~=preRate)
    [~,pval,~,st] = ttest(traceRate,preRate);
    tstat = st.tstat;
else
    tstat = NaN; pval = NaN;
end
end


function selected = chooseDiverseCandidates(C, n, opt)
selected = C([]);
used = false(1,numel(C));
inside = [C.primaryInsideFraction];
score = robustScale([C.rankScore]);
fieldScore = robustScale([C.fieldClarity]);
traceMod = [C.traceModulated];
peakInside = [C.tracePeakInsidePrimary];
fieldPeakRatio = [C.primaryMapPeakRatio];
mapMass = [C.primaryMapMassFraction];
primaryNonTaskPeakRatio = [C.primaryNonTaskPeakRatio];
if numel(C)<n
    error('taskPlaceExamples:InsufficientArchetypes', ...
        ['Not enough cells passed screening for the requested number of examples. ' ...
         'Relax MinFieldBins, request fewer examples, or use SelectedNeurons.']);
end

% Three spatial relationships. With the default n=3 this selects one cell
% whose activity is outside, one inside, and one split across the field.
labels = {'Outside primary field', ...
    'Inside primary field', ...
    'Mixed inside and outside'};
counts = categoryCounts(n);
targets = { ...
    score + 2*fieldScore + 3*(1-inside), ...                    % outside
    score + 2*fieldScore + 3*inside + 3*mapMass, ...            % inside
    score + 2*fieldScore - 2*abs(inside-0.5) ...                % mixed
        - 2*abs(mapMass-0.5) + 3*primaryNonTaskPeakRatio};
eligible = { ...
    traceMod & ~peakInside & isfinite(inside) & inside<=opt.OutsideMaxOverlap & ...
        isfinite(fieldPeakRatio) & fieldPeakRatio<=opt.OutsideMaxFieldPeak, ...
    traceMod & peakInside & isfinite(inside) & isfinite(mapMass) & ...
        (inside>=opt.OverlapMinOverlap | mapMass>=opt.InsideMinMapMass) & ...
        fieldPeakRatio>=opt.OverlapMinFieldPeak, ...
    traceMod & isfinite(inside) & inside>=opt.MixedMinOverlap & ...
        inside<=opt.MixedMaxOverlap & isfinite(mapMass) & ...
        mapMass>=opt.MixedMinOverlap & mapMass<=opt.MixedMaxOverlap & ...
        primaryNonTaskPeakRatio>=0.75};
for category=1:3
    for exampleNo=1:counts(category)
        q=targets{category};
        q(used | ~eligible{category} | ~isfinite(q))=-inf;
        [best,idx]=max(q);
        if ~isfinite(best) && category==2
            % Visual inside-field fallback: require the displayed map peak
            % to lie inside the displayed primary-field boundary.
            relaxed=traceMod & peakInside & ~used;
            q=targets{category}; q(~relaxed | ~isfinite(q))=-inf;
            [best,idx]=max(q);
        elseif ~isfinite(best) && category==3
            % Mixed fallback retains nontrivial activity on both sides.
            relaxed=traceMod & isfinite(inside) & ~used & ...
                inside>0.10 & inside<0.90 & primaryNonTaskPeakRatio>=0.50;
            q=targets{category}; q(~relaxed | ~isfinite(q))=-inf;
            [best,idx]=max(q);
            if ~isfinite(best)
                relaxed=traceMod & ~used & isfinite(inside) & inside>0 & inside<1 & ...
                    isfinite(mapMass) & mapMass>0 & mapMass<1 & ...
                    isfinite(primaryNonTaskPeakRatio);
                q=targets{category}; q(~relaxed | ~isfinite(q))=-inf;
                [best,idx]=max(q);
            end
        end
        if ~isfinite(best)
            error('taskPlaceExamples:MissingArchetype', ...
                ['Not enough eligible "%s" examples passed screening. ' ...
                 'Adjust thresholds, request fewer examples, or use SelectedNeurons.'], ...
                labels{category});
        end
        selected(end+1)=C(idx); %#ok<AGROW>
        selected(end).archetype=labels{category};
        used(idx)=true;
    end
end
end


function C = labelArchetypesByRow(C)
labels={'Outside primary field','Inside primary field','Mixed inside and outside'};
counts=categoryCounts(numel(C));
idx=0;
for category=1:3
    for exampleNo=1:counts(category)
        idx=idx+1;
        C(idx).archetype=labels{category};
    end
end
end


function counts = categoryCounts(n)
counts=floor(n/3)*ones(1,3);
counts(1:mod(n,3))=counts(1:mod(n,3))+1;
end


function z = robustScale(x)
x = double(x);
lo = min(x,[],'omitnan'); hi = max(x,[],'omitnan');
if ~isfinite(lo) || ~isfinite(hi) || hi<=lo, z = zeros(size(x)); else, z=(x-lo)/(hi-lo); end
end


function [fig, T] = plotSelectedCandidates(C, opt)
n = numel(C);
fig = figure('Color','w','Visible',opt.Visible,'Position',[50 40 1500 max(700,245*n)]);
tl = tiledlayout(fig,n,3,'TileSpacing','compact','Padding','loose');
titles = {'A. Non-tEBC spatial rate map', ...
    'B. Task-period spatial map (0–2 s)', ...
    'C. CS-aligned PETH'};
if opt.TraceOnly
    titles{2}='B. Trace interval spatial map (0.25–0.75 s)';
end
for row = 1:n
    c = C(row);
    % Display normalization is condition-specific; raw maps and the mask
    % remain untouched for rate summaries and spatial classification.
    nonMapDisplay = peakNormalizeMap(c.nonMap,opt.MapNormalizationPercentile);
    traceMapDisplay = peakNormalizeMap(c.traceMap,opt.MapNormalizationPercentile);
    nonColorLim = normalizedColorLimits(nonMapDisplay);
    traceColorLim = normalizedColorLimits(traceMapDisplay);

    ax = nexttile(tl,(row-1)*3+1);
    plotRateMap(ax,nonMapDisplay,c.nonEdges,nonColorLim,[]);
    if row==1, title(ax,titles{1},'FontWeight','bold','FontSize',14); end
    text(ax,0.02,0.98,sprintf('%s | %s | cell %d | peaks %.2f / %.2f Hz',c.animal, ...
        strrep(c.session,'_','-'),c.cellID,c.nonTaskPeakRate,c.tracePeakRate), ...
        'Units','normalized', ...
        'VerticalAlignment','top','FontSize',8,'Color',[0.15 0.15 0.15], ...
        'BackgroundColor',[1 1 1],'Margin',2,'Interpreter','none');

    ax = nexttile(tl,(row-1)*3+2);
    plotRateMap(ax,traceMapDisplay,c.traceEdges,traceColorLim,[]);
    overlayMaskOutline(ax,c.primaryFieldMask,c.nonEdges,1.5);
    if row==1, title(ax,titles{2},'FontWeight','bold','FontSize',14); end

    ax = nexttile(tl,(row-1)*3+3); hold(ax,'on');
    plotPETH(ax,c,opt);
    if row==1, title(ax,titles{3},'FontWeight','bold','FontSize',14); end
end
title(tl,'Representative task-modulated place cells','FontSize',17,'FontWeight','bold');
questionText=['Where does task/trace activity occur relative to the neuron''s ' ...
    'primary established place field?'];
subtitle(tl,questionText,'FontSize',12);
if opt.TraceOnly
    intervalCaption=['trace-period activity can occur outside, remain inside, or span ' ...
        'both regions relative to established non-tEBC spatial firing patterns. '];
else
    intervalCaption=['task-period activity can occur outside, remain inside, or span ' ...
        'both regions relative to established non-tEBC spatial firing patterns. '];
end
caption = ['Examples show that ' intervalCaption ...
    'Spatial maps are peak-normalized within each condition to visualize the location ' ...
    'of activity; absolute event rates are quantified separately in Fig. Sx.'];
annotation(fig,'textbox',[0.08 0.002 0.84 0.045],'String',caption, ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',10, ...
    'FontAngle','italic','Interpreter','none');
T = candidateTable(C);
end


function rateNorm = peakNormalizeMap(rate,upperPercentile)
positiveValues=rate(isfinite(rate) & rate>0);
if ~isempty(positiveValues)
    displayCeiling=prctile(positiveValues,upperPercentile);
    if ~isfinite(displayCeiling) || displayCeiling<=0
        displayCeiling=max(positiveValues);
    end
    rateNorm=min(rate/displayCeiling,1);
else
    rateNorm=zeros(size(rate));
end
end


function colorLim = normalizedColorLimits(rateNorm)
positiveValues=rateNorm(isfinite(rateNorm) & rateNorm>0);
if isempty(positiveValues)
    lowerLimit=0.01;
else
    lowerLimit=prctile(positiveValues,10);
    lowerLimit=max(min(lowerLimit,0.95),0.05);
end
colorLim=[lowerLimit 1];
end


function plotRateMap(ax,rate,edges,colorLim,mask)
imagesc(ax,edges.xCenters,edges.yCenters,rate);
set(ax,'YDir','normal','Box','off','TickDir','out');
axis(ax,'image');
colormap(ax,parula);
set(ax,'CLim',colorLim);
if ~isempty(mask), overlayMaskOutline(ax,mask,edges,4); end
cb=colorbar(ax); cb.Label.String='Normalized event rate'; cb.FontSize=11; cb.Label.FontSize=12;
set(ax,'FontSize',12,'LineWidth',1);
end


function overlayMaskOutline(ax,mask,edges,lineWidth)
hold(ax,'on');
[X,Y]=meshgrid(edges.xCenters,edges.yCenters);
if any(mask(:)) && any(~mask(:))
    contour(ax,X,Y,double(mask),[0.5 0.5],'Color','w', ...
        'LineStyle','-','LineWidth',lineWidth);
end
end


function plotPETH(ax,c,opt)
t=c.pethTime; m=c.pethMean; s=c.pethSEM;
yl=[min(m-s,[],'omitnan'),max(m+s,[],'omitnan')];
if ~all(isfinite(yl)) || yl(2)<=yl(1), yl=[-1 1]; end
pad=0.08*diff(yl); yl=yl+[-pad pad];
patch(ax,[opt.TraceWin(1) opt.TraceWin(2) opt.TraceWin(2) opt.TraceWin(1)], ...
    [yl(1) yl(1) yl(2) yl(2)],[0.25 0.65 0.35], ...
    'FaceAlpha',0.27,'EdgeColor','none');
fill(ax,[t fliplr(t)],[m+s fliplr(m-s)],[0.16 0.16 0.16], ...
    'FaceAlpha',0.18,'EdgeColor','none');
plot(ax,t,m,'k','LineWidth',1.5);
xline(ax,0,'k--','LineWidth',1.4);
xline(ax,opt.TraceWin(2),'--','Color',[0.85 0.10 0.10],'LineWidth',1.4);
text(ax,opt.TraceWin(2),yl(2),' US','Color',[0.65 0.15 0.10], ...
    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',8);
xlim(ax,opt.PETHWindow); ylim(ax,yl);
yline(ax,0,':','Color',[0.65 0.65 0.65]);
xlabel(ax,'Time from CS onset (s)');
ylabel(ax,'Event rate (Hz)');
set(ax,'Box','off','TickDir','out','FontSize',12,'LineWidth',1);
end


function [x,y,source] = arenaBoundary(rat,day,pos)
x=[]; y=[]; source='';
names={'arenaBoundary','arena_boundary','arena','boundaries'};
for k=1:numel(names)
    if ~isfield(rat,names{k}), continue, end
    v=rat.(names{k});
    if isstruct(v)
        fields={day,['boundary_' day],['arena_' day]};
        for j=1:numel(fields)
            if isfield(v,fields{j}), v=v.(fields{j}); break, end
        end
    end
    if isnumeric(v) && size(v,2)>=2 && size(v,1)>=3
        if size(v,2)>=3, v=v(:,end-1:end); else, v=v(:,1:2); end
        ok=all(isfinite(v),2); v=v(ok,:);
        if size(v,1)>=3, x=[v(:,1);v(1,1)]; y=[v(:,2);v(1,2)]; source=names{k}; return, end
    end
end
xy=pos(:,2:3); ok=all(isfinite(xy),2); xy=xy(ok,:);
if size(xy,1)>=3
    try
        idx=boundary(xy(:,1),xy(:,2),0.9);
    catch
        idx=convhull(xy(:,1),xy(:,2));
    end
    x=xy(idx,1); y=xy(idx,2); source='position-derived';
end
end


function mi = fetchSpatialMI(rat,day,cellID)
mi=NaN; fld=['MI_' day];
if isfield(rat,'MI_noCSUS15_shuff') && isfield(rat.MI_noCSUS15_shuff,fld)
    v=rat.MI_noCSUS15_shuff.(fld);
    if size(v,1)>=cellID && size(v,2)>=3, mi=v(cellID,3); end
end
end


function flag = fetchTraceFlag(rat,day,cellID)
flag=NaN; fld=['tn_' day];
if isfield(rat,'traceneurons') && isfield(rat.traceneurons,fld)
    v=rat.traceneurons.(fld);
    if numel(v)>=cellID, flag=v(cellID); end
end
end


function keep = fetchRateMask(rat,day,cellID)
keep=true; fld=['ratemask_' day];
if isfield(rat,'ratemask') && isfield(rat.ratemask,fld)
    v=rat.ratemask.(fld);
    if numel(v)>=cellID, keep=isfinite(v(cellID)) && logical(v(cellID)); end
end
end


function tf = hasDayInputs(rat,day)
tf=isfield(rat,'pos') && isfield(rat.pos,['pos_' day]) && ...
   isfield(rat,'Ca_peaks') && isfield(rat.Ca_peaks,['CA_peaks_' day]) && ...
   isfield(rat,'CS_times') && isfield(rat.CS_times,['CS_' day]);
end


function days = analysisDays(rat,requested)
if ~isempty(requested), days=cellstr(requested); return, end
if exist('autoDateList','file')==2
    days=autoDateList(rat);
else
    fields=fieldnames(rat.Ca_peaks);
    fields=fields(startsWith(fields,'CA_peaks_'));
    days=cellfun(@(s) erase(s,'CA_peaks_'),fields,'UniformOutput',false);
    days=sort(days);
end
if isfield(rat,'An'), idx=find(strcmp(days,rat.An),1); else, idx=[]; end
if isempty(idx), idx=numel(days); end
days=days(max(1,idx-2):idx);
end


function S = normalizeSelection(x)
if isempty(x), S=table; return, end
if istable(x)
    vars=lower(string(x.Properties.VariableNames));
    ia=find(ismember(vars,["animal","rat"]),1);
    id=find(ismember(vars,["session","day","date"]),1);
    ic=find(ismember(vars,["cellid","cell","neuron"]),1);
    if isempty(ia)||isempty(id)||isempty(ic)
        error('SelectedNeurons table needs Animal, Session, and CellID columns.');
    end
    S=table(string(x{:,ia}),string(x{:,id}),double(x{:,ic}), ...
        'VariableNames',{'Animal','Session','CellID'});
elseif isstruct(x)
    S=struct2table(x);
    S=normalizeSelection(S);
elseif iscell(x) && size(x,2)==3
    S=table(string(x(:,1)),string(x(:,2)),cell2mat(x(:,3)), ...
        'VariableNames',{'Animal','Session','CellID'});
else
    error('SelectedNeurons must be a table, struct, or N-by-3 cell array.');
end
end


function tf = isExcludedCell(S,animal,session,cellID)
if isempty(S)
    tf=false;
    return
end
tf=any(S.Animal==string(animal) & S.Session==string(session) & ...
    S.CellID==cellID);
end


function T = candidateTable(C)
if isempty(C)
    T=table; return
end
T=table(string({C.archetype})',string({C.animal})',string({C.session})',[C.cellID]',[C.spatialMI]', ...
    [C.traceModulated]',[C.traceStatistic]',[C.traceP]',[C.csUsEventFraction]', ...
    [C.nonTaskPeakRate]',[C.tracePeakRate]', ...
    [C.traceInsideFraction]',[C.traceOutsideFraction]', ...
    [C.primaryInsideFraction]',[C.primaryOutsideFraction]', ...
    [C.tracePeakInsidePrimary]',[C.primaryMapMassFraction]', ...
    [C.primaryMapPeakRatio]',[C.primaryFieldContainsPeak]', ...
    [C.primaryNonTaskPeakRatio]',[C.nonTaskEventCount]', ...
    [C.traceEventCount]',[C.placeFieldCount]',[C.placeFieldMass]', ...
    [C.peakContrast]',[C.fieldClarity]',[C.traceFieldCount]',[C.rankScore]', ...
    'VariableNames',{'Archetype','Animal','Session','CellID','SpatialMI','TraceModulated', ...
    'TraceStatistic','TraceP','CSUSEventFraction', ...
    'NonTaskPeakRate','TracePeakRate','TraceInsideMask','TraceOutsideMask', ...
    'PrimaryFieldInside','PrimaryFieldOutside','TracePeakInsidePrimary', ...
    'PrimaryMapMassFraction','PrimaryMapPeakRatio','PrimaryFieldContainsPeak', ...
    'PrimaryNonTaskPeakRatio','NonTaskEvents','TraceEvents', ...
    'PlaceFieldCount','PlaceFieldMass', ...
    'PeakContrast','FieldClarity','TraceFieldCount','RankScore'});
end


function saveReviewerFigure(fig,path)
[folder,~,ext]=fileparts(path);
if ~isempty(folder) && ~exist(folder,'dir'), mkdir(folder), end
switch lower(ext)
    case '.fig', savefig(fig,path);
    case {'.pdf','.svg'}, exportgraphics(fig,path,'ContentType','vector');
    otherwise, exportgraphics(fig,path,'Resolution',600);
end
end


function C = emptyCandidateStruct()
C=struct('animal',{},'session',{},'cellID',{},'spatialMI',{}, ...
    'placeSignificant',{},'traceFlag',{},'traceStatistic',{},'traceP',{}, ...
    'traceModulated',{},'traceRateDelta',{},'csUsEventFraction',{}, ...
    'nonTaskPeakRate',{},'tracePeakRate',{}, ...
    'traceInsideFraction',{},'traceOutsideFraction',{}, ...
    'primaryInsideFraction',{},'primaryOutsideFraction',{}, ...
    'tracePeakInsidePrimary',{},'primaryMapMassFraction',{}, ...
    'primaryMapPeakRatio',{},'primaryFieldContainsPeak',{}, ...
    'primaryNonTaskPeakRatio',{},'nonTaskEventCount',{}, ...
    'traceEventCount',{},'placeFieldCount',{},'placeFieldMass',{}, ...
    'peakContrast',{},'fieldClarity',{},'traceFieldCount',{},'archetype',{},'rankScore',{}, ...
    'nonMap',{},'traceMap',{},'nonMask',{},'primaryFieldMask',{}, ...
    'nonEdges',{},'traceEdges',{},'pethTime',{},'trialPETH',{},'pethMean',{}, ...
    'pethSEM',{},'arenaX',{},'arenaY',{},'arenaSource',{},'maskCutoff',{});
end
