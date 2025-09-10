function R = populationSpatialOrthogonality(ratName, varargin)
% populationSpatialOrthogonality
% Build K×N population vectors over K linear spatial bins (PC1; v>=MinSpeed),
% compute similarities (cosine OR Pearson correlation), and PLOT
% (heatmaps + pooled MDS + grouped bars).
%
% Usage:
%   R = populationSpatialOrthogonality('rat0314');
%   R = populationSpatialOrthogonality('rat0314','NBins',6,'BinMode','equal_size',...
%         'MinSpeed',4,'Mode','zscore','Similarity','corr');
%   R = populationSpatialOrthogonality('rat0314','Days','2022_08_13');               % single day
%   R = populationSpatialOrthogonality('rat0314','Days',{'2022_08_13','2022_08_15'});% custom list
%
% Name-Value:
%   'NBins'            – # spatial bins (default 6)
%   'BinMode'          – 'equal_occ' | 'equal_size' (default 'equal_size')
%   'MinSpeed'         – cm/s threshold on velocity timebase (default 4)
%   'Mode'             – 'raw' | 'demean' | 'zscore' (default 'raw')
%   'Similarity'       – 'cosine' | 'corr' (Pearson) (default 'cosine')
%   'Days'             – char/string (one date) or cellstr of dates; if empty, uses last up-to-3 up to An
%   'MICutoff'         – [] (off) or scalar in [0,1]; INCLUDE if MI(:,3) >= cutoff
%   'MIExcludeCutoff'  – [] (off) or scalar in [0,1]; EXCLUDE if MI(:,3) >= cutoff (from MI_CSUS15_shuff)
%   'DayWindowMode'    – 'fraction' (0..1) or 'seconds' (relative to vt start) (default 'fraction')
%   'DayWindow'        – [] or [a b] with b>a (default [])
%   'SaveFig'          – true/false (default false)
%   'FigName'          – filename when SaveFig=true

% ---- args ----
p = inputParser;
addParameter(p,'NBins',6, @(x) isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'BinMode','equal_size',@(s) any(strcmpi(s,{'equal_occ','equal_size'})));
addParameter(p,'MinSpeed',4,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Mode','raw',@(s) any(strcmpi(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','cosine',@(s) any(strcmpi(s,{'cosine','corr','pearson'})));
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty(d));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));  % NEW
addParameter(p,'DayWindowMode','fraction',@(s) any(strcmpi(s,{'fraction','seconds'})));
addParameter(p,'DayWindow',[],@(v) isempty(v) || (isnumeric(v)&&numel(v)==2&&v(2)>v(1)));
addParameter(p,'SaveFig',false,@islogical);
addParameter(p,'FigName','',@(s) ischar(s) || isstring(s));
parse(p,varargin{:});

K            = p.Results.NBins;
binMode      = lower(p.Results.BinMode);
vMin         = p.Results.MinSpeed;
modeStr      = lower(p.Results.Mode);
simStr       = lower(p.Results.Similarity);
daysArg      = p.Results.Days;
miCutoff     = p.Results.MICutoff;
miExCutoff   = p.Results.MIExcludeCutoff;  % NEW
doSave       = p.Results.SaveFig;
figName      = string(p.Results.FigName);
dayWinMode   = lower(p.Results.DayWindowMode);
dayWin       = p.Results.DayWindow;

% ---- compute (helper below) ----
R = populationLinearSpatialOrthogonality_byRat( ...
        ratName, 'NBins',K, 'BinMode',binMode, 'MinSpeed',vMin, ...
        'Mode',modeStr, 'Similarity',simStr, 'Days',daysArg, ...
        'MICutoff',miCutoff, 'MIExcludeCutoff',miExCutoff, ... % NEW
        'DayWindowMode',dayWinMode, 'DayWindow',dayWin);

% ---- collect per-day panels (similarity & norms) ----
nDays  = numel(R.perDay);
days   = R.days;

have       = false(1,nDays);
meanOff    = nan(1,nDays);
meanAng    = nan(1,nDays);
rowNormsAll= cell(1,nDays);

for d = 1:nDays
  if d>numel(R.perDay) || isempty(R.perDay(d).cosSimK), continue; end
  C = R.perDay(d).cosSimK;
  have(d) = true;
  offMask = ~eye(K);
  meanOff(d) = mean(C(offMask),'omitnan');
  if strcmpi(R.params.Similarity,'cosine')
      meanAng(d) = real(acos(max(-1,min(1,meanOff(d)))));
  else
      meanAng(d) = NaN;
  end

  if ~isempty(R.perDay(d).rateKxN_raw)
      Xmode = normalizeRatesPerCell(R.perDay(d).rateKxN_raw, R.params.Mode);
      rowNormsAll{d} = vecnorm(Xmode, 2, 2);
  else
      rowNormsAll{d} = nan(K,1);
  end
end

% pooled summaries
Cpool = R.pooled.cosSimK;
offMask = ~eye(K);
meanOff_pool = mean(Cpool(offMask),'omitnan');
if strcmpi(R.params.Similarity,'cosine')
    meanAng_pool = real(acos(max(-1,min(1,meanOff_pool))));
else
    meanAng_pool = NaN;
end

% plotting ranges/labels by metric
isCorr   = any(strcmpi(R.params.Similarity, {'corr','pearson'}));
labShort = tern(isCorr,"corr","cos");
clims    = tern(isCorr,[-1 1],[0 1]);

% ---- figure layout ----
nHeat = sum(have) + 1;
nCols = max(3, nHeat);
nRows = 2;

fig = figure('Color','w','Position',[100 100 320*nCols 800]);
tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');

% ---- heatmaps (per-day + pooled) ----
plotIdx = 0;
for d = 1:nDays
  if ~have(d), continue; end
  plotIdx = plotIdx + 1;
  nexttile(plotIdx);
  imagesc(R.perDay(d).cosSimK, clims); axis square; colormap(gca, parula); colorbar;
  xticks(1:K); yticks(1:K); xlabel('Bin'); ylabel('Bin');
  if isCorr
      title(sprintf('%s  |  mean off-diag %s = %.2f', days{d}, labShort, meanOff(d)));
  else
      title(sprintf('%s  |  mean off-diag %s = %.2f (%.0f°)', ...
            days{d}, labShort, meanOff(d), rad2deg(meanAng(d))));
  end
end
% pooled
plotIdx = plotIdx + 1;
nexttile(plotIdx);
imagesc(Cpool, clims); axis square; colormap(gca, parula); colorbar;
xticks(1:K); yticks(1:K); xlabel('Bin'); ylabel('Bin');
if isCorr
    title(sprintf('Pooled  |  mean off-diag %s = %.2f', labShort, meanOff_pool));
else
    title(sprintf('Pooled  |  mean off-diag %s = %.2f (%.0f°)', ...
          labShort, meanOff_pool, rad2deg(meanAng_pool)));
end

% ---- MDS of pooled (distance = 1 - similarity) ----
nexttile([1 max(1, floor(nCols/2))]);
D = 1 - max(-1,min(1,Cpool));  D = (D + D')/2; D(1:K+1:end) = 0;
try
  Y = mdscale(D, 2, 'Criterion','stress');
catch
  [V,~] = eigs((Cpool + Cpool')/2, 2);
  Y = V;
end
plot(Y(:,1), Y(:,2), '-o','LineWidth',1.5,'MarkerFaceColor',[.1 .1 .1]); hold on
text(Y(:,1), Y(:,2), compose('%d',1:K), 'VerticalAlignment','bottom','FontSize',10);
axis equal; grid on; xlabel('MDS-1'); ylabel('MDS-2');
title(sprintf('Pooled: MDS of %d spatial bins (dist = 1 - %s)', K, labShort));

% ---- grouped bars: vector norms per bin (per day + per-bin mean) ----
nexttile([1 max(1, ceil(nCols/2))]); cla; hold on
validIdx = find(have);
nValid   = numel(validIdx);

if nValid == 0
    text(0.5,0.5,'No valid days after filtering.', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
    axis off;
else
    Rnorms = nan(K, nValid); lbls = strings(1, nValid);
    for j = 1:nValid
        d = validIdx(j);
        Rnorms(:, j) = rowNormsAll{d}(:);
        lbls(j)      = string(days{d});
    end
    rnMean = mean(Rnorms, 2, 'omitnan');
    Ybars  = [Rnorms, rnMean];

    ymax = max(Ybars(:), [], 'omitnan');
    if isfinite(ymax) && ymax > 0
        b = bar(Ybars, 'grouped'); box on
        for j = 1:nValid, b(j).DisplayName = lbls(j); end
        b(nValid+1).DisplayName = 'Mean(non-NaN)';
        b(nValid+1).FaceAlpha   = 0.45;

        xticks(1:K); xlim([0.5 K+0.5])
        xlabel('Spatial bin (PC1 order)');
        ylabel(sprintf('||population||_2 (%s)', upper(R.params.Mode)))
        title('Population vector norms per spatial bin');
        legend('Location','bestoutside')
        ylim([0, ymax * 1.1])
    else
        text(0.5,0.5,'No finite bar values (after DayWindow / filters).', ...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
        axis off;
    end
end

sgtitle(sprintf('%s | K=%d bins | %s | v\\geq%.0f cm/s | %s | %s', ...
        ratName, K, strrep(binMode,'_','-'), vMin, upper(R.params.Mode), upper(labShort)), ...
        'FontWeight','bold');

% ---- save if requested ----
if doSave
  if strlength(figName)==0
    figName = sprintf('%s_linearSpatial_K%d_%s_vmin%.0f_%s_%s.png', ...
              ratName, K, binMode, vMin, modeStr, labShort);
  end
  exportgraphics(fig, figName, 'Resolution', 300);
end

% summaries
R.summary.meanOff_perDay  = meanOff;
R.summary.meanAng_perDay  = meanAng;      % NaN for corr
R.summary.meanOff_pooled  = meanOff_pool;
R.summary.meanAng_pooled  = meanAng_pool; % NaN for corr
end

% =======================================================================
function R = populationLinearSpatialOrthogonality_byRat(ratName, varargin)
% Running-only spatial population code across K linear bins (PC1-projected),
% with LOCKED PC1 and CONSISTENT bin edges across days.
% Excludes [CS, CS+2] on the velocity timebase; only keeps v >= MinSpeed.

% ---------- args ----------
p = inputParser;
addParameter(p,'MinSpeed',4,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'NBins',6,@(x) isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'BinMode','equal_size',@(s) any(strcmpi(s,{'equal_occ','equal_size'})));
addParameter(p,'LockAxis','x',@(s) any(strcmpi(s,{'x','y'})));
addParameter(p,'Mode','raw',@(s) any(strcmpi(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','cosine',@(s) any(strcmpi(s,{'cosine','corr','pearson'})));
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty(d));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x))); % NEW
addParameter(p,'DayWindowMode','fraction',@(s) any(strcmpi(s,{'fraction','seconds'})));
addParameter(p,'DayWindow',[],@(v) isempty(v) || (isnumeric(v)&&numel(v)==2&&v(2)>v(1)));
parse(p,varargin{:});

vMin        = p.Results.MinSpeed;
K           = p.Results.NBins;
binMode     = lower(p.Results.BinMode);
lockAx      = lower(p.Results.LockAxis);
modeStr     = lower(p.Results.Mode);
simStr      = lower(p.Results.Similarity);
daysArg     = p.Results.Days;
miCutoff    = p.Results.MICutoff;
miExCutoff  = p.Results.MIExcludeCutoff;   % NEW
dayWinMode  = lower(p.Results.DayWindowMode);
dayWin      = p.Results.DayWindow;

% ---------- resolve days ----------
rat   = evalin('base', ratName);
dates = autoDateList(rat);
iAn   = find(strcmp(dates, rat.An), 1);

if ~isempty(daysArg)
    userDays = cellstr(daysArg);
    days     = cellfun(@(s) resolveDayToken(rat, s), userDays, 'UniformOutput', false);
    days     = unique(days, 'stable');
else
    if isempty(iAn), error('An date not found for %s.', ratName); end
    days = dates(max(1,iAn-2):iAn);
end
days = unique(days, 'stable');
fprintf('Using days: %s\n', strjoin(days, ', '));

% ---------- PASS 1: gather ALL running, NON-TRIAL samples ----------
all_x = []; all_y = [];
perDayRun = struct('vt',[],'xv',[],'yv',[],'dt',[],'csTimes',[],'winUsed',[]);
goodDay = false(1,numel(days));

for d = 1:numel(days)
    D    = days{d};
    Pfld = sprintf('pos_%s', D);
    if ~isfield(rat.pos, Pfld), continue; end

    pos = rat.pos.(Pfld);                    % [t x y]
    vel  = ca_velocity(pos);                 % [speed; time]
    vt   = vel(2,:)';  vmag = vel(1,:)';
    xv   = interp1(pos(:,1), pos(:,2), vt, 'linear', NaN);
    yv   = interp1(pos(:,1), pos(:,3), vt, 'linear', NaN);

    Cfld = sprintf('CS_%s', D);
    if isfield(rat.CS_times,Cfld), csTimes = rat.CS_times.(Cfld); else, csTimes = []; end
    inTrialVel = false(size(vt));
    for iCS = 1:numel(csTimes)
        inTrialVel = inTrialVel | (vt >= csTimes(iCS) & vt < csTimes(iCS)+2);
    end

    keep = (vmag >= vMin) & ~inTrialVel & isfinite(xv) & isfinite(yv);
    vt = vt(keep); xv = xv(keep); yv = yv(keep);

    % optional day window
    if ~isempty(dayWin)
        tStart = vt(1); tEnd = vt(end);
        switch dayWinMode
            case 'fraction'
                f0 = max(0, min(1, dayWin(1)));
                f1 = max(0, min(1, dayWin(2)));
                w0 = tStart + f0*(tEnd - tStart);
                w1 = tStart + f1*(tEnd - tStart);
            case 'seconds'
                w0 = tStart + dayWin(1);
                w1 = tStart + dayWin(2);
        end
        if w1 <= w0, warning('DayWindow invalid for %s; skipping day.', D); continue; end
        inWin = (vt >= w0) & (vt < w1);
        vt = vt(inWin); xv = xv(inWin); yv = yv(inWin);
    else
        w0 = vt(1); w1 = vt(end);
    end

    if numel(vt) < 10, continue; end

    dt = diff(vt); dt = [dt; median(dt,'omitnan')];

    perDayRun(d).vt      = vt;
    perDayRun(d).xv      = xv;
    perDayRun(d).yv      = yv;
    perDayRun(d).dt      = dt;
    perDayRun(d).csTimes = csTimes;
    perDayRun(d).winUsed = [w0 w1];
    goodDay(d) = true;

    all_x = [all_x; xv]; %#ok<AGROW>
    all_y = [all_y; yv]; %#ok<AGROW>
end
if ~any(goodDay), error('No days with sufficient running/non-trial data.'); end

% ---------- global PC1 and sign lock ----------
Xall = [all_x, all_y];
Xall = Xall - mean(Xall,1,'omitnan');
[coeff, score] = pca(Xall, 'NumComponents',1, 'Centered',false);
s_all = score(:,1);
if strcmp(lockAx,'x'), rho = corr(s_all, all_x, 'rows','complete');
else,                  rho = corr(s_all, all_y, 'rows','complete'); end
if rho < 0, coeff(:,1) = -coeff(:,1); s_all = -s_all; end

% ---------- shared edges along global PC1 ----------
switch binMode
    case 'equal_occ',  edges = quantile(s_all, linspace(0,1,K+1));
    case 'equal_size', edges = linspace(min(s_all), max(s_all), K+1);
end
edges = unique(edges); if numel(edges) < K+1, edges = linspace(min(s_all), max(s_all), K+1); end

% ---------- PASS 2: per-day rates ----------
perDay = struct('rateKxN_raw',[],'cosSimK',[],'angleK',[],'dotK',[],'edges',edges);
goodIdxList = find(goodDay);
cosStack = NaN(K,K,numel(goodIdxList));

for gi = 1:numel(goodIdxList)
    dOrig = goodIdxList(gi);
    D     = days{dOrig};
    Sfld  = sprintf('CA_peaks_%s', D);
    Mfld  = sprintf('ratemask_%s', D);

    if ~isfield(rat.Ca_peaks,Sfld), warning('Missing %s, skip.', Sfld); continue; end
    S = rat.Ca_peaks.(Sfld);

    % base mask from ratemask (or all-true if missing)
    if ~isfield(rat.ratemask,Mfld)
        if iscell(S), nCells = numel(S); else, nCells = size(S,1); end
        goodMask = true(nCells,1);
    else
        goodMask = rat.ratemask.(Mfld) == 1;
    end

    % ---- MI include (from MI_noCSUS15_shuff) ----
    if ~isempty(miCutoff)
        MIfld  = 'MI_noCSUS15_shuff';
        MIroot = sprintf('MI_%s', D);
        if isfield(evalin('base',ratName), MIfld)
            MIstruct = evalin('base',ratName);
            if isfield(MIstruct.(MIfld), MIroot)
                MIvals = MIstruct.(MIfld).(MIroot);
                keepMI = MIvals(:,3) >= miCutoff;
                goodMask = goodMask(:) & logical(padToLength(keepMI, numel(goodMask)));
            end
        end
    end

    % ---- MI EXCLUDE (from MI_CSUS15_shuff)  NEW ----
    if ~isempty(miExCutoff)
        EXfld  = 'MI_CSUS15_shuff';
        MIroot = sprintf('MI_%s', D);
        if isfield(evalin('base',ratName), EXfld)
            MIstruct = evalin('base',ratName);
            if isfield(MIstruct.(EXfld), MIroot)
                MIvals = MIstruct.(EXfld).(MIroot);
                exMask = MIvals(:,3) >= miExCutoff;           % mark for exclusion
                goodMask = goodMask(:) & ~logical(padToLength(exMask, numel(goodMask)));
            end
        end
    end

    goodIdx = find(goodMask(:));
    if isempty(goodIdx), continue; end

    vt      = perDayRun(dOrig).vt;
    xv      = perDayRun(dOrig).xv;
    yv      = perDayRun(dOrig).yv;
    dt      = perDayRun(dOrig).dt;
    csTimes = perDayRun(dOrig).csTimes;
    winUsed = perDayRun(dOrig).winUsed;

    % project running onto global PC1
    XY  = [xv, yv]; XYc = XY - mean([all_x, all_y],1,'omitnan');
    s   = XYc * coeff(:,1);

    % occupancy per bin
    binIdx = discretize(s, edges);
    valid  = ~isnan(binIdx);
    occT   = accumarray(binIdx(valid), dt(valid), [K,1], @sum, NaN);
    occT(occT<=0) = NaN;

    % rates per neuron
    N = numel(goodIdx);
    rateKxN_raw = nan(K, N);
    for j = 1:N
        c  = goodIdx(j);
        st = getCellSpikes_local(S, c);
        if numel(st) < 3, continue; end
        % drop trial spikes
        inTrialSpk = false(size(st));
        for iCS = 1:numel(csTimes)
            inTrialSpk = inTrialSpk | (st >= csTimes(iCS) & st < csTimes(iCS)+2);
        end
        st = st(~inTrialSpk);
        % restrict to window
        st = st(st >= winUsed(1) & st < winUsed(2));
        if numel(st) < 3, continue; end

        % map spikes and bin
        xs = interp1(vt, xv, st, 'linear', NaN);
        ys = interp1(vt, yv, st, 'linear', NaN);
        XYs = [xs, ys]; XYs = XYs(all(isfinite(XYs),2),:);
        if isempty(XYs), continue; end

        XYsC = XYs - mean([all_x, all_y],1,'omitnan');
        ss   = XYsC * coeff(:,1);
        bSpk = discretize(ss, edges);

        for b = 1:K
            nSpk = sum(bSpk == b);
            rateKxN_raw(b,j) = nSpk ./ occT(b);
        end
    end

    % transform + similarity
    rateForSim = normalizeRatesPerCell(rateKxN_raw, modeStr);
    [dotK, S, angK] = computeSimilarityMatrix(rateForSim, simStr);

    perDay(gi).rateKxN_raw = rateKxN_raw;
    perDay(gi).cosSimK     = S;
    perDay(gi).angleK      = angK;
    perDay(gi).dotK        = dotK;
    perDay(gi).edges       = edges;

    % stack once per valid day
    offMask = ~eye(size(S,1));
    if any(isfinite(S(offMask)), 'all')
        cosStack(:,:,gi) = S;
    end
end

% ---------- pooled ----------
if ~isempty(cosStack)
    validSlices = squeeze(~all(all(isnan(cosStack),1),2));
    cosStack = cosStack(:,:,validSlices);
    if isempty(cosStack)
        pooledCos = nan(K,K); pooledAng = nan(K,K);
    elseif size(cosStack,3) == 1
        pooledCos = cosStack(:,:,1);
        if strcmpi(simStr,'cosine')
            pooledAng = real(acos(max(-1,min(1,pooledCos))));
        else
            pooledAng = nan(size(pooledCos));
        end
    else
        pooledCos = nanmean(cosStack,3);
        if strcmpi(simStr,'cosine')
            pooledAng = real(acos(max(-1,min(1,pooledCos))));
        else
            pooledAng = nan(size(pooledCos));
        end
    end
else
    pooledCos = nan(K,K); pooledAng = nan(K,K);
end

daysKept = days(goodIdxList);
R = struct();
R.ratVar = ratName;
R.days   = daysKept;
R.perDay = perDay;
R.pooled = struct('cosSimK', pooledCos, 'angleK', pooledAng, 'edges', edges);
R.params = struct('MinSpeed', vMin, 'NBins', K, 'BinMode', binMode, ...
                  'LockAxis', lockAx, 'Mode', modeStr, 'Similarity', simStr, ...
                  'MICutoff', miCutoff, 'MIExcludeCutoff', miExCutoff, 'Days', {days});
end

% ================= helpers =================
function X = normalizeRatesPerCell(R, modeStr)
switch lower(modeStr)
    case 'raw'
        X = R;
    case 'demean'
        mu = mean(R, 1, 'omitnan');
        X  = R - mu;
    case 'zscore'
        mu = mean(R, 1, 'omitnan');
        sd = std(R, 0, 1, 'omitnan');
        sd(~isfinite(sd) | sd==0) = 1;
        X  = (R - mu) ./ sd;
    otherwise
        error('Unknown Mode: %s', modeStr);
end
end

function [dotM, S, angM] = computeSimilarityMatrix(X, simStr)
dotM = X * X.';                         % K×K
switch lower(simStr)
    case {'cosine'}
        rowNorms = vecnorm(X,2,2);
        den      = rowNorms * rowNorms.';
        S        = dotM ./ max(den, eps);
        S        = max(-1, min(1, S));
        angM     = real(acos(S));
    case {'corr','pearson'}
        validPerRow = sum(isfinite(X), 2) >= 2;
        if ~any(validPerRow)
            S    = nan(size(X,1));
            angM = nan(size(S));
            return
        end
        C = corrcoef(X', 'Rows', 'pairwise');   % corr of rows across cells
        S = C(1:size(X,1), 1:size(X,1));
        angM = nan(size(S));
    otherwise
        error('Unknown Similarity: %s', simStr);
end
end

function st = getCellSpikes_local(S, c)
  if iscell(S), st = S{c}(:); else, st = S(c,:).'; end
  st = st(~isnan(st) & st>0);
end

function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end

function v = padToLength(v, L)
v = v(:);
if numel(v) < L, v(end+1:L) = false; end
if numel(v) > L, v = v(1:L); end
end

function tok = resolveDayToken(rat, s)
s = strrep(strrep(s,'-','_'),'/','_');
m = regexp(s, '^(\d{4})_(\d{1,2})_(\d{1,2})$', 'tokens', 'once');
if isempty(m), tok = s; return; end
yy = m{1}; mm = str2double(m{2}); dd = str2double(m{3});
padded   = sprintf('%s_%02d_%02d', yy, mm, dd);
unpadded = sprintf('%s_%d_%d',    yy, mm, dd);
forms = {padded, unpadded};
for k = 1:numel(forms)
    cand = forms{k};
    if (isfield(rat.pos,      ['pos_'      cand])) || ...
       (isfield(rat.Ca_peaks, ['CA_peaks_' cand])) || ...
       (isfield(rat.CS_times, ['CS_'       cand])) || ...
       (isfield(rat.ratemask, ['ratemask_' cand]))
        tok = cand; return
    end
end
tok = padded;
end
