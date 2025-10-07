function R = populationSpatialOrthogonality(ratName, varargin)
% populationSpatialOrthogonality
% (top-level): multi-rat pooled page, then returns struct.
% Single-rat path delegates to populationSpatialOrthogonality_singleRat.

% ---------------- parse shared args ----------------
p = inputParser;
addParameter(p,'NBins',16, @(x) isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'GridRC',[],@(v) (isempty(v) || (isnumeric(v) && numel(v)==2 && all(v>=1))));
addParameter(p,'BinMode','equal_size',@(s) any(strcmpi(s,{'equal_occ','equal_size'})));
addParameter(p,'MinSpeed',4,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Mode','raw',@(s) any(strcmpi(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','corr',@(s) any(strcmpi(s,{'cosine','corr','pearson'})));
addParameter(p,'Adjacency','rook',@(s) any(strcmpi(s,{'moore','rook'})));
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty(d));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'DayWindowMode','fraction',@(s) any(strcmpi(s,{'fraction','seconds'})));
addParameter(p,'DayWindow',[],@(v) isempty(v) || (isnumeric(v)&&numel(v)==2&&v(2)>v(1)));
addParameter(p,'SaveFig',false,@islogical);
addParameter(p,'FigName','',@(s) ischar(s) || isstring(s));
addParameter(p,'ShowPerRatHeatmaps',true,@islogical);  % (kept for compatibility)
addParameter(p,'PerRatFigures',true,@islogical);  % NEW: make a page per rat

parse(p,varargin{:});

K            = p.Results.NBins;
GridRC       = p.Results.GridRC;
binMode      = lower(p.Results.BinMode);
vMin         = p.Results.MinSpeed;
modeStr      = lower(p.Results.Mode);
simStr       = lower(p.Results.Similarity);
adjStr       = lower(p.Results.Adjacency);
daysArg      = p.Results.Days;
miCutoff     = p.Results.MICutoff;
miExCutoff   = p.Results.MIExcludeCutoff;
dayWinMode   = lower(p.Results.DayWindowMode);
dayWin       = p.Results.DayWindow;
doSave       = p.Results.SaveFig;
figName      = string(p.Results.FigName);
doPerRat    = p.Results.PerRatFigures;  % NEW


isCorr   = any(strcmpi(simStr, {'corr','pearson'}));
labShort = tern(isCorr,"corr","cos");
clims    = tern(isCorr,[-1 1],[0 1]);

if ~isempty(GridRC)
    if prod(GridRC) ~= K
        error('GridRC (%d x %d) must multiply to NBins=%d.', GridRC(1), GridRC(2), K);
    end
end

% ===================== MODE SWITCH =====================
if ~iscell(ratName)
    % -------- SINGLE-RAT (plots the 2x4 per-animal figure) --------
    R = populationSpatialOrthogonality_singleRat( ...
            ratName, 'NBins',K,'GridRC',GridRC,'BinMode',binMode,'MinSpeed',vMin, ...
            'Mode',modeStr,'Similarity',simStr,'Adjacency',adjStr,'Days',daysArg, ...
            'MICutoff',miCutoff,'MIExcludeCutoff',miExCutoff, ...
            'DayWindowMode',dayWinMode,'DayWindow',dayWin, ...
            'SaveFig',doSave,'FigName',figName);
    return
end

% -------- MULTI-RAT MODE --------
ratNames = ratName(:)';
nRats    = numel(ratNames);

% compute per-rat pooled results
Rrats = cell(1,nRats);
for r = 1:nRats
    Rrats{r} = populationLinearSpatialOrthogonality_byRat( ...
                    ratNames{r}, 'NBins',K,'GridRC',GridRC,'BinMode',binMode,'MinSpeed',vMin, ...
                    'Mode',modeStr,'Similarity',simStr,'Days',daysArg, ...
                    'MICutoff',miCutoff,'MIExcludeCutoff',miExCutoff, ...
                    'DayWindowMode',dayWinMode,'DayWindow',dayWin);
end

% stack per-rat pooled similarity matrices
cosStack_rats = NaN(K,K,nRats);
haveMeta      = false(1,nRats);
metaList      = cell(1,nRats);
for r = 1:nRats
    C = Rrats{r}.pooled.cosSimK;
    if ~isempty(C)
        cosStack_rats(:,:,r) = C;
        if isfield(Rrats{r}.pooled,'binMeta')
            metaList{r} = Rrats{r}.pooled.binMeta;
            haveMeta(r) = true;
        end
    end
end

validSlices    = squeeze(~all(all(isnan(cosStack_rats),1),2));
cosStack_rats  = cosStack_rats(:,:,validSlices);
Cgrand         = tern(isempty(cosStack_rats), nan(K,K), nanmean(cosStack_rats,3));

% choose a representative binMeta (first available)
grandMeta = [];
ixMeta = find(haveMeta,1,'first');
if ~isempty(ixMeta), grandMeta = metaList{ixMeta}; end

% summary stats
offMask        = ~eye(K);
meanOff_grand  = mean(Cgrand(offMask),'omitnan');
meanAng_grand  = tern(~isCorr, real(acos(max(-1,min(1,meanOff_grand)))) , NaN);

% --------- NEW: per-rat 2x4 figures (optional) ---------
if p.Results.PerRatFigures
    for r = 1:nRats
        populationSpatialOrthogonality_singleRat( ...
            ratNames{r}, ...
            'NBins',K, 'GridRC',GridRC, 'BinMode',binMode, 'MinSpeed',vMin, ...
            'Mode',modeStr, 'Similarity',simStr, 'Adjacency',adjStr, ...
            'Days',daysArg, 'MICutoff',miCutoff, 'MIExcludeCutoff',miExCutoff, ...
            'DayWindowMode',dayWinMode, 'DayWindow',dayWin, ...
            'SaveFig',doSave, 'FigName',"");
        drawnow;
    end
end

% GRAND spatial stats
grandStats = [];
if ~isempty(grandMeta)
    grandStats = computeSpatialGridStats(Cgrand, grandMeta, ...
                    'Adjacency', adjStr, 'MantelPerms', 500, ...
                    'IsCorr', isCorr);
end

% Per-rat distance curves for bottom panel
ratDist = struct('centers',[],'mean',[]);
ratHas  = false(1,nRats);
for r = 1:nRats
    if isfield(Rrats{r}.pooled,'binMeta') && ~isempty(Rrats{r}.pooled.binMeta) ...
            && ~isempty(Rrats{r}.pooled.cosSimK)
        sr = computeSpatialGridStats(Rrats{r}.pooled.cosSimK, Rrats{r}.pooled.binMeta, ...
                'Adjacency', adjStr, 'MantelPerms', 500, 'IsCorr', isCorr);
        ratDist(r).centers = sr.dist_centers;
        ratDist(r).mean    = sr.dist_mean;
        ratHas(r) = any(isfinite(sr.dist_mean));
    end
end

% ------------------- POOLED (ALL RATS) FIGURE -------------------
figAll = figure('Color','w','Position',[80 80 1200 780]);
tiledlayout(figAll,2,3,'TileSpacing','compact','Padding','compact');


% (1,1) GRAND pooled heatmap
nexttile(1);
imagesc(Cgrand, clims); axis square; colormap(gca, parula); colorbar;
xticks(1:K); yticks(1:K); xlabel('Bin'); ylabel('Bin');
if ~isempty(grandStats)
    if isCorr
        title(sprintf('GRAND | mean off %s=%.2f | adj(%s)=%.2f | Mantel r=%.2f (p=%.3g)', ...
              labShort, meanOff_grand, adjStr, grandStats.mean_adj, grandStats.mantel_r, grandStats.mantel_p));
    else
        title(sprintf('GRAND | mean off %s=%.2f (%.0f°) | adj(%s)=%.2f | Mantel r=%.2f (p=%.3g)', ...
              labShort, meanOff_grand, rad2deg(meanAng_grand), adjStr, grandStats.mean_adj, grandStats.mantel_r, grandStats.mantel_p));
    end
else
    title(sprintf('GRAND | mean off-diag %s=%.2f', labShort, meanOff_grand));
end

% (1,2:3) Adjacent vs Non-adjacent bars
nexttile([1 2]); axBars = gca; cla(axBars);
[adjVals, nonadjVals] = deal(nan(1,nRats));
for r = 1:nRats
    C = Rrats{r}.pooled.cosSimK;
    if isempty(C), continue; end
    bm = Rrats{r}.pooled.binMeta;
    st = computeSpatialGridStats(C, bm, 'Adjacency', adjStr, 'MantelPerms', 500, 'IsCorr', isCorr);
    adjVals(r)    = st.mean_adj;
    nonadjVals(r) = st.mean_nonadj;
end
ok = isfinite(adjVals) & isfinite(nonadjVals);
adjVals = adjVals(ok); nonadjVals = nonadjVals(ok);
nUse = numel(adjVals);
mAdj = mean(adjVals); mNon = mean(nonadjVals);
seAdj = std(adjVals,0,2)/sqrt(max(1,nUse));
seNon = std(nonadjVals,0,2)/sqrt(max(1,nUse));
bar(axBars,[1 2],[mAdj mNon],0.65,'FaceColor',[0.80 0.85 1.00],'EdgeColor','k'); hold on
errorbar(axBars,[1 2],[mAdj mNon],[seAdj seNon],'k','LineStyle','none','LineWidth',1.3)
for i = 1:nUse, plot(axBars,[1 2],[adjVals(i) nonadjVals(i)],'-o','LineWidth',1.4); end
xticks(axBars,[1 2]); xticklabels(axBars,{'Adjacent','Non-adjacent'});
ylabel(axBars, tern(isCorr,'Mean Pearson similarity','Mean cosine similarity'));
title(axBars, sprintf('Adjacent > Non-adjacent  (n=%d rats)', nUse));
box(axBars,'on');

% (2,1:3) Lag curves (per rat + grand)
nexttile([1 3]); axLag = gca; hold(axLag,'on');
if ~isempty(grandStats) && any(isfinite(grandStats.dist_mean))
    plot(axLag, grandStats.dist_centers, grandStats.dist_mean, '-', 'LineWidth', 2.2);
end
for r = 1:nRats
    if ratHas(r)
        if ~isempty(grandStats) && any(grandStats.dist_centers)
            yy = interp1(ratDist(r).centers, ratDist(r).mean, grandStats.dist_centers, 'linear','extrap');
            plot(axLag, grandStats.dist_centers, yy, '-', 'LineWidth', 0.9);
        else
            plot(axLag, ratDist(r).centers, ratDist(r).mean, '-', 'LineWidth', 0.9);
        end
    end
end
xlabel(axLag,'Bin center distance'); ylabel(axLag,sprintf('Mean %s similarity', labShort));
title(axLag,'Similarity vs physical distance (per rat + grand)'); box(axLag,'off');

% Pretty title
if ~isempty(GridRC)
    sgtitle(figAll, sprintf('ALL RATS | %dx%d grid (%d bins) | %s | v>=%.0f | %s | adj=%s', ...
        GridRC(1),GridRC(2),K, strrep(binMode,'_','-'), vMin, upper(modeStr), upper(adjStr)), 'FontWeight','bold');
else
    sgtitle(figAll, sprintf('ALL RATS | K=%d (PC1) | %s | v>=%.0f | %s', ...
        K, strrep(binMode,'_','-'), vMin, upper(modeStr)), 'FontWeight','bold');
end

if doSave
    if strlength(figName)==0
        if ~isempty(GridRC)
            figName = sprintf('ALLRATS_grid_%dx%d_K%d_%s_vmin%.0f_%s_%s.png', ...
                        GridRC(1),GridRC(2),K, binMode, vMin, modeStr, labShort);
        else
            figName = sprintf('ALLRATS_linearSpatial_K%d_%s_vmin%.0f_%s_%s.png', ...
                        K, binMode, vMin, modeStr, labShort);
        end
    end
    exportgraphics(figAll, figName, 'Resolution', 300);
end

% -------- return struct --------
R = struct();
R.mode = 'multi-rat';
R.params = struct('NBins',K,'GridRC',GridRC,'BinMode',binMode,'MinSpeed',vMin,'Mode',modeStr, ...
                  'Similarity',simStr,'Adjacency',adjStr,'Days',{daysArg});
R.byRat = Rrats;
R.pooledAcrossRats.cosSimK        = Cgrand;
R.pooledAcrossRats.meanOff        = meanOff_grand;
R.pooledAcrossRats.meanAng        = meanAng_grand;
R.pooledAcrossRats.binMeta        = grandMeta;
R.pooledAcrossRats.spatialStats   = grandStats;
end


% ===================== SINGLE-RAT VIEW (grid-aware) =====================
function R = populationSpatialOrthogonality_singleRat(ratName, varargin)
% Same computations as before, but now emits a 2x4 figure per-animal:
% top row: day1, day2, day3, pooled; bottom row: bars (span 2), lag curve (span 2)

p = inputParser;
addParameter(p,'NBins',16, @(x) isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'GridRC',[],@(v) (isempty(v) || (isnumeric(v)&&numel(v)==2&&all(v>=1))));
addParameter(p,'BinMode','equal_occ',@(s) any(strcmpi(s,{'equal_occ','equal_size'})));
addParameter(p,'MinSpeed',4,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Mode','raw',@(s) any(strcmpi(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','cosine',@(s) any(strcmpi(s,{'cosine','corr','pearson'})));
addParameter(p,'Adjacency','moore',@(s) any(strcmpi(s,{'moore','rook'})));
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty(d));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'DayWindowMode','fraction',@(s) any(strcmpi(s,{'fraction','seconds'})));
addParameter(p,'DayWindow',[],@(v) isempty(v) || (isnumeric(v)&&numel(v)==2&&v(2)>v(1)));
addParameter(p,'SaveFig',false,@islogical);
addParameter(p,'FigName','',@(s) ischar(s) || isstring(s));
parse(p,varargin{:});

K            = p.Results.NBins;
GridRC       = p.Results.GridRC;
binMode      = lower(p.Results.BinMode);
vMin         = p.Results.MinSpeed;
modeStr      = lower(p.Results.Mode);
simStr       = lower(p.Results.Similarity);
adjStr       = lower(p.Results.Adjacency);
daysArg      = p.Results.Days;
miCutoff     = p.Results.MICutoff;
miExCutoff   = p.Results.MIExcludeCutoff;
doSave       = p.Results.SaveFig;
figName      = string(p.Results.FigName);
dayWinMode   = lower(p.Results.DayWindowMode);
dayWin       = p.Results.DayWindow;

% ---- compute (helper) ----
R = populationLinearSpatialOrthogonality_byRat( ...
        ratName, 'NBins',K, 'GridRC',GridRC, 'BinMode',binMode, 'MinSpeed',vMin, ...
        'Mode',modeStr, 'Similarity',simStr, 'Days',daysArg, ...
        'MICutoff',miCutoff, 'MIExcludeCutoff',miExCutoff, ...
        'DayWindowMode',dayWinMode,'DayWindow',dayWin);

% ---- collect per-day panels ----
K = size(R.pooled.cosSimK,1);
nDays = numel(R.perDay);
have  = false(1,nDays);
for d = 1:nDays
    have(d) = (d<=numel(R.perDay)) && ~isempty(R.perDay(d).cosSimK);
end

Cpool  = R.pooled.cosSimK;
offMask= ~eye(K);
isCorr = any(strcmpi(R.params.Similarity, {'corr','pearson'}));
labShort = tern(isCorr,"corr","cos");
clims    = tern(isCorr,[-1 1],[0 1]);

meanOff_pool = mean(Cpool(offMask),'omitnan');
meanAng_pool = tern(~isCorr, real(acos(max(-1,min(1,meanOff_pool)))) , NaN);

% pooled spatial stats (used for bars + lag)
poolStats = [];
if isfield(R.pooled,'binMeta')
    poolStats = computeSpatialGridStats(Cpool, R.pooled.binMeta, ...
                 'Adjacency', adjStr, 'MantelPerms', 500, 'IsCorr', isCorr);
end

% ---------- observed distance-curve stats (no perms) ----------
distStats = [];
if isfield(R.pooled,'binMeta') && ~isempty(R.pooled.binMeta) && ~isempty(Cpool)
    distStats = computeSpatialGridStats(Cpool, R.pooled.binMeta, ...
                   'Adjacency', adjStr, 'MantelPerms', 0, 'IsCorr', isCorr);
end

% ---------- shuffled null for distance curve (NaN-safe & edge-consistent) ----------
nullLo = []; nullHi = []; nullMean = []; goodCols = [];

if ~isempty(distStats) && isfield(distStats,'dist_edges') && ~isempty(distStats.dist_edges)
    edges    = distStats.dist_edges;
    dCenters = distStats.dist_centers;
    Kloc     = size(Cpool,1);
    UT       = triu(true(Kloc),1);

    % Build UT masks using THE SAME edges as the observed curve
    % (works for grid or PC1)
    Dfull = nan(Kloc);
    if isfield(R.pooled.binMeta,'mode') && strcmpi(R.pooled.binMeta.mode,'grid') && ...
       isfield(R.pooled.binMeta,'centers2D') && size(R.pooled.binMeta.centers2D,2)==2
        XY = R.pooled.binMeta.centers2D;
        [I,J] = find(UT);
        dvec = sqrt(sum((XY(I,:) - XY(J,:)).^2,2));
    else
        % 1D fallback
        if isfield(R.pooled.binMeta,'centers1D')
            c = R.pooled.binMeta.centers1D(:);
        else
            c = (1:Kloc)'; % last resort
        end
        [I,J] = find(UT);
        dvec = abs(c(I) - c(J));
    end

    dMasks = cell(1, numel(edges)-1);
    for k = 1:numel(edges)-1
        dMasks{k} = (dvec >= edges(k) & dvec < edges(k+1));
    end

    % pooled split×cell matrix (normalized same as observed)
    Mpool = [];
    for d = 1:numel(R.perDay)
        rr = R.perDay(d).rateKxN_raw;
        if ~isempty(rr)
            Mpool = [Mpool, normalizeRatesPerCell(rr, R.params.Mode)]; %#ok<AGROW>
        end
    end

    nCells = size(Mpool,2);
    nSh    = 500;                               % 200–1000 typical
    useCos = any(strcmpi(R.params.Similarity, {'cosine'}));

    if Kloc >= 2 && nCells >= 3 && ~isempty(dMasks)
        % gridRC hint for 2D shuffle (falls back to 1D inside if mismatched)
        gridRC = [];
        if isfield(R.pooled,'binMeta') && isfield(R.pooled.binMeta,'gridRC') && ~isempty(R.pooled.binMeta.gridRC)
            gridRC = R.pooled.binMeta.gridRC;
        end

        lagNull = nan(nSh, numel(dCenters));
        for b = 1:nSh
          Cb = shuffle_full_per_cell(Mpool, useCos);
            yb = Cb(UT);
            for k = 1:numel(dCenters)
                mk = dMasks{k};
                if any(mk)
                    lagNull(b,k) = mean(yb(mk), 'omitnan');
                end
            end
        end

        % Keep only columns that had any finite null values (some bins can be empty)
        goodCols = find(any(isfinite(lagNull),1));
        if ~isempty(goodCols)
            lagNull  = lagNull(:, goodCols);
            nullLo   = prctile(lagNull,  2.5, 1);
            nullHi   = prctile(lagNull, 97.5, 1);
            nullMean = mean(lagNull, 1, 'omitnan');
        end
    end
end

% ---- (optional) light bootstrap CI for observed curve (off by default) ----
distLag = struct('centers',[],'mean',[],'ci_low',[],'ci_high',[],'rho',NaN,'p_rho',NaN);
if isfield(R.pooled,'binMeta') && ~isempty(R.pooled.binMeta)
    [~, edgesB, ctrsB, masksB] = buildDistanceBins(R.pooled.binMeta, size(Cpool,1), 8);
    UTloc = triu(true(size(Cpool,1)),1);
    yObs  = Cpool(UTloc);
    means = nan(1,numel(masksB));
    for k = 1:numel(masksB)
        idxk = masksB{k}(UTloc); means(k) = mean(yObs(idxk),'omitnan');
    end
    distLag.centers = ctrsB; distLag.mean = means;
    [distLag.rho, distLag.p_rho] = corr(ctrsB(:), means(:), 'Type','Spearman','Rows','complete');

    % (disabled by default; set doCI=true to enable)
    doCI = false;
    if doCI
        Mpool = [];
        for d = 1:numel(R.perDay)
            rr = R.perDay(d).rateKxN_raw;
            if isempty(rr), continue; end
            Xr = normalizeRatesPerCell(rr, R.params.Mode);
            Mpool = [Mpool, Xr]; %#ok<AGROW>
        end
        nCells = size(Mpool,2);
        if nCells >= 5
            nBoot = 200;
            lagBoot = nan(nBoot, numel(masksB));
            useCosB  = any(strcmpi(R.params.Similarity, {'cosine'}));
            for b = 1:nBoot
                cols = randsample(nCells, nCells, true);
                Cb   = pvSimFromRates_local(Mpool(:, cols), useCosB);
                yb   = Cb(UTloc);
                for k = 1:numel(masksB)
                    idxk = masksB{k}(UTloc);
                    lagBoot(b,k) = mean(yb(idxk),'omitnan');
                end
            end
            distLag.ci_low  = prctile(lagBoot, 2.5, 1);
            distLag.ci_high = prctile(lagBoot,97.5, 1);
        end
    end
end

% ======== FIXED 2×4 LAYOUT (PER-ANIMAL FIGURE) ========
dayIdx = find(have);
if numel(dayIdx) > 3, dayIdx = dayIdx(1:3); end

figPer = figure('Color','w','Position',[100 100 1400 700]);
tiledlayout(figPer,2,4,'TileSpacing','compact','Padding','compact');

% Top row: up to 3 day heatmaps
for j = 1:3
    nexttile(j); cla
    if j <= numel(dayIdx)
        d = dayIdx(j);
        C = R.perDay(d).cosSimK;
        imagesc(C, clims); axis square; colormap(gca, parula); colorbar;
        xticks(1:K); yticks(1:K); xlabel('Bin'); ylabel('Bin');
        mOff = mean(C(offMask),'omitnan');
        if isCorr
            title(sprintf('%s | mean off %s = %.2f', R.days{d}, labShort, mOff));
        else
            title(sprintf('%s | mean off %s = %.2f (%.0f°)', R.days{d}, labShort, mOff, ...
                 rad2deg(real(acos(max(-1,min(1,mOff)))))));
        end
    else
        axis off
        text(0.5,0.5,'(no day)', 'HorizontalAlignment','center', ...
             'VerticalAlignment','middle','FontSize',12);
    end
end

% Top row: pooled heatmap
nexttile(4); cla
imagesc(Cpool, clims); axis square; colormap(gca, parula); colorbar;
xticks(1:K); yticks(1:K); xlabel('Bin'); ylabel('Bin');
if ~isempty(poolStats)
    if isCorr
        title(sprintf('Pooled | mean off %s=%.2f | adj(%s)=%.2f', ...
              labShort, meanOff_pool, adjStr, poolStats.mean_adj));
    else
        title(sprintf('Pooled | mean off %s=%.2f', labShort, meanOff_pool));
    end
else
    if isCorr
        title(sprintf('Pooled | mean off %s=%.2f', labShort, meanOff_pool));
    else
        title(sprintf('Pooled | mean off %s=%.2f (%.0f°)', labShort, meanOff_pool, rad2deg(meanAng_pool)));
    end
end

% Bottom-left: Adjacent vs Non-adjacent bars (span 2 cols)
nexttile(5, [1 2]); cla; hold on
if ~isempty(poolStats) && ~isempty(Cpool)
    % --- build adjacency mask for this rat's geometry ---
    K  = size(Cpool,1);
    UT = triu(true(K),1);
    sUT = Cpool(UT);

    if isfield(R.pooled,'binMeta') && strcmpi(R.pooled.binMeta.mode,'grid')
        Rg = R.pooled.binMeta.gridRC(1); Cg = R.pooled.binMeta.gridRC(2);
        [ri,cj] = ind2sub([Rg Cg], (1:K)');
        [I,J] = find(UT);
        dCheb = max(abs(ri(I)-ri(J)), abs(cj(I)-cj(J)));
        dMan  = abs(ri(I)-ri(J)) + abs(cj(I)-cj(J));
        switch lower(adjStr)
            case 'moore',  adjMaskVec = (dCheb == 1);
            case 'rook',   adjMaskVec = (dMan  == 1);
        end
    else
        % 1D / PC1 fallback
        [I,J] = find(UT);
        dIdx = abs(I-J);
        adjMaskVec = (dIdx == 1);
    end

    sAdj    = sUT(adjMaskVec);    sAdj    = sAdj(isfinite(sAdj));
    sNonAdj = sUT(~adjMaskVec);   sNonAdj = sNonAdj(isfinite(sNonAdj));

    mAdj = mean(sAdj,'omitnan');
    mNon = mean(sNonAdj,'omitnan');

    % --- inferential stats (Welch t-test + Wilcoxon) ---
    p_t  = NaN; p_w = NaN;
    size(sNonAdj)
    size(sAdj)
    try
        [~,p_t] = ttest2(sAdj, sNonAdj);  % unequal variances
        p_w     = ranksum(sAdj, sNonAdj);
    catch
        % leave NaNs if any edge case
    end

    % --- plot bars + error bars (SEM of pair dists) ---
    sem = @(v) std(v,0,'omitnan')/sqrt(max(1,nnz(isfinite(v))));
    bh = bar([1 2],[mAdj mNon], 0.65, 'FaceColor',[0.80 0.85 1.00], 'EdgeColor','k'); %#ok<NASGU>
    errorbar([1 2],[mAdj mNon],[sem(sAdj) sem(sNonAdj)],'k','LineStyle','none','LineWidth',1.3)

    % scatter all pair values (jitter for visibility)
    rng(0);
    jadj  = 1 + 0.12*(rand(size(sAdj))-0.5);
    jnon  = 2 + 0.12*(rand(size(sNonAdj))-0.5);
    plot(jadj, sAdj,  '.', 'MarkerSize',10, 'Color',[0.10 0.35 0.80])
    plot(jnon, sNonAdj,'.', 'MarkerSize',10, 'Color',[0.55 0.15 0.15])

    xticks([1 2]); xticklabels({'Adjacent','Non-adjacent'});
    ylabel( tern(isCorr,'Pearson similarity','Cosine similarity') );
    box on

    yl = ylim;
    yb = yl(2) - 0.06*range(yl);
    plot([1 1 2 2],[yb-0.01 yb yb yb-0.01],'k-','LineWidth',1.25)
    txt = sprintf('%s (t p=%.2g; rank p=%.2g)', sigStars(p_t), p_t, p_w);
    text(1.5, yb+0.01*range(yl), txt, 'HorizontalAlignment','center', 'FontSize',12, 'FontWeight','bold');

    title('Adjacent vs Non-adjacent (pooled pairs within rat)');
    ylim([min([sAdj; sNonAdj],[],'omitnan')-0.05*range([sAdj; sNonAdj])  yl(2)*1.05]);
else
    axis off
    text(0.5,0.5,'No adjacency stats available','HorizontalAlignment','center','VerticalAlignment','middle');
end


% -------- Bottom-right: similarity vs physical distance (tiles 7–8) --------
nexttile(7, [1 2]); cla; hold on
hasStats = ~isempty(distStats) && isfield(distStats,'dist_centers') && isfield(distStats,'dist_mean');
if hasStats
    x = distStats.dist_centers(:);
    y = distStats.dist_mean(:);
    ok = isfinite(x) & isfinite(y);

    % --- plot shuffled null band + dashed mean (if available) ---
    if ~isempty(nullLo) && ~isempty(nullHi) && ~isempty(goodCols)
        xNull = x(goodCols);
        xx = [xNull; flipud(xNull)];
        yy = [nullLo(:); flipud(nullHi(:))];
        fill(xx, yy, [0.8 0.85 1], 'EdgeColor','none', 'FaceAlpha',0.45);
    end
    if ~isempty(nullMean) && ~isempty(goodCols)
        plot(x(goodCols), nullMean(:), '--', 'LineWidth', 1.4);
    end

    % --- observed curve ---
    if nnz(ok) >= 2
        plot(x(ok), y(ok), 'k-', 'LineWidth', 1.8);
    end

    xlabel('Bin center distance');
    ylabel(sprintf('Mean %s similarity', labShort));
    if isfield(distStats,'mantel_r') && isfinite(distStats.mantel_r)
        if isfield(distStats,'mantel_p') && isfinite(distStats.mantel_p)
            title(sprintf('Similarity vs distance'));
        else
            title(sprintf('Similarity vs distance'));
        end
    else
        title('Similarity vs distance');
    end
    box off

    % Optional legend
    lg = {}; lh = [];
    if ~isempty(nullLo) && ~isempty(nullHi) && ~isempty(goodCols)
        lg{end+1} = 'Shuffled 95% band';  %#ok<AGROW>
        lh(end+1) = plot(nan,nan,'-','Color',[0.8 0.85 1],'LineWidth',8); %#ok<AGROW>
    end
    if ~isempty(nullMean) && ~isempty(goodCols)
        lg{end+1} = 'Shuffled mean';      %#ok<AGROW>
        lh(end+1) = plot(nan,nan,'--k','LineWidth',1.4); %#ok<AGROW>
    end
    lg{end+1} = 'Observed'; lh(end+1) = plot(nan,nan,'k-','LineWidth',1.8); %#ok<AGROW>
    if ~isempty(lh), legend(lh, lg, 'Location','northeast'); end
else
    axis off
    text(0.5,0.5,'No distance metadata','HorizontalAlignment','center','VerticalAlignment','middle');
end

% Page title
if ~isempty(GridRC)
    sgtitle(figPer, sprintf('%s | %dx%d grid (%d bins) | %s | v>=%.0f | %s | adj=%s', ...
        ratName, GridRC(1),GridRC(2),K, strrep(binMode,'_','-'), vMin, upper(R.params.Similarity), upper(adjStr)), ...
        'FontWeight','bold');
else
    sgtitle(figPer, sprintf('%s | K=%d (PC1) | %s | v>=%.0f | %s', ...
        ratName, K, strrep(binMode,'_','-'), vMin, upper(R.params.Similarity)), ...
        'FontWeight','bold');
end

% Optional save
if doSave
  if strlength(figName)==0
    if ~isempty(GridRC)
        figName = sprintf('%s_grid_%dx%d_K%d_%s_vmin%.0f_%s_%s.png', ...
                  ratName, GridRC(1),GridRC(2),K, binMode, vMin, lower(R.params.Mode), lower(R.params.Similarity));
    else
        figName = sprintf('%s_linearSpatial_K%d_%s_vmin%.0f_%s_%s.png', ...
                  ratName, K, binMode, vMin, lower(R.params.Mode), lower(R.params.Similarity));
    end
  end
  exportgraphics(figPer, figName, 'Resolution', 300);
end

% summaries (unchanged bits)
R.summary.meanOff_pooled  = meanOff_pool;
R.summary.meanAng_pooled  = meanAng_pool;
R.spatialStats_pooled     = poolStats;
R.params.Adjacency        = adjStr;
end


% =======================================================================
function R = populationLinearSpatialOrthogonality_byRat(ratName, varargin)
% Running-only spatial population code across bins
% Either: 1D along global PC1 (legacy), or 2D grid over physical space.
% Excludes [CS, CS+2] on the velocity timebase; only keeps v >= MinSpeed.

% ---------- args ----------
p = inputParser;
addParameter(p,'MinSpeed',4,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'NBins',6,@(x) isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'GridRC',[],@(v) (isempty(v) || (isnumeric(v)&&numel(v)==2&&all(v>=1))));
addParameter(p,'BinMode','equal_size',@(s) any(strcmpi(s,{'equal_occ','equal_size'})));
addParameter(p,'LockAxis','x',@(s) any(strcmpi(s,{'x','y'}))); % only used in 1D mode
addParameter(p,'Mode','raw',@(s) any(strcmpi(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','cosine',@(s) any(strcmpi(s,{'cosine','corr','pearson'})));
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty(d));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'DayWindowMode','fraction',@(s) any(strcmpi(s,{'fraction','seconds'})));
addParameter(p,'DayWindow',[],@(v) isempty(v) || (isnumeric(v)&&numel(v)==2&&v(2)>v(1)));
parse(p,varargin{:});

vMin        = p.Results.MinSpeed;
K           = p.Results.NBins;
GridRC      = p.Results.GridRC;
binMode     = lower(p.Results.BinMode);
lockAx      = lower(p.Results.LockAxis);
modeStr     = lower(p.Results.Mode);
simStr      = lower(p.Results.Similarity);
daysArg     = p.Results.Days;
miCutoff    = p.Results.MICutoff;
miExCutoff  = p.Results.MIExcludeCutoff;
dayWinMode  = lower(p.Results.DayWindowMode);
dayWin      = p.Results.DayWindow;

useGrid = ~isempty(GridRC);
if useGrid && prod(GridRC)~=K
    error('GridRC must multiply to NBins.');
end

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
    pos = smoothpos(pos);
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

% ---------- DEFINE BINS (grid or PC1) ----------
binMeta = struct('mode','', 'gridRC',[], 'edges1D',[], 'centers1D',[], ...
                 'edgesX',[], 'edgesY',[], 'centers2D',[], 'xyMin',[], 'xyMax',[]);
if useGrid
    Rg = GridRC(1); Cg = GridRC(2);
    % edges along x/y
    switch binMode
        case 'equal_size'
            ex = linspace(min(all_x), max(all_x), Cg+1);
            ey = linspace(min(all_y), max(all_y), Rg+1);
        case 'equal_occ'
            ex = quantile(all_x, linspace(0,1,Cg+1));
            ey = quantile(all_y, linspace(0,1,Rg+1));
    end
    ex = unique(ex); if numel(ex)<Cg+1, ex = linspace(min(all_x), max(all_x), Cg+1); end
    ey = unique(ey); if numel(ey)<Rg+1, ey = linspace(min(all_y), max(all_y), Rg+1); end
    cx = (ex(1:end-1)+ex(2:end))/2;
    cy = (ey(1:end-1)+ey(2:end))/2;
    [CX,CY] = meshgrid(cx, cy);             % rows=Y (1..Rg), cols=X (1..Cg)
    binMeta.mode      = 'grid';
    binMeta.gridRC    = [Rg Cg];
    binMeta.edgesX    = ex;
    binMeta.edgesY    = ey;
    binMeta.centers2D = [CX(:), CY(:)];
    binMeta.xyMin     = [min(all_x) min(all_y)];
    binMeta.xyMax     = [max(all_x) max(all_y)];
else
    % global PC1 and sign lock
    Xall = [all_x, all_y];
    Xall = Xall - mean(Xall,1,'omitnan');
    [coeff, score] = pca(Xall, 'NumComponents',1, 'Centered',false);
    s_all = score(:,1);
    if strcmp(lockAx,'x'), rho = corr(s_all, all_x, 'rows','complete');
    else,                  rho = corr(s_all, all_y, 'rows','complete'); end
    if rho < 0, coeff(:,1) = -coeff(:,1); s_all = -s_all; end

    switch binMode
        case 'equal_occ',  edges = quantile(s_all, linspace(0,1,K+1));
        case 'equal_size', edges = linspace(min(s_all), max(s_all), K+1);
    end
    edges = unique(edges); if numel(edges) < K+1, edges = linspace(min(s_all), max(s_all), K+1); end
    centers = (edges(1:end-1)+edges(2:end))/2;

    binMeta.mode      = 'pc1';
    binMeta.gridRC    = [1 K];
    binMeta.edges1D   = edges;
    binMeta.centers1D = centers;
    binMeta.xyMin     = [min(all_x) min(all_y)];
    binMeta.xyMax     = [max(all_x) max(all_y)];
    binMeta.coeffPC1  = coeff(:,1);
    binMeta.muXY      = mean([all_x, all_y],1,'omitnan');
end

% ---------- PASS 2: per-day rates ----------
perDay = struct('rateKxN_raw',[],'cosSimK',[],'angleK',[],'dotK',[],'binMeta',binMeta);
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

    % ---- MI include ----
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

    % ---- MI EXCLUDE (task MI) ----
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

    % map running samples to bins
    if useGrid
        bx = discretize(xv, binMeta.edgesX);
        by = discretize(yv, binMeta.edgesY);
        valid  = isfinite(bx) & isfinite(by);
        binIdx = sub2ind([GridRC(1), GridRC(2)], by(valid), bx(valid)); % row-major
        occT   = accumarray(binIdx, dt(valid), [K,1], @sum, NaN);
    else
        XY  = [xv, yv]; XYc = XY - binMeta.muXY;
        s   = XYc * binMeta.coeffPC1(:);
        b1  = discretize(s, binMeta.edges1D);
        valid = ~isnan(b1);
        occT  = accumarray(b1(valid), dt(valid), [K,1], @sum, NaN);
    end
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

        % map spikes to bins
        xs = interp1(vt, xv, st, 'linear', NaN);
        ys = interp1(vt, yv, st, 'linear', NaN);
        XYs = [xs, ys]; XYs = XYs(all(isfinite(XYs),2),:);
        if isempty(XYs), continue; end

        if useGrid
            bxs = discretize(XYs(:,1), binMeta.edgesX);
            bys = discretize(XYs(:,2), binMeta.edgesY);
            ok  = isfinite(bxs) & isfinite(bys);
            bSpk = sub2ind([GridRC(1), GridRC(2)], bys(ok), bxs(ok));
        else
            XYsC = XYs - binMeta.muXY;
            ss   = XYsC * binMeta.coeffPC1(:);
            bSpk = discretize(ss, binMeta.edges1D);
            bSpk = bSpk(isfinite(bSpk));
        end

        if isempty(bSpk), continue; end
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
    perDay(gi).binMeta     = binMeta;

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
R.pooled = struct('cosSimK', pooledCos, 'angleK', pooledAng, 'binMeta', binMeta);
R.params = struct('MinSpeed', vMin, 'NBins', K, 'GridRC', GridRC, 'BinMode', binMode, ...
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



function [p_perm, diff_perm] = adjacencyPermutationP(S, adjMask, nonAdjMask, nPerm)
% Permute bin labels; recompute (mean_adj - mean_nonadj).
K  = size(S,1);
UT = triu(true(K),1);
validBins = find(~all(isnan(S),2));

% observed
sAdj    = S(adjMask);    sAdj    = sAdj(isfinite(sAdj));
sNonAdj = S(nonAdjMask); sNonAdj = sNonAdj(isfinite(sNonAdj));
diff_obs = mean(sAdj) - mean(sNonAdj);

diff_perm = nan(nPerm,1);
for t = 1:nPerm
    perm = 1:K;
    perm(validBins) = perm(validBins(randperm(numel(validBins))));
    Sp = S(perm,perm);
    sA = Sp(adjMask);     sA = sA(isfinite(sA));
    sN = Sp(nonAdjMask);  sN = sN(isfinite(sN));
    if isempty(sA) || isempty(sN)
        diff_perm(t) = NaN;
        continue
    end
    diff_perm(t) = mean(sA) - mean(sN);
end
diff_perm = diff_perm(isfinite(diff_perm));
if isempty(diff_perm)
    p_perm = NaN;
else
    % one-sided: is observed adj>nonadj larger than permuted?
    p_perm = (1 + sum(diff_perm >= diff_obs)) / (numel(diff_perm) + 1);
end
end



function Sstats = computeSpatialGridStats(S, binMeta, varargin)


p = inputParser;
addParameter(p,'Adjacency','moore',@(s) any(strcmpi(s,{'moore','rook'})));
addParameter(p,'MantelPerms',10000,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'IsCorr',true,@islogical);
parse(p,varargin{:});
adjStr = lower(p.Results.Adjacency);
nPerm  = p.Results.MantelPerms;

K = size(S,1);
S = (S + S')./2; S(1:K+1:end) = NaN;

UT = triu(true(K),1);
[idxI, idxJ] = find(UT);
sUT = S(UT);

% distances + adjacency
if isfield(binMeta,'mode') && strcmpi(binMeta.mode,'grid') && ...
   isfield(binMeta,'centers2D') && size(binMeta.centers2D,2)==2
    Rg = binMeta.gridRC(1); Cg = binMeta.gridRC(2);
    [rows, cols] = ind2sub([Rg Cg], (1:K)');
    XY = binMeta.centers2D;
    DX = XY(idxI,1) - XY(idxJ,1);
    DY = XY(idxI,2) - XY(idxJ,2);
    De = sqrt(DX.^2 + DY.^2);

    Dcheb = max(abs(rows(idxI)-rows(idxJ)), abs(cols(idxI)-cols(idxJ)));
    Dman  = abs(rows(idxI)-rows(idxJ)) + abs(cols(idxI)-cols(idxJ));
    switch adjStr
        case 'moore', adjMaskVec = (Dcheb == 1);
        case 'rook',  adjMaskVec = (Dman  == 1);
    end

    % --- GRID branch ---
    kmax = max(Dcheb(:));
    lag_k1 = nan(1,kmax);
    for k = 1:kmax
        lag_k1(k) = mean(sUT(Dcheb==k), 'omitnan');
    end
    lag_mean = [1, lag_k1];       % <-- prepend zero-lag = 1

    % --- 1-D fallback branch ---
    kmax = max(Didx(:));
    lag_k1 = nan(1,kmax);
    for k = 1:kmax
        lag_k1(k) = mean(sUT(Didx==k), 'omitnan');
    end
    lag_mean = [1, lag_k1];       % <-- prepend zero-lag = 1

    for k = 1:kmax, lag_mean(k) = mean(sUT(Dcheb==k), 'omitnan'); end

    Dfull = nan(K); Dfull(UT) = De; Dfull = (Dfull + Dfull')./2; Dfull(1:K+1:end) = 0;
else
    % 1D fallback
    if isfield(binMeta,'centers1D') && ~isempty(binMeta.centers1D)
        centers = binMeta.centers1D(:);
    else
        centers = (1:K).';
    end
    De   = abs(centers(idxI) - centers(idxJ));
    Didx = abs(idxI - idxJ);
    adjMaskVec = (Didx == 1);
    kmax = max(Didx(:));
    lag_mean = nan(1,kmax);
    for k = 1:kmax, lag_mean(k) = mean(sUT(Didx==k), 'omitnan'); end

    Dfull = nan(K); Dfull(UT) = De; Dfull = (Dfull + Dfull')./2; Dfull(1:K+1:end) = 0;
end

% summary
mean_off    = mean(sUT, 'omitnan');
sAdj        = sUT(adjMaskVec);
sNonAdj     = sUT(~adjMaskVec);
mean_adj    = mean(sAdj,    'omitnan');
mean_nonadj = mean(sNonAdj, 'omitnan');
adj_diff    = mean_adj - mean_nonadj;

% --- Mantel (robust) ---
% If distances are degenerate, return r=0, p=1 (instead of NaN)
dOK = De(isfinite(De));
if numel(dOK) < 3 || numel(unique(dOK)) < 2
    mantel_r = 0;
    mantel_p = 1;
else
    [mantel_r, mantel_p] = mantelTest_fromDistanceMatrix(S, Dfull, nPerm);
    % Last-resort guard: if anything still NaN, emit null-like values
    if ~isfinite(mantel_r)
        UT = triu(true(K),1);
        s = S(UT); d = Dfull(UT);
        ok = isfinite(s) & isfinite(d);
        if nnz(ok) >= 3 && std(s(ok))>0 && std(d(ok))>0
            mantel_r = corr(s(ok), -d(ok), 'type','Pearson');
            mantel_p = NaN;  % no perms → no p-value
        else
            mantel_r = 0; mantel_p = 1;
        end
    end
end



% linear slope/R2
ok = isfinite(sUT) & isfinite(De(:));
if nnz(ok) >= 3 && std(De(ok))>0 && std(sUT(ok))>0
    X = [ones(nnz(ok),1) De(ok)];
    b = X \ sUT(ok);
    yhat = X*b;
    ssr = sum((sUT(ok)-yhat).^2);
    sst = sum((sUT(ok)-mean(sUT(ok))).^2);
    slope = b(2);
    r2    = max(0, 1 - ssr/sst);
else
    slope = NaN; r2 = NaN;
end

% distance bins (guard degenerate case)
if all(~isfinite(De))
    edges = [0 1]; centers = 0.5; dist_mean = NaN;
else
    nb = min(10, max(4, round(sqrt(nnz(isfinite(De))))));
    dOK = De(isfinite(De));
    if isempty(dOK)
        edges = [0 1]; centers = 0.5; dist_mean = NaN;
    else
        edges = unique([min(dOK), quantile(dOK, linspace(0,1,nb)), max(dOK)]);
        edges(1) = edges(1)-eps;
        if numel(edges) < nb+1
            edges = linspace(min(dOK), max(dOK), nb+1);
            edges(1) = edges(1)-eps;
        end
        centers = (edges(1:end-1)+edges(2:end))/2;
        dist_mean = nan(1,numel(centers));
        for k = 1:numel(centers)
            dist_mean(k) = mean(sUT(De>=edges(k) & De<edges(k+1)), 'omitnan');
        end
    end
end

% -------- prepend zero-distance point (lag 0) --------
centers   = centers(:);
dist_mean = dist_mean(:);
centers   = [0; centers];
dist_mean = [1; dist_mean];
% -----------------------------------------------------

% then later keep your existing assignments:
Sstats.dist_edges   = edges;
Sstats.dist_centers = centers;


Sstats = struct();
Sstats.mean_off       = mean_off;
Sstats.mean_adj       = mean_adj;
Sstats.mean_nonadj    = mean_nonadj;
Sstats.adj_diff       = adj_diff;
Sstats.adj_perm_p     = NaN;          % (unchanged here)
Sstats.n_adj_pairs    = numel(sAdj);
Sstats.n_nonadj_pairs = numel(sNonAdj);
Sstats.lag_mean       = lag_mean;
Sstats.mantel_r       = mantel_r;
Sstats.mantel_p       = mantel_p;
Sstats.slope          = slope;
Sstats.r2             = r2;
Sstats.adjType        = adjStr;
Sstats.diff_perm      = [];           % (unchanged here)
Sstats.dist_edges     = edges;
Sstats.dist_centers   = centers;
Sstats.dist_mean      = dist_mean;
end



function [p_perm, diff_perm] = adjacencyPermutationP_vec(S, adjMaskVec, nPerm)
% Permute bin labels; recompute (mean_adj - mean_nonadj) in UT vector space.
K  = size(S,1);
UT = triu(true(K),1);
validBins = find(~all(isnan(S),2));

sUT = S(UT);
adjMaskVec = adjMaskVec(:) ~= 0;       % ensure logical column
nonAdjMaskVec = ~adjMaskVec;

% observed
sAdj    = sUT(adjMaskVec);    sAdj    = sAdj(isfinite(sAdj));
sNonAdj = sUT(nonAdjMaskVec); sNonAdj = sNonAdj(isfinite(sNonAdj));
diff_obs = mean(sAdj) - mean(sNonAdj);

diff_perm = nan(nPerm,1);
for t = 1:nPerm
    perm = 1:K;
    perm(validBins) = perm(validBins(randperm(numel(validBins))));
    Sp  = S(perm,perm);
    sUTp = Sp(UT);
    sAp  = sUTp(adjMaskVec);    sAp  = sAp(isfinite(sAp));
    sNp  = sUTp(nonAdjMaskVec); sNp  = sNp(isfinite(sNp));
    if isempty(sAp) || isempty(sNp)
        diff_perm(t) = NaN;
        continue
    end
    diff_perm(t) = mean(sAp) - mean(sNp);
end
diff_perm = diff_perm(isfinite(diff_perm));
if isempty(diff_perm)
    p_perm = NaN;
else
    % one-sided: is observed adj>nonadj larger than permuted?
    p_perm = (1 + sum(diff_perm >= diff_obs)) / (numel(diff_perm) + 1);
end
end

function [F, OUT] = plotAdjVsNonadjBars(R, varargin)
% Plot mean similarity for ADJACENT vs NON-ADJACENT spatial bins.
% Overlays per-rat lines and reports group stats.
%
% Inputs:
%   R : output struct from populationSpatialOrthogonality (multi-rat)
%
% Options:
%   'Adjacency'   : 'moore' (8-neigh, default) or 'rook' (4-neigh)
%   'MantelPerms' : permutations for per-rat stats (only used if you re-compute here)
%   'Colors'      : nRats×3 colormap for the lines (default: lines(nRats))
%   'YAxis'       : label override (default auto from R.params.Similarity)
%   'YLim'        : [ymin ymax] (optional)
%
% Outputs:
%   F   : figure handle
%   OUT : struct with per-rat values and stats

p = inputParser;
addParameter(p,'Adjacency','moore',@(s) any(strcmpi(s,{'moore','rook'})));
addParameter(p,'MantelPerms',500,@(x) isnumeric(x)&&isscalar(x));  % 0 = don't recompute
addParameter(p,'Colors',[],@(c) isempty(c) || (isnumeric(c)&&size(c,2)==3));
addParameter(p,'YAxis','',@(s) ischar(s) || isstring(s));
addParameter(p,'YLim',[],@(v) isempty(v) || (isnumeric(v)&&numel(v)==2));
parse(p,varargin{:});
adjStr   = lower(p.Results.Adjacency);
nPerm    = p.Results.MantelPerms;
clr      = p.Results.Colors;
ylabIn   = string(p.Results.YAxis);
ylimsIn  = p.Results.YLim;

assert(isfield(R,'byRat') && iscell(R.byRat), 'R must be multi-rat output with R.byRat cell array.');

nRats = numel(R.byRat);
adjVals    = nan(1,nRats);
nonadjVals = nan(1,nRats);
permP      = nan(1,nRats);

% determine label
if isfield(R,'params') && isfield(R.params,'Similarity')
    isCorr = any(strcmpi(R.params.Similarity,{'corr','pearson'}));
else
    isCorr = true;
end
ylab = ylabIn;
if strlength(ylab)==0
    ylab = tern(isCorr,'Mean Pearson similarity','Mean cosine similarity');
end

% collect per-rat adjacent/non-adjacent means
for r = 1:nRats
    Ri = R.byRat{r};
    if ~isfield(Ri,'pooled') || ~isfield(Ri.pooled,'cosSimK') || isempty(Ri.pooled.cosSimK)
        continue
    end
    C = Ri.pooled.cosSimK;
    K = size(C,1);

    % geometry
    if isfield(Ri.pooled,'binMeta') && ~isempty(Ri.pooled.binMeta)
        bm = Ri.pooled.binMeta;
    else
        % fallback: 1-D bins along index
        centers = (1:K);
        bm = struct('mode','pc1','centers1D',centers(:));
    end

    % compute stats (use your fast vectorized version)
    if exist('computeSpatialGridStats','file')==2
        stats = computeSpatialGridStats(C, bm, ...
                    'Adjacency', adjStr, ...
                    'MantelPerms', max(500,nPerm), ...
                    'IsCorr', isCorr);
        adjVals(r)    = stats.mean_adj;
        nonadjVals(r) = stats.mean_nonadj;
        if isfield(stats,'adj_perm_p'), permP(r) = stats.adj_perm_p; end
    else
        error('computeSpatialGridStats.m not found on path.');
    end
end

% group stats
ok = isfinite(adjVals) & isfinite(nonadjVals);
adjVals    = adjVals(ok);
nonadjVals = nonadjVals(ok);
ratIdx     = find(ok);
nUse = numel(adjVals);
assert(nUse>=2,'Need at least 2 rats with data.');

[~,p_t,~,tt] = ttest(adjVals, nonadjVals);           % paired t-test
p_w = signrank(adjVals, nonadjVals);                 % Wilcoxon (non-param)
d_cohen = (mean(adjVals - nonadjVals)) / std(adjVals - nonadjVals, 1);

% bar values & SEM
mAdj = mean(adjVals);  mNon = mean(nonadjVals);
seAdj = std(adjVals,0,2)/sqrt(nUse);
seNon = std(nonadjVals,0,2)/sqrt(nUse);

% colors
if isempty(clr), clr = lines(nUse); end

% figure
F = figure('Color','w','Position',[200 200 540 360]); hold on
barX = [1 2];
barY = [mAdj mNon];
bh = bar(barX, barY, 0.65, 'FaceColor',[0.80 0.85 1.00], 'EdgeColor','k'); %#ok<NASGU>
eh = errorbar(barX, barY, [seAdj seNon], 'k','LineStyle','none','LineWidth',1.5);

% per-rat lines
for i = 1:nUse
    plot(barX, [adjVals(i) nonadjVals(i)], '-o', ...
        'Color', clr(i,:), 'MarkerFaceColor', clr(i,:), 'LineWidth',1.8, 'MarkerSize',5);
end

% axes & labels
xticks(barX); xticklabels({'Adjacent','Non-adjacent'});
ylabel(ylab);
box on
if ~isempty(ylimsIn)
    ylim(ylimsIn);
else
    ymax = max([barY + [seAdj seNon], adjVals, nonadjVals]);
    ylim([0, max(0.05, ymax)*1.15]);
end

% significance bracket
yl = ylim; yb = yl(2) - 0.05*range(yl);
plot([1 1 2 2],[yb-0.01 yb yb yb-0.01],'k-','LineWidth',1.25);
txt = sigStars(p_t);
text(1.5, yb+0.01*range(yl), txt, 'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');

% legend (rats)
lgd = legend(arrayfun(@(k) sprintf('Rat %d', k), 1:nUse, 'Uni',0), ...
    'Location','northeastoutside'); %#ok<NASGU>

title(sprintf('Adjacent > Non-adjacent  (paired t: p=%.3g; Wilcoxon: p=%.3g; d=%.2f)', ...
      p_t, p_w, d_cohen));

% outputs
OUT = struct();
OUT.rats_used   = ratIdx;
OUT.adjacent    = adjVals;
OUT.nonadjacent = nonadjVals;
OUT.mean_adj    = mAdj;
OUT.mean_nonadj = mNon;
OUT.sem_adj     = seAdj;
OUT.sem_nonadj  = seNon;
OUT.paired_t    = struct('p',p_t,'tstat',tt.tstat,'df',tt.df);
OUT.signrank_p  = p_w;
OUT.cohens_d    = d_cohen;
OUT.per_rat_perm_p = permP(ok);

end

% ---------- tiny helpers ----------

  function s = sigStars(p)
  if ~isfinite(p), s = 'n.s.'; return; end
  if p < 1e-4, s='****';
  elseif p < 1e-3, s='***';
  elseif p < 1e-2, s='**';
  elseif p < 0.05, s='*';
  else, s='n.s.'; end
  end

  function [Dfull, edges, centers, masks] = buildDistanceBins(binMeta, K, nBins)
  % Build a symmetric K×K Euclidean distance matrix across bins and
  % return nBins roughly-even bins with masks for UT pairs.
  UT = triu(true(K),1);
  Dfull = nan(K);   % distances
  if isfield(binMeta,'mode') && strcmpi(binMeta.mode,'grid') && ...
     isfield(binMeta,'centers2D') && size(binMeta.centers2D,2)==2
      XY = binMeta.centers2D;         % K×2
      [I,J] = find(UT);
      dvec = sqrt(sum((XY(I,:) - XY(J,:)).^2,2));
      Dfull(UT) = dvec; Dfull = (Dfull + Dfull')./2; Dfull(1:K+1:end) = 0;
  else
      % 1D PC1 fallback
      if isfield(binMeta,'centers1D')
          c = binMeta.centers1D(:);
      else
          c = (1:K)'; % last-resort fallback
      end
      [I,J] = find(UT);
      dvec = abs(c(I) - c(J));
      Dfull(UT) = dvec; Dfull = (Dfull + Dfull')./2; Dfull(1:K+1:end) = 0;
  end
  % Choose bin edges (quantiles give ~balanced counts)
  dUT = Dfull(UT);
  dUT = dUT(isfinite(dUT));
  if isempty(dUT)
      edges = linspace(0,1,nBins+1); centers = movmean(edges,2,'Endpoints','discard'); masks = repmat({false(K)},1,nBins);
      return
  end
  edges = unique([min(dUT), quantile(dUT, linspace(0,1,nBins)), max(dUT)]);
  edges(1) = edges(1)-eps;  % include min
  if numel(edges) < nBins+1
      edges = linspace(min(dUT), max(dUT), nBins+1);
      edges(1) = edges(1)-eps;
  end
  centers = (edges(1:end-1)+edges(2:end))/2;
  masks = cell(1,numel(centers));
  for k = 1:numel(centers)
      Mk = (Dfull >= edges(k) & Dfull < edges(k+1));
      % Exclude diagonal
      Mk(1:K+1:end) = false;
      masks{k} = Mk;
  end
  end

  function C = pvSimFromRates_local(M, useCosine)
  % M: K×N (split/bin × cells)
  if useCosine
      rn = vecnorm(M,2,2);
      C  = (M*M') ./ max(rn*rn', eps);
      C  = max(-1,min(1,C));
  else
      C  = corr(M','rows','pairwise');  % Pearson across cells
  end
  C = (C + C')./2;
  C(1:size(C,1)+1:end) = 1;
  end





function dMasks = distMasksFromEdges(binMeta, K, edges)
% Return cell array of logical masks over the UT vector for each distance bin.
UT = triu(true(K),1);

% Make a full distance matrix consistent with binMeta
Dfull = nan(K);
if isfield(binMeta,'mode') && strcmpi(binMeta.mode,'grid') && ...
   isfield(binMeta,'centers2D') && size(binMeta.centers2D,2)==2
    XY = binMeta.centers2D;
    [I,J] = find(UT);
    dvec = sqrt(sum((XY(I,:) - XY(J,:)).^2,2));
else
    % 1D fallback
    if isfield(binMeta,'centers1D'), c = binMeta.centers1D(:);
    else,                            c = (1:K)';               end
    [I,J] = find(UT);
    dvec = abs(c(I) - c(J));
end

% Bin the UT distances using provided edges
dMasks = cell(1, numel(edges)-1);
for k = 1:numel(edges)-1
    dMasks{k} = (dvec >= edges(k) & dvec < edges(k+1));
end
end

function Cb = shuffle_full_per_cell(M, useCosine)
% M: K×N (K spatial bins × N cells), already normalized like the observed.
% For each cell, randomly permute its K-bin rate vector (independently).
% This destroys spatial adjacency while preserving each cell's overall rate distribution.

[K, N] = size(M);
Ms = zeros(size(M));
for j = 1:N
    v = M(:,j);
    idx = randperm(K);
    Ms(:,j) = v(idx);
end

% Return similarity in same family (cosine or Pearson)
if useCosine
    rn = vecnorm(Ms,2,2);
    Cb = (Ms*Ms') ./ max(rn*rn', eps);
    Cb = max(-1,min(1,Cb));
else
    Cb = corr(Ms','rows','pairwise');  % Pearson across cells
end
Cb = (Cb+Cb')./2;
Cb(1:K+1:end) = 1;
end

function [rMantel, pMantel] = mantelTest_fromDistanceMatrix(S, D, nPerm)
% Robust Mantel: Pearson on UT vectors, with safe fallbacks.
% If vectors are constant / degenerate, returns r=0, p=1.
% If nPerm==0 -> p=NaN (effect size only).

if nargin < 3 || isempty(nPerm), nPerm = 0; end
K = size(S,1);

% Symmetrize; blank diagonals
S = (S + S.')/2; S(1:K+1:end) = NaN;
D = (D + D.')/2; D(1:K+1:end) = NaN;

UT = triu(true(K),1);
s = S(UT);
d = D(UT);

% Finite pairs only
ok = isfinite(s) & isfinite(d);
s = s(ok); d = d(ok);

% Not enough data
if numel(s) < 3
    rMantel = NaN; pMantel = NaN; return
end

% Zero-variance / degenerate distances guard
if std(s)==0 || std(d)==0 || numel(unique(d))<2
    rMantel = 0; pMantel = 1; return
end

% Observed Pearson (proximity = -distance)
rObs = corr(s, -d, 'type','Pearson', 'rows','complete');
rMantel = rObs;

% Permutations
if nPerm <= 0
    pMantel = NaN;  % effect size only
    return
end

perm_r = nan(nPerm,1);
validBins = find(~all(isnan(S),2));  % rows/cols that have any data
for t = 1:nPerm
    perm = 1:K;
    if ~isempty(validBins)
        perm(validBins) = perm(validBins(randperm(numel(validBins))));
    else
        perm = randperm(K);
    end
    Sp = S(perm,perm);
    Dp = D(perm,perm);

    sp = Sp(UT); dp = Dp(UT);
    okp = isfinite(sp) & isfinite(dp);
    sp  = sp(okp); dp = dp(okp);

    if numel(sp) < 3 || std(sp)==0 || std(dp)==0 || numel(unique(dp))<2
        perm_r(t) = NaN; continue
    end
    perm_r(t) = corr(sp, -dp, 'type','Pearson', 'rows','complete');
end

perm_r = perm_r(isfinite(perm_r));
if isempty(perm_r)
    pMantel = NaN;
else
    % one-sided: is observed association stronger than label-shuffled?
    pMantel = (1 + sum(perm_r >= rObs)) / (numel(perm_r) + 1);
end
end


function debugMantelInputs(S, binMeta)
K = size(S,1);
S = (S + S.')/2; S(1:K+1:end) = NaN;
UT = triu(true(K),1);
s = S(UT);

% Make D the same way computeSpatialGridStats does
if isfield(binMeta,'mode') && strcmpi(binMeta.mode,'grid') && ...
   isfield(binMeta,'centers2D') && size(binMeta.centers2D,2)==2
    XY = binMeta.centers2D;
    Dfull = squareform(pdist(XY)); % Euclidean
else
    if isfield(binMeta,'centers1D') && ~isempty(binMeta.centers1D)
        c = binMeta.centers1D(:);
    else
        c = (1:K).';
    end
    Dfull = squareform(pdist(c));
end
Dfull(1:K+1:end) = NaN;
d = Dfull(UT);

ok = isfinite(s) & isfinite(d);
fprintf('[Mantel DEBUG] K=%d | UT pairs=%d | finite=%d | std(s)=%.3g | std(d)=%.3g | unique(d)=%d\n', ...
    K, nnz(UT), nnz(ok), std(s(ok)), std(d(ok)), numel(unique(d(ok))));
end

function stars = sigStar(p)
if     p<1e-3, stars='***';
elseif p<1e-2, stars='**';
elseif p<0.05, stars='*';
else,          stars='n.s.';
end
end
