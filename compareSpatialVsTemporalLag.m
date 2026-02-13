function OUT = compareSpatialVsTemporalLag(ratNames, varargin)
% Overlays spatial (populationSpatialOrthogonality) and temporal
% (populationOrthogonality_group) lag curves, with shape stats and optional
% shuffle/null reference bands.
%
% Example:
% OUT = compareSpatialVsTemporalLag( ...
%   {'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%   'NBins',[6 8 11 15 20], 'Similarity','corr', 'BinMode','equal_occ', 'Mode','demean', ...
%   'ShowShuffled',true, 'NShuff',500);

% ---------- parse args ----------
p = inputParser;  p.KeepUnmatched = true;
addParameter(p,'Mode','raw',@(s) any(strcmpi(s,{'raw','demean','zscore','meanrate'})));
addParameter(p,'Similarity','corr',@(s) any(strcmpi(s,{'corr','pearson','cosine','spearman'})));
addParameter(p,'NBins',[5,8,12,15],@(x) isnumeric(x) && all(x>=2));               % vector OK
addParameter(p,'GridRC',[],@(v) isempty(v) || (isnumeric(v)&&numel(v)==2&&all(v>=1)));
addParameter(p,'BinMode','equal_occ',@(s) any(strcmpi(s,{'equal_occ','equal_size'})));
addParameter(p,'MinSpeed',4,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Adjacency','moore',@(s) any(strcmpi(s,{'moore','rook'})));
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty(d));
addParameter(p,'NSplits',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=2));  % default later to NBins
addParameter(p,'WinSecs',[0 2],@(v) isnumeric(v)&&numel(v)==2&&v(2)>v(1));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'Spatial',struct(),@(s) isstruct(s));
addParameter(p,'Temporal',struct(),@(s) isstruct(s));
addParameter(p,'ShowShuffled',true,@(x) islogical(x) || isnumeric(x));
addParameter(p,'NShuff',500,@(x) isnumeric(x) && x>=10);
addParameter(p,'SpatialX','distance',@(s) any(strcmpi(s,{'step','distance'})));
addParameter(p,'NormalizeX',true,@(x) islogical(x) || isnumeric(x));
parse(p,varargin{:});
C = p.Results;

binsList = C.NBins(:).';

% ---------- MULTI-BIN MODE ----------
if numel(binsList) > 1
    colors = lines(numel(binsList));
    F = figure('Color','w','Position',[120 120 780 460]);
    ax = axes('Parent',F); hold(ax,'on');

    legendStr = {};
    allY = []; allX = [];
    OUT.multi = cell(1, numel(binsList));

    for ii = 1:numel(binsList)
        Ci = C; Ci.NBins = binsList(ii);
        OUTi = compareSpatialVsTemporalLag_one(ratNames, Ci);

        % Optional: show one shuffle/null reference from the first run (scaled x)
        if C.ShowShuffled && ii==1
            if ~isempty(OUTi.curves.temporal_shuff.x)
              xt = OUTi.curves.temporal_shuff.x;

                muT = OUTi.curves.temporal_shuff.mu; seT = OUTi.curves.temporal_shuff.se;
                fill(ax, [xt; flipud(xt)], [muT-seT; flipud(muT+seT)], [0 0 0], ...
                    'FaceAlpha',0.05, 'EdgeColor','none');
                plot(ax, xt, muT, '-', 'LineWidth',1.0,'Color',[0.4 0.4 0.4]);
            end
            if ~isempty(OUTi.curves.spatial_shuff.x)
                xs = OUTi.curves.spatial_shuff.x;
                muS = OUTi.curves.spatial_shuff.mu; seS = OUTi.curves.spatial_shuff.se;
                fill(ax, [xs; flipud(xs)], [muS-seS; flipud(muS+seS)], [0 0 0], ...
                    'FaceAlpha',0.05, 'EdgeColor','none');
                plot(ax, xs, muS, '--', 'LineWidth',1.0,'Color',[0.4 0.4 0.4]);
            end
        end

        % temporal solid, spatial dashed (same color); scale x to [0,1] if desired
        plot(ax, OUTi.curves.temporal.x, OUTi.curves.temporal.y, '-',  'LineWidth',2,'Color',colors(ii,:));
        plot(ax, OUTi.curves.spatial.x, OUTi.curves.spatial.y,  '--', 'LineWidth',2,'Color',colors(ii,:));

        legendStr{end+1} = sprintf('Temporal (%d bins)', binsList(ii));
        legendStr{end+1} = sprintf('Spatial  (%d bins)', binsList(ii));


        allX = [allX; OUTi.curves.temporal.x(:); OUTi.curves.spatial.x(:)];

        OUT.multi{ii} = OUTi;
        print_shape_stats(binsList(ii), OUTi);
    end

    xlabel(ax,'Lag (fraction of window)'); ylabel(ax,'Mean similarity'); grid(ax,'on'); box(ax,'on');
    xlim(ax,[0, max(allX)]); ylim(ax,[min(allY)-0.05, max(allY)+0.05]);
    title(ax,'Lag curves across binning schemes (temporal solid, spatial dashed)');
    legend(ax, legendStr, 'Location','northeastoutside');
    OUT.figure_multi = F;
    return
end

% ---------- SINGLE-BIN MODE ----------
OUT = compareSpatialVsTemporalLag_one(ratNames, C);
print_shape_stats(C.NBins, OUT);

% quick plot (single overlay + shuffle bands)
xT = OUT.curves.temporal.x;  yT = OUT.curves.temporal.y;  yTse = OUT.curves.temporal.y_sem;
xS = OUT.curves.spatial.x;   yS = OUT.curves.spatial.y;

F = figure('Color','w','Position',[120 120 720 440]); hold on

% shuffled bands (plot first so solids draw on top)
if C.ShowShuffled
    % temporal null band + dotted mean
    if isfield(OUT.curves,'temporal_shuff') && ~isempty(OUT.curves.temporal_shuff.x)
        xt = OUT.curves.temporal_shuff.x; muT = OUT.curves.temporal_shuff.mu; seT = OUT.curves.temporal_shuff.se;
        fill([xt; flipud(xt)], [muT-seT; flipud(muT+seT)], [0 0 0], 'FaceAlpha',0.06, 'EdgeColor','none');
        plot(xt, muT, '-', 'LineWidth', 1.2, 'Color', [0.2 0.2 0.2]);
    end
    % spatial null band + dotted mean
    if isfield(OUT.curves,'spatial_shuff') && ~isempty(OUT.curves.spatial_shuff.x)
        xs = OUT.curves.spatial_shuff.x; muS = OUT.curves.spatial_shuff.mu; seS = OUT.curves.spatial_shuff.se;
        fill([xs; flipud(xs)], [muS-seS; flipud(muS+seS)], [0 0 0], 'FaceAlpha',0.06, 'EdgeColor','none');
        plot(xs, muS, '--', 'LineWidth', 1.2, 'Color', [0.2 0.2 0.2]);
    end
end

% temporal SEM band (if available) + solid lines
if any(isfinite(yTse))
    xx = [xT; flipud(xT)];
    yy = [yT - yTse; flipud(yT + yTse)];
    fill(xx, yy, [0 0 0], 'FaceAlpha', 0.10, 'EdgeColor','none');
end
p1 = plot(xT, yT, '-', 'LineWidth', 2.4);              % temporal
p2 = plot(xS, yS, '-', 'LineWidth', 2.4);              % spatial

if strcmpi(pick(C,'SpatialX','distance'),'distance') && pick(C,'NormalizeX',true)
    xlabel('Separation (fraction of max distance / max lag)');
else
    xlabel('Lag (bins)');
end
if any(strcmpi(C.Similarity,{'corr','pearson'})) || strcmpi(mapTemporalSimilarity(pick(C.Temporal,'Similarity',C.Similarity)),'pearson')
    ylabel('Mean Pearson similarity');
else
    ylabel(sprintf('Mean %s similarity', lower(mapTemporalSimilarity(pick(C.Temporal,'Similarity',C.Similarity)))));
end
legend([p1 p2], {'Temporal (populationOrthogonality\_group)', 'Spatial (populationSpatialOrthogonality)'}, 'Location','northeast');
title('Lag curves: temporal vs spatial (group summaries)');
grid on; box on
xlim([0, max([xT(:); xS(:)])]);
OUT.figure = F;
end

% ================= helper that does one NBins value =====================
function OUT = compareSpatialVsTemporalLag_one(ratNames, C)

% ---- build SPATIAL cfg ----
sp = struct( ...
   'NBins',C.NBins,'GridRC',C.GridRC,'BinMode',lower(C.BinMode), ...
   'MinSpeed',C.MinSpeed,'Mode',lower(C.Mode), ...
   'Similarity',lower(C.Similarity),'Adjacency',lower(C.Adjacency), ...
   'Days',C.Days,'MICutoff',C.MICutoff,'MIExcludeCutoff',C.MIExcludeCutoff, ...
   'PerRatFigures',false,'SaveFig',false);
if ~isempty(C.Spatial), sp = mergeStructs(sp, C.Spatial); end

% ---- build TEMPORAL cfg ----
temSim    = mapTemporalSimilarity(lower( pick(C.Temporal,'Similarity', C.Similarity) ));
temSplits = pick(C.Temporal,'NSplits', C.NSplits);
if isempty(temSplits), temSplits = sp.NBins; end

tp = struct( ...
   'NSplits',temSplits,'WinSecs',pick(C.Temporal,'WinSecs', C.WinSecs), ...
   'Mode',lower(pick(C.Temporal,'Mode', C.Mode)), ...
   'Similarity',temSim, ...
   'MICutoff',pick(C.Temporal,'MICutoff', C.MICutoff), ...
   'MIExcludeCutoff',pick(C.Temporal,'MIExcludeCutoff', C.MIExcludeCutoff), ...
   'DoStats',true,'GroupSaveFig',false);
if ~isempty(C.Temporal)
    tp = mergeStructs(tp, C.Temporal);
    tp.Similarity = mapTemporalSimilarity(tp.Similarity);
end

% ---------- run SPATIAL ----------
Rsp = populationSpatialOrthogonality(ratNames, ...
    'NBins',sp.NBins,'GridRC',sp.GridRC,'BinMode',sp.BinMode, ...
    'MinSpeed',sp.MinSpeed,'Mode',sp.Mode,'Similarity',sp.Similarity, ...
    'Adjacency',sp.Adjacency,'Days',sp.Days, ...
    'MICutoff',sp.MICutoff,'MIExcludeCutoff',sp.MIExcludeCutoff, ...
    'PerRatFigures',sp.PerRatFigures,'SaveFig',false);

% --- choose spatial x-axis: step-lag vs Euclidean distance ---
useDist = strcmpi(pick(C,'SpatialX','distance'),'distance');

% ---------- build SPATIAL curve (xS, sLag) ----------
sY = []; xS = [];
if isfield(Rsp,'pooledAcrossRats') && isfield(Rsp.pooledAcrossRats,'spatialStats') ...
        && ~isempty(Rsp.pooledAcrossRats.spatialStats)

    st = Rsp.pooledAcrossRats.spatialStats;

    if useDist && isfield(st,'dist_centers') && isfield(st,'dist_mean') ...
            && ~isempty(st.dist_centers) && ~isempty(st.dist_mean)
        xS = st.dist_centers(:);
        sY = st.dist_mean(:);
        if isempty(xS) || xS(1) ~= 0
            xS = [0; xS];
            sY = [1; sY];
        elseif ~isfinite(sY(1))
            sY(1) = 1;
        end
    elseif isfield(st,'lag_mean') && ~isempty(st.lag_mean)
        sY = st.lag_mean(:);
        xS = (0:numel(sY)-1).';
    end
end

% fallback if stats weren't present
if isempty(sY)
    Cpool_for_curve = getPooledSpatialMatrix(Rsp, any(strcmpi(sp.Similarity,{'corr','pearson'})));
    if ~isempty(Cpool_for_curve)
        if useDist && isfield(Rsp,'pooledAcrossRats') && isfield(Rsp.pooledAcrossRats,'binMeta') && ~isempty(Rsp.pooledAcrossRats.binMeta)
            st2 = computeSpatialGridStats(Cpool_for_curve, Rsp.pooledAcrossRats.binMeta, ...
                'Adjacency', sp.Adjacency, 'MantelPerms', 0, ...
                'IsCorr', any(strcmpi(sp.Similarity,{'corr','pearson'})));
            if isfield(st2,'dist_centers') && ~isempty(st2.dist_centers)
                xS = st2.dist_centers(:);
                sY = st2.dist_mean(:);
                if xS(1) ~= 0, xS = [0; xS]; sY = [1; sY]; end
            end
        end
        if isempty(sY)
            sY = lagMean_from_matrix(Cpool_for_curve);
            sY = [1; sY(:)];
            xS = (0:numel(sY)-1).';
        end
    end
end

% optional x normalization (makes overlay across NBins sane)
if pick(C,'NormalizeX',true)
    if ~isempty(xS) && max(xS) > 0, xS = xS ./ max(xS); end
end
sLag = sY;

% ---------- run TEMPORAL ----------
Gtp = populationOrthogonality_group(ratNames, ...
    'NSplits',tp.NSplits,'WinSecs',tp.WinSecs, ...
    'Mode',tp.Mode,'Similarity',tp.Similarity, ...
    'MICutoff',tp.MICutoff,'MIExcludeCutoff',tp.MIExcludeCutoff, ...
    'DoStats',true,'NBoot',500,'NPerm',500, ...
    'GroupSaveFig',false);

% temporal lag (group mean ± SEM)
yT = []; yTse = [];
if isfield(Gtp,'lagMat') && ~isempty(Gtp.lagMat)
    tLag = Gtp.lagMat;                   % nRats × NSplits
    yT   = mean(tLag,1,'omitnan').';
    yTse = std(tLag,0,1,'omitnan').' ./ sqrt(max(1,sum(isfinite(tLag),1)).');
elseif isfield(Gtp,'simMean') && ~isempty(Gtp.simMean)
    yT   = lagMean_from_matrix(Gtp.simMean);
    yTse = nan(size(yT));
else
    yT   = nan(sp.NBins,1);
    yTse = nan(sp.NBins,1);
end

xT = (0:numel(yT)-1).';
if pick(C,'NormalizeX',true) && max(xT)>0
    xT = xT ./ max(xT);
end

% --- shuffled/null lag curves (temporal + spatial) ---
nSh = pick(C,'NShuff',500);

% temporal null
temNull = struct('x',[],'mu',[],'se',[]);
Ctemp = getPooledTemporalMatrix(Gtp);
if ~isempty(Ctemp)
    [muT,seT] = lag_null_from_matrix(Ctemp, nSh);
    temNull.x  = (0:numel(muT)-1).';
    temNull.mu = muT;
    temNull.se = seT;
    if pick(C,'NormalizeX',true) && max(temNull.x)>0
        temNull.x = temNull.x ./ max(temNull.x);
    end
end

% spatial null
spaNull = struct('x',[],'mu',[],'se',[]);
Csp = getPooledSpatialMatrix(Rsp, any(strcmpi(sp.Similarity,{'corr','pearson'})));
if ~isempty(Csp)
    didDistanceNull = false;

    if strcmpi(pick(C,'SpatialX','distance'),'distance') && ...
       isfield(Rsp,'pooledAcrossRats') && isfield(Rsp.pooledAcrossRats,'binMeta') && ~isempty(Rsp.pooledAcrossRats.binMeta)

        bm = Rsp.pooledAcrossRats.binMeta;

        % edges: prefer what spatialStats computed; otherwise recompute
        edges = [];
        if isfield(Rsp.pooledAcrossRats,'spatialStats') && ~isempty(Rsp.pooledAcrossRats.spatialStats) && ...
           isfield(Rsp.pooledAcrossRats.spatialStats,'dist_edges') && ~isempty(Rsp.pooledAcrossRats.spatialStats.dist_edges)
            edges = Rsp.pooledAcrossRats.spatialStats.dist_edges;
        end
        if isempty(edges)
            stTmp = computeSpatialGridStats(Csp, bm, ...
                'Adjacency', sp.Adjacency, 'MantelPerms', 0, ...
                'IsCorr', any(strcmpi(sp.Similarity,{'corr','pearson'})));
            if isfield(stTmp,'dist_edges'), edges = stTmp.dist_edges; end
        end

        if ~isempty(edges) && numel(edges) >= 2
            [xSnull, muS, seS] = distance_null_from_matrix(Csp, bm, edges, nSh);

            if pick(C,'NormalizeX',true) && max(xSnull)>0
                xSnull = xSnull ./ max(xSnull);
            end

            spaNull.x  = xSnull(:);
            spaNull.mu = muS(:);
            spaNull.se = seS(:);
            didDistanceNull = true;
        end
    end

    % fallback: step-lag null
    if ~didDistanceNull || isempty(spaNull.x)
        [muS,seS] = lag_null_from_matrix(Csp, nSh);
        spaNull.x  = (0:numel(muS)-1).';
        spaNull.mu = muS;
        spaNull.se = seS;

        if pick(C,'NormalizeX',true) && max(spaNull.x)>0
            spaNull.x = spaNull.x ./ max(spaNull.x);
        end
    end
end

% ---- SHAPE STATS on pooled curves ----
[shapeStats, expFits] = compare_curve_shapes(yT, sLag, xT, xS);

% (Optional) per-rat tau paired test
tPerRat = {}; sPerRat = {};
if isfield(Gtp,'lagMat') && ~isempty(Gtp.lagMat)
    for k = 1:size(Gtp.lagMat,1), tPerRat{k,1} = Gtp.lagMat(k,:).'; end
end
if isfield(Rsp,'perRat')
    PR = Rsp.perRat;
    for k = 1:numel(PR)
        if isfield(PR(k),'spatialStats') && isfield(PR(k).spatialStats,'lag_mean') && ~isempty(PR(k).spatialStats.lag_mean)
            sPerRat{k,1} = PR(k).spatialStats.lag_mean(:);
        elseif isfield(PR(k),'cosSimK') && ~isempty(PR(k).cosSimK)
            sPerRat{k,1} = lagMean_from_matrix(PR(k).cosSimK);
        end
    end
end
pairedTau = struct();
if ~isempty(tPerRat) && ~isempty(sPerRat) && numel(tPerRat)==numel(sPerRat)
    pairedTau = paired_tau_stats(tPerRat, sPerRat);
end

% ---------- outputs ----------
OUT = struct();
OUT.spatial  = struct('R',Rsp);
OUT.temporal = struct('G',Gtp);
OUT.curves   = struct('spatial',struct('x',xS,'y',sLag), ...
                      'temporal',struct('x',xT,'y',yT,'y_sem',yTse), ...
                      'spatial_shuff', spaNull, ...
                      'temporal_shuff', temNull);
OUT.shapeStats = shapeStats;
OUT.expFits    = expFits;
OUT.pairedTau  = pairedTau;
end

% ===== helpers =====
function S = mergeStructs(A,B)
S = A; if isempty(B), return; end
f = fieldnames(B); for i=1:numel(f), S.(f{i}) = B.(f{i}); end
end

function v = pick(S, field, defaultVal)
if isfield(S,field) && ~isempty(S.(field)), v = S.(field); else, v = defaultVal; end
end

function s = mapTemporalSimilarity(s)
if any(strcmpi(s,{'corr','pearson'})), s = 'pearson'; end
if any(strcmpi(s,{'cos'})),            s = 'cosine';  end
end

function y = lagMean_from_matrix(C)
K = size(C,1); y = nan(K,1);
for L = 0:(K-1)
    idx = diag(true(K-L,1), L) | diag(true(K-L,1), -L);
    y(L+1) = mean(C(idx), 'omitnan');
end
end

function [S, fits] = compare_curve_shapes(yT, yS, ~, ~)
L = min(numel(yT), numel(yS)); yT = yT(1:L); yS = yS(1:L); x = (0:L-1)';
zT = zscore(yT); zS = zscore(yS);
S.pearson  = corr(zT, zS, 'rows','pairwise', 'type','Pearson');
S.spearman = corr(zT, zS, 'rows','pairwise', 'type','Spearman');
S.rmse     = sqrt(mean((zT - zS).^2, 'omitnan'));
aucT = trapz(x, zT); aucS = trapz(x, zS);
S.auc_diff = abs(aucT - aucS);
nPerm = 500; diffs = zeros(nPerm,1);
for k = 1:nPerm
    flip = (rand(L,1) > 0.5)*2 - 1;
    diffs(k) = abs(trapz(x, zT.*flip) - trapz(x, zS.*flip));
end
S.auc_p_perm = mean(diffs >= S.auc_diff);
fitfun = @(p,xx) p(1)*exp(-xx/max(eps,p(2))) + p(3);
opts   = optimset('Display','off');
pT = fminsearch(@(p) nansum((yT - fitfun(p,x)).^2), [1,2,0], opts);
pS = fminsearch(@(p) nansum((yS - fitfun(p,x)).^2), [1,2,0], opts);
fits.temporal = pT;   fits.spatial = pS;
S.tau_temporal = max(eps, pT(2));
S.tau_spatial  = max(eps, pS(2));
S.tau_ratio    = S.tau_temporal / S.tau_spatial;
end

function P = paired_tau_stats(tLagPerRat, sLagPerRat)
n = numel(tLagPerRat);
tauT = nan(n,1); tauS = nan(n,1);
for k = 1:n
    [S,~]   = compare_curve_shapes(tLagPerRat{k}, sLagPerRat{k}, [], []);
    tauT(k) = S.tau_temporal;
    tauS(k) = S.tau_spatial;
end
[~,P.ttest,~,stats] = ttest(tauT, tauS);
P.tauT = tauT; P.tauS = tauS; P.tstat = stats.tstat;
end

function print_shape_stats(nbins, OUTi)
S = OUTi.shapeStats;
fprintf(['NBins=%d | r=%.3f, ρ=%.3f, RMSE_z=%.3f, ΔAUC=%.3f (perm p=%.4g), ' ...
         'τ_time=%.3f, τ_space=%.3f, τ_ratio=%.3f'], ...
        nbins, S.pearson, S.spearman, S.rmse, S.auc_diff, S.auc_p_perm, ...
        S.tau_temporal, S.tau_spatial, S.tau_ratio);
if isfield(OUTi,'pairedTau') && ~isempty(OUTi.pairedTau) && isfield(OUTi.pairedTau,'ttest')
    fprintf(', paired-τ t=%.2f, p=%.3g\n', OUTi.pairedTau.tstat, OUTi.pairedTau.ttest);
else
    fprintf('\n');
end
end

function [mu,se,allLags] = lag_null_from_matrix(C, nShuff)
K = size(C,1);
allLags = nan(K,nShuff);
for ii = 1:nShuff
    p = randperm(K);
    Cp = C(p,p);                     % permute bin order (destroys locality)
    allLags(:,ii) = lagMean_from_matrix(Cp);
end
mu = mean(allLags,2,'omitnan');
se = std(allLags,0,2,'omitnan') ./ sqrt(nShuff);
end

function [x, mu, se] = distance_null_from_matrix(C, binMeta, edges, nShuff)
% Null for similarity-vs-distance curve:
% permute bin labels in C (Cp=C(p,p)) while keeping geometry fixed,
% then average Cp(UT) within each Euclidean-distance bin.

K  = size(C,1);
UT = triu(true(K),1);
[I,J] = find(UT);

% distances for UT pairs from fixed geometry
if isfield(binMeta,'mode') && strcmpi(binMeta.mode,'grid') && ...
   isfield(binMeta,'centers2D') && size(binMeta.centers2D,2)==2
    XY = binMeta.centers2D;
    dvec = sqrt(sum((XY(I,:) - XY(J,:)).^2, 2));
else
    if isfield(binMeta,'centers1D') && ~isempty(binMeta.centers1D)
        c = binMeta.centers1D(:);
    else
        c = (1:K)';
    end
    dvec = abs(c(I) - c(J));
end

% distance-bin masks in UT-vector space
nBins = numel(edges) - 1;
masks = cell(1,nBins);
for b = 1:nBins
    masks{b} = (dvec >= edges(b) & dvec < edges(b+1));
end

centers = (edges(1:end-1) + edges(2:end)) / 2;
centers = centers(:);

allVals = nan(nBins, nShuff);
for s = 1:nShuff
    p  = randperm(K);
    Cp = C(p,p);
    v  = Cp(UT);
    for b = 1:nBins
        mk = masks{b};
        if any(mk)
            allVals(b,s) = mean(v(mk), 'omitnan');
        end
    end
end

mu = mean(allVals, 2, 'omitnan');
se = std(allVals, 0, 2, 'omitnan') ./ sqrt(nShuff);

% prepend distance=0 anchor (match your observed curve convention)
x  = [0; centers];
mu = [1; mu];
se = [0; se];
end

function C = getPooledTemporalMatrix(Gtp)
C = [];
if isfield(Gtp,'simMean') && ~isempty(Gtp.simMean), C = Gtp.simMean; return; end
if isfield(Gtp,'simMat')  && ~isempty(Gtp.simMat),  C = Gtp.simMat;  return; end
if isfield(Gtp,'C')       && ~isempty(Gtp.C),       C = Gtp.C;       return; end
end

function Cpool = getPooledSpatialMatrix(Rsp, isCorr)
Cpool = [];

% pooledAcrossRats candidates
if isfield(Rsp,'pooledAcrossRats') && ~isempty(Rsp.pooledAcrossRats)
    P = Rsp.pooledAcrossRats;

    cand = {'cosSimK','simMat','simMean','C','corrMat'};
    for i = 1:numel(cand)
        f = cand{i};
        if isfield(P,f) && ~isempty(P.(f))
            Cpool = P.(f);
            return
        end
    end
end

% perRat fallback: pool matrices across rats
if isfield(Rsp,'perRat') && ~isempty(Rsp.perRat)
    PR = Rsp.perRat;
    mats = {};
    for k = 1:numel(PR)
        if isfield(PR(k),'cosSimK') && ~isempty(PR(k).cosSimK)
            mats{end+1} = PR(k).cosSimK; %#ok<AGROW>
        elseif isfield(PR(k),'simMat') && ~isempty(PR(k).simMat)
            mats{end+1} = PR(k).simMat; %#ok<AGROW>
        end
    end
    if ~isempty(mats)
        Cpool = pool_mats(mats, isCorr);
    end
end
end


function C = pool_mats(mats, isCorr)
% mats: cell array of KxK matrices (same K)
K = size(mats{1},1);
stack = nan(K,K,numel(mats));
for i = 1:numel(mats)
    M = mats{i};
    if ~isequal(size(M),[K K]), continue; end
    stack(:,:,i) = M;
end

if isCorr
    stack = max(min(stack, 0.999999), -0.999999);
    Z = atanh(stack);
    Zm = mean(Z, 3, 'omitnan');
    C = tanh(Zm);
else
    C = mean(stack, 3, 'omitnan');
end
end
