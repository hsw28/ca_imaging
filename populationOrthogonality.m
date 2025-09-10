function R = populationOrthogonality(ratName, varargin)
% populationOrthogonality
% Build S×N population vectors (0–2 s default), compute similarities, and PLOT.
%
% Usage:
%   R = populationOrthogonality('rat0314');
%   R = populationOrthogonality('rat0314','WinSecs',[0 2],'NSplits',6,'Mode','zscore','Similarity','pearson');
%
% Name-Value:
%   'WinSecs'          – [t0 t1], default [0 2]
%   'NSplits'          – default 6
%   'Mode'             – 'raw' | 'demean' | 'zscore' (default 'raw')
%   'Similarity'       – 'cosine' | 'pearson' | 'spearman' (default 'pearson')
%   'MICutoff'         – [] (off) or scalar in [0,1]; INCLUDE if MI(:,3) >= cutoff (from MI_CSUS15_shuff)
%   'MIExcludeCutoff'  – [] (off) or scalar in [0,1]; EXCLUDE if MI(:,3) >= cutoff (from MI_noCSUS15_shuff)
%   'SaveFig'          – true/false, default false
%   'FigName'          – filename if SaveFig=true (default auto)
%   'DoStats'    : true/false (default true)
%   'NBoot'      : 1000   (cell bootstrap for CIs)
%   'NPerm'      : 5000   (permutations for Mantel & epoch test)
%   'EpochEdges' : [0 0.25 0.75 0.85 2]  (relative to WinSecs)

p = inputParser;
addParameter(p,'WinSecs',[0 2], @(v) isnumeric(v) && numel(v)==2 && v(2)>v(1));
addParameter(p,'NSplits',16, @(x) isnumeric(x) && isscalar(x) && x>=2);
addParameter(p,'Mode','raw',@(s) any(validatestring(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','pearson',@(s) any(validatestring(lower(s),{'cosine','pearson','spearman'})));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'SaveFig',false,@islogical);
addParameter(p,'FigName','',@(s) ischar(s) || isstring(s));
% NEW
addParameter(p,'DoStats',true,@islogical);
addParameter(p,'NBoot',500,@isscalar);
addParameter(p,'NPerm',500,@isscalar);
addParameter(p,'EpochEdges',[0 0.25 0.75 0.85 2],@(v) isnumeric(v) && numel(v)==5);
parse(p,varargin{:});
winSecs    = p.Results.WinSecs;
nSplits    = p.Results.NSplits;
modeStr    = lower(p.Results.Mode);
simStr     = lower(p.Results.Similarity);
miCutoff   = p.Results.MICutoff;
miXCutoff  = p.Results.MIExcludeCutoff;
doSave     = p.Results.SaveFig;
figName    = string(p.Results.FigName);
doStats    = p.Results.DoStats;
nBoot      = p.Results.NBoot;
nPerm      = p.Results.NPerm;
epochEdges = p.Results.EpochEdges;

% ---- compute (reuses the “byRat” helper) ----
R = populationOrthogonality_byRat(ratName, ...
      'WinSecs',winSecs,'NSplits',nSplits,'Mode',modeStr,'Similarity',simStr, ...
      'MICutoff',miCutoff,'MIExcludeCutoff',miXCutoff);

% ---- collect per-day panels (similarity & norms) ----
days = R.days; nDays = numel(days);
haveDay = false(1,nDays);
meanOff = nan(1,nDays);
meanAng = nan(1,nDays);
rowNormsAll = cell(1,nDays);

for d = 1:nDays
  if isempty(R.perDay(d).simMat), continue; end
  C = R.perDay(d).simMat;
  haveDay(d) = true;
  offMask = ~eye(nSplits);
  meanOff(d) = mean(C(offMask),'omitnan');
  if strcmp(simStr,'cosine')
      meanAng(d) = real(acos(max(-1,min(1,meanOff(d))))); % radians
  else
      meanAng(d) = NaN;
  end
  if ~isempty(R.perDay(d).rate5xN)
    rowNormsAll{d} = vecnorm(R.perDay(d).rate5xN,2,2);
  else
    rowNormsAll{d} = nan(nSplits,1);
  end
end

% pooled summaries
Cpool = R.pooled.simMat;
offMask = ~eye(nSplits);
meanOff_pool = mean(Cpool(offMask),'omitnan');
if strcmp(simStr,'cosine')
    meanAng_pool = real(acos(max(-1,min(1,meanOff_pool))));
else
    meanAng_pool = NaN;
end

% =====================  STATS on pooled matrix  =========================
stats = struct();   % will attach to R.stats
if doStats && all(isfinite(Cpool(:)))
    % gather pooled split×cell matrix (concat cells across days)
    Mpool = [];
    for d = 1:nDays
        M = R.perDay(d).rate5xN;   % (splits × cells)
        if ~isempty(M), Mpool = [Mpool, M]; end %#ok<AGROW>
    end
    nCells = size(Mpool,2);


    % A) lag–similarity curve + bootstrap CI
    lags = 0:(nSplits-1);  utLag = cell(1,nSplits);
    lagMean = nan(1,nSplits);
    for L = lags
        idx = diag(true(nSplits-L,1), L) | diag(true(nSplits-L,1), -L);
        utLag{L+1} = find(idx);                         %#ok<AGROW>
        lagMean(L+1) = mean(Cpool(idx),'omitnan');
    end
    lagBoot = nan(nBoot, nSplits);
    if nCells >= 5
        for b = 1:nBoot
            cols = randsample(nCells, nCells, true);
            Cb   = simFromM(Mpool(:, cols));
            for L = lags
                lagBoot(b, L+1) = mean(Cb(utLag{L+1}), 'omitnan');
            end
        end
        stats.lag.ci_low  = prctile(lagBoot, 2.5, 1);
        stats.lag.ci_high = prctile(lagBoot,97.5, 1);
    else
        stats.lag.ci_low = nan(1,nSplits); stats.lag.ci_high = stats.lag.ci_low;
    end
    stats.lag.mean = lagMean;
    [stats.lag.rho, stats.lag.p_rho] = corr((lags(:)), lagMean(:), 'Type','Spearman','Rows','complete');
    halfTarget = lagMean(1) * 0.5;
    ix = find(lagMean <= halfTarget, 1, 'first');
    if isempty(ix)
        stats.lag.halfwidth_bins = NaN;
    else
        stats.lag.halfwidth_bins = ix - 1;
    end

    % B) Mantel test vs temporal proximity model
    Dmodel = -abs((1:nSplits) - (1:nSplits)'); ut = triu(true(nSplits),1);
    r_true = corr( Cpool(ut), Dmodel(ut), 'Type','Spearman','Rows','complete');
    r_perm = nan(nPerm,1);
    ord = 1:nSplits;
    for ppp = 1:nPerm
        perm = ord(randperm(nSplits));
        Xp = Cpool(perm, perm);
        r_perm(ppp) = corr( Xp(ut), Dmodel(ut), 'Type','Spearman','Rows','complete');
    end
    stats.mantel.r = r_true;
    stats.mantel.p = mean(r_perm >= r_true);

    % C) RSA: similarity ~ SameEpoch + Proximity
    tEdges = linspace(winSecs(1), winSecs(2), nSplits+1);
    tCtr   = (tEdges(1:end-1)+tEdges(2:end))/2;
    epID   = zeros(1,nSplits);
    E = epochEdges(:)';   % relative to WinSecs
    for k = 1:nSplits
        t = tCtr(k);
        if     t>=E(1) && t<E(2), epID(k)=1;
        elseif t>=E(2) && t<E(3), epID(k)=2;
        elseif t>=E(3) && t<E(4), epID(k)=3;
        else,  epID(k)=4;
        end
    end
    SameEpoch = double(epID(:)==epID(:)');
    Proximity = Dmodel;

    y  = Cpool(ut);
    Xr = [ SameEpoch(ut) Proximity(ut) ];
    Xr = zscore(Xr);                   % scale predictors
    b  = Xr \ y;
    stats.rsa.beta_same = b(1);
    stats.rsa.beta_prox = b(2);

    % bootstrap CIs for betas
    betab = nan(nBoot,2);
    if nCells >= 5
        for bidx = 1:nBoot
            cols = randsample(nCells, nCells, true);
            Cb   = simFromM(Mpool(:, cols));
            yb   = Cb(ut);
            betab(bidx,:) = (Xr \ yb).';
        end
        stats.rsa.beta_same_CI = prctile(betab,[2.5 97.5]);
        stats.rsa.beta_prox_CI = prctile(betab,[2.5 97.5]);
    else
        stats.rsa.beta_same_CI = [NaN NaN]; stats.rsa.beta_prox_CI = [NaN NaN];
    end

    % D) Within vs between epoch similarity (permutation p)
    W = []; B = [];
    for i=1:nSplits
        for j=i+1:nSplits
            if epID(i)==epID(j), W(end+1)=Cpool(i,j); %#ok<AGROW>
            else,                 B(end+1)=Cpool(i,j); %#ok<AGROW>
            end
        end
    end
    stats.epoch.within_mean  = mean(W,'omitnan');
    stats.epoch.between_mean = mean(B,'omitnan');
    diff_true = stats.epoch.within_mean - stats.epoch.between_mean;
    diff_perm = nan(nPerm,1);
    for ppp = 1:nPerm
        ep = epID(randperm(nSplits));
        Wp = []; Bp = [];
        for i=1:nSplits
            for j=i+1:nSplits
                if ep(i)==ep(j), Wp(end+1)=Cpool(i,j); else, Bp(end+1)=Cpool(i,j); end %#ok<AGROW>
            end
        end
        diff_perm(ppp) = mean(Wp,'omitnan') - mean(Bp,'omitnan');
    end
    stats.epoch.diff   = diff_true;
    stats.epoch.p_perm = mean(diff_perm >= diff_true);
end

% ---- figure layout ----
% ---- figure layout ----
nHeat = sum(haveDay) + 1;
nCols = max(3, nHeat);
nRows = 2 + double(doStats);   % <— add one extra row for stats panels
fig = figure('Color','w','Position',[100 100 320*nCols 820]);
tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');




% ---- heatmaps (per-day + pooled) ----
plotIdx = 0;
for d = 1:nDays
  if ~haveDay(d), continue; end
  plotIdx = plotIdx + 1;
  nexttile(plotIdx);
  clim = strcmp(simStr,'cosine')*[0 1] + ~strcmp(simStr,'cosine')*[-1 1];
  imagesc(R.perDay(d).simMat, clim);
  axis square; colormap(gca, parula); colorbar;
  xticks(1:nSplits); yticks(1:nSplits);
  xlabel('Split'); ylabel('Split');
  if strcmp(simStr,'cosine')
      ttlExtra = sprintf('mean off-diag cos = %.2f (%.0f°)', meanOff(d), rad2deg(meanAng(d)));
  else
      ttlExtra = sprintf('mean off-diag %s = %.2f', simStr, meanOff(d));
  end
  title(sprintf('%s  |  %s', days{d}, ttlExtra));
end
% pooled
plotIdx = plotIdx + 1;
nexttile(plotIdx);
clim = strcmp(simStr,'cosine')*[0 1] + ~strcmp(simStr,'cosine')*[-1 1];
imagesc(Cpool, clim);
axis square; colormap(gca, parula); colorbar;
xticks(1:nSplits); yticks(1:nSplits);
xlabel('Split'); ylabel('Split');
if strcmp(simStr,'cosine')
    ttlExtra = sprintf('mean off-diag cos = %.2f (%.0f°)', meanOff_pool, rad2deg(meanAng_pool));
else
    ttlExtra = sprintf('mean off-diag %s = %.2f', simStr, meanOff_pool);
end
title(['Pooled  |  ' ttlExtra]);

% ---- MDS of pooled (distance = 1 - similarity) ----
%{
nexttile([1 max(1, floor(nCols/2))]);
D = 1 - max(-1,min(1,Cpool));
D = (D + D')/2; D(1:nSplits+1:end) = 0;
try, Y = mdscale(D, 2, 'Criterion','stress'); catch, [V,~]=eigs((Cpool+Cpool')/2,2); Y=V; end
plot(Y(:,1), Y(:,2), '-o','LineWidth',1.5,'MarkerFaceColor',[.1 .1 .1]); hold on
text(Y(:,1), Y(:,2), compose('%d',1:nSplits), 'VerticalAlignment','bottom','FontSize',10);
axis equal; grid on
xlabel('MDS-1'); ylabel('MDS-2');
title(sprintf('Pooled: MDS of splits (dist = 1 - %s)', simStr));


% ---- grouped bars: vector norms per split (per day + pooled mean) ----
nexttile([1 max(1, ceil(nCols/2))]); cla; hold on
validIdx = find(haveDay); nValid = numel(validIdx);
Rnorms = nan(nSplits, nValid); lbls = strings(1, nValid);
for j = 1:nValid, d = validIdx(j); Rnorms(:, j) = rowNormsAll{d}(:); lbls(j)=string(days{d}); end
rnMean = mean(Rnorms, 2, 'omitnan'); Ybars = [Rnorms, rnMean];
b = bar(Ybars, 'grouped'); box on
for j = 1:nValid, b(j).DisplayName = lbls(j); end
b(nValid+1).DisplayName = 'Mean(non-NaN)'; b(nValid+1).FaceAlpha = 0.45;
xticks(1:nSplits); xlim([0.5 nSplits+0.5])
xlabel('Split'); ylabel(sprintf('||population||_2 (%s)', upper(R.params.Mode)))
title('Population vector norms per split')
legend('Location','bestoutside')
ylim([0, max(Ybars(:), [], 'omitnan') * 1.1])

%}

% ---- NEW PANELS: Lag curve + Within-vs-Between ----
if doStats && isfield(stats,'lag')

  % Within vs between
  nexttile([1 max(1, ceil(nCols/3))]); cla; hold on
  bar([1 2], [stats.epoch.within_mean, stats.epoch.between_mean], 0.6);
  set(gca,'XTick',[1 2],'XTickLabel',{'Within','Between'});
  ylabel(sprintf('Mean %s similarity', simStr));
  title(sprintf('Within > Between?  \\Delta=%.3f, p_{perm}=%.3g', ...
        stats.epoch.diff, stats.epoch.p_perm));
  box off

  % Text box: Mantel + RSA
  annotation('textbox',[0.78 0.02 0.20 0.12],'String', sprintf( ...
    'Mantel r=%.2f, p=%.1g\nRSA  β_same=%.3f [%0.3f %0.3f]\nRSA  β_prox=%.3f [%0.3f %0.3f]', ...
    stats.mantel.r, stats.mantel.p, ...
    stats.rsa.beta_same, stats.rsa.beta_same_CI(1), stats.rsa.beta_same_CI(2), ...
    stats.rsa.beta_prox, stats.rsa.beta_prox_CI(1), stats.rsa.beta_prox_CI(2) ), ...
    'EdgeColor',[.8 .8 .8],'BackgroundColor',[1 1 1], 'FontSize',10);

    % Lag curve
    nexttile([1 max(1, 2)]); hold on
    x = 0:(nSplits-1);
    plot(x, stats.lag.mean, 'k-', 'LineWidth',1.8);
    if all(isfinite(stats.lag.ci_low))
        xx = [x fliplr(x)];
        yy = [stats.lag.ci_low fliplr(stats.lag.ci_high)];
        fill(xx, yy, [0.7 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.5);
        plot(x, stats.lag.mean, 'k-', 'LineWidth',1.8);
    end
    xlabel('Lag (bins)'); ylabel(sprintf('Mean %s similarity', simStr));
    ttl = sprintf('Lag curve (\\rho=%.2f, p=%.1g; half-width=%s bins)', ...
                  stats.lag.rho, stats.lag.p_rho, num2str(stats.lag.halfwidth_bins));
    title(ttl); box off


end

% ---- save if requested ----
if doSave
  if strlength(figName)==0
    figName = sprintf('%s_orthogonality_%dsplits_%0.1f-%0.1fs_%s_%s.png', ...
              ratName, nSplits, winSecs(1), winSecs(2), modeStr, simStr);
  end
  exportgraphics(fig, figName, 'Resolution', 300);
end

% summaries + stats back out
R.summary.meanOff_perDay = meanOff;
R.summary.meanAng_perDay = meanAng;
R.summary.meanOff_pooled = meanOff_pool;
R.summary.meanAng_pooled = meanAng_pool;
R.stats = stats;


    function C = simFromM(Min)
        if strcmp(simStr,'cosine')
            rn = vecnorm(Min,2,2);
            C  = (Min*Min') ./ max(rn*rn', eps);
            C  = max(-1,min(1,C));
        else
            C  = corr(Min','rows','pairwise','Type',upper(simStr));
        end
        C = (C + C')/2;
        C(1:nSplits+1:end) = 1;
    end
end  % <- end of populationOrthogonality






function R = populationOrthogonality_byRat(ratName, varargin)
% Build S×N population vectors of mean trial-window rates,
% and compute similarity matrices according to 'Similarity'.

p = inputParser;
addParameter(p,'WinSecs',[0 2], @(v) isnumeric(v) && numel(v)==2 && v(2)>v(1));
addParameter(p,'NSplits',6, @(x) isnumeric(x) && isscalar(x) && x>=2);
addParameter(p,'Mode','raw',@(s) any(validatestring(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','cosine',@(s) any(validatestring(lower(s),{'cosine','pearson','spearman'})));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));          % INCLUDE if >=
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));   % EXCLUDE if >=
parse(p,varargin{:});
winSecs   = p.Results.WinSecs;
nSplits   = p.Results.NSplits;
modeStr   = lower(p.Results.Mode);
simStr    = lower(p.Results.Similarity);
miCutoff  = p.Results.MICutoff;
miXCutoff = p.Results.MIExcludeCutoff;

rat   = evalin('base', ratName);
dates = autoDateList(rat);
iAn   = find(strcmp(dates, rat.An), 1);
if isempty(iAn), error('Could not locate An date in autoDateList for %s.', ratName); end
days = dates(max(1,iAn-2):iAn);

perDay = struct('rate5xN_raw',[],'rate5xN',[],'simMat',[],'angle5',[],'dot5',[]);
simStack = nan(nSplits,nSplits,numel(days));

for d = 1:numel(days)
  D = days{d};
  Sfld = sprintf('CA_peaks_%s', D);
  Cfld = sprintf('CS_%s',       D);
  Mfld = sprintf('ratemask_%s', D);
  if ~isfield(rat.Ca_peaks,Sfld) || ~isfield(rat.CS_times,Cfld)
    continue
  end

  S        = rat.Ca_peaks.(Sfld);
  csTimes  = rat.CS_times.(Cfld)(:);

  % ---- baseline "good" mask from ratemask (or keep all if missing) ----
  if isfield(rat,'ratemask') && isfield(rat.ratemask,Mfld)
      baseMask = rat.ratemask.(Mfld) == 1;
  else
      % derive nCells from S
      if iscell(S), nCells = numel(S); else, nCells = size(S,1); end
      baseMask = true(nCells,1);
  end

  % ---- INCLUDE filter: MI_CSUS15_shuff >= MICutoff ----
  if ~isempty(miCutoff)
      incOK = false(size(baseMask));
      if isfield(rat, 'MI_CSUS15_shuff')
          MIroot = sprintf('MI_%s', D);
          if isfield(rat.MI_CSUS15_shuff, MIroot)
              incVals = rat.MI_CSUS15_shuff.(MIroot); % n×4
              keepInc = incVals(:,3) >= miCutoff;
              incOK   = logical(padToLength(keepInc, numel(baseMask)));
          end
      end
      baseMask = baseMask & incOK;
  end

  % ---- EXCLUDE filter: MI_noCSUS15_shuff >= MIExcludeCutoff ----
  if ~isempty(miXCutoff)
      excMask = false(size(baseMask));
      if isfield(rat, 'MI_noCSUS15_shuff')
          MIroot2 = sprintf('MI_%s', D);
          if isfield(rat.MI_noCSUS15_shuff, MIroot2)
              excVals = rat.MI_noCSUS15_shuff.(MIroot2); % n×4
              bad     = excVals(:,3) >= miXCutoff;
              excMask = logical(padToLength(bad, numel(baseMask)));
          end
      end
      baseMask = baseMask & ~excMask;
  end

  good = find(baseMask(:));
  if isempty(csTimes) || isempty(good), continue, end

  [rates, ok] = perTrialRates_CS(S, csTimes, good, winSecs);
  if ~ok, continue, end

  % split trials into contiguous groups, avg per neuron
  T = size(rates,1);
  edges = round(linspace(0, T, nSplits+1));
  rate5xN_raw = nan(nSplits, size(rates,2));
  for k = 1:nSplits
    idx = (edges(k)+1):edges(k+1);
    if ~isempty(idx)
      rate5xN_raw(k,:) = mean(rates(idx,:), 1, 'omitnan');
    end
  end

  % per-cell normalization across splits
  rate5xN = normalizeRatesPerCell(rate5xN_raw, modeStr);

  % similarities
  [dot5, simMat, angle5] = pvSimilaritiesFlexible(rate5xN, simStr);

  perDay(d).rate5xN_raw = rate5xN_raw;
  perDay(d).rate5xN     = rate5xN;
  perDay(d).simMat      = simMat;
  perDay(d).angle5      = angle5;
  perDay(d).dot5        = dot5;

  simStack(:,:,d) = simMat;
end

validD = ~isnan(squeeze(simStack(1,1,:)));
if any(validD)
  pooledSim = nanmean(simStack(:,:,validD), 3);
else
  pooledSim = nan(nSplits,nSplits);
end

R = struct();
R.ratVar = ratName;
R.days   = days;
R.perDay = perDay;
R.pooled = struct('simMat', pooledSim);
R.params = struct('WinSecs', winSecs, 'NSplits', nSplits, 'Mode', modeStr, ...
                  'Similarity', simStr, 'MICutoff', miCutoff, ...
                  'MIExcludeCutoff', miXCutoff);
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
    X = (R - mu) ./ sd;
  otherwise
    error('Unknown Mode: %s', modeStr);
end
end

function [rates, ok] = perTrialRates_CS(S, csTimes, good, winSecs)
t0 = winSecs(1); t1 = winSecs(2); dur = t1 - t0;
T = numel(csTimes); N = numel(good);
rates = nan(T, N); ok = true;
for j = 1:N
  c  = good(j);
  st = getCellSpikes(S, c);
  if isempty(st),  rates(:,j) = NaN; continue; end
  for tr = 1:T
    w0 = csTimes(tr) + t0; w1 = csTimes(tr) + t1;
    rates(tr,j) = sum(st >= w0 & st < w1) / dur;
  end
end
if all(isnan(rates(:))), ok = false; end
end

function [dotM, simM, angM] = pvSimilaritiesFlexible(X, simStr)
% X is (splits × cells). Return dot-products; similarity per 'simStr'.
rowNorms = vecnorm(X, 2, 2);
dotM     = X * X.';
switch simStr
  case 'cosine'
    den   = rowNorms * rowNorms.';
    simM  = dotM ./ max(den, eps);
    simM  = max(-1, min(1, simM));
    angM  = real(acos(simM));        % radians
  case 'pearson'
    simM  = pairwiseCorr(X, 'Pearson');
    angM  = NaN(size(simM));
  case 'spearman'
    simM  = pairwiseCorr(X, 'Spearman');
    angM  = NaN(size(simM));
  otherwise
    error('Unknown Similarity: %s', simStr);
end
end

function C = pairwiseCorr(X, typeStr)
S = size(X,1);
C = nan(S,S);
for i = 1:S
  for j = i:S
    xi = X(i,:); xj = X(j,:);
    good = isfinite(xi) & isfinite(xj);
    if nnz(good) >= 2
      C(i,j) = corr(xi(good).', xj(good).', 'Type', typeStr);
      C(j,i) = C(i,j);
    end
  end
end
end

function st = getCellSpikes(S, c)
  if iscell(S), st = S{c}(:);
  else,         st = S(c,:).';
  end
  st = st(~isnan(st) & st>0);
end

function v = padToLength(v, L)
% pad/truncate logical or numeric vector v to length L (pad with false/0)
v = v(:);
if numel(v) < L, v(end+1:L) = 0; end
if numel(v) > L, v = v(1:L);   end
end


function S = analyzePVSimStats(R, epochEdges, nBoot, nPerm)
% epochEdges: [0 0.25 0.75 0.85 2] (relative to your WinSecs), used to label bins
% nBoot: #bootstraps over cells (default 1000); nPerm: #Mantel perms (default 5000)

if nargin<2 || isempty(epochEdges), epochEdges = [0 0.25 0.75 0.85 2]; end
if nargin<3, nBoot = 1000; end
if nargin<4, nPerm = 5000; end

S = struct; nSplits = size(R.pooled.simMat,1);
lags = 0:(nSplits-1);
S.lags = lags;

% ---------- A) Lag curves (pooled) with cell-bootstrap ----------
X = R.pooled.simMat;         % similarity across splits (pooled)
lagMean = nan(1,nSplits);
for L = lags
    idx = diag(true(nSplits-L,1), L) | diag(true(nSplits-L,1), -L);
    lagMean(L+1) = mean(X(idx), 'omitnan');
end
S.lag.mean = lagMean;

% bootstrap using stored per-day split×cell matrices (faster than recomputing)
useCosine = strcmpi(R.params.Similarity,'cosine');
getSim = @(M) pvSimFromRates(M, useCosine);   % nested below

% stack all days' split×cell matrices (concatenate cells)
Mpool = [];
for d = 1:numel(R.perDay)
    M = R.perDay(d).rate5xN;
    if ~isempty(M), Mpool = [Mpool, M]; end %#ok<AGROW>
end
nCells = size(Mpool,2);

lagBoot = nan(nBoot, nSplits);
for b = 1:nBoot
    cols = randsample(nCells, nCells, true);  % resample cells
    Cb   = getSim(Mpool(:, cols));
    for L = lags
        idx = diag(true(nSplits-L,1), L) | diag(true(nSplits-L,1), -L);
        lagBoot(b, L+1) = mean(Cb(idx), 'omitnan');
    end
end
S.lag.ci_low  = prctile(lagBoot, 2.5, 1);
S.lag.ci_high = prctile(lagBoot,97.5, 1);

% monotonic decay & half-width
[ S.lag.rho, S.lag.rho_p ] = corr((lags(:)), lagMean(:), 'Type','Spearman','Rows','complete');
target = lagMean(1) * 0.5;     % half of lag-0
ix = find(lagMean <= target, 1, 'first');
S.lag.halfwidth = (ix-1);      % in lags (bins); empty if never reaches half

% ---------- B) Mantel test (temporal proximity) ----------
Dmodel = -abs((1:nSplits) - (1:nSplits)');   % larger (closer to 0 lag) => higher similarity
ut = triu(true(nSplits),1);
r_true = corr( X(ut), Dmodel(ut), 'Type','Spearman','Rows','complete');
r_perm = nan(nPerm,1);
ord = 1:nSplits;
for p = 1:nPerm
    perm = ord(randperm(nSplits));
    Xp = X(perm, perm);
    r_perm(p) = corr( Xp(ut), Dmodel(ut), 'Type','Spearman','Rows','complete');
end
S.mantel.r = r_true;
S.mantel.p = mean(r_perm >= r_true);    % one-sided: more positive than null

% ---------- C) RSA: proximity vs same-epoch ----------
binTimes = linspace(R.params.WinSecs(1), R.params.WinSecs(2), nSplits+1);
binCtrs  = (binTimes(1:end-1)+binTimes(2:end))/2;
epochID = zeros(1,nSplits); E = epochEdges(:)';
for k = 1:nSplits
    t = binCtrs(k);
    if     t>=E(1) && t<E(2), epochID(k)=1;  % CS
    elseif t>=E(2) && t<E(3), epochID(k)=2;  % Trace
    elseif t>=E(3) && t<E(4), epochID(k)=3;  % US
    else,  epochID(k)=4;                      % Post
    end
end
SameEpoch = double( epochID(:)==epochID(:)' );
Proximity = Dmodel;                     % −|i−j|
y  = X(ut);
Xr = [ SameEpoch(ut) Proximity(ut) ];
Xr = zscore(Xr);                         % scale predictors
b  = Xr \ y;                             % OLS on upper triangle
S.rsa.beta_same   = b(1);
S.rsa.beta_prox   = b(2);

% bootstrap CIs for betas (resample cells)
betab = nan(nBoot,2);
for bidx=1:nBoot
    cols = randsample(nCells, nCells, true);
    Cb   = getSim(Mpool(:, cols));
    yb   = Cb(ut);
    betab(bidx,:) = (Xr \ yb).';
end
S.rsa.beta_same_CI = prctile(betab,[2.5 97.5]);
S.rsa.beta_prox_CI = prctile(betab,[2.5 97.5]);

% ---------- D) Within/between epoch contrasts ----------
W = []; B = [];
for i=1:nSplits
    for j=i+1:nSplits
        if epochID(i)==epochID(j), W(end+1)=X(i,j); %#ok<AGROW>
        else,                       B(end+1)=X(i,j); %#ok<AGROW>
        end
    end
end
S.epoch.within_mean  = mean(W,'omitnan');
S.epoch.between_mean = mean(B,'omitnan');
% permutation test: shuffle epoch labels
diff_true = S.epoch.within_mean - S.epoch.between_mean;
diff_perm = nan(nPerm,1);
for p = 1:nPerm
    ep = epochID(randperm(nSplits));
    Wp = []; Bp = [];
    for i=1:nSplits
        for j=i+1:nSplits
            if ep(i)==ep(j), Wp(end+1)=X(i,j); else, Bp(end+1)=X(i,j); end %#ok<AGROW>
        end
    end
    diff_perm(p) = mean(Wp) - mean(Bp);
end
S.epoch.diff   = diff_true;
S.epoch.p_perm = mean(diff_perm >= diff_true);
end

% ---- helper to recompute similarity from split×cell matrix M (S×N) ----
function C = pvSimFromRates(M, useCosine)
    if useCosine
        rn = vecnorm(M,2,2); C = (M*M') ./ max(rn*rn', eps);
        C = max(-1,min(1,C));
    else
        C = corr(M','rows','pairwise');  % Pearson
    end
end
