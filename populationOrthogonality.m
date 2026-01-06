function R = populationOrthogonality(ratName, varargin)
% populationOrthogonality
% Build S×N population vectors (0–2 s default), compute similarities, and PLOT.
%
% Usage:
%   R = populationOrthogonality('rat0314');
%   R = populationOrthogonality('rat0314','WinSecs',[0 2],'NSplits',16,'Mode','zscore','Similarity','pearson');
%
% Name-Value:
%   'WinSecs'          – [t0 t1], default [0 2]
%   'NSplits'          – default 16
%   'Mode'             – 'raw' | 'demean' | 'zscore' (default 'raw')
%   'Similarity'       – 'cosine' | 'pearson' | 'spearman' (default 'pearson')
%   'MICutoff'         – [] (off) or scalar in [0,1]; INCLUDE if MI(:,3) >= cutoff (from MI_CSUS15_shuff)
%   'MIExcludeCutoff'  – [] (off) or scalar in [0,1]; EXCLUDE if MI(:,3) >= cutoff (from MI_noCSUS15_shuff)
%   'SaveFig'          – true/false, default false
%   'FigName'          – filename if SaveFig=true (default auto)
%   'DoStats'          – true/false (default true)
%   'NBoot'            – 500 (cell bootstrap for CI of lag curve)
%   'NPerm'            – 500 (permutation test for within>between)
%   'EpochEdges'       – [0 0.25 0.75 0.85 2] relative to WinSecs (CS, trace, US, post)
%
% Adds per-cell shuffled null for lag curve (95% band + dashed mean),
% plus Within vs Between epoch bars with paired t, Wilcoxon, and perm p.
% Layout mirrors populationSpatialOrthogonality_singleRat: top heatmaps, bottom bars + lag.

p = inputParser;
addParameter(p,'WinSecs',[0 2], @(v) isnumeric(v) && numel(v)==2 && v(2)>v(1));
addParameter(p,'NSplits',15, @(x) isnumeric(x) && isscalar(x) && x>=2);
addParameter(p,'Mode','demean',@(s) any(validatestring(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','pearson',@(s) any(validatestring(lower(s),{'cosine','pearson','spearman'})));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'SaveFig',false,@islogical);
addParameter(p,'FigName','',@(s) ischar(s) || isstring(s));
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

% ---- collect per-day panels (similarity) ----
days = R.days; nDays = numel(days);
haveDay = false(1,nDays);
meanOff = nan(1,nDays);
meanAng = nan(1,nDays);

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

if doStats && any(isfinite(Cpool(:)))
    % stack pooled split×cell matrix across days
    Mpool = [];
    for d = 1:nDays
        M = R.perDay(d).rate5xN;   % (splits × cells)
        if ~isempty(M), Mpool = [Mpool, M]; end %#ok<AGROW>
    end
    nCells = size(Mpool,2);

    % ---------- Lag curve (observed) ----------
    lags = 0:(nSplits-1);
    utLag = cell(1,nSplits);
    lagMean = nan(1,nSplits);
    for L = lags
        idx = diag(true(nSplits-L,1), L) | diag(true(nSplits-L,1), -L);
        utLag{L+1} = find(idx);
        lagMean(L+1) = mean(Cpool(idx),'omitnan');
    end
    stats.lag.mean = lagMean;
    [stats.lag.rho, stats.lag.p_rho] = corr((lags(:)), lagMean(:), 'Type','Spearman','Rows','complete');
    halfTarget = lagMean(1)*0.5;
    ix = find(lagMean <= halfTarget, 1, 'first');
    stats.lag.halfwidth_bins = tern(isempty(ix), NaN, ix-1);

    % ---------- Lag null via per-cell spatially-permuted splits ----------
    % Shuffle each cell's split vector independently (destroys adjacency, preserves marginals)
    useCos = strcmpi(simStr,'cosine');
    if nCells >= 3
        nSh = 500;                       % mirrors spatial code
        lagNull = nan(nSh, nSplits);
        for b = 1:nSh
            Cb = simFromM( shuffle_full_per_cell(Mpool), useCos, simStr );
            for L = lags
                lagNull(b, L+1) = mean(Cb(utLag{L+1}), 'omitnan');
            end
        end
        stats.lag.null_mean = mean(lagNull,1,'omitnan');
        stats.lag.null_lo   = prctile(lagNull,  2.5, 1);
        stats.lag.null_hi   = prctile(lagNull, 97.5, 1);
    else
        stats.lag.null_mean = nan(1,nSplits);
        stats.lag.null_lo   = nan(1,nSplits);
        stats.lag.null_hi   = nan(1,nSplits);
    end

    % ---------- Within vs Between epoch ----------
    % Build epoch labels from centers of time bins
    tEdges = linspace(winSecs(1), winSecs(2), nSplits+1);
    tCtr   = (tEdges(1:end-1)+tEdges(2:end))/2;
    epID   = zeros(1,nSplits);
    E = epochEdges(:)';   % relative to WinSecs
    for k = 1:nSplits
        t = tCtr(k) - winSecs(1);  % convert to 0..(t1-t0)
        if     t>=E(1) && t<E(2), epID(k)=1;  % CS
        elseif t>=E(2) && t<E(3), epID(k)=2;  % Trace
        elseif t>=E(3) && t<E(4), epID(k)=3;  % US
        else,  epID(k)=4;                      % Post
        end
    end
    ut = triu(true(nSplits),1);
    W = []; B = [];
    for i=1:nSplits
        for j=i+1:nSplits
            if epID(i)==epID(j), W(end+1)=Cpool(i,j); %#ok<AGROW>
            else,                 B(end+1)=Cpool(i,j); %#ok<AGROW>
            end
        end
    end
    W = W(isfinite(W)); B = B(isfinite(B));
    stats.epoch.within_mean  = mean(W,'omitnan');
    stats.epoch.between_mean = mean(B,'omitnan');
    stats.epoch.diff         = stats.epoch.within_mean - stats.epoch.between_mean;

    % unpaired (Welch) + Wilcoxon on pair distributions
    try
        [~,p_t] = ttest2(W, B);
    catch, p_t = NaN; end
    try
        p_w = ranksum(W, B);
    catch, p_w = NaN; end
    stats.epoch.p_ttest   = p_t;
    stats.epoch.p_wilcox  = p_w;

    % permutation p (shuffle epoch labels)
    diff_perm = nan(nPerm,1);
    ord = 1:nSplits;
    for ppp = 1:nPerm
        ep = epID(ord(randperm(nSplits)));
        Wp = []; Bp = [];
        for i=1:nSplits
            for j=i+1:nSplits
                if ep(i)==ep(j), Wp(end+1)=Cpool(i,j); else, Bp(end+1)=Cpool(i,j); end %#ok<AGROW>
            end
        end
        diff_perm(ppp) = mean(Wp,'omitnan') - mean(Bp,'omitnan');
    end
    stats.epoch.p_perm = (1 + sum(diff_perm >= stats.epoch.diff)) / (nnz(isfinite(diff_perm)) + 1);
end

% ---------- RSA + Mantel (add back) ----------
try
    Smore = analyzePVSimStats(R, epochEdges, nBoot, nPerm);
    if isfield(Smore,'mantel'), stats.mantel = Smore.mantel; end
    if isfield(Smore,'rsa'),    stats.rsa    = Smore.rsa;    end
catch ME
    warning('analyzePVSimStats failed: %s', ME.message);
end



% =====================  FIGURE (match reference layout)  =================
% Top row: up to 3 per-day heatmaps + pooled; bottom: bars (span 2) + lag (span 2)
dayIdx = find(haveDay);
if numel(dayIdx) > 3, dayIdx = dayIdx(1:3); end

fig = figure('Color','w','Position',[100 100 1400 700]);
tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');

% --- top: per-day heatmaps (up to 3) ---
for j = 1:3
    nexttile(j); cla
    if j <= numel(dayIdx)
        d = dayIdx(j);
        C = R.perDay(d).simMat;
        %clim = strcmp(simStr,'cosine')*[0 1] + ~strcmp(simStr,'cosine')*[-1 1];
        clim = [-.7, 1];
        imagesc(C, clim); axis square; colormap(gca, parula); colorbar;
        xticks(1:nSplits); yticks(1:nSplits); xlabel('Split'); ylabel('Split');
        if strcmp(simStr,'cosine')
            ttlExtra = sprintf('mean off cos = %.2f (%.0f°)', meanOff(d), rad2deg(real(acos(max(-1,min(1,meanOff(d)))))));
        else
            ttlExtra = sprintf('mean off %s = %.2f', simStr, meanOff(d));
        end
        title(sprintf('%s  |  %s', R.days{d}, ttlExtra));
    else
        axis off
        text(0.5,0.5,'(no day)', 'HorizontalAlignment','center', ...
             'VerticalAlignment','middle','FontSize',12);
    end
end

% --- top: pooled heatmap ---
nexttile(4); cla
clim = strcmp(simStr,'cosine')*[0 1] + ~strcmp(simStr,'cosine')*[-1 1];
imagesc(Cpool, clim); axis square; colormap(gca, parula); colorbar;
xticks(1:nSplits); yticks(1:nSplits); xlabel('Split'); ylabel('Split');
if strcmp(simStr,'cosine')
    ttlExtra = sprintf('Pooled | mean off cos=%.2f (%.0f°)', meanOff_pool, rad2deg(meanAng_pool));
else
    ttlExtra = sprintf('Pooled | mean off %s=%.2f', simStr, meanOff_pool);
end
title(ttlExtra);

% --- bottom-left: Within vs Between bars (span 2) ---
nexttile(5,[1 2]); cla; hold on
if doStats && isfield(stats,'epoch') && isfinite(stats.epoch.within_mean) && isfinite(stats.epoch.between_mean)
    mW = stats.epoch.within_mean; mB = stats.epoch.between_mean;

    % compute SEM from distributions if available
    sem = @(v) std(v,0,'omitnan')/sqrt(max(1,nnz(isfinite(v))));
    if exist('W','var') && exist('B','var') && ~isempty(W) && ~isempty(B)
        seW = sem(W); seB = sem(B);
    else
        seW = NaN; seB = NaN;
    end

    bar([1 2],[mW mB],0.65,'FaceColor',[0.80 0.85 1.00],'EdgeColor','k');
    errorbar([1 2],[mW mB],[seW seB],'k','LineStyle','none','LineWidth',1.3);

    xticks([1 2]); xticklabels({'Within','Between'});
    ylabel(sprintf('Mean %s similarity', simStr)); box on

    yl = ylim; yb = yl(2) - 0.06*range(yl);
    plot([1 1 2 2],[yb-0.01 yb yb yb-0.01],'k-','LineWidth',1.25)
    txt = sprintf('%s (t p=%.2g; rank p=%.2g; perm p=%.2g)', ...
        sigStars(stats.epoch.p_ttest), stats.epoch.p_ttest, stats.epoch.p_wilcox, stats.epoch.p_perm);
    text(1.5, yb+0.01*range(yl), txt, 'HorizontalAlignment','center', 'FontSize',12, 'FontWeight','bold');

    title('Within vs Between epoch (pooled pairs)');
else
    axis off
    text(0.5,0.5,'No epoch stats available','HorizontalAlignment','center','VerticalAlignment','middle');
end

% --- bottom-right: lag curve with shuffled band (span 2) ---
nexttile(7,[1 2]); cla; hold on
if doStats && isfield(stats,'lag') && any(isfinite(stats.lag.mean))
    x = 0:(nSplits-1);
    % null band
    if isfield(stats.lag,'null_lo') && any(isfinite(stats.lag.null_lo))
        xx = [x fliplr(x)];
        yy = [stats.lag.null_lo fliplr(stats.lag.null_hi)];
        fill(xx, yy, [0.8 0.85 1], 'EdgeColor','none', 'FaceAlpha',0.45);
    end
    if isfield(stats.lag,'null_mean') && any(isfinite(stats.lag.null_mean))
        plot(x, stats.lag.null_mean, '--', 'LineWidth', 1.4);
    end
    % observed
    plot(x, stats.lag.mean, 'k-', 'LineWidth', 1.8);
    xlabel('Lag (bins)'); ylabel(sprintf('Mean %s similarity', simStr));
    title(sprintf('Lag curve (\\rho=%.2f, p=%.1g; half-width=%s bins)', ...
          stats.lag.rho, stats.lag.p_rho, num2str(stats.lag.halfwidth_bins)));
    box off
    % legend
    lg = {}; lh = [];
    if any(isfinite(stats.lag.null_lo))
        lg{end+1} = 'Shuffled 95% band';  lh(end+1)=plot(nan,nan,'-','Color',[0.8 0.85 1],'LineWidth',8); %#ok<AGROW>
    end
    if any(isfinite(stats.lag.null_mean))
        lg{end+1} = 'Shuffled mean';      lh(end+1)=plot(nan,nan,'--k','LineWidth',1.4); %#ok<AGROW>
    end
    lg{end+1} = 'Observed';               lh(end+1)=plot(nan,nan,'k-','LineWidth',1.8); %#ok<AGROW>
    if ~isempty(lh), legend(lh, lg, 'Location','northeast'); end
else
    axis off
    text(0.5,0.5,'No lag stats available','HorizontalAlignment','center','VerticalAlignment','middle');
end

% Page title (match style)
sgtitle(fig, sprintf('%s | %dsplits | %0.1f–%0.1fs | %s | %s', ...
       ratName, nSplits, winSecs(1), winSecs(2), upper(modeStr), upper(simStr)), 'FontWeight','bold');

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

% ---------------- local helpers ----------------
function C = simFromM(Min, useCos, simStrLocal)
    if useCos
        rn = vecnorm(Min,2,2);
        C  = (Min*Min') ./ max(rn*rn', eps);
        C  = max(-1,min(1,C));
    else
        C  = corr(Min','rows','pairwise','Type',upper(simStrLocal));
    end
    C = (C + C')/2; C(1:size(C,1)+1:end) = 1;
end

function Ms = shuffle_full_per_cell(M)
    % Independently permute the S-bin vector for each cell (columns).
    [S,N] = size(M);
    Ms = zeros(S,N);
    for j=1:N
        idx = randperm(S);
        Ms(:,j) = M(idx,j);
    end
end

function out = tern(cond, a, b)
    if cond, out=a; else, out=b; end
end

function s = sigStars(pval)
  if ~isfinite(pval), s = 'n.s.'; return; end
  if pval < 1e-4, s='****';
  elseif pval < 1e-3, s='***';
  elseif pval < 1e-2, s='**';
  elseif pval < 0.05, s='*';
  else, s='n.s.'; end
end

end

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
getSim = @(M) simFromM_local(M, useCosine, lower(R.params.Similarity));


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

function C = simFromM_local(Min, useCos, simStrLocal)
% Build split×split similarity from split×cell matrix Min
if useCos
    rn = vecnorm(Min,2,2);
    C  = (Min*Min') ./ max(rn*rn', eps);
    C  = max(-1,min(1,C));
else
    C  = corr(Min','rows','pairwise','Type',upper(simStrLocal));  % 'pearson' or 'spearman'
end
C = (C + C')/2; C(1:size(C,1)+1:end) = 1;
end
