function Rall = spatialPV_timeblocks_evenOdd(ratNames, varargin)
% spatialPV_timeblocks_evenOdd
% Split-half spatial PV correlation for running, excluding [CS, CS+2].
% New: SplitMode 'perbin' to avoid striping from uneven A/B coverage.
%
% Example:
% R = spatialPV_timeblocks_evenOdd({'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%        'BlockLen',2,'NBins',16,'GridRC',[4 4],'BinMode','equal_occ', ...
%        'Mode','demean','Similarity','corr','MinSpeed',4, ...
%        'SplitMode','perbin','PhaseRepeats',2,'DoPlot',true);

% ---------- args ----------
p = inputParser;
addParameter(p,'Days','An-2:An');
addParameter(p,'NBins',15,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'GridRC',[],@(v) isempty(v) || (isnumeric(v)&&numel(v)==2&&all(v>=1)));
addParameter(p,'BinMode','equal_occ',@(s) any(strcmpi(s,{'equal_occ','equal_size'})));
addParameter(p,'MinSpeed',4,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Mode','demean',@(s) any(strcmpi(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','corr',@(s) any(strcmpi(s,{'corr','pearson','cosine'})));
addParameter(p,'BlockLen',2,@(x) isnumeric(x)&&isscalar(x)&&x>0);   % s (only for SplitMode=blocks)
addParameter(p,'PhaseRepeats',1,@(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'Pairwise',true,@islogical);                        % corr Rows='pairwise' vs 'complete'
addParameter(p,'NShuff',100,@(x) isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'MinOccPerBin',2,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'MinCellsPerCorr',10,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'CellMinBinsFrac',0.30,@(x) isnumeric(x)&&isscalar(x)&&x>0&&x<=1); % per-cell stability
addParameter(p,'MinCellsBinFrac',0.10,@(x) isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'DoPlot',true,@islogical);
addParameter(p,'Debug',false,@islogical);
addParameter(p,'SaveFig',false,@islogical);
addParameter(p,'FigName','',@(s) ischar(s) || isstring(s));
addParameter(p,'SubsampleN',0,@(x) isnumeric(x)&&isscalar(x)&&x>=0);  % 0 = off
addParameter(p,'SplitMode','withinbin_random', ...   % 'withinbin_evenodd' | 'withinbin_random' | 'blocks'
    @(s) any(strcmpi(s,{'withinbin_evenodd','withinbin_random','blocks'})));
addParameter(p,'NormAcrossCells',false,@islogical);   % extra row-wise (bin-wise) zscore

parse(p,varargin{:});
opt = p.Results;

% normalize ratNames -> cellstr
if ischar(ratNames) || isstring(ratNames)
    ratNames = {char(ratNames)};
elseif ~iscell(ratNames)
    error('ratNames must be char, string, or cell array of these.');
end

% ---------- per-rat loop ----------
perRat = cell(1,numel(ratNames));
pooled_Z_stack = [];
pooled_wts     = [];
binMeta_master = [];

if opt.DoPlot && numel(ratNames) > 1
    figure('Color','w','Position',[100 100 1200 max(420, 260*numel(ratNames))]);
end

for r = 1:numel(ratNames)
    rn = ratNames{r};
    Ri = run_one_rat(rn, opt);
    perRat{r} = Ri;

    if ~isempty(Ri)
        if opt.DoPlot && numel(ratNames) > 1
            row = r;
            subplot(numel(ratNames),2,2*(row-1)+1);
            clim = chooseClim(opt.Similarity);
            imagesc(Ri.pooled.C, clim); axis image; colormap(gca, parula); colorbar;
            set(gca,'XTick',1:size(Ri.pooled.C,1),'YTick',1:size(Ri.pooled.C,1));
            xlabel('B bins'); ylabel('A bins');
            title(sprintf('%s: PV corr (A vs B)', rn), 'Interpreter','none');

            subplot(numel(ratNames),2,2*(row-1)+2); hold on
            plot(Ri.pooled.lag.centers, Ri.pooled.lag.mean, 'k-', 'LineWidth', 1.8);
            xlabel('Bin distance'); ylabel(yLabel(opt.Similarity));
            title(sprintf('%s: similarity vs distance', rn), 'Interpreter','none');
            grid on
        end

        C = Ri.pooled.C;
        Z = atanh(max(min(C, 0.999999), -0.999999));
        pooled_Z_stack = cat(3, pooled_Z_stack, Z);
        pooled_wts     = [pooled_wts; sqrt(max(eps, Ri.pooled.totalIncludedSec))]; %#ok<AGROW>
        if isempty(binMeta_master) && isfield(Ri.pooled,'binMeta')
            binMeta_master = Ri.pooled.binMeta;
        end
    end
end

% ---------- pooled-across-rats ----------
if ~isempty(pooled_Z_stack)
    W  = reshape(pooled_wts,1,1,[]);
    Zg = sum(pooled_Z_stack .* W, 3, 'omitnan') ./ sum(W,3,'omitnan');
    Cg = tanh(Zg);
    lagG = spatialLagCurve(Cg, binMeta_master);   % includes distance=0

    Rall.pooledAcrossRats.C        = Cg;
    Rall.pooledAcrossRats.lag      = lagG;
    Rall.pooledAcrossRats.binMeta  = binMeta_master;

    if opt.DoPlot
        figure('Color','w','Position',[100 100 1050 430]);
        subplot(1,2,1);
        clim = chooseClim(opt.Similarity);
        clim = [-.5,.8];
        imagesc(Cg, clim); axis image; colormap(gca, parula); colorbar;
        xlabel('B bins'); ylabel('A bins');
        title('POOLED across rats: PV corr (A vs B)');
        subplot(1,2,2);
        plot(lagG.centers, lagG.mean, 'k-', 'LineWidth', 2);
        xlabel('Bin distance'); ylabel(yLabel(opt.Similarity));
        title('POOLED: similarity vs distance (starts at 0)'); grid on
    end
end

Rall.perRat = perRat;
Rall.params = opt;

if opt.SaveFig
    fn = string(opt.FigName);
    if strlength(fn)==0, fn = "spatialPV_timeblocks_AB.png"; end
    exportgraphics(gcf, fn, 'Resolution', 300);
end
end % ========================= main =========================


% ========================= per-rat core =========================
function R = run_one_rat(ratName, opt)
rat = evalin('base', ratName);

% ---- resolve days ----
allDays = autoDateList(rat);
if ischar(opt.Days) || isstring(opt.Days)
    if strcmpi(opt.Days,'An-2:An')
        idx = find(strcmp(allDays, rat.An));
        assert(~isempty(idx) && idx>=3, 'Not enough days before An in %s.', ratName);
        daysToUse = allDays(idx-2:idx);
    else
        daysToUse = {char(opt.Days)};
    end
else
    daysToUse = opt.Days;
end
daysToUse = daysToUse(:)';

% ---- PASS 1: running timebase (exclude trials) & bins ----
[perDayRun, all_x, all_y] = buildRunningTimebase(rat, daysToUse, opt.MinSpeed);
assert(~isempty(perDayRun) && any([perDayRun.valid]), 'No valid running segments.');
binMeta = defineBins(all_x, all_y, opt.NBins, opt.GridRC, lower(opt.BinMode));

% ---- per-day: compute split-half C ----
Cs   = [];                 % Fisher-z stack across days
wts  = [];                 % per-day weights (sqrt seconds)
diagPs_all = [];           % per-bin p's pooled later

for d = 1:numel(perDayRun)
    if ~perDayRun(d).valid, continue; end

    day = perDayRun(d).day;
    vt  = perDayRun(d).vt;
    xv  = perDayRun(d).xv;
    yv  = perDayRun(d).yv;
    dt  = perDayRun(d).dt;

    % Map samples→bins ONCE for THIS day
    binIdxVT_all = mapSamplesToBins(vt, xv, yv, binMeta, opt.GridRC);

    Z_ph = [];     % per-phase Fisher-z
    for ph = 1:opt.PhaseRepeats
        % ---- A/B labels (within-bin to kill striping) ----
        switch lower(opt.SplitMode)
            case 'withinbin_evenodd'
                labVT = zeros(size(vt));
                for b = 1:prod(binMeta.gridRC)
                    idx = find(binIdxVT_all==b);
                    if numel(idx) < 2, continue; end
                    labVT(idx(1:2:end)) = 1;
                    labVT(idx(2:2:end)) = 2;
                end
              case 'withinbin_random'
                  rng(1000+ph);
                  labVT = zeros(size(vt));
                  for b = 1:prod(binMeta.gridRC)
                      idx = find(binIdxVT_all==b);
                      if numel(idx) < 2, continue; end
                      % weights = dwell per sample
                      w = dt(idx);
                      % greedily assign samples to equalize total dwell in A and B
                      [ia, ib] = split_balanced_by_dt(idx, w);
                      labVT(ia) = 1;          % A
                      labVT(ib) = 2;          % B
                  end

            otherwise % 'blocks' (legacy)
                phaseShift = (ph-1)/max(1,opt.PhaseRepeats)*opt.BlockLen;
                blocks = makeAlternatingBlocks(vt, perDayRun(d).segments, opt.BlockLen, phaseShift);
                labVT = labelTimeSamples(vt, blocks);
        end

        % Occupancy for this labeling
        [occA, occB] = occupancyByLabels(binIdxVT_all, dt, labVT, prod(binMeta.gridRC));

        if opt.Debug
            figure('Color','w','Position',[300 300 900 300]);
            subplot(1,3,1); imagesc(reshape(occA, binMeta.gridRC)); axis image; colorbar; title('A dwell (s)');
            subplot(1,3,2); imagesc(reshape(occB, binMeta.gridRC)); axis image; colorbar; title('B dwell (s)');
            subplot(1,3,3); imagesc(reshape(occA-occB, binMeta.gridRC)); axis image; colorbar; title('A - B dwell');
        end


        % Spikes & rates
        Sfld  = sprintf('CA_peaks_%s', day);
        if ~isfield(rat.Ca_peaks, Sfld), warning('Missing %s; skip day.', Sfld); continue; end
        Spk = rat.Ca_peaks.(Sfld);

        %%% NEW: apply ratemask==1 (filter cells before rate matrices) %%%
        Mfld = sprintf('ratemask_%s', day);
        if isfield(rat,'ratemask') && isstruct(rat.ratemask) && isfield(rat.ratemask, Mfld)
            rm = rat.ratemask.(Mfld) == 1;
            rm = rm(:);
            if iscell(Spk)
                if numel(rm) == numel(Spk)
                    Spk = Spk(rm);
                end
            else
                if size(Spk,1) == numel(rm)
                    Spk = Spk(rm,:);
                end
            end
        end
        %%% END NEW %%%

        csTimes = getCSTimes(rat, day);
        [rateA, rateB] = rateMatrices_AB(Spk, vt, xv, yv, labVT, csTimes, binMeta, opt.GridRC, occA, occB);

        % Gating: occupancy
        badA = ~isfinite(occA) | occA < opt.MinOccPerBin;
        badB = ~isfinite(occB) | occB < opt.MinOccPerBin;
        rateA(badA,:) = NaN;  rateB(badB,:) = NaN;

        % Gating: per-cell stability across BOTH splits
        K = size(rateA,1);
        minBinsPerCell = max(3, ceil(opt.CellMinBinsFrac * K));
        keep = (sum(isfinite(rateA),1) >= minBinsPerCell) & (sum(isfinite(rateB),1) >= minBinsPerCell);
        rateA = rateA(:, keep);  rateB = rateB(:, keep);
        Nc = size(rateA,2); if Nc < opt.MinCellsPerCorr, continue; end

        % Gating: per-bin min cells
        minCellsBin = max(opt.MinCellsPerCorr, ceil(opt.MinCellsBinFrac * Nc));
        rowOK = sum(isfinite(rateA),2) >= minCellsBin;
        colOK = sum(isfinite(rateB),2) >= minCellsBin;
        rateA(~rowOK,:) = NaN;  rateB(~colOK,:) = NaN;

        % Normalize & correlate
        % per-cell normalization across bins (you already have this)
        rateA = normalizeRatesPerCell(rateA, lower(opt.Mode));
        rateB = normalizeRatesPerCell(rateB, lower(opt.Mode));

        % NEW: per-bin normalization across cells (row-wise)
        if opt.NormAcrossCells
            rateA = row_zscore(rateA);   % zscore each bin across cells
            rateB = row_zscore(rateB);
        end

        if ~opt.Pairwise
            keepTest = all(isfinite([rateA; rateB]),1);
            if nnz(keepTest) < max(3, minCellsBin)
                fprintf('Listwise impossible this day (kept cells = %d); using pairwise.\n', nnz(keepTest));
            end
        end

        % --- Freeze one complete-case cell list for the whole matrix
        keepList = all(isfinite([rateA; rateB]), 1);
        nKeep    = nnz(keepList);
        if nKeep < max(3, opt.SubsampleN)
            % too few cells this day; skip it
            continue
        end
        rateA = rateA(:, keepList);
        rateB = rateB(:, keepList);

        % (now Pairwise vs. Complete is irrelevant; every (i,j) sees the same cells)

        % similarity
        C = crossSimilarity_AB(rateA, rateB, lower(opt.Similarity), ...
                               minCellsBin, opt.Pairwise, opt.SubsampleN);

        Z_ph = cat(3, Z_ph, atanh(max(min(C,0.999999), -0.999999)));

        % Shuffle p’s once per day (after last phase)
        if ph==opt.PhaseRepeats
            [p_diag, ~] = diagP_byShuffle(rateA, rateB, binMeta, opt, ...
                vt, xv, yv, dt, labVT, Spk, csTimes, occA, occB, ...
                minCellsBin, opt.SubsampleN, binIdxVT_all);
        end
    end

    if isempty(Z_ph), continue; end
    Z_day = mean(Z_ph, 3, 'omitnan');
    C_day = tanh(Z_day);

    Cs   = cat(3, Cs, atanh(max(min(C_day,0.999999),-0.999999)));
    wts  = [wts; sqrt(sum(dt,'omitnan'))];
    diagPs_all = [diagPs_all; p_diag(:)'];
end

if isempty(Cs)
    R = [];
    return
end

% ---- pooled within rat ----
W  = reshape(wts,1,1,[]);
Zp = sum(Cs .* W, 3, 'omitnan') ./ sum(W,3,'omitnan');
Cp = tanh(Zp);

diag_r = diag(Cp).';

% combine per-bin p across days (Stouffer)
pMat = diagPs_all;  pMat(isnan(pMat)) = 1;
Zs   = -norminv(pMat/2);
Wrep = repmat(wts,1,size(pMat,2));
Zcomb= sum(Wrep .* Zs, 1) ./ sqrt(sum(Wrep.^2,1));
p_pooled = 2*normcdf(-abs(Zcomb));

% FDR across bins & lag curve
mask_fdr = fdr_mask_local(p_pooled, 0.05);
lag = spatialLagCurve(Cp, binMeta);

R.perDay = [];
R.pooled.C                = Cp;
R.pooled.diag_r           = diag_r;
R.pooled.p_diag           = p_pooled;
R.pooled.mask_fdr         = mask_fdr;
R.pooled.lag              = lag;
R.pooled.binMeta          = binMeta;
R.pooled.totalIncludedSec = sum(wts.^2);
R.params = opt;  R.params.Days = daysToUse;
end


% ========================= helpers =========================
function [perDayRun, all_x, all_y] = buildRunningTimebase(rat, days, vMin)
perDayRun = struct('day',[],'vt',[],'xv',[],'yv',[],'dt',[],'segments',[],'valid',false);
all_x=[]; all_y=[];
for d=1:numel(days)
    D = days{d};
    Pf = sprintf('pos_%s', D);
    if ~isfield(rat.pos, Pf), continue; end
    pos = rat.pos.(Pf);
    vel = ca_velocity(pos);
    pos = local_smoothpos(pos);
    vt  = vel(2,:)';  vmag = vel(1,:)';
    xv  = interp1(pos(:,1), pos(:,2), vt, 'linear', NaN);
    yv  = interp1(pos(:,1), pos(:,3), vt, 'linear', NaN);

    csTimes = getCSTimes(rat, D);
    inTrial = false(size(vt));
    for i=1:numel(csTimes)
        inTrial = inTrial | (vt >= csTimes(i) & vt < csTimes(i)+2);
    end

    keep = (vmag >= vMin) & ~inTrial & isfinite(xv) & isfinite(yv);
    vt = vt(keep); xv = xv(keep); yv = yv(keep);
    if numel(vt) < 10, continue; end

    dt = diff(vt); dt = [dt; median(dt,'omitnan')];

    brk = find(diff(vt) > 0.5);
    segs = [[vt([1; brk+1]) vt([brk; end])]];

    perDayRun(d).day      = D;
    perDayRun(d).vt       = vt;
    perDayRun(d).xv       = xv;
    perDayRun(d).yv       = yv;
    perDayRun(d).dt       = dt;
    perDayRun(d).segments = segs;
    perDayRun(d).valid    = true;

    all_x = [all_x; xv]; %#ok<AGROW>
    all_y = [all_y; yv]; %#ok<AGROW>
end
end

function binMeta = defineBins(all_x, all_y, K, GridRC, binMode)
binMeta = struct('mode','', 'gridRC',[], 'edges1D',[], 'centers1D',[], ...
                 'edgesX',[], 'edgesY',[], 'centers2D',[], ...
                 'xyMin',[],'xyMax',[], 'coeffPC1',[],'muXY',[]);
if ~isempty(GridRC)
    Rg = GridRC(1); Cg = GridRC(2);
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
    [CX,CY] = meshgrid(cx, cy);
    binMeta.mode='grid'; binMeta.gridRC=[Rg Cg];
    binMeta.edgesX=ex; binMeta.edgesY=ey; binMeta.centers2D=[CX(:), CY(:)];
    binMeta.xyMin=[min(all_x) min(all_y)]; binMeta.xyMax=[max(all_x) max(all_y)];
else
    X = [all_x, all_y];
    X = X - mean(X,1,'omitnan');
    [coeff, score] = pca(X, 'NumComponents',1, 'Centered',false);
    s  = score(:,1);
    ex = linspace(min(s), max(s), K+1);
    cx = (ex(1:end-1)+ex(2:end))/2;
    binMeta.mode='pc1'; binMeta.gridRC=[1 K];
    binMeta.edges1D=ex; binMeta.centers1D=cx;
    binMeta.xyMin=[min(all_x) min(all_y)]; binMeta.xyMax=[max(all_x) max(all_y)];
    binMeta.coeffPC1=coeff(:,1); binMeta.muXY=mean([all_x, all_y],1,'omitnan');
end
end

function blocks = makeAlternatingBlocks(vt, segs, L, phaseShift)
% Nx3: [start end label] with label=1(A) or 2(B)
blocks = [];
for s = 1:size(segs,1)
    a = segs(s,1) + phaseShift; b = segs(s,2);
    if a>=b, continue; end
    k = 0; t0 = a;
    while t0 < b
        t1 = min(t0 + L, b);
        lab = 1 + mod(k,2);   % 1,2,1,2,...
        blocks = [blocks; [t0 t1 lab]]; %#ok<AGROW>
        k = k+1; t0 = t1;
    end
end
end

function labVT = labelTimeSamples(vt, blocks)
labVT = zeros(size(vt));
for i=1:size(blocks,1)
    m = vt>=blocks(i,1) & vt<blocks(i,2);
    labVT(m) = blocks(i,3);
end
end

function [occA, occB] = occupancyByLabels(binIdxVT, dt, labVT, K)
occA = accumarray_fill(binIdxVT(labVT==1), dt(labVT==1), [K 1], @sum);
occB = accumarray_fill(binIdxVT(labVT==2), dt(labVT==2), [K 1], @sum);
occA(occA<=0) = NaN;  occB(occB<=0) = NaN;
end


function [occA, occB, binIdxVT] = occupancyBySplit(vt, xv, yv, dt, labVT, binMeta, GridRC)
K = prod(binMeta.gridRC);
if strcmpi(binMeta.mode,'grid')
    bx = discretize(xv, binMeta.edgesX);
    by = discretize(yv, binMeta.edgesY);
    valid = isfinite(bx) & isfinite(by);
    binIdxVT = nan(size(vt));
    binIdxVT(valid) = sub2ind([GridRC(1), GridRC(2)], by(valid), bx(valid));
else
    XY = [xv, yv]; XYc = XY - binMeta.muXY;
    s  = XYc * binMeta.coeffPC1(:);
    b1 = discretize(s, binMeta.edges1D);
    binIdxVT = b1;  K = binMeta.gridRC(2);
end
occA = accumarray_fill(binIdxVT(labVT==1), dt(labVT==1), [K 1], @sum);
occB = accumarray_fill(binIdxVT(labVT==2), dt(labVT==2), [K 1], @sum);
occA(occA<=0) = NaN; occB(occB<=0) = NaN;
end

function [rateA, rateB] = rateMatrices_AB(Spk, vt, xv, yv, labVT, csTimes, binMeta, GridRC, occA, occB)
% K×N rate matrices for A/B splits
if iscell(Spk), N = numel(Spk); else, N = size(Spk,1); end
K = prod(binMeta.gridRC);
rateA = nan(K, N); rateB = nan(K, N);

if strcmpi(binMeta.mode,'grid')
    mapXYtoBin = @(xs,ys) sub2ind(GridRC, ...
        discretize(ys, binMeta.edgesY), discretize(xs, binMeta.edgesX));
else
    mapXYtoBin = @(xs,ys) discretize(([xs ys]-binMeta.muXY) * binMeta.coeffPC1(:), binMeta.edges1D);
end

for c = 1:N
    st = getCellSpikes_local(Spk, c);
    if isempty(st), continue; end
    for i=1:numel(csTimes)    % drop trial spikes
        st = st( ~(st>=csTimes(i) & st<csTimes(i)+2) );
    end
    if isempty(st), continue; end

    xs = interp1(vt, xv, st, 'linear', NaN);
    ys = interp1(vt, yv, st, 'linear', NaN);
    labS = interp1(vt, labVT, st, 'nearest','extrap');
    bSpk = mapXYtoBin(xs, ys);

    okA = isfinite(bSpk) & labS==1;
    okB = isfinite(bSpk) & labS==2;

    if any(okA)
        nA = accumarray_fill(bSpk(okA), 1, [K 1], @sum);
        rateA(:,c) = nA ./ occA;
    end
    if any(okB)
        nB = accumarray_fill(bSpk(okB), 1, [K 1], @sum);
        rateB(:,c) = nB ./ occB;
    end
end
end

function [C, nShared] = crossSimilarity_AB(XA, XB, simStr, minShared, pairwise, subsampleN)
% XA, XB: K×N (bins×cells) with NaNs for missing bin/cell entries
% simStr : 'corr'/'pearson' or 'cosine'
% minShared : minimum # of shared cells required for a pair (i,j)
% pairwise  : true=per-pair intersection, false=complete-case listwise
% subsampleN: if >0, compute each (i,j) using exactly subsampleN cells from A∩B

if nargin < 5 || isempty(pairwise),   pairwise   = true; end
if nargin < 6 || isempty(subsampleN), subsampleN = 0;    end

K   = size(XA,1);
Aok = isfinite(XA); Bok = isfinite(XB);
nShared = double(Aok) * double(Bok)';

C = nan(K,K,'like',XA);

% ---- Equal-N subsampling (removes residual striping) ----
if subsampleN > 0
    need   = max(minShared, subsampleN);
    useCos = any(strcmpi(simStr,{'cosine'}));
    for i = 1:K
        ai = isfinite(XA(i,:));
        aiVals = XA(i,:);
        for j = 1:K
            bj  = isfinite(XB(j,:));
            m   = ai & bj;
            ns  = nnz(m);
            if ns < need, C(i,j) = NaN; continue; end
            idx  = find(m);
            pick = idx(randperm(ns, subsampleN));      % fixed N without replacement
            a = aiVals(pick);
            b = XB(j,pick);
            if useCos
                C(i,j) = (a*b') / max(norm(a)*norm(b), eps);
                C(i,j) = max(-1,min(1,C(i,j)));
            else
                C(i,j) = corr(a', b', 'Rows','complete', 'Type','Pearson');
            end
        end
    end
    C(nShared < minShared) = NaN;
    return
end

% ---- Standard paths (no subsampling) ----
if any(strcmpi(simStr,{'cosine'}))
    rnA = vecnorm(XA,2,2); rnB = vecnorm(XB,2,2);
    C   = (XA*XB.') ./ max(rnA,eps) ./ max(rnB',eps);
    C   = max(-1,min(1,C));
else
    if pairwise
        C = corr(XA', XB', 'Rows','pairwise');   % per-pair intersection of cells
    else
        % listwise: single complete cell set for whole matrix (also removes stripes)
        keep = all(isfinite([XA; XB]), 1);
        if nnz(keep) >= max(3, minShared)
            C = corr(XA(:,keep)', XB(:,keep)', 'Rows','complete');
        else
            C = corr(XA', XB', 'Rows','pairwise'); % fallback if too few cells
        end
    end
end

C(nShared < minShared) = NaN;
end


function [p_diag, diag_null] = diagP_byShuffle(rateA, rateB, binMeta, opt, ...
    vt, xv, yv, dt, labVT, Spk, csTimes, occA0, occB0, ...
    minCellsBin, subsampleN, binIdxVT_all)

K = size(rateA,1);
diag_obs  = diag(crossSimilarity_AB(rateA, rateB, lower(opt.Similarity), minCellsBin, true)).';
diag_null = nan(opt.NShuff, K, 'single');

for s = 1:opt.NShuff
    % --- within-bin permutation that preserves A/B counts in each bin ---
    labSh = labVT;
    for b = 1:prod(binMeta.gridRC)
        ix = find(binIdxVT_all==b & labVT>0);
        if numel(ix) < 2, continue; end
        ia = ix(labVT(ix)==1); ib = ix(labVT(ix)==2);
        allix = [ia(:); ib(:)];
        perm  = allix(randperm(numel(allix)));
        labSh(perm(1:numel(ia)))         = 1;
        labSh(perm(numel(ia)+1:end))     = 2;
    end

    [occA_s, occB_s] = occupancyByLabels(binIdxVT_all, dt, labSh, K);
    [rateA_s, rateB_s] = rateMatrices_AB(Spk, vt, xv, yv, labSh, csTimes, binMeta, binMeta.gridRC, occA_s, occB_s);
    rateA_s = normalizeRatesPerCell(rateA_s, lower(opt.Mode));
    rateB_s = normalizeRatesPerCell(rateB_s, lower(opt.Mode));

    diag_null(s,:) = diag( ...
        crossSimilarity_AB(rateA_s, rateB_s, lower(opt.Similarity), ...
                           minCellsBin, opt.Pairwise, subsampleN) ).';
end


% p-values
p_diag = nan(1,K);
for k = 1:K
    r0 = diag_obs(k); nullk = diag_null(:,k); nullk = nullk(isfinite(nullk));
    if isnan(r0) || isempty(nullk), p_diag(k) = NaN;
    else, p_diag(k) = max(1/numel(nullk), mean(abs(nullk) >= abs(r0)));
    end
end
end

function out = accumarray_fill(idx, val, sz, fun, fill)
if nargin < 5, fill = NaN; end
idx = idx(:);
if isscalar(val), val = repmat(val, size(idx)); else, val = val(:); end
if isempty(idx), out = repmat(fill, sz); return; end
imax = prod(sz);
ok = isfinite(idx) & idx >= 1 & idx <= imax & isfinite(val);
if ~any(ok), out = repmat(fill, sz); return; end
out = accumarray(idx(ok), val(ok), sz, fun, fill);
end

function labBlk = makeBlockIDs(labVT)
labBlk = zeros(size(labVT)); id=0;
for i=1:numel(labVT)
    if labVT(i)==0, continue; end
    if i==1 || labVT(i-1)==0 || labVT(i-1)~=labVT(i), id = id+1; end
    labBlk(i) = id;
end
end

function st = getCellSpikes_local(S, c)
if iscell(S), st = S{c}(:); else, st = S(c,:).'; end
st = st(~isnan(st) & st>0);
end

function cs = getCSTimes(rat, day)
Cfld = sprintf('CS_%s', day);
if isfield(rat.CS_times, Cfld), cs = rat.CS_times.(Cfld)(:);
else, cs = []; end
end

function X = normalizeRatesPerCell(R, modeStr)
switch lower(modeStr)
    case 'raw',     X = R;
    case 'demean',  mu = mean(R, 1, 'omitnan'); X = R - mu;
    case 'zscore'
        mu = mean(R, 1, 'omitnan');
        sd = std(R, 0, 1, 'omitnan'); sd(~isfinite(sd) | sd==0) = 1;
        X = (R - mu) ./ sd;
    otherwise, error('Unknown Mode: %s', modeStr);
end
end

function lag = spatialLagCurve(C, binMeta)
K = size(C,1); UT = triu(true(K),1);
if strcmpi(binMeta.mode,'grid')
    XY = binMeta.centers2D; D = squareform(pdist(XY));
else
    if isfield(binMeta,'centers1D'), c = binMeta.centers1D(:); else, c = (1:K)'; end
    D = squareform(pdist(c));
end
d = D(UT); s = C(UT);
ok = isfinite(d) & isfinite(s);
d = d(ok); s = s(ok);

% distance=0 (diagonal)
lag.centers = 0;  lag.mean = mean(diag(C),'omitnan');

% 8 quantile bins for the rest
if ~isempty(d)
    edges = unique(quantile(d, linspace(0,1,9)));
    edges(1) = edges(1)-eps;
    for k = 1:numel(edges)-1
        mk = d>=edges(k) & d<edges(k+1);
        lag.centers(end+1) = mean([edges(k) edges(k+1)]); %#ok<AGROW>
        lag.mean(end+1)    = mean(s(mk),'omitnan');       %#ok<AGROW>
    end
end
lag.centers = lag.centers(:)'; lag.mean = lag.mean(:)';
end

function m = fdr_mask_local(p, q)
p = p(:)'; p(isnan(p))=1;
[ps,~] = sort(p); m = false(size(p));
k = find(ps <= (1:numel(ps))/numel(ps)*q, 1, 'last');
if ~isempty(k), thr = ps(k); m = p<=thr; end
end

function lbl = yLabel(simStr)
if any(strcmpi(simStr,{'corr','pearson'})), lbl = 'Mean Pearson similarity';
else, lbl = 'Mean cosine similarity'; end
end

function clim = chooseClim(simStr)
if any(strcmpi(simStr,{'corr','pearson'})), clim = [-1 1]; else, clim = [0 1]; end
end

function pos_s = local_smoothpos(pos)
try, pos_s = smoothpos(pos);
catch
    if size(pos,2) >= 3
        win = 5; pos_s = pos;
        pos_s(:,2) = movmean(pos(:,2),win,'omitnan');
        pos_s(:,3) = movmean(pos(:,3),win,'omitnan');
    else, pos_s = pos;
    end
end
end

function out = tern(cond, a, b)
if cond, out=a; else, out=b; end
end

function binIdxVT = mapSamplesToBins(vt, xv, yv, binMeta, GridRC)
% Map each running-time sample to a spatial bin index (NaN if out)
if strcmpi(binMeta.mode,'grid')
    bx = discretize(xv, binMeta.edgesX);
    by = discretize(yv, binMeta.edgesY);
    ok = isfinite(bx) & isfinite(by);
    binIdxVT = nan(size(vt));
    binIdxVT(ok) = sub2ind(GridRC, by(ok), bx(ok));
else
    XY  = [xv, yv] - binMeta.muXY;
    s   = XY * binMeta.coeffPC1(:);
    binIdxVT = discretize(s, binMeta.edges1D);
end
end

function X = row_zscore(X)
mu = mean(X,2,'omitnan');
sd = std(X,0,2,'omitnan');
sd(~isfinite(sd) | sd==0) = 1;
X = (X - mu) ./ sd;
end

function [Aidx, Bidx] = split_balanced_by_dt(idx, w)
% Partition idx into two sets with ~equal sum(w)
perm = idx(randperm(numel(idx)));           % shuffle
wperm = w(randperm(numel(w)));              % keep same shuffle as idx if you prefer
sA = 0; sB = 0; Aidx = []; Bidx = [];
for k = 1:numel(perm)
    if sA <= sB
        Aidx(end+1) = perm(k); %#ok<AGROW>
        sA = sA + wperm(k);
    else
        Bidx(end+1) = perm(k); %#ok<AGROW>
        sB = sB + wperm(k);
    end
end
end
