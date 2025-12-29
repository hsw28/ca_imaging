function OUT = regressOutSpeed_taskEpochs(ratNames, varargin)
% regressOutSpeed_taskEpochs
%
% Bin-wise, per-cell regression of firing rate on running speed during the
% CS+trace window, grouped by cell "epoch" labels (tone/trace/shock/post/...).
% Uses the same binning code as xcorr_speed_MUA (task mode).
%
% For each rat/day/cell:
%   FR(bin) = beta0 + beta1 * speed(bin) + residual(bin)
%
% OUTPUT (struct OUT):
%   .ratNames          : cell array of rat names
%   .epochNames        : cellstr of epoch labels (E1..Enepoch)
%   .perRat(r) with fields:
%       .days          : date strings used
%       .binT{d}       : [nBins x 1] bin center times
%       .binSpeed{d}   : [nBins x 1] bin speed
%       .rateRaw{d}    : [nBins x nCells] raw/normalized FR (Hz)
%       .rateResid{d}  : [nBins x nCells] FR with speed regressed out
%       .epochID{d}    : [nCells x 1] discrete epoch label (1..nEpoch)
%       .cellR{d}      : [nCells x 1] corr(speed, FR) per cell
%       .cellSlope{d}  : [nCells x 1] beta1 per cell
%       .cellResStd{d} : [nCells x 1] std of residuals per cell
%       .cellResAbsMean{d} : [nCells x 1] mean(|residual|) per cell
%
%   .epochSummary with fields:
%       .meanR_raw(e)      : mean r(speed,FR) across all cells in epoch e
%       .semR_raw(e)
%       .meanSlope(e)      : mean beta1 across cells in epoch e
%       .semSlope(e)
%       .meanResAbs(e)     : mean mean(|residual FR|) across cells in epoch e
%       .semResAbs(e)
%
% Plus plots:
%   1) Histogram of per-cell residual STD (all cells).
%   2) Bar plot of per-epoch mean |residual FR| ± SEM.
%
% Example:
%   OUT = regressOutSpeed_taskEpochs({'rat0222','rat0307',...}, ...
%        'Win',[0 2], 'binSize',1/7.5, ...
%        'EpochVar','epoch','EpochNum',5);

% ---------- options ----------
p = inputParser;
p.addParameter('Win',[0 2]);               % CS window used for task bins
p.addParameter('binSize',1/7.5);           % bin length (s)
p.addParameter('SpeedThresh',[]);          % not used for task here
p.addParameter('PopNorm','none',@(s)any(strcmpi(s,{'none','mean','demean'})));
p.addParameter('EpochVar','epoch',@(s)ischar(s) || isstring(s));
p.addParameter('EpochNum',4,@(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('DoPlot',true,@islogical);  % summary plots

p.parse(varargin{:});
Win       = p.Results.Win;
binSize   = p.Results.binSize;
vThresh   = p.Results.SpeedThresh; %#ok<NASGU> % kept for symmetry; unused
popNorm   = lower(p.Results.PopNorm);
EpochVar  = char(p.Results.EpochVar);
nEpoch    = round(p.Results.EpochNum);
DoPlot    = p.Results.DoPlot;

if nargin < 1 || isempty(ratNames)
    error('Specify ratNames as in xcorr_speed_MUA.');
end
nR = numel(ratNames);

% ---------- top-level OUT ----------
OUT = struct();
OUT.ratNames   = ratNames(:).';
OUT.epochNames = arrayfun(@(k)sprintf('E%d',k),1:nEpoch,'uni',false);

% we build per-rat in a cell array and convert to struct at the end
perRatCell = {};

% Collect across rats/days per epoch
allR          = cell(nEpoch,1);   % r(speed,FR)
allSlope      = cell(nEpoch,1);   % beta1
allResStd     = cell(nEpoch,1);   % std(residual)
allResAbsMean = cell(nEpoch,1);   % mean(|residual|)
allResStd_allCells = [];          % for global histogram (residual)
allRawStd_allCells = [];          % NEW: for global histogram (raw)


% ---------- main loop over rats ----------
for r = 1:nR
    ratName = ratNames{r};
    rat     = evalin('base', ratName);

    dates = autoDateList(rat);
    iAn   = find(strcmp(dates, rat.An), 1);
    if isempty(iAn) || iAn < 3
        warning('%s: not enough days before An; skipping.', ratName);
        continue;
    end
    days = dates(iAn-2:iAn);   % last 3 days up to An

    per = struct();
    per.days          = days;
    per.binT          = cell(numel(days),1);
    per.binSpeed      = cell(numel(days),1);
    per.rateRaw       = cell(numel(days),1);
    per.rateResid     = cell(numel(days),1);
    per.epochID       = cell(numel(days),1);
    per.cellR         = cell(numel(days),1);
    per.cellSlope     = cell(numel(days),1);
    per.cellResStd    = cell(numel(days),1);
    per.cellResAbsMean= cell(numel(days),1);
    per.cellResStd    = cell(numel(days),1);
    per.cellResAbsMean= cell(numel(days),1);
    per.cellRawStd    = cell(numel(days),1);   % NEW


    for d = 1:numel(days)
        D = days{d};

        % ----- spikes as cell array -----
        Sraw   = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));
        S      = normalizeCApeaks_local(Sraw);
        nCells = numel(S);

        % ----- position & speed -----
        pos          = smoothpos(rat.pos.(sprintf('pos_%s',D)));
        [vt, vv]     = velocityFromPos_local(pos);
        tPos         = pos(:,1);
        CS           = rat.CS_times.(sprintf('CS_%s',D));

        % ----- task bins: same as xcorr_speed_MUA, Mode='task' -----
        [binSpeed, rateMat, effDur] = taskBinsSpeedAndRates_withDur_local( ...
            S, CS, tPos, vt, vv, Win, binSize, []);  % no speed gate

        if isempty(binSpeed)
            warning('%s %s: no task bins; skipping day.', ratName, D);
            continue;
        end

        % Bin times ~ midpoints (approx; mainly for bookkeeping)
        nB              = numel(binSpeed);
        per.binT{d}     = (0:nB-1)'*binSize + (CS(1)+Win(1));
        per.binSpeed{d} = binSpeed(:);

        % ----- optional per-cell normalization (mean/demean) -----
        rateMat_use      = normalizePerCell_local(rateMat, popNorm);
        per.rateRaw{d}   = rateMat_use;

        % --- per-cell raw FR std across bins (before regression) ---
        cellRawStd = std(rateMat_use, 0, 1, 'omitnan')';  % [nCells x 1]
        per.cellRawStd{d} = cellRawStd;

        % ----- epoch labels per cell -----
        if ~isfield(rat, EpochVar)
            error('Rat %s is missing field %s for epoch grouping.', ratName, EpochVar);
        end
        GVstruct = rat.(EpochVar);
        if ~isstruct(GVstruct)
            error('%s.%s must be a struct with per-day fields.', ratName, EpochVar);
        end

        fns   = fieldnames(GVstruct);
        idxFn = contains(fns, D);
        if ~any(idxFn)
            error('%s: no field in %s contains date %s.', ratName, EpochVar, D);
        end
        if sum(idxFn) > 1
            warning('%s: multiple %s fields match %s; using first.', ratName, EpochVar, D);
        end
        fname   = fns{find(idxFn,1,'first')};
        epochID = GVstruct.(fname);
        epochID = epochID(:);
        if numel(epochID) < nCells
            error('%s: %s.%s has only %d labels for %d cells.', ...
                  ratName, EpochVar, fname, numel(epochID), nCells);
        end
        epochID = epochID(1:nCells);
        per.epochID{d} = epochID;

        % ----- per-cell regression FR ~ speed -----
        nBins         = numel(binSpeed);
        rateResid     = nan(nBins, nCells);
        cellR         = nan(nCells,1);
        cellSlope     = nan(nCells,1);
        cellResStd    = nan(nCells,1);
        cellResAbsMean= nan(nCells,1);

        x = binSpeed(:);
        for c = 1:nCells
            y = rateMat_use(:,c);
            mask = isfinite(x) & isfinite(y) & effDur>0;
            if nnz(mask) < 5
                continue;
            end
            xx = x(mask);
            yy = y(mask);

            % linear regression: y = beta0 + beta1*x
            X    = [ones(numel(xx),1), xx];
            beta = X \ yy;
            yhat = X * beta;
            resid = yy - yhat;

            rateResid(mask,c)     = resid;
            cellSlope(c)          = beta(2);
            cellR(c)              = corr(xx, yy, 'rows','complete');
            cellResStd(c)         = std(resid, 0, 'omitnan');
            cellResAbsMean(c)     = mean(abs(resid), 'omitnan');
        end

        per.rateResid{d}      = rateResid;
        per.cellR{d}          = cellR;
        per.cellSlope{d}      = cellSlope;
        per.cellResStd{d}     = cellResStd;
        per.cellResAbsMean{d} = cellResAbsMean;

        % collect across rats/days per epoch
        for e = 1:nEpoch
            cells_e = (epochID == e);
            if any(cells_e)
                allR{e}          = [allR{e};          cellR(cells_e)];
                allSlope{e}      = [allSlope{e};      cellSlope(cells_e)];
                allResStd{e}     = [allResStd{e};     cellResStd(cells_e)];
                allResAbsMean{e} = [allResAbsMean{e}; cellResAbsMean(cells_e)];
            end
        end

        % global collection for histogram of residual std
        allResStd_allCells = [allResStd_allCells; cellResStd(:)];
        allRawStd_allCells = [allRawStd_allCells; cellRawStd(:)];

        fprintf('%s %s: %d bins, %d cells; median r(speed,FR)=%.3f\n', ...
            ratName, D, nBins, nCells, median(cellR,'omitnan'));
    end

    % only append if we actually processed days
    perRatCell{end+1} = per; %#ok<AGROW>
end

% Convert per-rat cell array → struct array (or leave empty)
if ~isempty(perRatCell)
    OUT.perRat = [perRatCell{:}];
else
    OUT.perRat = struct([]);
end

% ---------- epoch-wise summary across all rats/days ----------
epochSummary = struct();
epochSummary.meanR_raw   = nan(1,nEpoch);
epochSummary.semR_raw    = nan(1,nEpoch);
epochSummary.meanSlope   = nan(1,nEpoch);
epochSummary.semSlope    = nan(1,nEpoch);
epochSummary.meanResAbs  = nan(1,nEpoch);
epochSummary.semResAbs   = nan(1,nEpoch);

for e = 1:nEpoch
    rvec  = allR{e};
    svec  = allSlope{e};
    rstd  = allResStd{e};      %#ok<NASGU>  % available if needed later
    rabs  = allResAbsMean{e};

    rvec  = rvec(isfinite(rvec));
    svec  = svec(isfinite(svec));
    rabs  = rabs(isfinite(rabs));

    if ~isempty(rvec)
        epochSummary.meanR_raw(e) = mean(rvec);
        epochSummary.semR_raw(e)  = std(rvec) ./ sqrt(numel(rvec));
    end
    if ~isempty(svec)
        epochSummary.meanSlope(e) = mean(svec);
        epochSummary.semSlope(e)  = std(svec) ./ sqrt(numel(svec));
    end
    if ~isempty(rabs)
        epochSummary.meanResAbs(e) = mean(rabs);
        epochSummary.semResAbs(e)  = std(rabs) ./ sqrt(numel(rabs));
    end
end

OUT.epochSummary = epochSummary;

% ---------- plots for residuals ----------
if DoPlot
  % 1) Histogram: raw FR std vs residual FR std (per cell)
  rawStdAll = allRawStd_allCells(isfinite(allRawStd_allCells));
  resStdAll = allResStd_allCells(isfinite(allResStd_allCells));

  if ~isempty(rawStdAll) && ~isempty(resStdAll)
      figure('Color','w','Position',[300 300 520 380]); hold on;

      edges = linspace(0, max([rawStdAll; resStdAll]), 50);

      histogram(rawStdAll,'Normalization','probability')
      hold on
      histogram(resStdAll,'Normalization','probability')
      % Raw FR variability
      %histogram(rawStdAll, edges, ...
      %    'Normalization','probability', ...
      %    'FaceColor',[0.1 0.4 0.8], ...
      %    'FaceAlpha',0.4, ...
      %    'EdgeColor','none');

      % Residual FR variability (after speed regression)
      %histogram(resStdAll, edges, ...
      %    'Normalization','probability', ...
      %    'FaceColor',[0.9 0.3 0.1], ...
      %    'FaceAlpha',0.5, ...
      %    'EdgeColor','none');

      xlabel('Per-cell FR std (Hz)');
      ylabel('Probability');
      title('Raw vs speed-residual FR variability across cells');
      legend({'Raw FR std','Residual FR std'}, 'Location','northeast');
      box on;
  else
      warning('Not enough data to plot raw vs residual std.');
  end

    % 2) Per-epoch mean |residual FR| ± SEM
    m = epochSummary.meanResAbs;
    s = epochSummary.semResAbs;

    figure('Color','w','Position',[380 380 450 350]); hold on;
    b = bar(1:nEpoch, m, 'FaceColor',[0.2 0.5 0.9]);
    errorbar(1:nEpoch, m, s, 'k.', 'LineWidth',1.5);
    set(gca,'XTick',1:nEpoch,'XTickLabel',OUT.epochNames);
    ylabel('Mean |residual FR| (Hz)');
    title('Epoch-wise residual firing (after speed regression)');
    box on;
end
end

% ===== helper functions (same behavior as in xcorr_speed_MUA) =====

function Scell = normalizeCApeaks_local(Sraw)
if iscell(Sraw)
    if isvector(Sraw)
        Scell = Sraw(:);
        return;
    end
    K = size(Sraw,2);
    scores = zeros(1,K);
    for j = 1:K
        col = Sraw(:,j);
        scores(j) = mean(cellfun(@(x) isnumeric(x) && isvector(x), col));
    end
    [bestScore, jbest] = max(scores);
    if bestScore == 0
        error('normalizeCApeaks_local:NoNumericColumn');
    end
    Scell = Sraw(:, jbest);
    Scell = cellfun(@(v) v(:), Scell, 'uni', false);
    return;
end
if isnumeric(Sraw)
    n = size(Sraw,1);
    Scell = cell(n,1);
    for i = 1:n
        v = Sraw(i,:).';
        v = v(isfinite(v) & v>0);
        Scell{i,1} = v;
    end
    return;
end
error('normalizeCApeaks_local:BadType %s', class(Sraw));
end

function [vt, vv] = velocityFromPos_local(pos)
t  = pos(:,1);
xy = pos(:,2:3);
dt = diff(t);
dt(dt==0) = NaN;
v  = [0; sqrt(sum(diff(xy).^2,2)) ./ dt];
v(~isfinite(v)) = 0;
vt = t;
vv = v;
end

function [binSpeed, rateMat, effDur] = taskBinsSpeedAndRates_withDur_local( ...
        S, CS, tPos, vt, vv, Win, binSec, vThr) %#ok<INUSD>
% vThr is ignored here (pass []).

nCells  = numel(S);
dtPos   = median(diff(tPos),'omitnan');
vOnPos  = interp1(vt, vv, tPos, 'linear','extrap');

binSpeed = [];
effDur   = [];
rateMat  = [];

for k = 1:numel(CS)
    t0 = CS(k)+Win(1);
    t1 = CS(k)+Win(2);
    if ~isfinite(t0) || ~isfinite(t1) || t1<=t0
        continue;
    end

    edges = t0:binSec:t1;
    if edges(end) < t1
        edges = [edges t1];
    end

    for b = 1:numel(edges)-1
        b0 = edges(b);
        b1 = edges(b+1);
        inBin = (tPos>=b0) & (tPos<b1);

        mask = inBin;  % no speed gate
        dur  = sum(mask) * dtPos;
        if dur <= 0
            continue;
        end

        sbin      = mean(vOnPos(mask), 'omitnan');
        intervals = [b0 b1];

        row = NaN(1,nCells);
        for c = 1:nCells
            sp = S{c};
            if isempty(sp)
                row(c) = 0;
                continue;
            end
            nSp    = sum(sp>=intervals(1) & sp<intervals(2));
            row(c) = nSp / dur;
        end

        binSpeed = [binSpeed; sbin];
        effDur   = [effDur;   dur];
        rateMat  = [rateMat;  row];
    end
end

ok = isfinite(binSpeed) & (effDur>0) & any(isfinite(rateMat),2);
binSpeed = binSpeed(ok);
effDur   = effDur(ok);
rateMat  = rateMat(ok,:);
end

function Y = normalizePerCell_local(X, mode)
switch lower(mode)
    case 'none'
        Y = X;
    case 'mean'
        mu = mean(X,1,'omitnan');
        mu(mu==0 | ~isfinite(mu)) = NaN;
        Y = X ./ mu;
    case 'demean'
        mu = mean(X,1,'omitnan');
        Y  = X - mu;
    otherwise
        error('PopNorm must be ''none'', ''mean'', or ''demean''.');
end
end
