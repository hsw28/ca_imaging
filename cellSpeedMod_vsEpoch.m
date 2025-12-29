function OUT = cellSpeedMod_vsEpoch(ratNames, varargin)
% cellSpeedMod_vsEpoch
% For each cell, compute r(speed, FR) across task-time bins and relate it
% to the cell's preferred task epoch (from rat.epoch).
%
% Usage:
%   OUT = cellSpeedMod_vsEpoch({'rat0222','rat0307',...}, ...
%          'Win',[0 2], 'binSize',1/7.5, 'MinBins',20, ...
%          'EpochVar','epoch', 'EpochNum',5, ...
%          'TrialDemean',false, 'DoPlot',true);
%
% Required rat fields (per day):
%   Ca_peaks.CA_peaks_<date>
%   pos.pos_<date>            % [t x y]
%   CS_times.CS_<date>        % CS onsets
%   epoch.epoch_<date>        % per-cell integer labels 1..EpochNum
%
% OUTPUT OUT fields:
%   .ratName   [N x 1]  string
%   .day       [N x 1]  string
%   .cellIdx   [N x 1]  double
%   .epoch     [N x 1]  double
%   .r_speed   [N x 1]  double   (Pearson r(speed, FR) across bins)
%   .perEpochMean   [EpochNum x 1]
%   .perEpochSEM    [EpochNum x 1]

% ---------- options ----------
p = inputParser;
p.addParameter('Win',[0 2]);          % task window relative to CS
p.addParameter('binSize',1/7.5);      % ~frame size
p.addParameter('MinBins',20);         % min # bins to compute r
p.addParameter('EpochVar','epoch');   % rat.(EpochVar) struct with labels
p.addParameter('EpochNum',4);         % # possible labels (1..EpochNum)
p.addParameter('TrialDemean',false,@islogical);  % demean within each trial??? no, should be across. whatever
p.addParameter('DoPlot',true,@islogical);

p.parse(varargin{:});
Win        = p.Results.Win;
binSize    = p.Results.binSize;
MinBins    = p.Results.MinBins;
EpochVar   = char(p.Results.EpochVar);
EpochNum   = p.Results.EpochNum;
TrialDemean= p.Results.TrialDemean;
DoPlot     = p.Results.DoPlot;

if nargin<1 || isempty(ratNames)
    error('Provide ratNames cell array.');
end

nRats = numel(ratNames);

% storage (growable)
ratName_all = {};
day_all     = {};
cellIdx_all = [];
epoch_all   = [];
r_all       = [];

for r = 1:nRats
    ratName = ratNames{r};
    rat     = evalin('base', ratName);

    % autoDateList should already exist in your codebase
    dates = autoDateList(rat);
    iAn   = find(strcmp(dates, rat.An),1);
    if isempty(iAn) || iAn < 3
        warning('%s: not enough days before An; skipping.', ratName);
        continue;
    end
    days = dates(iAn-2:iAn);   % last 3 days up to An

    if ~isfield(rat, EpochVar)
        error('%s is missing field %s for epoch labels.', ratName, EpochVar);
    end
    Estruct = rat.(EpochVar);

    for d = 1:numel(days)
        D = days{d};

        % ----- spikes -----
        Sraw   = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));
        S      = normalizeCApeaks_local(Sraw);
        nCells = numel(S);

        % ----- position & velocity -----
        pos       = smoothpos(rat.pos.(sprintf('pos_%s',D)));
        [vt, vv]  = velocityFromPos_local(pos);
        CS        = rat.CS_times.(sprintf('CS_%s',D));

        % ----- epoch labels for this day -----
        fns    = fieldnames(Estruct);
        idxFld = contains(fns, D);
        if ~any(idxFld)
            error('%s: no field in %s contains date %s.', ratName, EpochVar, D);
        end
        if sum(idxFld)>1
            warning('%s: multiple fields in %s match %s; using first.', ...
                    ratName, EpochVar, D);
        end
        fname    = fns{find(idxFld,1,'first')};
        epochDay = Estruct.(fname);
        epochDay = epochDay(:);
        if numel(epochDay) < nCells
            warning('%s: %s.%s has only %d labels for %d cells; truncating.', ...
                    ratName, EpochVar, fname, numel(epochDay), nCells);
            epochDay(end+1:nCells) = NaN;
        end
        epochDay = epochDay(1:nCells);

        % ----- task bins (fine timescale) -----
        [binSpeed, rateMat, effDur] = taskBinsSpeedAndRates_withDur_local( ...
            S, CS, pos(:,1), vt, vv, Win, binSize, []);

        if isempty(binSpeed) || isempty(rateMat)
            continue;
        end

        % optional within-trial demeaning of speed and FR
        if TrialDemean
            [binSpeed, rateMat] = trial_demean_bins_local(binSpeed, rateMat, Win, binSize);
        end

        % per-cell speed–rate r
        for c = 1:nCells
            labs = epochDay(c);
            if ~isfinite(labs) || labs<1 || labs>EpochNum
                continue;   % unlabeled / out of range
            end
            y = rateMat(:,c);
            x = binSpeed;

            mask = isfinite(x) & isfinite(y) & effDur>0;
            if nnz(mask) < MinBins
                continue;
            end

            rc = corr(x(mask), y(mask), 'rows','complete');

            ratName_all{end+1,1} = ratName; %#ok<AGROW>
            day_all{end+1,1}     = D;       %#ok<AGROW>
            cellIdx_all(end+1,1) = c;       %#ok<AGROW>
            epoch_all(end+1,1)   = labs;    %#ok<AGROW>
            r_all(end+1,1)       = rc;      %#ok<AGROW>
        end
    end
end

% ---------- pack output ----------
OUT = struct();
OUT.ratName    = string(ratName_all);
OUT.day        = string(day_all);
OUT.cellIdx    = cellIdx_all;
OUT.epoch      = epoch_all;
OUT.r_speed    = r_all;
OUT.TrialDemean= TrialDemean;
OUT.Win        = Win;
OUT.binSize    = binSize;

% per-epoch summary
perEpochMean = nan(EpochNum,1);
perEpochSEM  = nan(EpochNum,1);
for e = 1:EpochNum
    vals = r_all(epoch_all==e & isfinite(r_all));
    if isempty(vals), continue; end
    perEpochMean(e) = mean(vals);
    perEpochSEM(e)  = std(vals) / sqrt(numel(vals));
end
OUT.perEpochMean = perEpochMean;
OUT.perEpochSEM  = perEpochSEM;

% ---------- optional plot ----------
if DoPlot
    figure('Color','w','Position',[200 200 800 350]);
    subplot(1,2,1); hold on;
    boxplot(r_all, epoch_all, ...
        'Labels', arrayfun(@(k)sprintf('E%d',k),1:EpochNum,'uni',0));
    ylabel('r(speed, FR)');
    xlabel('Epoch label');
    title(sprintf('Per-cell speed modulation by epoch (TrialDemean=%d)',TrialDemean));

    subplot(1,2,2); hold on;
    bar(1:EpochNum, perEpochMean);
    errorbar(1:EpochNum, perEpochMean, perEpochSEM, 'k.', 'LineWidth',1.5);
    xlim([0.5 EpochNum+0.5]);
    set(gca,'XTick',1:EpochNum,...
        'XTickLabel',arrayfun(@(k)sprintf('E%d',k),1:EpochNum,'uni',0));
    ylabel('Mean r(speed, FR)');
    title('Epoch-wise mean speed modulation');
end
end

%% ===== Helpers (same logic as in xcorr_speed_MUA) =====

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
error('normalizeCApeaks_local:BadType');
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
        S, CS, tPos, vt, vv, Win, binSec, vThr)
% Split each trial window [CS+Win(1), CS+Win(2)) into binSec chunks.
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

        if isempty(vThr)
            mask = inBin;
        else
            mask = inBin & (vOnPos>=vThr);
        end

        dur = sum(mask) * dtPos;
        if dur <= 0
            continue;
        end

        sbin = mean(vOnPos(mask), 'omitnan');

        if isempty(vThr)
            intervals = [b0 b1];
        else
            intervals = maskToIntervals_local(tPos, mask, dtPos);
        end

        row = NaN(1,nCells);
        for c = 1:nCells
            sp = S{c};
            if isempty(sp)
                row(c) = 0;
                continue;
            end
            if isempty(vThr)
                nSp = sum(sp>=b0 & sp<b1);
            else
                nSp = countInIntervals_local(sp, intervals);
            end
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

function [binSpeedOut, rateMatOut] = trial_demean_bins_local(binSpeed, rateMat, Win, binSec)
% Demean speed and FR within each trial window.
nBinsPerTrial = round(diff(Win) / binSec);
if nBinsPerTrial <= 0
    binSpeedOut = binSpeed;
    rateMatOut  = rateMat;
    return;
end

nTotal = numel(binSpeed);
nFull  = floor(nTotal / nBinsPerTrial);
if nFull == 0
    binSpeedOut = binSpeed;
    rateMatOut  = rateMat;
    return;
end

idxKeep = 1:(nFull*nBinsPerTrial);

% reshape to [nTrials x nBins]
Smat = reshape(binSpeed(idxKeep),  nBinsPerTrial, nFull).';
nCells = size(rateMat,2);
Rmat = reshape(rateMat(idxKeep,:), nBinsPerTrial, nFull, nCells);
Rmat = permute(Rmat,[2 1 3]);   % [nTrials x nBins x nCells]

% demean per trial (row) for speed
Smu  = mean(Smat,2,'omitnan');
Smat = Smat - Smu;

% demean per trial for each cell
for c = 1:nCells
    Rc   = squeeze(Rmat(:,:,c));          % [nTrials x nBins]
    Rmu  = mean(Rc,2,'omitnan');
    Rc   = Rc - Rmu;
    Rmat(:,:,c) = Rc;
end

% reshape back
binSpeedOut = Smat.'; binSpeedOut = binSpeedOut(:);
Rmat2       = permute(Rmat,[2 1 3]);      % [nBins x nTrials x nCells]
rateMatOut  = reshape(Rmat2, nFull*nBinsPerTrial, nCells);

% append leftover bins unchanged
if nFull*nBinsPerTrial < nTotal
    tailIdx     = (nFull*nBinsPerTrial+1):nTotal;
    binSpeedOut = [binSpeedOut; binSpeed(tailIdx)];
    rateMatOut  = [rateMatOut;  rateMat(tailIdx,:)];
end
end

function intervals = maskToIntervals_local(t, mask, dt)
if nargin<3 || isempty(dt)
    dt = median(diff(t),'omitnan');
end
mask = mask(:);
t    = t(:);
dm     = diff([false; mask; false]);
starts = find(dm==1);
ends   = find(dm==-1)-1;
intervals = [t(starts)  t(ends)+dt];
end

function n = countInIntervals_local(spikes, intervals)
if isempty(spikes) || isempty(intervals)
    n = 0;
    return;
end
spikes = spikes(:);
inside = false(size(spikes));
for k = 1:size(intervals,1)
    inside = inside | (spikes>=intervals(k,1) & spikes<intervals(k,2));
end
n = sum(inside);
end
