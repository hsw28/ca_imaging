function OUT = PV_speedRegress_temporal(ratNames, varargin)
% PV_speedRegress_temporal
%
% For each rat, within the CS+Win window:
%   * For each day (last 3 pre-An):
%       1) Bin trials at binSize resolution.
%       2) Regress per-cell FR on running speed across all bins & trials:
%              FR_c(t) = b0_c + b1_c * speed(t) + eps_c(t)
%       3) Build per-bin population vectors (mean FR across trials)
%          using raw FR and residual FR.
%       4) Compute time×time PV similarity matrices (raw vs residual).
%          If SplitHalf=true, similarity is even-vs-odd trials.
%   * Then average similarity matrices across days for that rat.
%   * Then average across rats.
%
% Usage:
%   OUT = PV_speedRegress_temporal({'rat0222','rat0307',...}, ...
%           'Win',[0 2], 'binSize',1/7.5, ...
%           'PopNorm','none', ...       % 'none' | 'mean' | 'demean'
%           'Similarity','corr', ...    % 'corr' | 'cosine'
%           'MinTrials',10, ...
%           'SplitHalf',false, ...
%           'DoPlot',true);
%
% OUTPUT OUT:
%   .meta.Win, .meta.binSize, .meta.Similarity, .meta.PopNorm,
%   .meta.SplitHalf, .meta.ratNames
%   .perRat(r) with fields:
%       .ratName
%       .days             : dates used
%       .nBinsPerTrial
%       .PV_raw_days      : cell{nDays} of [nBins x nCells_day] (all trials)
%       .PV_resid_days    : same
%       .S_raw_days       : cell{nDays} of [nBins x nBins]
%       .S_resid_days     : same
%       .S_raw            : mean across days [nBins x nBins]
%       .S_resid          : mean across days [nBins x nBins]
%   .S_raw_mean           : mean S_raw across rats
%   .S_resid_mean         : mean S_resid across rats

% ---------- options ----------
p = inputParser;
p.addParameter('Win',[0 2]);             % trial window relative to CS
p.addParameter('binSize',1/7.5);         % bin length (s)
p.addParameter('PopNorm','demean', ...
    @(s)any(strcmpi(s,{'none','mean','demean'})));
p.addParameter('Similarity','corr', ...
    @(s)any(strcmpi(s,{'corr','cosine'})));
p.addParameter('MinTrials',5,@(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('SplitHalf',true,@islogical);   % NEW
p.addParameter('DoPlot',true,@islogical);

p.parse(varargin{:});
Win        = p.Results.Win;
binSize    = p.Results.binSize;
popNorm    = lower(p.Results.PopNorm);
SimType    = lower(p.Results.Similarity);
MinTrials  = p.Results.MinTrials;
SplitHalf  = p.Results.SplitHalf;
DoPlot     = p.Results.DoPlot;

if nargin<1 || isempty(ratNames)
    error('Provide ratNames cell array.');
end
nR = numel(ratNames);

% ---------- prepare OUT ----------
OUT = struct();
OUT.meta = struct('Win',Win,'binSize',binSize,'PopNorm',popNorm, ...
                  'Similarity',SimType,'SplitHalf',SplitHalf, ...
                  'ratNames',{ratNames});

perRat = struct([]);

S_raw_all   = {};
S_resid_all = {};

% ---------- main loop over rats ----------
for r = 1:nR
    ratName = ratNames{r};
    rat     = evalin('base', ratName);

    dates = autoDateList(rat);
    iAn   = find(strcmp(dates, rat.An),1);
    if isempty(iAn) || iAn < 3
        warning('%s: not enough days before An; skipping rat.', ratName);
        continue;
    end
    days = dates(iAn-2:iAn);   % last 3 days up to An

    PV_raw_days   = {};
    PV_resid_days = {};
    S_raw_days    = {};
    S_resid_days  = {};
    nBinsPerTrial = [];

    for d = 1:numel(days)
        D = days{d};

        % ----- spikes -----
        Sraw   = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));
        S      = normalizeCApeaks_local(Sraw);
        nCells_day = numel(S);

        % ----- position & speed -----
        pos       = smoothpos(rat.pos.(sprintf('pos_%s',D)));
        [vt, vv]  = velocityFromPos_local(pos);
        tPos      = pos(:,1);
        CS        = rat.CS_times.(sprintf('CS_%s',D));

        % ----- flattened task bins for this day -----
        [binSpeed_flat, rateMat_flat, effDur_flat] = ...
            taskBinsSpeedAndRates_withDur_local( ...
                S, CS, tPos, vt, vv, Win, binSize, []);

        if isempty(binSpeed_flat)
            continue;
        end

        % Determine bins-per-trial and reshape into [nBins x nTrials x nCells_day]
        nBins = round(diff(Win) / binSize);
        if nBins < 1
            error('Win/binSize gives <1 bin per trial.');
        end
        nTotal   = numel(binSpeed_flat);
        nTrials  = floor(nTotal / nBins);
        if nTrials < MinTrials
            fprintf('%s %s: only %d trials after binning; skipping this day.\n', ...
                ratName, D, nTrials);
            continue;
        end
        idxKeep  = 1:(nBins*nTrials);

        binSpeed_flat = binSpeed_flat(idxKeep);
        effDur_flat   = effDur_flat(idxKeep); %#ok<NASGU>
        rateMat_flat  = rateMat_flat(idxKeep, :);  % [nBins*nTrials x nCells_day]

        % reshape speed: [nBins x nTrials]
        speed_trials = reshape(binSpeed_flat, nBins, nTrials);

        % reshape FR: [nBins x nTrials x nCells_day]
        tmp = reshape(rateMat_flat.', nCells_day, nBins, nTrials);
        rate_trials = permute(tmp, [2 3 1]);  % [nBins x nTrials x nCells_day]

        if isempty(nBinsPerTrial)
            nBinsPerTrial = nBins;
        else
            if nBins ~= nBinsPerTrial
                error('Day %s in rat %s has different nBinsPerTrial.', D, ratName);
            end
        end

        % ---------- optional per-cell normalization (within this day) ----------
        rate_reshaped = reshape(rate_trials, nBins*nTrials, nCells_day);
        rate_reshaped = normalizePerCell_local(rate_reshaped, popNorm);
        rate_trials   = reshape(rate_reshaped, nBins, nTrials, nCells_day);

        % ---------- regress out speed per cell for this day ----------
        x = speed_trials(:);     % [nBins*nTrials x 1]
        rate_resid_trials = nan(size(rate_trials));

        for c = 1:nCells_day
            y = rate_trials(:,:,c);
            y = y(:);

            mask = isfinite(x) & isfinite(y);
            if nnz(mask) < 5
                continue;
            end
            xx = x(mask);
            yy = y(mask);

            Xmat = [ones(numel(xx),1), xx];
            beta = Xmat \ yy;
            yhat = Xmat * beta;
            resid = yy - yhat;

            y_full = nan(size(y));
            y_full(mask) = resid;
            rate_resid_trials(:,:,c) = reshape(y_full, nBins, nTrials);
        end

        % ---------- per-bin PVs (mean across ALL trials) ----------
        PV_raw_day   = squeeze(mean(rate_trials,       2, 'omitnan')); % [nBins x nCells_day]
        PV_resid_day = squeeze(mean(rate_resid_trials, 2, 'omitnan'));

        % ---------- PV similarity matrices for this day ----------
        if ~SplitHalf
            % classic within-all-trials similarity
            S_raw_day   = compute_PV_sim(PV_raw_day,   SimType);
            S_resid_day = compute_PV_sim(PV_resid_day, SimType);
        else
            % split-half: even vs odd trials
            idxOdd  = 1:2:nTrials;
            idxEven = 2:2:nTrials;

            if isempty(idxOdd) || isempty(idxEven)
                warning('%s %s: not enough trials for split-half; using all-trials.', ...
                        ratName, D);
                S_raw_day   = compute_PV_sim(PV_raw_day,   SimType);
                S_resid_day = compute_PV_sim(PV_resid_day, SimType);
            else
                % PVs in each half
                PV_raw_even   = squeeze(mean(rate_trials(:,idxEven,:),       2,'omitnan'));
                PV_raw_odd    = squeeze(mean(rate_trials(:,idxOdd,:),        2,'omitnan'));
                PV_res_even   = squeeze(mean(rate_resid_trials(:,idxEven,:), 2,'omitnan'));
                PV_res_odd    = squeeze(mean(rate_resid_trials(:,idxOdd,:),  2,'omitnan'));

                S_raw_day   = compute_PV_sim_cross(PV_raw_even,  PV_raw_odd,  SimType);
                S_resid_day = compute_PV_sim_cross(PV_res_even,  PV_res_odd,  SimType);
            end
        end

        fprintf('%s %s: %d trials, %d bins/trial, %d cells.\n', ...
            ratName, D, nTrials, nBins, nCells_day);

        PV_raw_days{end+1}   = PV_raw_day;   %#ok<AGROW>
        PV_resid_days{end+1} = PV_resid_day; %#ok<AGROW>
        S_raw_days{end+1}    = S_raw_day;    %#ok<AGROW>
        S_resid_days{end+1}  = S_resid_day;  %#ok<AGROW>
    end

    if isempty(S_raw_days)
        warning('%s: no usable days/trials; skipping rat.', ratName);
        continue;
    end

    % ---------- average S across days for this rat ----------
    nDays_use = numel(S_raw_days);
    nBins     = size(S_raw_days{1},1);
    Sraw_stack   = nan(nBins,nBins,nDays_use);
    Sresid_stack = nan(nBins,nBins,nDays_use);
    for k = 1:nDays_use
        Sraw_stack(:,:,k)   = S_raw_days{k};
        Sresid_stack(:,:,k) = S_resid_days{k};
    end
    S_raw_rat   = mean(Sraw_stack,   3, 'omitnan');
    S_resid_rat = mean(Sresid_stack, 3, 'omitnan');

    % store per rat
    pr = struct();
    pr.ratName       = ratName;
    pr.days          = days;
    pr.nBinsPerTrial = nBinsPerTrial;
    pr.PV_raw_days   = PV_raw_days;
    pr.PV_resid_days = PV_resid_days;
    pr.S_raw_days    = S_raw_days;
    pr.S_resid_days  = S_resid_days;
    pr.S_raw         = S_raw_rat;
    pr.S_resid       = S_resid_rat;

    perRat = [perRat; pr]; %#ok<AGROW>
    S_raw_all{end+1}   = S_raw_rat;   %#ok<AGROW>
    S_resid_all{end+1} = S_resid_rat; %#ok<AGROW>
end

OUT.perRat = perRat;

% ---------- group-average similarity across rats ----------
if ~isempty(S_raw_all)
    nBins = size(S_raw_all{1},1);
    nUse  = numel(S_raw_all);
    Sraw_stack   = nan(nBins,nBins,nUse);
    Sresid_stack = nan(nBins,nBins,nUse);
    for k = 1:nUse
        Sraw_stack(:,:,k)   = S_raw_all{k};
        Sresid_stack(:,:,k) = S_resid_all{k};
    end
    OUT.S_raw_mean   = mean(Sraw_stack,   3, 'omitnan');
    OUT.S_resid_mean = mean(Sresid_stack, 3, 'omitnan');
else
    OUT.S_raw_mean   = [];
    OUT.S_resid_mean = [];
end

% ---------- plotting ----------
if DoPlot && ~isempty(OUT.S_raw_mean)
    figure('Color','w','Position',[200 200 900 350]);

    subplot(1,2,1);
    imagesc(OUT.S_raw_mean, [-1 1]);
    axis square; colorbar;
    if SplitHalf
        title(sprintf('PV similarity (raw, split-half) – %s', SimType));
    else
        title(sprintf('PV similarity (raw) – %s', SimType));
    end
    xlabel('Time bin'); ylabel('Time bin');

    subplot(1,2,2);
    imagesc(OUT.S_resid_mean, [-1 1]);
    axis square; colorbar;
    if SplitHalf
        title(sprintf('PV similarity (speed-residual, split-half) – %s', SimType));
    else
        title(sprintf('PV similarity (speed-residual) – %s', SimType));
    end
    xlabel('Time bin'); ylabel('Time bin');
end
end

% =========================================================
% ==================== LOCAL HELPERS ======================
% =========================================================

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
        S, CS, tPos, vt, vv, Win, binSec, vThr)
% Flattened task bins across all trials for a single day.
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

        sbin      = mean(vOnPos(mask), 'omitnan');
        intervals = [b0 b1];

        row = NaN(1,nCells);
        for c = 1:nCells
            sp = S{c};
            if isempty(sp)
                row(c) = 0;
                continue;
            end
            if isempty(vThr)
                nSp = sum(sp>=intervals(1) & sp<intervals(2));
            else
                nSp = sum(sp>=intervals(1) & sp<intervals(2));
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

function Y = normalizePerCell_local(X, mode)
switch lower(mode)
    case 'none'
        Y = X;
    case 'mean'
        mu = mean(X,1,'omitnan');
        mu(mu==0 | ~isfinite(mu)) = NaN;
        Y = bsxfun(@rdivide, X, mu);
    case 'demean'
        mu = mean(X,1,'omitnan');
        Y  = bsxfun(@minus, X, mu);
    otherwise
        error('PopNorm must be ''none'' or ''mean'' or ''demean''.');
end
end

function S = compute_PV_sim(PV, SimType)
[nBins, ~] = size(PV);
S = nan(nBins, nBins);

switch lower(SimType)
    case 'corr'
        S = corrcoef(PV.');  % [nBins x nBins]
    case 'cosine'
        for i = 1:nBins
            vi = PV(i,:);
            for j = 1:nBins
                vj = PV(j,:);
                mask = isfinite(vi) & isfinite(vj);
                vi2 = vi(mask);
                vj2 = vj(mask);
                if isempty(vi2) || isempty(vj2)
                    S(i,j) = NaN;
                    continue;
                end
                num = sum(vi2 .* vj2);
                den = sqrt(sum(vi2.^2) * sum(vj2.^2));
                if den==0
                    S(i,j) = NaN;
                else
                    S(i,j) = num / den;
                end
            end
        end
    otherwise
        error('Unknown Similarity type: %s', SimType);
end
end

function S = compute_PV_sim_cross(PV_A, PV_B, SimType)
% PV_A, PV_B: [nBins x nCells] (e.g. even vs odd trials)
nBins = size(PV_A,1);
S     = nan(nBins, nBins);

switch lower(SimType)
    case 'corr'
        for i = 1:nBins
            for j = 1:nBins
                vi = PV_A(i,:);
                vj = PV_B(j,:);
                mask = isfinite(vi) & isfinite(vj);
                if ~any(mask), continue; end
                S(i,j) = corr(vi(mask).', vj(mask).');
            end
        end
    case 'cosine'
        for i = 1:nBins
            for j = 1:nBins
                vi = PV_A(i,:);
                vj = PV_B(j,:);
                mask = isfinite(vi) & isfinite(vj);
                vi2 = vi(mask); vj2 = vj(mask);
                if isempty(vi2) || isempty(vj2), continue; end
                num = sum(vi2 .* vj2);
                den = sqrt(sum(vi2.^2) * sum(vj2.^2));
                if den==0, continue; end
                S(i,j) = num / den;
            end
        end
    otherwise
        error('Unknown Similarity type: %s', SimType);
end
end
