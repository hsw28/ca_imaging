function R = taskFiringVariability(ratName, varargin)
% taskFiringVariability
% Resolves days (default An-2:An), runs the per-day helper, pools, and plots.
%
% Usage:
%   R = taskFiringVariability('rat0222');                     % uses An-2:An
%   R = taskFiringVariability('rat0222','Days','2023_05_04'); % single day
%   R = taskFiringVariability('rat0222','Days',{'2023-05-04','2023/05/09'})
%
% Pass-through options (same defaults as inner function):
%   'WinSecs',[0 2], 'SpeedThresh',4, 'NSD',1, 'MinTrials',10, 'MinSpikes',3
%   'CentroidMethod','mode' (mode|medoid|mean|speed), 'CentroidBinCm',2.5
% Extra for wrapper:
%   'Plot' (true): pooled figure across all days

% ---------- args ----------
p = inputParser;
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty(d));
addParameter(p,'WinSecs',[0 2], @(v) isnumeric(v)&&numel(v)==2&&v(2)>v(1));
addParameter(p,'SpeedThresh',4, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'NSD',1, @(x) isnumeric(x)&&isscalar(x)&&isfinite(x));
addParameter(p,'MinTrials',10, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'MinSpikes',0, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'CentroidMethod','mode', ...
    @(s) any(validatestring(lower(s),{'mode','medoid','mean','speed'})));
addParameter(p,'CentroidBinCm',2.5, @(x) isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Plot',true, @islogical);
parse(p,varargin{:});

daysArg        = p.Results.Days;
winSecs        = p.Results.WinSecs;
vMin           = p.Results.SpeedThresh;
Nsd            = p.Results.NSD;
minT           = p.Results.MinTrials;
minSpk         = p.Results.MinSpikes;
centroidMethod = lower(p.Results.CentroidMethod);
centroidBinCm  = p.Results.CentroidBinCm;
doPlot         = p.Results.Plot;

% ---------- resolve days ----------
rat   = evalin('base', ratName);
dates = autoDateList(rat);
iAn   = find(strcmp(dates, rat.An), 1);

if ~isempty(daysArg)
    userDays = cellstr(daysArg);
    days     = cellfun(@(s) resolveDayToken(rat, s), userDays, 'UniformOutput', false);
else
    if isempty(iAn), error('An date not found for %s.', ratName); end
    days = dates(max(1,iAn-2):iAn);
end
days = unique(days,'stable');
fprintf('Using days: %s\n', strjoin(days, ', '));

% ---------- run per-day helper ----------
perDay = struct('date',[],'perCell',[],'pairs',[]);
all_trial_rates_inPF  = [];   % per-trial across ALL cells/days
all_trial_rates_outPF = [];
perCell_mean_in  = [];        % per-cell means for paired stats
perCell_mean_out = [];

pair_dist_all = [];           % all pairwise distances (cm), PF ignored
pair_ndiff_all = [];          % normalized |Δ|/mean per pair (PF ignored)

for d = 1:numel(days)
    Rd = taskFiringSpatialVariability_oneDay( ...
            ratName, days{d}, ...
            'WinSecs',winSecs, 'SpeedThresh',vMin, 'NSD',Nsd, ...
            'MinTrials',minT, 'MinSpikes',minSpk, ...
            'CentroidMethod',centroidMethod, 'CentroidBinCm',centroidBinCm, ...
            'Plot',false);

    perDay(d).date   = Rd.date;
    perDay(d).perCell= Rd.perCell;
    perDay(d).pairs  = Rd.pairs;

    % ---- pools for Fig: PF-vs (left panel) ----
    all_trial_rates_inPF  = [all_trial_rates_inPF;  Rd.extra.trialRates.inPF];
    all_trial_rates_outPF = [all_trial_rates_outPF; Rd.extra.trialRates.outPF];
    perCell_mean_in  = [perCell_mean_in;  Rd.extra.perCellMean.inPF(:)];
    perCell_mean_out = [perCell_mean_out; Rd.extra.perCellMean.outPF(:)];

    % ---- pools for Fig: distance vs variability (right panel) ----
    pair_dist_all  = [pair_dist_all;  Rd.pairs.dist_cm];
    pair_ndiff_all = [pair_ndiff_all; Rd.pairs.ndiff];
end

% ---------- assemble output ----------
R = struct();
R.ratVar = ratName;
R.days   = days;
R.params = struct('WinSecs',winSecs,'SpeedThresh',vMin,'NSD',Nsd, ...
                  'MinTrials',minT,'MinSpikes',minSpk, ...
                  'CentroidMethod',centroidMethod,'CentroidBinCm',centroidBinCm);
R.perDay = perDay;

R.extra.trialRates.inPF   = all_trial_rates_inPF;
R.extra.trialRates.outPF  = all_trial_rates_outPF;
R.extra.perCellMean.inPF  = perCell_mean_in;
R.extra.perCellMean.outPF = perCell_mean_out;

R.pairs.dist_cm           = pair_dist_all;    % PF ignored
R.pairs.ndiff             = pair_ndiff_all;   % PF ignored

% ---------- pooled figure (optional quick-check) ----------
if doPlot
    figure('Color','w','Position',[160 120 1150 480]);

    % LEFT: In-field vs Out-of-field (per-cell paired means)
    subplot(1,2,1); hold on; box on
    mi = R.extra.perCellMean.inPF(:);
    mo = R.extra.perCellMean.outPF(:);
    goodPair = isfinite(mi) & isfinite(mo);
    mi = mi(goodPair); mo = mo(goodPair);
    K  = numel(mi);
    for k=1:K
        plot([1 2],[mi(k) mo(k)],'-','Color',[0.7 0.7 0.7]);
    end
    scatter(ones(K,1), mi, 28, 'filled');
    scatter(2*ones(K,1), mo, 28, 'filled');
    xlim([0.5 2.5]); xticks([1 2]); xticklabels({'In PF (mean/cell)','Out PF (mean/cell)'});
    ylabel('Mean task rate (Hz)');
    title('Per-cell means (paired)');
    try
        [~,p_t,~,st] = ttest(mi, mo);
        text(0.55, 0.92, sprintf('paired t-test: p=%.3g', p_t), 'Units','normalized');
    catch, end

    % RIGHT: distance vs normalized variability (PF ignored)
    subplot(1,2,2); hold on; box on
    d  = R.pairs.dist_cm(:);
    nd = R.pairs.ndiff(:);
    good = isfinite(d) & isfinite(nd);
    d = d(good); nd = nd(good);
    scatter(d, nd, 8, 'filled', 'MarkerFaceAlpha',0.25);
    xlabel('Centroid distance between trials (cm)');
    ylabel('Normalized |Δrate| / mean');
    title('Pairwise variability vs distance');

    if numel(d) >= 3
        X = [ones(numel(d),1), d];
        b = X \ nd;
        xfit = linspace(min(d), max(d), 100)';
        yfit = [ones(numel(xfit),1) xfit]*b;
        plot(xfit, yfit, 'LineWidth', 1.8);
        [rval, pval] = corr(d, nd, 'rows','complete', 'type','Pearson');
        text(0.02, 0.95, sprintf('r = %.2f, p = %.3g', rval, pval), ...
             'Units','normalized', 'VerticalAlignment','top');
        legend({'pairs','linear fit'},'Location','best');
    end
end
end

% ================= per-day helper =================
function R = taskFiringSpatialVariability_oneDay(ratName, dateStr, varargin)
% Defines PF from NON-TASK running (speed>=SpeedThresh; CS→CS+2 excluded),
% computes per-trial task rates, tags trials In/Out PF (by centroid),
% builds (distance, |Δ|/mean) pairs ignoring PF, and returns pools.
%
% Returns:
%   R.perCell(c).rates         [Tkeep x 1]
%   R.perCell(c).centroids     [Tkeep x 2]
%   R.perCell(c).inPF_trial    [Tkeep x 1] logical
%   R.extra.trialRates.(inPF/outPF)  pooled per-trial
%   R.extra.perCellMean.(inPF/outPF) paired means per cell
%   R.pairs.dist_cm, R.pairs.ndiff   pooled pairs across cells (PF ignored)

p = inputParser;
addParameter(p,'WinSecs',[0 2], @(v) isnumeric(v)&&numel(v)==2&&v(2)>v(1));
addParameter(p,'SpeedThresh',4, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'NSD',1, @(x) isnumeric(x)&&isscalar(x)&&isfinite(x));
addParameter(p,'MinTrials',10, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'MinSpikes',3, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'CentroidMethod','mode', ...
    @(s) any(validatestring(lower(s),{'mode','medoid','mean','speed'})));
addParameter(p,'CentroidBinCm',2.5, @(x) isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Plot',false, @islogical);
parse(p,varargin{:});
winSecs        = p.Results.WinSecs;
vMin           = p.Results.SpeedThresh;
Nsd            = p.Results.NSD;
minT           = p.Results.MinTrials;
minSpk         = p.Results.MinSpikes;
centroidMethod = lower(p.Results.CentroidMethod);
centroidBinCm  = p.Results.CentroidBinCm;

% -------- fetch & sanitize day token --------
rat = evalin('base', ratName);
D   = strrep(strrep(dateStr,'-','_'),'/','_');
forms = {D};
m = regexp(D, '^(\d{4})_(\d{1,2})_(\d{1,2})$','tokens','once');
if ~isempty(m)
    yy = m{1}; mm = str2double(m{2}); dd = str2double(m{3});
    forms = {sprintf('%s_%02d_%02d',yy,mm,dd), sprintf('%s_%d_%d',yy,mm,dd)};
end
found = '';
for f = 1:numel(forms)
    cand = forms{f};
    if isfield(rat.Ca_peaks,['CA_peaks_' cand]) && isfield(rat.pos,['pos_' cand]) ...
       && isfield(rat.CS_times,['CS_' cand])
        found = cand; break
    end
end
if isempty(found), error('Could not find day %s in rat.', dateStr); end
D = found;

% -------- pull data --------
S   = rat.Ca_peaks.(['CA_peaks_' D]);      % spikes
mask = rat.ratemask.(['ratemask_' D]);
pos = rat.pos.(['pos_' D]);                % [t x y]
CS  = rat.CS_times.(['CS_' D])(:);         % trials
if isempty(CS), error('No CS times for %s.', D); end

if iscell(S), nCells = numel(S); else, nCells = size(S,1); end
ts   = pos(:,1); xy = pos(:,2:3);

% -------- non-task running stream (for PF masks) --------
vel      = ca_velocity(pos);              % [speed; time]
vt       = vel(2,:)'; vmag = vel(1,:)';
xv       = interp1(ts, xy(:,1), vt, 'linear', NaN);
yv       = interp1(ts, xy(:,2), vt, 'linear', NaN);
inTrialV = false(size(vt));
for k = 1:numel(CS), inTrialV = inTrialV | (vt >= CS(k) & vt < CS(k)+2); end
keepV    = (vmag >= vMin) & ~inTrialV & isfinite(xv) & isfinite(yv);
runPos   = [vt(keepV), xv(keepV), yv(keepV)];

% -------- trial centroids in task window (configurable) --------
t0 = winSecs(1); t1 = winSecs(2); dur = t1 - t0;
T  = numel(CS);
trialWindows = [CS+t0, CS+t1];

% Build grid for 'mode' method (centers/edges aligned to visited runPos)
[Xedges, Yedges, xCtr, yCtr] = make_xy_edges(runPos(:,2), runPos(:,3), centroidBinCm);

% Sample durations for dwell-time weights
dtAll = [diff(ts); median(diff(ts))];

% Compute centroids by selected method
trialCentroids = compute_trial_centroids(ts, xy, trialWindows, centroidMethod, ...
                                         vt, vmag, Xedges, Yedges, xCtr, yCtr, dtAll);

validTrials = all(isfinite(trialCentroids),2);

% -------- per-cell loop --------
perCell = struct('rates',[],'centroids',trialCentroids(validTrials,:), ...
                 'inPF_trial',[]);
% pooled per-trial & per-cell-means
trialRates_in  = []; trialRates_out = [];
mean_in = []; mean_out = [];

% pooled pairs for distance scatter (PF ignored)
pool_dist = []; pool_ndiff = [];

for c = 1:nCells
    st = getCellSpikes_local(S, c);
    if mask(c)==0, continue; end
    if isempty(st), continue; end

    % remove trial spikes for PF map
    inTrialS = false(size(st));
    for k = 1:numel(CS), inTrialS = inTrialS | (st >= CS(k) & st < CS(k)+2); end
    st_non = st(~inTrialS);

    % PF mask from non-task running
    rate_non = CA_normalizePosData(st_non(:), runPos, 2.5, 1);  % (y,x)
    if isempty(rate_non) || all(rate_non(:)==0), continue, end
    rn = rate_non; rn(rn==0) = NaN;
    mu = nanmean(rn(:)); sd = nanstd(rn(:));
    mask_non = rn > (mu + Nsd*sd);
    if ~any(mask_non(:)), mask_non = rn > 0; if ~any(mask_non(:)), continue, end, end


    % per-trial task rate + inPF flag by centroid
    % ---- PF center (weighted COM of suprathreshold non-task map) ----
    Xedges_non = linspace(min(runPos(:,2)), max(runPos(:,2)), size(rate_non,2)+1);
    Yedges_non = linspace(min(runPos(:,3)), max(runPos(:,3)), size(rate_non,1)+1);
    xCtr_non = (Xedges_non(1:end-1)+Xedges_non(2:end))/2;
    yCtr_non = (Yedges_non(1:end-1)+Yedges_non(2:end))/2;

    [yy,xx] = find(mask_non);
    if isempty(xx)
        continue
    end
    w_pf = rn(mask_non);                   % weight by non-task rate
    cx_pf = sum(xCtr_non(xx).*w_pf) / max(eps,sum(w_pf));
    cy_pf = sum(yCtr_non(yy).*w_pf) / max(eps,sum(w_pf));

    % PF center as rate-weighted COM within the PF mask
    W = rn;                 % rate map with zeros -> NaN already above
    W(~mask_non) = 0;       % keep only PF bins
    [XX,YY] = meshgrid(xCtr_non, yCtr_non);
    totW = nansum(W(:));
    if totW > 0
        pf_center = [ nansum(W(:).*XX(:))/totW , nansum(W(:).*YY(:))/totW ];
    else
        % fallback: geometric center of PF bins
        [by,bx] = find(mask_non);
        if isempty(by)
            pf_center = [mean(xCtr_non,'omitnan') mean(yCtr_non,'omitnan')];
        else
            pf_center = [ mean(xCtr_non(unique(bx)),'omitnan') , mean(yCtr_non(unique(by)),'omitnan') ];
        end
    end

    % ----- per-trial task rate + inPF flag by centroid -----
    rates = nan(T,1);
    inPF_trial = false(T,1);

    for tr = 1:T
        w0 = trialWindows(tr,1); w1 = trialWindows(tr,2);
        if ~isfinite(trialCentroids(tr,1)), continue, end
        nSpk = sum(st >= w0 & st < w1);
        rates(tr) = nSpk / dur;

        cx = trialCentroids(tr,1); cy = trialCentroids(tr,2);
        bx = discretize(cx, Xedges_non);
        by = discretize(cy, Yedges_non);
        if bx>=1 && by>=1 && bx<=size(mask_non,2) && by<=size(mask_non,1)
            inPF_trial(tr) = mask_non(by,bx);
        end
    end

    keep = validTrials & isfinite(rates);
    rates = rates(keep);
    cents = trialCentroids(keep,:);
    inPF  = inPF_trial(keep);
    % distance from each kept trial centroid to PF center
dist_to_pf = sqrt( (cents(:,1)-cx_pf).^2 + (cents(:,2)-cy_pf).^2 );




    if numel(rates) < minT || sum(rates) < minSpk, continue, end


    % store per-trial pools
    trialRates_in  = [trialRates_in;  rates(inPF)];
    trialRates_out = [trialRates_out; rates(~inPF)];

    % per-cell paired means (require both)
    mi = mean(rates(inPF),'omitnan');
    mo = mean(rates(~inPF),'omitnan');
    if isfinite(mi) && isfinite(mo)
        mean_in(end+1)  = mi; %#ok<AGROW>
        mean_out(end+1) = mo; %#ok<AGROW>
    end

    % all pairwise distances & normalized diffs (PF ignored)
    if numel(rates) >= 2
        Dxy = squareform(pdist(cents,'euclidean'));
        [I,J] = find(triu(true(size(Dxy)),1));
        dij   = Dxy(sub2ind(size(Dxy), I, J));
        ri    = rates(I);
        rj    = rates(J);
        good  = isfinite(ri) & isfinite(rj) & isfinite(dij);
        dij   = dij(good);
        ndiff = abs(ri(good) - rj(good)) ./ max(eps, 0.5*(ri(good) + rj(good)));
        pool_dist  = [pool_dist;  dij];   %#ok<AGROW>
        pool_ndiff = [pool_ndiff; ndiff]; %#ok<AGROW>
    end

    % per-cell record
    perCell(c).rates       = rates;
    perCell(c).centroids   = cents;
    perCell(c).inPF_trial  = inPF;
    perCell(c).pf_center   = pf_center;
    perCell(c).dist_to_pf  = dist_to_pf;


end

% assemble
R = struct();
R.ratVar = ratName;
R.date   = D;
R.params = struct('WinSecs',winSecs,'SpeedThresh',vMin,'NSD',Nsd, ...
                  'MinTrials',minT,'MinSpikes',minSpk, ...
                  'CentroidMethod',centroidMethod,'CentroidBinCm',centroidBinCm);
R.perCell = perCell;

R.extra.trialRates.inPF   = trialRates_in;
R.extra.trialRates.outPF  = trialRates_out;
R.extra.perCellMean.inPF  = mean_in(:);
R.extra.perCellMean.outPF = mean_out(:);

R.pairs.dist_cm = pool_dist;
R.pairs.ndiff   = pool_ndiff;
end

% ================= helpers =================
function tok = resolveDayToken(rat, s)
% Accepts 'YYYY-MM-DD', 'YYYY_MM_DD', or 'YYYY/MM/DD'; returns whichever exists.
s = strrep(strrep(char(s),'-','_'),'/','_');
m = regexp(s,'^(\d{4})_(\d{1,2})_(\d{1,2})$','tokens','once');
if isempty(m), tok = s; return, end
yy = m{1}; mm = str2double(m{2}); dd = str2double(m{3});
padded   = sprintf('%s_%02d_%02d',yy,mm,dd);
unpadded = sprintf('%s_%d_%d',yy,mm,dd);
cands = {padded, unpadded};
tok = padded;
for k = 1:numel(cands)
    cand = cands{k};
    if (isfield(rat.Ca_peaks, ['CA_peaks_' cand])) || ...
       (isfield(rat.pos,      ['pos_'      cand])) || ...
       (isfield(rat.CS_times, ['CS_'       cand]))
        tok = cand; return
    end
end
end

function st = getCellSpikes_local(S, c)
  if iscell(S), st = S{c}(:); else, st = S(c,:).'; end
  st = st(~isnan(st) & st>0);
end

function [Xedges,Yedges,xCtr,yCtr] = make_xy_edges(x, y, binCm)
% Build edges covering visited range with ~binCm spacing; return centers.
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
Xedges = xmin:binCm:(xmin + ceil((xmax-xmin)/binCm)*binCm);
Yedges = ymin:binCm:(ymin + ceil((ymax-ymin)/binCm)*binCm);
if numel(Xedges)<2, Xedges = [xmin xmin+binCm]; end
if numel(Yedges)<2, Yedges = [ymin ymin+binCm]; end
xCtr = (Xedges(1:end-1)+Xedges(2:end))/2;
yCtr = (Yedges(1:end-1)+Yedges(2:end))/2;
end

function C = compute_trial_centroids(ts, xy, trialWindows, method, vt, vmag, ...
                                     Xedges, Yedges, xCtr, yCtr, dtAll)
% Returns [T x 2] representative positions per trial window.
T = size(trialWindows,1);
C = nan(T,2);
for tr = 1:T
    w0 = trialWindows(tr,1); w1 = trialWindows(tr,2);
    maskT = ts>=w0 & ts<w1;
    if nnz(maskT) < 3, continue, end

    switch method
        case 'mean'
            C(tr,:) = mean(xy(maskT,:),1,'omitnan');

        case 'speed'
            v_on_ts = interp1(vt, vmag, ts(maskT), 'linear', 'extrap');
            w = v_on_ts(:); w(~isfinite(w)) = 0;
            X = xy(maskT,1); Y = xy(maskT,2);
            C(tr,:) = [sum(X.*w)/max(eps,sum(w)), sum(Y.*w)/max(eps,sum(w))];

        case 'medoid'
            P = xy(maskT,:);
            if size(P,1) < 2, C(tr,:) = P(1,:); continue, end
            D = squareform(pdist(P,'euclidean'));
            [~,imin] = min(sum(D,2));
            C(tr,:) = P(imin,:);

          case 'mode'
              % dwell-time weights for irregular sampling
              wT = dtAll(maskT);
              counts = histcounts2_weighted(xy(maskT,1), xy(maskT,2), Xedges, Yedges, wT);
              if ~any(counts(:))  % rare fallback to unweighted
                  counts = histcounts2_weighted(xy(maskT,1), xy(maskT,2), Xedges, Yedges, ones(nnz(maskT),1));
              end
              [~,imax] = max(counts(:));
              [by,bx] = ind2sub(size(counts), imax);
              C(tr,:) = [xCtr(bx), yCtr(by)];


        otherwise
            error('Unknown CentroidMethod: %s', method);
    end
end
end

function C = histcounts2_weighted(x, y, Xedges, Yedges, w)
% Version-safe weighted 2D histogram.
% C is [numYbins x numXbins], matching histcounts2's counts orientation.
bx = discretize(x, Xedges);
by = discretize(y, Yedges);
w  = w(:);
good = isfinite(bx) & isfinite(by) & isfinite(w);
nx = numel(Xedges) - 1;
ny = numel(Yedges) - 1;
C = accumarray([by(good), bx(good)], w(good), [ny, nx], @sum, 0);
end
