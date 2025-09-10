function R = epochPSTH_ZscoreHists(ratNames, varargin)
% epochPSTH_ZscoreHists  Per-epoch z-scored PSTH histograms (no clustering)
%
% Z-score logic (per cell):
%   1) Build trial-aligned PSTH (133 ms bins) from -5 to +2 s.
%   2) Baseline = all bins in [-5,0) across all trials.
%   3) z = (binCount - mu_base) / (sigma_base + eps0).
%   4) For each epoch (CS/Trace/US/Post), take the cell's mean z across bins.
%
% Usage:
%   R = epochPSTH_ZscoreHists({'rat0222','rat0307','rat0313','rat0314','rat0816'});
%
% Options:
%   'Bin'         : 1/7.5           % 133 ms imaging frame (sec)
%   'Window'      : [-15 2]          % peri-CS window (sec)
%   'EpochEdges'  : [-15 0 0.25 0.75 0.85 2.00]  % baseline + 4 epochs
%   'MinTrials'   : 20              % require at least this many trials
%   'MinBaseSpk'  : 0               % min baseline spikes across all trials
%   'Eps'         : 1e-6            % variance stabilizer for z
%   'Plot'        : true
%   'BinsHist'    : -3:0.25:8       % histogram x-axis for z
%
% Returns:
%   R.zEpoch   : [nCells x 4] mean z per epoch (pooled)
%   R.labels   : {'CS','Trace','US','Post'}
%   R.fracPos  : fraction > 0 per epoch (rat-weighted)
%   R.nCells   : total cells included; R.nRats
%
% Notes:
%   - Counts per bin (not rates) are z-scored; this matches Poisson-ish
%     intuition and avoids tiny-window rate artifacts.
%   - Cells with near-zero baseline variance get eps in denominator.

p = inputParser;
addParameter(p,'Bin',1/7.5,@isscalar);
addParameter(p,'Window',[-5 2],@(v) numel(v)==2);
addParameter(p,'EpochEdges',[-2 0 0.25 0.75 0.85 2.00]);
addParameter(p,'MinTrials',20,@isscalar);
addParameter(p,'MinBaseSpk',5,@isscalar);
addParameter(p,'Eps',1e-5,@isscalar);
addParameter(p,'Plot',true,@islogical);
addParameter(p,'BinsHist',-.5:0.25:4);
parse(p,varargin{:});
bin   = p.Results.Bin;
win   = p.Results.Window;
E     = p.Results.EpochEdges;
minTr = p.Results.MinTrials;
minBs = p.Results.MinBaseSpk;
eps0  = p.Results.Eps;
dop   = p.Results.Plot;
xbins = p.Results.BinsHist;

labels = {'CS','Trace','US','Post'};
nE = 4;

% build bin edges for the PSTH
edgesP = win(1):bin:win(2);
nBins  = numel(edgesP)-1;

% helper to map epoch -> bin indices
epochBinIdx = cell(1,nE);
for e = 1:nE
    [~,i0] = min(abs(edgesP - E(e+0)));
    [~,i1] = min(abs(edgesP - E(e+1)));
    epochBinIdx{e} = max(1,i0):min(nBins,i1-1);
end
% baseline bins:
[~,ib0] = min(abs(edgesP - E(1)));
[~,ib1] = min(abs(edgesP - E(2)));
baseBins = max(1,ib0):min(nBins,ib1-1);

Z_by_rat = cell(numel(ratNames),1);

for r = 1:numel(ratNames)
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    iAn   = find(strcmp(dates,rat.An),1);
    days  = dates(max(1,iAn-2):iAn);

    Zep = [];  % rows=cells, cols=4 epochs (mean z per epoch)

    for d = 1:numel(days)
        D   = days{d};
        S   = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));
        CS  = rat.CS_times.(sprintf('CS_%s',D));
        M   = rat.ratemask.(sprintf('ratemask_%s',D))==1;

        if numel(CS) < minTr, continue, end

        % precompute trial-relative bin edges for speed
        T0  = CS(:) + win(1);
        T1  = CS(:) + win(2);
        nTr = numel(CS);

        nC  = size(S,1);
        for c = 1:nC
            if ~M(c), continue, end
            st = S(c,:); st = st(~isnan(st) & st>0);
            if isempty(st), continue, end

            % build trial-aligned spike-count matrix: [nTr x nBins]
            Cmat = zeros(nTr, nBins, 'double');
            for t = 1:nTr
                % histogram counts per bin for this trial
                % shift spikes into trial-relative coordinates
                rel = st(st>=T0(t) & st<T1(t)) - CS(t);
                if ~isempty(rel)
                    Cmat(t,:) = histcounts(rel, edgesP);
                end
            end

            % baseline mean & std from baseline bins across trials
            baseCounts = Cmat(:, baseBins);
            sumBase = sum(baseCounts,'all');
            if sumBase < minBs, continue, end   % too few baseline events

            muB = mean(baseCounts,'all');
            sdB = std(baseCounts(:), 0, 'omitnan');

            % guard against near-zero variance
            denom = sdB + eps0;

            % z-score each bin of the PSTH (per cell)
            Z = (Cmat - muB) / denom;    % nTr x nBins

            % mean z per epoch (averaged across bins & trials)
            zEpoch = zeros(1,nE);
            for e = 1:nE
                ze = Z(:, epochBinIdx{e});
                zEpoch(e) = mean(ze(:), 'omitnan');
            end
            Zep = [Zep; zEpoch]; %#ok<AGROW>
        end
    end

    Z_by_rat{r} = Zep;
end

% pooled & rat-weighted summaries
Z_all = cat(1, Z_by_rat{:});

fracPos = zeros(1,nE);
for e = 1:nE
    fr = nan(numel(ratNames),1);
    for r = 1:numel(ratNames)
        v = Z_by_rat{r};
        if isempty(v), continue, end
        fr(r) = mean(v(:,e) > 0, 'omitnan');
    end
    fracPos(e) = mean(fr,'omitnan');
end

if dop
    figure('Color','w','Position',[100 100 1200 420]);
    tl = tiledlayout(1,5,'Padding','compact','TileSpacing','compact');

    % per-epoch z histograms
    for e = 1:nE
        nexttile; hold on
        histogram(Z_all(:,e), xbins, 'Normalization','probability');
        xline(0,'k-'); yline(0,'k-');
        xlabel(sprintf('%s z',labels{e})); ylabel('Fraction of cells');
        title(labels{e}); box off
    end

    % summary bar: fraction positive z per epoch (rat-weighted)
    nexttile; hold on
    bar(1:nE, fracPos, 0.6, 'FaceColor',[0.2 0.5 0.9], 'EdgeColor','none');
    ylim([0 1]); xlim([0.5 4.5]); yline(0.5,'k:');
    xticks(1:nE); xticklabels(labels);
    ylabel('Frac(z > 0)'); title('Rat-weighted positive fraction');
    box off

    title(tl,'Per-epoch PSTH z-scores: distributions + summary');
end

% return
R.labels  = labels;
R.zEpoch  = Z_all;
R.fracPos = fracPos;
R.nCells  = size(Z_all,1);
R.nRats   = numel(ratNames);
R.params  = p.Results;
end
