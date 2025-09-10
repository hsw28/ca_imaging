function plotTaskVar_InOut_vs_Distance(ratNames, varargin)
% plotTaskVar_InOut_vs_Distance
% Figure with three columns and (#rats + 1 pooled) rows:
%   - Col 1: Paired per-cell means (in-PF vs out-PF) [cells with both]
%   - Col 2: Bar + error bars (mean ± SEM) of ALL trial rates (in-PF vs out-PF)
%            [includes cells that only fired in one condition]; Welch t-test
%   - Col 3: Distance between trial centroids vs variability (|Δ|/mean), all pairs
%
% Usage:
%   plotTaskVar_InOut_vs_Distance({'rat0222','rat0314'}, 'Days',[], 'WinSecs',[0 2], ...)
%
% Pass-through options: 'Days', 'WinSecs',[0 2], 'SpeedThresh',4, 'NSD',1,
% 'MinTrials',10, 'MinSpikes',3

% ---------- args ----------
p = inputParser;
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty(d));
addParameter(p,'WinSecs',[0 2], @(v) isnumeric(v)&&numel(v)==2&&v(2)>v(1));
addParameter(p,'SpeedThresh',4, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'NSD',1, @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)); % number of SDs for PF mask
addParameter(p,'MinTrials',10, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'MinSpikes',0, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'CentroidMethod','mode', ...
    @(s) any(validatestring(lower(s),{'mode','medoid','mean','speed'})));
addParameter(p,'CentroidBinCm',2.5, @(x) isnumeric(x)&&isscalar(x)&&x>0);
parse(p,varargin{:});
daysArg        = p.Results.Days;
winSecs        = p.Results.WinSecs;
vMin           = p.Results.SpeedThresh;
Nsd            = p.Results.NSD;
minT           = p.Results.MinTrials;
minSpk         = p.Results.MinSpikes;
centroidMethod = lower(p.Results.CentroidMethod);
centroidBinCm  = p.Results.CentroidBinCm;

if ischar(ratNames) || isstring(ratNames)
    ratNames = cellstr(ratNames);
end
nRats = numel(ratNames);

% Containers for pooled row
pooled_in_means   = [];
pooled_out_means  = [];
pooled_in_trials  = [];
pooled_out_trials = [];
pooled_dist       = [];
pooled_cv         = [];
pooled_dist_pf    = [];
pooled_rates_all  = [];

% Precompute results per rat
ratResults = cell(1, nRats);
for r = 1:nRats
    rat = ratNames{r};
    R = getRatR(rat, daysArg, winSecs, vMin, Nsd, minT, minSpk, centroidMethod, centroidBinCm);
    ratResults{r} = R;
end

% ---------- build figure ----------
nRows = nRats + 1;
nCols = 1;  % two active tiles per row
figure('Color','w','Position',[100 50 1500 320 + 220*nRows]);
t = tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');

for r = 1:nRats
    R = ratResults{r};

    % --- Col 1: paired per-cell in-PF vs out-PF means ---
%    [inVec, outVec] = collectInOutPerCell(R);
%    nexttile; hold on; box on
%    if isempty(inVec)
%        text(0.5,0.5,'No cells with both in/out-PF','HorizontalAlignment','center');
%        axis off
%    else
%        pairedDotPlot(inVec, outVec);
%        [~,p_pair,~,stats] = ttest(inVec, outVec);
%        title(sprintf('%s | paired in/out  (p=%.3g, t=%.2f, n=%d)', ...
%              ratNames{r}, p_pair, stats.tstat, numel(inVec)));
%    end
%    ylabel(ratNames{r});
%    if r == nRats
%        xlabel('Condition');
%    end

    % --- Col 2: bar graph of ALL trial firing rates (unpaired; includes one-sided cells) ---
%    [inTrials, outTrials] = collectAllTrialRates(R);  % pooled across days/cells
%    nexttile; hold on; box on
%    if isempty(inTrials) && isempty(outTrials)
%        text(0.5,0.5,'No in/out trial rates','HorizontalAlignment','center'); axis off
%    else
%        mu_in   = mean(inTrials,'omitnan');
%        mu_out  = mean(outTrials,'omitnan');
%        sem_in  = std(inTrials,'omitnan') / max(1,sqrt(nnz(isfinite(inTrials))));
%        sem_out = std(outTrials,'omitnan')/ max(1,sqrt(nnz(isfinite(outTrials))));
%        bar(1, mu_in,  'FaceColor',[0.7 0.7 0.9], 'EdgeColor','k'); hold on
%        bar(2, mu_out, 'FaceColor',[0.9 0.7 0.7], 'EdgeColor','k');
%        errorbar([1 2],[mu_in mu_out],[sem_in sem_out],'k','LineStyle','none','LineWidth',1.2);
%        xlim([0.5 2.5]); xticks([1 2]); xticklabels({'In-PF','Out-PF'});
%        ylabel('In-trial firing rate (Hz)');
%        if ~isempty(inTrials) && ~isempty(outTrials)
%            [~,p_unp] = ttest2(inTrials, outTrials, 'Vartype','unequal'); % Welch
%            title(sprintf('%s | ALL trials (unpaired): p=%.3g  (n_{in}=%d, n_{out}=%d)', ...
%                  ratNames{r}, p_unp, numel(inTrials), numel(outTrials)));
%        else
%            title(sprintf('%s | ALL trials (unpaired): insufficient data', ratNames{r}));
%        end
%    end

    % --- ACTIVE Col A: distance vs variability (PF ignored) ---
    [distAll, cvAll] = collectDistVsVarAllTrials(R);
%{
    axA = nexttile;  % <— capture axes
    if isempty(distAll)
        box(axA,'off'); axis(axA,'off');
        text(axA,0.5,0.5,'No trial pairs','HorizontalAlignment','center');
    else
        tallBinScatterFit(axA, distAll, cvAll, 60, 'Nshow', 5e4);
        xlabel(axA,'Centroid distance between trials (cm)');
        ylabel(axA,'|Δ rate| / mean');
        [m,b,rho,pv] = fitLineAndCorr(distAll, cvAll);
        xx = linspace(min(distAll), max(distAll), 200);
        hold(axA,'on'); plot(axA, xx, m*xx + b, 'k-', 'LineWidth',1.2);
        title(axA, sprintf('%s | r=%.2f, p=%.3g  (n=%d pairs)', ratNames{r}, rho, pv, numel(cvAll)));
    end
%}
    % --- ACTIVE Col B: rate vs distance to PF center ---
    [distPF, trialRates] = collectRateVsPFdistance(R);

    axB = nexttile;  % <— capture axes
    if isempty(distPF)
        box(axB,'off'); axis(axB,'off');
        text(axB,0.5,0.5,'No rate/dist-to-PF data','HorizontalAlignment','center');
    else
        tallBinScatterFit(axB, distPF, trialRates, 60, 'Nshow', 5e4);
        xlabel(axB,'Distance to PF center (cm)');
        ylabel(axB,'In-trial firing rate (Hz)');
        [m2,b2,r2,p2] = fitLineAndCorr(distPF, trialRates);
        xx2 = linspace(min(distPF), max(distPF), 200);
        hold(axB,'on'); plot(axB, xx2, m2*xx2 + b2, 'k-', 'LineWidth',1.2);
        title(axB, sprintf('%s | rate vs PF-dist  r=%.2f, p=%.3g  (n=%d trials)', ...
                ratNames{r}, r2, p2, numel(trialRates)));
    end

    % add to pooled (keep these even if Col 1/2 are commented)
    [inVec, outVec] = collectInOutPerCell(R);
    [inTrials, outTrials] = collectAllTrialRates(R);
    pooled_in_means   = [pooled_in_means;  inVec];          %#ok<AGROW>
    pooled_out_means  = [pooled_out_means; outVec];         %#ok<AGROW>
    pooled_in_trials  = [pooled_in_trials;  inTrials(:)];   %#ok<AGROW>
    pooled_out_trials = [pooled_out_trials; outTrials(:)];  %#ok<AGROW>
    pooled_dist       = [pooled_dist; distAll];             %#ok<AGROW>
    pooled_cv         = [pooled_cv;   cvAll];               %#ok<AGROW>
    pooled_dist_pf    = [pooled_dist_pf;   distPF(:)];      %#ok<AGROW>
    pooled_rates_all  = [pooled_rates_all; trialRates(:)];  %#ok<AGROW>
end

% ---------- pooled row ----------
% Col 1 pooled (paired means)
%nexttile; hold on; box on
%if isempty(pooled_in_means)
%    text(0.5,0.5,'No pooled cells with both in/out-PF','HorizontalAlignment','center'); axis off
%else
%    pairedDotPlot(pooled_in_means, pooled_out_means);
%    [~,p_pair,~,stats] = ttest(pooled_in_means, pooled_out_means);
%    title(sprintf('POOLED | paired in/out  (p=%.3g, t=%.2f, n=%d)', ...
%          p_pair, stats.tstat, numel(pooled_in_means)));
%end
%ylabel('POOLED'); xlabel('Condition');

% Col 2 pooled (ALL trials unpaired)
%nexttile; hold on; box on
%if isempty(pooled_in_trials) && isempty(pooled_out_trials)
%    text(0.5,0.5,'No pooled trial rates','HorizontalAlignment','center'); axis off
%else
%    mu_in   = mean(pooled_in_trials,'omitnan');
%    mu_out  = mean(pooled_out_trials,'omitnan');
%    sem_in  = std(pooled_in_trials,'omitnan') / max(1,sqrt(nnz(isfinite(pooled_in_trials))));
%    sem_out = std(pooled_out_trials,'omitnan')/ max(1,sqrt(nnz(isfinite(pooled_out_trials))));
%    bar(1, mu_in,  'FaceColor',[0.7 0.7 0.9], 'EdgeColor','k'); hold on
%    bar(2, mu_out, 'FaceColor',[0.9 0.7 0.7], 'EdgeColor','k');
%    errorbar([1 2],[mu_in mu_out],[sem_in sem_out],'k','LineStyle','none','LineWidth',1.2);
%    xlim([0.5 2.5]); xticks([1 2]); xticklabels({'In-PF','Out-PF'});
%    ylabel('In-trial firing rate (Hz)');
%    if ~isempty(pooled_in_trials) && ~isempty(pooled_out_trials)
%        [~,p_unp] = ttest2(pooled_in_trials, pooled_out_trials, 'Vartype','unequal');
%        title(sprintf('POOLED | ALL trials (unpaired): p=%.3g  (n_{in}=%d, n_{out}=%d)', ...
%              p_unp, numel(pooled_in_trials), numel(pooled_out_trials)));
%    else
%        title('POOLED | ALL trials (unpaired): insufficient data');
%    end
%end

% Col 3 pooled (distance vs variability)
%{
axP1 = nexttile;                           % << capture the axes
if isempty(pooled_dist)
    box(axP1,'off'); axis(axP1,'off');
    text(axP1,0.5,0.5,'No pooled trial pairs','HorizontalAlignment','center');
else
    tallBinScatterFit(axP1, pooled_dist, pooled_cv, 80, 'Nshow', 1e5);  % << pass ax
    xlabel(axP1,'Centroid distance between trials (cm)');
    ylabel(axP1,'|Δ rate| / mean');
    [m,b,rho,pv] = fitLineAndCorr(pooled_dist, pooled_cv);
    xx = linspace(min(pooled_dist), max(pooled_dist), 300);
    hold(axP1,'on'); plot(axP1, xx, m*xx + b, 'k-', 'LineWidth',1.2);
    title(axP1, sprintf('POOLED | r=%.2f, p=%.3g  (n=%d pairs)', rho, pv, numel(pooled_cv)));
end
%}

% NEW pooled: firing rate vs distance from PF center

axP2 = nexttile;
if isempty(pooled_dist_pf)
    box off; axis off
    text(0.5,0.5,'No pooled rate/dist-to-PF data','HorizontalAlignment','center');
else
    tallBinScatterFit(axP2, pooled_dist_pf, pooled_rates_all, 80, 'Nshow', 1e5);
    xlabel('Distance to PF center (cm)');
    ylabel('In-trial firing rate (Hz)');
    [m2,b2,r2,p2] = fitLineAndCorr(pooled_dist_pf, pooled_rates_all);
    xx2 = linspace(min(pooled_dist_pf), max(pooled_dist_pf), 300);
    hold on; plot(xx2, m2*xx2 + b2, 'k-', 'LineWidth',1.2);
    title(sprintf('POOLED | rate vs PF-dist  r=%.2f, p=%.3g  (n=%d trials)', ...
          r2, p2, numel(pooled_rates_all)));
end

% (keeping your original title line, but adding a clearer one below)
% title(t, 'Task firing variability: (1) Paired per-cell means  |  (2) ALL trials bar (unpaired)  |  (3) Distance vs Variability', ...
%       'FontWeight','bold');
title(t, 'Task firing variability: (A) Distance vs Variability  |  (B) Rate vs Distance-to-PF', ...
      'FontWeight','bold');

end % main function

% -------------------- helpers --------------------

function R = getRatR(ratName, daysArg, winSecs, vMin, Nsd, minT, minSpk, centroidMethod, centroidBinCm)
R = taskFiringVariability(ratName, ...
        'Days',daysArg, ...
        'WinSecs',winSecs, 'SpeedThresh',vMin, 'NSD',Nsd, ...
        'MinTrials',minT, 'MinSpikes',minSpk, ...
        'CentroidMethod',centroidMethod, 'CentroidBinCm',centroidBinCm, ...
        'Plot',false);
end

function [inVec, outVec] = collectInOutPerCell(R)
% Per-cell paired means (only cells with >=1 trial in both in/out)
inVec  = [];
outVec = [];
if ~isfield(R,'perDay') || isempty(R.perDay), return, end
for d = 1:numel(R.perDay)
    PC = R.perDay(d).perCell;
    if isempty(PC), continue, end
    for c = 1:numel(PC)
        if ~isfield(PC(c),'rates') || isempty(PC(c).rates) || ~isfield(PC(c),'inPF_trial')
            continue
        end
        rates = PC(c).rates(:);
        inpf  = PC(c).inPF_trial(:);
        if numel(inpf) ~= numel(rates), continue, end
        mi = mean(rates(inpf==1), 'omitnan');
        mo = mean(rates(inpf==0), 'omitnan');
        ni = nnz(inpf==1 & isfinite(rates));
        no = nnz(inpf==0 & isfinite(rates));
        if isfinite(mi) && isfinite(mo) && ni>=1 && no>=1
            inVec  = [inVec;  mi]; %#ok<AGROW>
            outVec = [outVec; mo]; %#ok<AGROW>
        end
    end
end
end

function [inTrials, outTrials] = collectAllTrialRates(R)
% ALL individual trial rates pooled (includes one-sided cells)
inTrials  = [];
outTrials = [];
if isfield(R,'extra') && isfield(R.extra,'trialRates')
    if isfield(R.extra.trialRates,'inPF'),  inTrials  = R.extra.trialRates.inPF;  end
    if isfield(R.extra.trialRates,'outPF'), outTrials = R.extra.trialRates.outPF; end
end
inTrials  = inTrials(isfinite(inTrials));
outTrials = outTrials(isfinite(outTrials));
end

function [distAll, cvAll] = collectDistVsVarAllTrials(R)
% For each perDay/perCell, form all trial pairs (i<j):
%   distance = euclidean distance between trial centroids
%   variability = |Δ rate| / mean(rate_i, rate_j)
% Ignores in/out PF.
distAll = [];
cvAll   = [];
if ~isfield(R,'perDay') || isempty(R.perDay), return, end
for d = 1:numel(R.perDay)
    PC = R.perDay(d).perCell;
    if isempty(PC), continue, end
    for c = 1:numel(PC)
        if ~isfield(PC(c),'rates') || isempty(PC(c).rates) || ...
           ~isfield(PC(c),'centroids') || isempty(PC(c).centroids)
            continue
        end
        rates = PC(c).rates(:);
        cents = PC(c).centroids;
        keep  = isfinite(rates) & all(isfinite(cents),2);
        rates = rates(keep);
        cents = cents(keep,:);
        if numel(rates) < 2, continue, end
        D = squareform(pdist(cents,'euclidean'));
        [I,J] = find(triu(true(size(D)),1));
        dvec  = D(sub2ind(size(D), I, J));
        diffA = abs(rates(I) - rates(J));
        cv    = diffA ./ max(eps, 0.5*(rates(I)+rates(J)));
        distAll = [distAll; dvec(:)]; %#ok<AGROW>
        cvAll   = [cvAll;   cv(:)];   %#ok<AGROW>
    end
end
end

function pairedDotPlot(inVec, outVec)
n = numel(inVec);
x1 = ones(n,1);
x2 = 2*ones(n,1);
for k = 1:n
    plot([x1(k) x2(k)], [inVec(k) outVec(k)], '-', 'Color',[0 0 0 0.2]); hold on
end
scatter(x1, inVec, 18, 'filled');
scatter(x2, outVec,18, 'filled');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'In-PF','Out-PF'});
ylabel('In-trial firing rate (Hz)');
end

function [m,b,rho,pv] = fitLineAndCorr(x,y)
x = x(:); y = y(:);
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);
if numel(x) < 3
    m = NaN; b = NaN; rho = NaN; pv = NaN; return
end
P = polyfit(x, y, 1);
m = P(1); b = P(2);
[rho,pv] = corr(x,y,'Rows','complete','Type','Pearson');
end

function [distPF, rates] = collectRateVsPFdistance(R)
% Concatenate per-trial distance-to-PF-center and trial firing rate.
% Robust to length mismatches by recomputing from centroids/pf_center when available,
% otherwise truncates to the common minimum length.

distPF = [];
rates  = [];
if ~isfield(R,'perDay') || isempty(R.perDay), return, end

for d = 1:numel(R.perDay)
    PC = R.perDay(d).perCell;
    if isempty(PC), continue, end

    for c = 1:numel(PC)
        if ~isfield(PC(c),'rates') || isempty(PC(c).rates)
            continue
        end
        r = PC(c).rates(:);

        dp = [];  % candidate distances
        have_dp = isfield(PC(c),'dist_to_pf') && ~isempty(PC(c).dist_to_pf);
        have_cent = isfield(PC(c),'centroids') && ~isempty(PC(c).centroids);
        have_center = isfield(PC(c),'pf_center') && numel(PC(c).pf_center)==2 && all(isfinite(PC(c).pf_center));

        if have_dp
            dp = PC(c).dist_to_pf(:);
        end

        % If lengths mismatch, try to recompute from centroids + pf_center
        if isempty(dp) || numel(dp) ~= numel(r)
            if have_cent && have_center
                cents = PC(c).centroids;
                % align lengths to be safe
                n = min(numel(r), size(cents,1));
                if n < 1, continue, end
                r = r(1:n);
                cents = cents(1:n,:);
                ctr = PC(c).pf_center(:).';  % [1x2]
                dp = sqrt(sum((cents - ctr).^2, 2));
            else
                % Last resort: truncate to common min length if we at least have dp
                if ~isempty(dp)
                    n = min(numel(r), numel(dp));
                    if n < 1, continue, end
                    r  = r(1:n);
                    dp = dp(1:n);
                else
                    % No way to form matched pairs
                    continue
                end
            end
        end

        % Final clean
        good = isfinite(r) & isfinite(dp);
        if any(good)
            distPF = [distPF; dp(good)]; %#ok<AGROW>
            rates  = [rates;  r(good)];  %#ok<AGROW>
        end
    end
end
end


function tallBinScatterFit(ax, x, y, nbins, varargin)
% Efficient scatter for large N with binned median ± IQR overlay.
ip = inputParser;
addParameter(ip,'Nshow',5e4,@(z)isnumeric(z)&&isscalar(z)&&z>=0);
parse(ip,varargin{:});
Nshow = ip.Results.Nshow;

hold(ax,'on'); box(ax,'on');

x = x(:); y = y(:);
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);
n = numel(x);

% backdrop subsample
if Nshow>0 && n>Nshow, idx = randperm(n, Nshow); else, idx = 1:n; end
scatter(ax, x(idx), y(idx), 10, 'filled', 'MarkerFaceAlpha',0.25);

% edges = linspace(min(x), max(x), nbins+1);
% cent  = (edges(1:end-1)+edges(2:end))/2;
% b = discretize(x, edges);
% med = nan(nbins,1); q25 = med; q75 = med;
% for k = 1:nbins
%     yk = y(b==k);
%     if ~isempty(yk)
%         med(k) = median(yk,'omitnan');
%         q = quantile(yk,[0.25 0.75]);
%         q25(k) = q(1); q75(k) = q(2);
%     end
% end
% plot(ax, cent, med, 'k-', 'LineWidth',2);
% plot(ax, cent, q25, ':', 'LineWidth',1, 'Color',[0.6 0.6 0.6]);
% plot(ax, cent, q75, ':', 'LineWidth',1, 'Color',[0.6 0.6 0.6]);
end
