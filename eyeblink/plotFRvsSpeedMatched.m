function plotFRvsSpeedMatched

%for each point in a trial:
% Look up its velocity
%Find all background bins whose velocity is the same (velBinCand == vb).
%Randomly sample the same number of background bins as there are trials (needN).
%Count spikes inside those sampled bins (using the patched loop sum(sTimes >= t0 & sTimes < t0+binSize)).
%Store those spike counts into FR_match.
%
%After we have looped over all velocity categories present in the trials,
%FR_trial (spikes per bin inside the CS window) and
%FR_match (spikes per velocity-matched non-trial bin) have the same length
%and very similar speed distribution.

%Statistics
%Paired comparison (default) – we use a paired-samples t-test
%between FR_trial and FR_match.
%If p < α (default 0.05) the neuron is deemed “velocity-modulated within
%trials”.

%We also compute ΔFR = mean(FR_trial) – mean(FR_match) so you can see
%the direction (positive = higher firing in trials).

%returns
%pctSigDiff:	percentage of neurons in that day whose paired test was significant (p < α)
%deltaFR:	one ΔFR value per velocity-modulated neuron
%pEach:	the t-test p-value for every neuron



% -------------------------------------------------------------------------
%  For each of the five rats:
%     • walks through the same 3 “training-window” days you already use
%     • calls trialVsSpeedMatched  (zero shuffling)
%     • stores the % of velocity-modulated cells whose firing during trials
%       STILL differs from velocity-matched non-trial bins
%
%  Generates three summary graphics:
%    • For each rat and each of the 3 “training-window” days
%      – counts how many velocity-modulated cells still differ in firing
%        when compared with speed-matched non-trial epochs
%    • Builds three graphics:
%        (1) bar-plot per day              (all rats)
%        (2) bar-plot of the rat means     (one bar / rat)
%        (3) swarm of ΔFR  (trial – speed matched) for every cell
%

ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};      % <-- change back to full list if desired
nRats    = numel(ratNames);

% analysis parameters -----------------------------------------------------
win      = [0 2];          % trial window [s] relative to CS
binSize  = 1/7.5;          % 133 ms
alpha    = 0.05;           % significance level

% collectors --------------------------------------------------------------
pctPerDay   = nan(nRats,3);      % rat × day(-2,-1,0)
pctPerRat   = nan(nRats,1);
allDeltaFR  = cell(nRats,1);

for r = 1:nRats
    rat = evalin('base',ratNames{r});      % assumes struct in base workspace
    dateList = autoDateList(rat);          % helper from your pipeline
    idx      = find(strcmp(dateList, rat.An));
    theseDays = dateList(idx-2:idx);       % 3-day training window

    deltaTmp = [];

    for d = 1:3
        dateStr   = theseDays{d};

        spikeMat  = rat.Ca_peaks.(['CA_peaks_' dateStr]);
        ts        = rat.pos.(['pos_' dateStr])(:,1);
        pos       = rat.pos.(['pos_' dateStr])(:,2:3);
        csTimes   = rat.CS_times.(['CS_' dateStr]);

        % ---------  single-day stats  ---------
        S = trialVsSpeedMatched(spikeMat, ts, pos, csTimes, ...
                'win',win,'binSize',binSize,'alpha',alpha);

        pctPerDay(r,d)  = S.pctSigDiff;
        deltaTmp        = [deltaTmp ; S.deltaFR(:)];   %#ok<AGROW>
    end

    pctPerRat(r)  = mean(pctPerDay(r, :), 'omitnan');
    allDeltaFR{r} = deltaTmp;
end

%% -----------  plotting  --------------------------------------------------
figure('Color','w','Position',[200 300 1400 400]);

% (1) every day
subplot(1,3,1); hold on;
bar(pctPerDay','grouped');
xticks(1:nRats); xticklabels(ratNames);
ylabel('% vel-mod cells diff');
title('% significant  –  every day');
legend({'day-2','day-1','day-0'},'Location','northwest');

% (2) rat means
subplot(1,3,2); hold on;
bar(pctPerRat,'FaceColor',[0.2 0.6 0.8]);
xticks(1:nRats); xticklabels(ratNames);
ylabel('mean % across 3 days');
title('Rat-wise mean');

% (3) ΔFR distribution
subplot(1,3,3); hold on;
for r = 1:nRats
    swarmchart(r*ones(size(allDeltaFR{r})), allDeltaFR{r}, ...
               5,'filled','MarkerFaceAlpha',0.4);
end
plot(xlim, [0 0],'k--');
xticks(1:nRats); xticklabels(ratNames);
ylabel('\DeltaFR  (trial – matched) [spk/bin]');
title('\DeltaFR distribution');

%% ------------ console summary -------------------------------------------
fprintf('\n===========  SUMMARY  ===========\n');
for r = 1:nRats
    fprintf('%s  –  %.1f ± %.1f %% (across 3 days)\n',...
        ratNames{r}, pctPerRat(r), std(pctPerDay(r,:),[],'omitnan'));
end
fprintf('Grand mean:  %.1f %% ± %.1f\n',...
        mean(pctPerRat,'omitnan'), std(pctPerRat,'omitnan'));
end
% ======================================================================
%                       Helper functions
% ======================================================================

function stats = trialVsSpeedMatched(spikeMat, ts, pos, csTimes, varargin)
% Compare firing in trial bins with velocity-matched non-trial bins
% (zero shuffling; purely deterministic matching).
% ----------------------------------------------------------------------

% ---- defaults --------------------------------------------------------
p = inputParser;
p.addParameter('win',      [0 2]);
p.addParameter('binSize',   1/7.5);
p.addParameter('nVelBins',  15);
p.addParameter('alpha',     0.05);
p.addParameter('test',     'signrank');
p.parse(varargin{:});
win      = p.Results.win;
binSize  = p.Results.binSize;
nVelBins = p.Results.nVelBins;
alpha    = p.Results.alpha;
testType = lower(p.Results.test);

% ---- instantaneous speed --------------------------------------------
dt    = diff(ts);
dx    = diff(pos);
speed = [0; hypot(dx(:,1),dx(:,2))./dt];
speed(~isfinite(speed)) = 0;

% ---- mask points NOT inside any trial -------------------------------
inTrial = false(size(ts));
for t = 1:numel(csTimes)
    inTrial = inTrial | (ts >= csTimes(t)+win(1) & ts <= csTimes(t)+win(2));
end
outMask = ~inTrial;

speedOut = speed(outMask);
timeOut  = ts(outMask);

velEdges = linspace(min(speed), max(speed), nVelBins+1);

% ---- pre-bin trial speed --------------------------------------------
nBins   = round(diff(win)/binSize);
nTrials = numel(csTimes);

binnedSpeed = nan(nTrials,nBins);
edgeMat     = nan(nTrials,nBins+1);
for t = 1:nTrials
    edges          = linspace(csTimes(t)+win(1), csTimes(t)+win(2), nBins+1);
    edgeMat(t,:)   = edges;
    tk             = ts>=edges(1) & ts<edges(end);
    spdT           = speed(tk);
    [~,~,binIdx]   = histcounts(ts(tk),edges);
    for b = 1:nBins
        binnedSpeed(t,b) = mean(spdT(binIdx==b));
    end
end

% ---- outputs ---------------------------------------------------------
nCells    = size(spikeMat,1);
deltaFR   = nan(nCells,1);
pVal      = nan(nCells,1);
isVelMod  = false(nCells,1);
sigDiff   = false(nCells,1);

% ---- main cell loop --------------------------------------------------
for c = 1:nCells
    sTimes = spikeMat(c,:);    sTimes = sTimes(~isnan(sTimes));
    if numel(sTimes)<5, continue; end

    % (1) velocity-modulation within trials
    [isMod,~,~] = findVelocityModCells( ...
        sTimes, ts, pos, csTimes, win, binSize, 0.05);

    if ~isMod, continue; end
    isVelMod(c) = true;

    % (2) FR in trials & matched bins
    FR_trial = [];
    FR_match = [];

    for t = 1:nTrials
        edges      = edgeMat(t,:);
        spdThis    = binnedSpeed(t,:);
        good       = ~isnan(spdThis);

        % spikes in real trial bins
        FR_trial = [FR_trial histcounts(sTimes, edges)]; %#ok<AGROW>

        % choose velocity-matched out-of-trial bins
        [~,~,binID_trial] = histcounts(spdThis, velEdges);

        badMask = timeOut>=edges(1) & timeOut<=edges(end);
        candIdx = find(~badMask);
        velBinCand = discretize(speedOut(candIdx), velEdges);

        for vb = unique(binID_trial(good))
            needN = sum(binID_trial==vb & good);
            pool  = candIdx(velBinCand==vb);

            if numel(pool)<needN         % sample with replacement
                pool = [pool ; randsample(pool, needN-numel(pool), true)];
            else                         % sample without replacement
                pool = randsample(pool, needN, false);
            end

            for k = 1:numel(pool)
                t0 = timeOut(pool(k));
                FR_match(end+1,1) = sum(sTimes>=t0 & sTimes<t0+binSize); %#ok<AGROW>
            end
        end
    end

    if numel(FR_trial)<10 || numel(FR_match)<10, continue; end

    m = min(numel(FR_trial), numel(FR_match));
    switch testType
        case 'ttest'
            [~,p_] = ttest(FR_trial(1:m), FR_match(1:m));
        otherwise
            p_ = signrank(FR_trial(1:m), FR_match(1:m));
    end

    pVal(c)    = p_;
    deltaFR(c) = mean(FR_trial) - mean(FR_match);
    sigDiff(c) = p_ < alpha;
end

% ---- package stats ---------------------------------------------------
stats.deltaFR     = deltaFR;
stats.pVal        = pVal;
stats.isVelMod    = isVelMod;
stats.sigDiff     = sigDiff;
stats.pctSigDiff  = 100*mean(sigDiff(isVelMod));

fprintf('Velocity-mod cells: %d / %d\n', sum(isVelMod), nCells);
fprintf('  …of those, %.1f %% differ (alpha %.2f)\n',...
        stats.pctSigDiff, alpha);
end
% ---------------------------------------------------------------------
function [isVelMod, meanR, pVal] = findVelocityModCells( ...
         spikeTimes, ts, pos, csTimes, win, binSize, pThresh)
% Test whether ONE neuron is speed-modulated inside trial epochs
% -------------------------------------------------------------------------

% defaults
if nargin<5, win     = [0 2];   end
if nargin<6, binSize = 1/7.5;   end
if nargin<7, pThresh = 0.05;    end

xy = pos;

% instantaneous speed ----------------------------------------------------
dt    = diff(ts);
dx    = diff(xy);
speed = [0; hypot(dx(:,1),dx(:,2))./dt];
speed(~isfinite(speed)) = 0;

% parameters -------------------------------------------------------------
nBins   = round(diff(win)/binSize);
nTrials = numel(csTimes);

rPerTrial = nan(nTrials,1);

for t = 1:nTrials
    edges = linspace(csTimes(t)+win(1), csTimes(t)+win(2), nBins+1);

    % assign every time-stamp to a bin (NaN if out of range)
    inBin = discretize(ts, edges);

    % ----------- speed per bin (exclude NaNs!) -----------
    valid   = ~isnan(inBin);
    spdMean = accumarray(inBin(valid), speed(valid), [nBins 1], @mean, NaN);

    % ----------- spikes per bin --------------------------
    spkCounts = histcounts(spikeTimes, edges);

    good = ~isnan(spdMean) & ~isnan(spkCounts(:));
    if sum(good) > 2
        rPerTrial(t) = corr(spdMean(good), spkCounts(good)', 'type','Pearson');
    end
end

rPerTrial = rPerTrial(~isnan(rPerTrial));
if numel(rPerTrial) < 3
    isVelMod = false; meanR = NaN; pVal = NaN;
    return
end

[~,p]   = ttest(rPerTrial);
meanR   = mean(rPerTrial);
pVal    = p;
isVelMod = p < pThresh;
end
