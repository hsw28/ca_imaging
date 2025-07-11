function plotFRvsSpeedWithinTrials(varargin)
% -------------------------------------------------------------------------
%  Within-trial FR-vs-speed correlations (Pearson ρ) with two layers
%  of significance testing:
%
%  (1) **Per neuron**
%      – Pearson r across ∼15 trace-window bins × (N trials).
%      – p-value from t-test of r against 0.
%      – Multiple-comparison control per rat (FDR or Bonferroni).
%      – Shuffle test: 500 speed-shuffles per neuron; neuron deemed
%        “shuffle-sig” if its mean(r_shuffles) exceeds observed |r|.
%
%  (2) **Population level (per rat)**
%      – Compare |mean(r)| to full null distribution of population-mean |r|
%        computed from the same shuffles (p(A-v-S) in console).
%
%  Two figure windows are produced:
%      • ALL neurons
%      • INCLUDED neurons      (≥ minSpk spikes in trace window)
%
%  Optional name-value pairs
%      'correction' : 'fdr' (default) | 'bonferroni'
%      'minSpk'     : minimum spikes in trace window for inclusion (default 5)
%
%  Requires: mafdr (Statistics & ML Toolbox),
%            autoDateList, ca_velocity, rat structs in base workspace.
% -------------------------------------------------------------------------

%% ---------------- USER CONSTANTS ---------------------------------------
ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
win      = [0 2];          % trace window relative to CS onset (s)
binSize  = 1/7.5;          % 133 ms Ca bins
nShuff   = 500;            % speed shuffles per neuron
alphaFW  = 0.05;           % family-wise α (per rat)
minSpk   = 10;
%% ------------------------------------------------------------------------

%% ---------- parse optional arguments -----------------------------------
p = inputParser;
p.addParameter('correction','fdr',@(s)ischar(s)&&ismember(lower(s),{'fdr','bonferroni'}));
p.addParameter('minSpk',5,@(x)isempty(x)||(isscalar(x)&&x>=0));
p.parse(varargin{:});
corrMethod = lower(p.Results.correction);
minSpk     = p.Results.minSpk;

fprintf('\n>>> Multiple-comparison: %s  (α = %.3g)\n',upper(corrMethod),alphaFW);
if isempty(minSpk)
    fprintf('>>> Inclusion: all neurons (no spike threshold).\n');
else
    fprintf('>>> Inclusion: ≥ %d spikes in trace window [%.1f %.1f] s.\n', ...
            minSpk,win);
end

%% -------- summary containers ------------------------------------------
nRats  = numel(ratNames);
blank  = struct('nTest',zeros(nRats,1),'nSig',zeros(nRats,1), ...
                'nSigShuf',zeros(nRats,1),'sigPct',zeros(nRats,1));
Sall = blank;  Sinc = blank;

corrAllPerRat = cell(nRats,1);   shMeanAllPerRat = cell(nRats,1);
corrIncPerRat = cell(nRats,1);   shMeanIncPerRat = cell(nRats,1);

pPopAll = NaN(nRats,1);  pPopInc = NaN(nRats,1);

%% ================= MAIN RAT LOOP ======================================
for r = 1:nRats
    rat      = evalin('base',ratNames{r});
    dList    = autoDateList(rat);
    idx      = find(strcmp(dList,rat.An));
    days     = dList(idx-2:idx);  % last 3 days

    % per-rat accumulators
    r_all    = [];  p_all    = [];  shMat_all = [];
    r_inc    = [];  p_inc    = [];  shMat_inc = [];

    %% ----- per-day loop -------------------------------------------
    for d = 1:3
        day    = days{d};
        spkMat = rat.Ca_peaks.(['CA_peaks_' day]);
        posRaw = rat.pos.(['pos_'       day]);  % [t x y]
        csT    = rat.CS_times.(['CS_'   day]);  % CS onsets

        % build speed trace
        vDat   = ca_velocity(posRaw');
        speed  = interp1(vDat(2,:),vDat(1,:),posRaw(:,1),'linear','extrap');
        speed(~isfinite(speed)) = 0;
        ts     = posRaw(:,1);

        nTr    = numel(csT);
        nBins  = round(diff(win)/binSize);

        % precompute bin edges & speedBins
        binEdges  = nan(nTr,nBins+1);
        speedBins = nan(nTr,nBins);
        for t = 1:nTr
            edges = linspace(csT(t)+win(1),csT(t)+win(2),nBins+1);
            binEdges(t,:)  = edges;
            for b = 1:nBins
                idxb = ts>=edges(b)&ts<edges(b+1);
                speedBins(t,b) = mean(speed(idxb));
            end
        end

        % full shuffle matrix (trials×bins×shuffles)
        shSpeed = nan(nTr,nBins,nShuff);
        for s = 1:nShuff
            for t = 1:nTr
                shSpeed(t,:,s) = speedBins(t,randperm(nBins));
            end
        end

        %% ----- per-neuron loop ----------------------------------
        for ni = 1:size(spkMat,1)
            spk = spkMat(ni,:); spk = spk(~isnan(spk));
            Sall.nTest(r) = Sall.nTest(r)+1;

            % build trial-wise r and shuffle-null
            rT   = nan(nTr,1);
            rSh  = nan(nTr,nShuff);
            nSpk = 0;

            for t = 1:nTr
                counts = histcounts(spk,binEdges(t,:));
                nSpk   = nSpk + sum(counts);
                valid  = ~isnan(counts) & ~isnan(speedBins(t,:));
                if nnz(valid)>=3
                    rT(t) = corr(counts(valid).', speedBins(t,valid).','type','Pearson');
                    for s = 1:nShuff
                        rSh(t,s) = corr(counts(valid).', shSpeed(t,valid,s).','type','Pearson');
                    end
                end
            end

            rT = rT(~isnan(rT));
            if numel(rT)<3, continue; end

            % observed mean and its t-test p
            rMu      = mean(rT);
            [~,pPear] = ttest(rT);

            % per-neuron null = mean over trials, permutation p, centered
                       shMean   = mean(rSh,1,'omitnan');         % 1×nShuff
                       mu0      = mean(shMean);
                      shZ      = shMean - mu0;                 % center null
                       obsZ     = rMu      - mu0;               % center obs
                       pPerm    = mean(abs(shZ) >= abs(obsZ));  % two-tailed permutation p

            % store ALL
            r_all(end+1)      = rMu;
            p_all(end+1)      = pPear;
            shMat_all(end+1,:) = shMean;

            % store INCLUDED
            if isempty(minSpk) || nSpk>=minSpk
                r_inc(end+1)      = rMu;
                p_inc(end+1)      = pPear;
                shMat_inc(end+1,:) = shMean;
                Sinc.nTest(r)     = Sinc.nTest(r)+1;
            end
        end
    end  % days

    % population-level shuffle p
    if ~isempty(r_all)
        nullAll    = mean(shMat_all,2);  % mean per neuron
        pPopAll(r) = mean(abs(nullAll) >= abs(mean(r_all)));
    end
    if ~isempty(r_inc)
        nullInc    = mean(shMat_inc,2);
        pPopInc(r) = mean(abs(nullInc) >= abs(mean(r_inc)));
    end

    % multiple-comparison correction
    sigP_all = applyMCC(p_all, corrMethod, alphaFW);
    sigP_inc = applyMCC(p_inc, corrMethod, alphaFW);

    sigS_all = applyMCC(shuffleP(p_all,shMat_all), corrMethod, alphaFW);
    sigS_inc = applyMCC(shuffleP(p_inc,shMat_inc), corrMethod, alphaFW);

    Sall.nSig(r)     = sum(sigP_all);
    Sall.nSigShuf(r) = sum(sigS_all);
    Sall.sigPct(r)   = 100*mean(sigP_all);

    Sinc.nSig(r)     = sum(sigP_inc);
    Sinc.nSigShuf(r) = sum(sigS_inc);
    Sinc.sigPct(r)   = 100*mean(sigP_inc);

    corrAllPerRat{r}   = r_all;
    shMeanAllPerRat{r} = mean(shMat_all,2);
    corrIncPerRat{r}   = r_inc;
    shMeanIncPerRat{r} = mean(shMat_inc,2);
end

%% -------- build summaries & visuals ------------------------------------
[sumAll,sumInc] = buildSummaries( Sall, Sinc, ...
                   corrAllPerRat, shMeanAllPerRat, ...
                   corrIncPerRat, shMeanIncPerRat );

printConsole(sumAll,sumInc, ratNames, corrMethod, pPopAll, pPopInc);
makeFigure(  sumAll, 'ALL neuron-days',     corrMethod, ratNames);
makeFigure(  sumInc, 'INCLUDED neuron-days',corrMethod, ratNames);
end

%% ======================= HELPERS =========================================

function pP = shuffleP(p,rSh)
% re-compute per-cell permutation p from rSh matrix
% rSh: [nNeurons×nShuffles]
obs   = mean(rSh,2);
pP    = mean(abs(rSh) >= abs(obs), 2);
end

function sig = applyMCC(p,method,alphaFW)
  % identical to sigTest + fdr_bh from your other code
  p(isnan(p)) = 1;
  switch method
    case 'bonferroni'
      thr = alphaFW/numel(p);
    otherwise  % 'fdr'
      thr = fdr_bh(p,alphaFW);
  end
  sig = (p < thr);
end

function thr = fdr_bh(p,q)
    p = sort(p(:));
    m = numel(p);
    if m == 0
        thr = 0;
        return;
    end
    k = find(p <= (1:m)'/m*q, 1, 'last');
    if isempty(k)
        thr = 0;
    else
        thr = p(k);
    end
end




function [Sall,Sinc] = buildSummaries(Sall,Sinc,cA,sA,cI,sI)
ms = @(C) cellfun(@(v)[mean(v,'omitnan'),std(v,'omitnan')/sqrt(numel(v))],C,'uni',false);
mA = ms(cA);  sA = ms(sA);
mI = ms(cI);  sI = ms(sI);

Sall.meanCorr   = cellfun(@(x)x(1), mA);
Sall.semCorr    = cellfun(@(x)x(2), mA);
Sall.meanShuff  = cellfun(@(x)x(1), sA);
Sall.semShuff   = cellfun(@(x)x(2), sA);
Sall.allCorrByRat   = cA;
Sall.sigNeuronPct = struct( ...
  'Pearson', Sall.sigPct, ...
  'Shuffle',100*Sall.nSigShuf./max(Sall.nTest,1) );

Sinc.meanCorr   = cellfun(@(x)x(1), mI);
Sinc.semCorr    = cellfun(@(x)x(2), mI);
Sinc.meanShuff  = cellfun(@(x)x(1), sI);
Sinc.semShuff   = cellfun(@(x)x(2), sI);
Sinc.allCorrByRat = cI;
Sinc.sigNeuronPct = struct( ...
  'Pearson', Sinc.sigPct, ...
  'Shuffle',100*Sinc.nSigShuf./max(Sinc.nTest,1) );
end

function printConsole(Sall,Sinc,ratNames,method,pPopAll,pPopInc)
fprintf('\n############################################################\n');
fprintf('###  FR-vs-Speed SUMMARY PER RAT  (MC: %s)\n',upper(method));
fprintf('############################################################\n');
hdr = @(lbl)fprintf(['\n---------------- %-17s----------------\n' ...
  '%-8s %-11s %-11s %-11s %-11s %-11s %-11s\n'],lbl, ...
  'Rat','nIncl/Total','p(A-v-S)','mean r','SD r','%Sig(P)','%Sig(sh)');
hdr('ALL NEURONS');
for k=1:numel(ratNames)
  v = Sall.allCorrByRat{k};
  fprintf('%-8s %4d/%-6d %11.3g %11.3f %11.3f %11.1f %11.1f\n',...
    ratNames{k},Sall.nTest(k),Sall.nTest(k),pPopAll(k), ...
    mean(v,'omitnan'),std(v,'omitnan'), ...
    Sall.sigNeuronPct.Pearson(k),Sall.sigNeuronPct.Shuffle(k));
end
hdr('INCLUDED NEURONS');
for k=1:numel(ratNames)
  v = Sinc.allCorrByRat{k};
  fprintf('%-8s %4d/%-6d %11.3g %11.3f %11.3f %11.1f %11.1f\n',...
    ratNames{k},Sinc.nTest(k),Sinc.nTest(k),pPopInc(k), ...
    mean(v,'omitnan'),std(v,'omitnan'), ...
    Sinc.sigNeuronPct.Pearson(k),Sinc.sigNeuronPct.Shuffle(k));
end
fprintf('############################################################\n\n');
end

function makeFigure(S,figTitle,method,ratNames)
nR = numel(ratNames);
figure('Color','w','Position',[100 300 1600 420],'Name',figTitle);
sgtitle(sprintf('%s  (%s correction)',figTitle,upper(method)));

% 1) grouped bar
subplot(1,4,1); hold on;
bh = bar([S.meanCorr S.meanShuff],'grouped');
errorbar(bh(1).XEndPoints,S.meanCorr,S.semCorr,'k.');
errorbar(bh(2).XEndPoints,S.meanShuff,S.semShuff,'k.');
legend({'Actual','Shuffle'},'Location','northwest');
xticks(1:nR); xticklabels(ratNames); ylabel('Mean r'); title('Mean ± SEM');
fprintf('REMEMBER FOR ALL CELLS THE FIRST GRAPH WILL BE IRRELEVANT BC YOU CANT AVERAGE FOR NON INCLUDED CELLS ANYWAY')


% 2) boxplot of neuron‐wise
subplot(1,4,2);
data = S.allCorrByRat(~cellfun(@isempty,S.allCorrByRat));
gdat=[]; lab=[];
for i=1:numel(data)
  gdat=[gdat;data{i}(:)]; lab=[lab;i*ones(numel(data{i}),1)];
end
if ~isempty(gdat), boxplot(gdat,lab,'Labels',ratNames); ylabel('r'); end
title('Neuron-wise');

% 3) %Sig Pearson
subplot(1,4,3); bar(S.sigNeuronPct.Pearson); ylim([0 100]);
xticks(1:nR); xticklabels(ratNames); title('%Sig (Pearson)');

% 4) %Sig Shuffle
subplot(1,4,4); bar(S.sigNeuronPct.Shuffle); ylim([0 100]);
xticks(1:nR); xticklabels(ratNames); title('%Sig (shuffle)');
end
