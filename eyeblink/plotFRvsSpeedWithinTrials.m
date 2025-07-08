function plotFRvsSpeedWithinTrials
% -------------------------------------------------------------------------
% Within–trial FR-vs-speed correlations with Bonferroni control.
% One Bonferroni correction per rat (α = .05 / #tested-neurons).
% -------------------------------------------------------------------------

ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
nRats    = numel(ratNames);

win      = [0 2];             % CS window (s)
binSize  = 1/7.5;             % sec
nBins    = round(diff(win)/binSize);
nShuff   = 500;               % shuffles / neuron
alphaFW  = 0.05;              % desired family-wise error rate

% collectors --------------------------------------------------------------
meanCorrs = nan(nRats,1);  semCorrs = nan(nRats,1);
meanShuff = nan(nRats,1);  semShuff = nan(nRats,1);
pvals     = nan(nRats,1);                     % session-level (no Bonf.)
sigPctPearson   = nan(nRats,1);               % Bonf. adj.
sigPctShuffle   = nan(nRats,1);               % Bonf. adj.
allCorrByRat    = cell(nRats,1);

for r = 1:nRats
    rat = evalin('base', ratNames{r});
    dateList = autoDateList(rat);
    idx      = find(strcmp(dateList,rat.An));
    theseDays = dateList(idx-2:idx);          % 3 days

    corrEach   = [];
    shuffAll   = [];
    pPearson   = [];      % raw p from t-test of per-trial r’s
    pShuffle   = [];

    for d = 1:3
        dateStr   = theseDays{d};
        spikeMat  = rat.Ca_peaks.(['CA_peaks_' dateStr]);
        ts        = rat.pos.(['pos_' dateStr])(:,1);
        pos       = rat.pos.(['pos_' dateStr])(:,2:3);
        csTimes   = rat.CS_times.(['CS_' dateStr]);

        % instantaneous speed --------------------------------------------
        dt    = diff(ts);
        dx    = diff(pos);
        speed = [0; hypot(dx(:,1),dx(:,2))./dt];
        speed(~isfinite(speed)) = 0;

        % pre-compute speed bins for every trial --------------------------
        nTrials = numel(csTimes);
        bSpd    = nan(nTrials,nBins);             % speed per bin
        binEdges = nan(nTrials,nBins+1);
        for t = 1:nTrials
            edges = linspace(csTimes(t)+win(1), csTimes(t)+win(2), nBins+1);
            binEdges(t,:) = edges;
            for b = 1:nBins
                idx              = ts>=edges(b) & ts<edges(b+1);
                bSpd(t,b)        = mean(speed(idx));
            end

        end

        % single shared shuffle table to speed things up ------------------
        shuffSpd = nan(nTrials,nBins,nShuff);
        for s = 1:nShuff
            for t = 1:nTrials
                shuffSpd(t,:,s) = bSpd(t,randperm(nBins));
            end
        end

        % ------------------ loop over neurons ----------------------------
        for ni = 1:size(spikeMat,1)
            spikes = spikeMat(ni,:);  spikes = spikes(~isnan(spikes));
            if numel(spikes)<5, continue; end

            perTrialR = nan(1,nTrials);
            for t = 1:nTrials
                counts = histcounts(spikes, binEdges(t,:));
                good   = ~isnan(counts) & ~isnan(bSpd(t,:));
                if sum(good)>=3
                    perTrialR(t) = corr(counts(good)', bSpd(t,good)','type','Pearson');
                end
            end
            perTrialR = perTrialR(~isnan(perTrialR));
            if numel(perTrialR)<3, continue; end

            rMean = mean(perTrialR);
            corrEach(end+1) = rMean;             % save

            % raw p from t-test of the trial correlations
            [~,pTmp] = ttest(perTrialR);
            pPearson(end+1) = pTmp;

            % shuffle ----------------------------------------------------
            shCorr = nan(nShuff,1);
            for s = 1:nShuff
                rSh = nan(1,nTrials);
                for t = 1:nTrials
                    counts = histcounts(spikes, binEdges(t,:));
                    ss    = shuffSpd(t,:,s);
                    good  = ~isnan(counts)&~isnan(ss);
                    if sum(good)>=3
                        rSh(t) = corr(counts(good)', ss(good)','type','Pearson');
                    end
                end
                shCorr(s) = mean(rSh(~isnan(rSh)));
            end
            shuffAll = [shuffAll ; shCorr];

            % shuffle p for this neuron (two-tailed)
            null = shCorr(~isnan(shCorr));
            pShuffle(end+1) = mean(abs(null-mean(null)) >= abs(rMean-mean(null)));
        end
    end % day loop

    % --------------- Bonferroni threshold for this rat -------------------
    nTests     = numel(pPearson);                  % neurons actually tested
    alphaBonf  = alphaFW / nTests;

    sigPctPearson(r) = 100*mean(pPearson  < alphaBonf);
    sigPctShuffle(r) = 100*mean(pShuffle  < alphaBonf);

    meanCorrs(r) = mean(corrEach,'omitnan');
    semCorrs(r)  = std(corrEach,'omitnan')/sqrt(numel(corrEach));
    meanShuff(r) = mean(shuffAll,'omitnan');
    semShuff(r)  = std(shuffAll,'omitnan')/sqrt(sum(~isnan(shuffAll)));

    % session-level comparison (not Bonferroni: one value per rat)
    pvals(r) = mean(mean(shuffAll) >= meanCorrs(r));

    allCorrByRat{r} = corrEach;
end % rat loop
% -------------------------------------------------------------------------
% ------------------------------ PLOTS ------------------------------------
figure('Color','w','Position',[200 300 1600 400]);
% 1. barplot --------------------------------------------------------------
subplot(1,4,1); hold on;
barD = [meanCorrs, meanShuff];  semD = [semCorrs, semShuff];
bh   = bar(barD,'grouped');
errorbar(bh(1).XEndPoints,meanCorrs,semCorrs,'k','linestyle','none','capsize',8)
errorbar(bh(2).XEndPoints,meanShuff,semShuff,'k','linestyle','none','capsize',8)
xticks(1:nRats); xticklabels(ratNames);
ylabel('Mean FR–speed r');
title('Mean ± SEM (Bonf.)');
legend({'Actual','Shuffle'});

for r = 1:nRats
    yMax = max(barD(r,:))+max(semD(r,:));
    if pvals(r)<0.001/nRats, txt='***';
    elseif pvals(r)<0.01/nRats, txt='**';
    elseif pvals(r)<0.05/nRats, txt='*';
    else, txt=''; end
    if ~isempty(txt)
        text(mean(bh(1).XEndPoints([r,r])), yMax*1.05, txt,...
             'horiz','center','fontsize',14)
    end
end

% 2. neuron r boxplot -----------------------------------------------------
subplot(1,4,2); hold on;
data = allCorrByRat(~cellfun(@isempty,allCorrByRat));
gdat = []; glab = [];
for i = 1:numel(data)
    gdat = [gdat; data{i}(:)];
    glab = [glab; i*ones(numel(data{i}),1)];
end
boxplot(gdat,glab,'labels',ratNames);
ylabel('Neuron Pearson r');
title('Within-trial correlations');

% 3. % sig Pearson --------------------------------------------------------
subplot(1,4,3);
bar(sigPctPearson); ylim([0 100]);
xticks(1:nRats); xticklabels(ratNames);
ylabel('% neurons sig. (Pearson)');
title(sprintf('p<%.3g Bonf.',alphaFW));

% 4. % sig shuffle --------------------------------------------------------
subplot(1,4,4);
bar(sigPctShuffle); ylim([0 100]);
xticks(1:nRats); xticklabels(ratNames);
ylabel('% neurons sig. (shuffle)');
title(sprintf('p<%.3g Bonf.',alphaFW));

% -------------- console summary -----------------------------------------
fprintf('\n==== Within-trial Bonferroni summary ====\n');
for r=1:nRats
    fprintf('%s  –  %.1f%% sig (Pearson)\n',ratNames{r},sigPctPearson(r));
end
fprintf('Grand mean (Pearson): %.1f ± %.1f %%\n',...
        mean(sigPctPearson,'omitnan'), std(sigPctPearson,'omitnan'));
end
