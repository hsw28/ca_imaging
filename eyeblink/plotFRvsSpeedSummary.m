function plotFRvsSpeedSummary
% -------------------------------------------------------------------------
% Same functionality as before, but every decision about “significant” now
% uses a Bonferroni‐corrected threshold: α / (# neurons tested in that rat).
% -------------------------------------------------------------------------

ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
nRats    = numel(ratNames);

win    = [0 2];       % CS window
nShuff = 500;         % shuffle repetitions
alpha  = 0.05;        % family-wise error rate

% collectors --------------------------------------------------------------
meanCorrs  = nan(nRats,1);   semCorrs  = nan(nRats,1);
meanShuff  = nan(nRats,1);   semShuff  = nan(nRats,1);
pvals      = nan(nRats,1);
sigNeuronPct        = nan(nRats,1);   % Bonferroni-based
sigNeuronPct_shuff  = nan(nRats,1);   % Bonferroni-based
allCorrByRat        = cell(nRats,1);

for r = 1:nRats
    rat = evalin('base',ratNames{r});           % assumes data in workspace
    dateList   = autoDateList(rat);
    idx        = find(strcmp(dateList,rat.An));
    theseDays  = dateList(idx-2:idx);

    allCorrs  = [];
    allShuffs = [];
    sigReal   = 0;
    sigShuff  = 0;
    nNeurons  = 0;                              % counts *tested* neurons

    % -------- loop over the 3 days --------------------------------------
    for d = 1:3
        dateStr  = theseDays{d};
        spikeMat = rat.Ca_peaks.(['CA_peaks_' dateStr]);
        ts       = rat.pos.(['pos_' dateStr])(:,1);
        pos      = rat.pos.(['pos_' dateStr])(:,2:3);
        csTimes  = rat.CS_times.(['CS_' dateStr]);

        dt    = diff(ts);
        dx    = diff(pos);
        speed = [0; hypot(dx(:,1),dx(:,2))./dt];
        speed(~isfinite(speed)) = 0;

        for ni = 1:size(spikeMat,1)
            spikes = spikeMat(ni,:);  spikes = spikes(~isnan(spikes));
            if numel(spikes)<5, continue; end

            % --- firing-rate & speed in CS window -----------------------
            fr  = nan(size(csTimes));
            spd = nan(size(csTimes));
            for t = 1:numel(csTimes)
                t1 = csTimes(t)+win(1);  t2 = csTimes(t)+win(2);
                fr(t)  = sum(spikes>=t1 & spikes<=t2) / diff(win);
                mask   = ts>=t1 & ts<=t2;
                spd(t) = mean(speed(mask));
            end

            valid = ~isnan(fr)&~isnan(spd);
            if sum(valid)<5, continue; end

            nNeurons = nNeurons + 1;   % one more test
            [rv,pReal] = corr(fr(valid)', spd(valid)','type','Pearson');

            allCorrs(end+1) = rv;

            % ---------- Bonferroni threshold for *this* rat -------------
            % (we don't know nNeurons yet, so we handle below once loop
            % is finished)

            % ---------- shuffle for that neuron -------------------------
            shCorr = nan(nShuff,1);
            for s = 1:nShuff
                spdSh = spd(randperm(numel(spd)));
                good  = ~isnan(fr)&~isnan(spdSh);
                if sum(good)>=5
                    shCorr(s) = corr(fr(good)',spdSh(good)','type','Pearson');
                end
            end
            allShuffs = [allShuffs; shCorr];

            % store p-value to evaluate later with Bonferroni
            corrPsTmp{nNeurons}     = pReal;            %#ok<*AGROW>
            shuffNullTmp{nNeurons}  = shCorr(~isnan(shCorr));
        end
    end % day loop

    % ---- BONFERRONI: corrected α for this rat --------------------------
    alphaBonf = alpha / nNeurons;

    % ---- go back through saved p-values to count significant -----------
    for k = 1:nNeurons
        if corrPsTmp{k} < alphaBonf
            sigReal = sigReal + 1;
        end
        % shuffled null (two-tailed)
        rv = allCorrs(k);
        null = shuffNullTmp{k};
        if ~isempty(null)
            pSh = mean(abs(null-mean(null)) >= abs(rv-mean(null)));
            if pSh < alphaBonf
                sigShuff = sigShuff + 1;
            end
        end
    end

    % ---------- grand stats for this rat --------------------------------
    meanCorrs(r) = mean(allCorrs,'omitnan');
    semCorrs(r)  = std(allCorrs,'omitnan')/sqrt(nNeurons);
    meanShuff(r) = mean(allShuffs,'omitnan');
    semShuff(r)  = std(allShuffs,'omitnan')/sqrt(sum(~isnan(allShuffs)));
    pvals(r)     = mean(abs(allShuffs) >= abs(meanCorrs(r)));  % global test

    sigNeuronPct(r)        = 100*sigReal  / nNeurons;
    sigNeuronPct_shuff(r)  = 100*sigShuff / nNeurons;
    allCorrByRat{r}        = allCorrs;
end % rat loop
% -------------------------------------------------------------------------
% -------------------------  PLOTS  ---------------------------------------
figure('Color','w','Position',[200 300 1600 400]);

% 1. Mean correlations ----------------------------------------------------
subplot(1,4,1); hold on;
barData = [meanCorrs, meanShuff];
semData = [semCorrs , semShuff ];
bh = bar(barData,'grouped');   % two colours
x1 = bh(1).XEndPoints;  x2 = bh(2).XEndPoints;
errorbar(x1,meanCorrs ,semCorrs ,'k','linestyle','none','capsize',8,'linewidth',1.2)
errorbar(x2,meanShuff ,semShuff ,'k','linestyle','none','capsize',8,'linewidth',1.2)
xticks(1:nRats); xticklabels(ratNames);
ylabel('Mean FR–speed r');
legend({'Actual','Shuffle'},'location','northwest');
title('Mean correlation (Bonferroni)');

% stars (Bonferroni)
for r = 1:nRats
    yMax = max(barData(r,:))+max(semData(r,:));
    if     pvals(r) < 0.001/ nRats          , txt = '***';
    elseif pvals(r) < 0.01 / nRats          , txt =  '**';
    elseif pvals(r) < 0.05 / nRats          , txt =   '*';
    else,  txt = '';  end
    if ~isempty(txt)
        text(mean([x1(r),x2(r)]), yMax*1.05, txt,...
             'FontSize',14,'HorizontalAlignment','center')
    end
end

% 2. Boxplot of neuron r --------------------------------------------------
subplot(1,4,2); cla; hold on;
boxData      = allCorrByRat(~cellfun(@isempty,allCorrByRat));
groupedData  = []; groupLabels = [];
for i = 1:numel(boxData)
    groupedData  = [groupedData ; boxData{i}(:)];
    groupLabels  = [groupLabels ; i*ones(numel(boxData{i}),1)];
end
boxplot(groupedData, groupLabels, 'Labels',ratNames);
ylabel('Pearson r (neuron)');
title('Neuron-wise correlation');

% 3. % significant neurons (raw p) ---------------------------------------
subplot(1,4,3); cla; hold on;
bar(sigNeuronPct);
xticks(1:nRats); xticklabels(ratNames);
ylabel('% significant neurons');
title(sprintf('Pos r, p<%.3g (Bonferroni)',alpha));

% 4. % significant neurons (shuffle) -------------------------------------
subplot(1,4,4); cla; hold on;
bar(sigNeuronPct_shuff);
xticks(1:nRats); xticklabels(ratNames);
ylabel('% neurons sig. vs shuffle');
title(sprintf('Shuffle test, p<%.3g (Bonf.)',alpha));

% ------------------ console summary -------------------------------------
fprintf('\n===========  SUMMARY (Bonferroni) ===========\n');
for r = 1:nRats
    fprintf('%s  –  %.1f %% of neurons significant (α=%.4g)\n',...
            ratNames{r}, sigNeuronPct(r), alpha/numel(allCorrByRat{r}));
end
fprintf('Grand mean (%% sig.): %.1f ± %.1f\n',...
        mean(sigNeuronPct,'omitnan'), std(sigNeuronPct,'omitnan'));
end
