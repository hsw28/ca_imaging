function [shuffMeans, shuffSEMs] = plotRateMapSimilaritySummary(nShuffle)
% ------------------------------------------------------------------------
% Task- vs non-task rate-map similarity for five rats.
% Computes BOTH
%   (1) population-mean test (rat-level)  and
%   (2) per-neuron significance (fraction sig. neurons).
% ------------------------------------------------------------------------
if nargin < 1, nShuffle = 1000; end      % >1 needed for population null

ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
metrics  = {'Pearson','Spearman','Cosine'};
nRats    = numel(ratNames);
nMetrics = numel(metrics);

% pre-allocate
actualMeans   = nan(nRats,nMetrics);
shuffMeans    = nan(nRats,nMetrics);
actualSEMs    = nan(nRats,nMetrics);
shuffSEMs     = nan(nRats,nMetrics);
pRat          = nan(nRats,nMetrics);     % rat-level p
fracSig       = nan(nRats,nMetrics);     % % neurons p<0.05

for r = 1:nRats
    rat   = evalin('base',ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates,rat.An));
    these = dates(idx-2:idx);

    actualCorrs  = [];                 % neuron × metric
    shuffCorrs   = [];                 % neuron × metric (mean over shuffles)
    pNeuron      = [];                 % neuron × metric (per-neuron p)

    % accumulate shuffle sums for rat-level null
    sumRep   = zeros(nShuffle,nMetrics);
    countRep = zeros(nShuffle,nMetrics);

    for d = 1:3
        pos  = rat.pos.(sprintf('pos_%s',these{d}));
        ts   = pos(:,1);
        cs   = rat.CS_times.(sprintf('CS_%s',these{d}));
        mask = rat.ratemask.(sprintf('ratemask_%s',these{d}));
        spkMat = rat.Ca_peaks.(sprintf('CA_peaks_%s',these{d}));

        % trial / non-trial time masks
        taskMask = false(size(ts));
        for i = 1:numel(cs)
            taskMask = taskMask | (ts>=cs(i) & ts<=cs(i)+2);
        end
        taskTS    = ts(taskMask);
        nonTaskTS = ts(~taskMask);

        for ni = find(mask(:)==1)'          % iterate kept neurons
            spikes = spkMat(ni,:);  spikes = spikes(~isnan(spikes));

            % split spikes
            taskSpk    = spikes(ismembertol(spikes,taskTS,1e-3));
            nonTaskSpk = setdiff(spikes,taskSpk,'stable');

            R1 = CA_normalizePosData(taskSpk(:),    pos, 2.5, 1);
            R2 = CA_normalizePosData(nonTaskSpk(:), pos, 2.5, 1);
            R1 = imgaussfilt(R1,0.75); R2 = imgaussfilt(R2,0.75);
            R1(R1==0)=NaN; R2(R2==0)=NaN;
            v1 = R1(:); v2 = R2(:);
            good = ~isnan(v1)&~isnan(v2);
            if nnz(good)<10, continue; end
            v1=v1(good); v2=v2(good);

            act = [corr(v1,v2,'type','Pearson'), ...
                   corr(v1,v2,'type','Spearman'), ...
                   dot(v1,v2)/(norm(v1)*norm(v2))];
            actualCorrs(end+1,:) = act;

            % ---------- shuffles for this neuron ----------
            shThis = nan(nShuffle,nMetrics);
            allSpk = [taskSpk(:); nonTaskSpk(:)];
            lbl    = [ones(numel(taskSpk),1); zeros(numel(nonTaskSpk),1)];

            parfor s = 1:nShuffle
                perm = lbl(randperm(numel(lbl)));
                s1   = allSpk(perm==1);     s2 = allSpk(perm==0);
                if numel(s1)<3 || numel(s2)<3, continue; end
                R1s = CA_normalizePosData(s1,pos,2.5,1);
                R2s = CA_normalizePosData(s2,pos,2.5,1);
                R1s = imgaussfilt(R1s,0.75); R2s = imgaussfilt(R2s,0.75);
                R1s(R1s==0)=NaN; R2s(R2s==0)=NaN;
                vs1 = R1s(:); vs2 = R2s(:);
                goodS = ~isnan(vs1)&~isnan(vs2);
                if nnz(goodS)<10, continue; end
                vs1=vs1(goodS); vs2=vs2(goodS);
                shThis(s,:) = [corr(vs1,vs2,'type','Pearson'), ...
                               corr(vs1,vs2,'type','Spearman'), ...
                               dot(vs1,vs2)/(norm(vs1)*norm(vs2))];
            end
            shuffCorrs(end+1,:) = nanmean(shThis,1);

            % per-neuron p
            pNeuron(end+1,:) = 1-mean(shThis >= act, 1,'omitnan');

            % accumulate for rat-level shuffle mean
            valid = ~isnan(shThis);
            sumRep(valid)   = sumRep(valid)   + shThis(valid);
            countRep(valid) = countRep(valid) + 1;
        end
    end

    if isempty(actualCorrs), continue; end
    actualMeans(r,:) = nanmean(actualCorrs,1);
    shuffMeans(r,:)  = nanmean(shuffCorrs,1);
    actualSEMs(r,:)  = nanstd(actualCorrs,0,1)/sqrt(size(actualCorrs,1));
    shuffSEMs(r,:)   = nanstd(shuffCorrs,0,1)/sqrt(size(shuffCorrs,1));

    % -------- (A) Rat-level p-values --------

    ratShuffMean = sumRep ./ countRep;                 % nShuffle × nMetrics
    for m = 1:nMetrics
      fprintf('shuff')
      sort(ratShuffMean(:,m))
      fprintf('actual')
      actualMeans(r,m)
        pRat(r,m) = 1-mean(ratShuffMean(:,m) >= actualMeans(r,m), 'omitnan')
    end

    % -------- (B) Fraction of significant neurons --------
    fracSig(r,:) = mean(pNeuron < 0.05, 1,'omitnan');  % % sig at α=0.05
end

%% =====================  FIGURE 1 : Rat-level means ======================
figure('Color','w','Position',[200 300 1200 400],'Name','Population means');
for m = 1:nMetrics
    subplot(1,3,m); cla; hold on;
    b = bar([actualMeans(:,m), shuffMeans(:,m)], 'grouped');
    b(1).FaceColor = [0.3 0.5 0.9]; b(2).FaceColor = [.7 .7 .7];
    e1 = errorbar(b(1).XEndPoints, actualMeans(:,m), actualSEMs(:,m),'k.');
    e2 = errorbar(b(2).XEndPoints, shuffMeans(:,m),  shuffSEMs(:,m),'k.');
    xticks(1:nRats); xticklabels(ratNames); ylabel('Similarity');
    title([metrics{m} ' (task vs non)']);
    legend({'Actual','Shuffle'},'location','southwest'); ylim([0 1]);
    % stars
    for r = 1:nRats
        if pRat(r,m)<0.001, star='***';
        elseif pRat(r,m)<0.01, star='**';
        elseif pRat(r,m)<0.05, star='*';
        else, star='';
        end
        if ~isempty(star)
            text(r, max(b.YData(r,:))*1.05, star,'HorizontalAlignment','center');
        end
    end
end

%% ===== FIGURE 2 : Fraction of significant neurons (α = 0.05) ============
figure('Color','w','Position',[250 320 1200 350],'Name','Fraction sig. neurons');
for m = 1:nMetrics
    subplot(1,3,m); bar(fracSig(:,m)*100,'FaceColor',[0.4 0.7 0.4]);
    ylim([0 100]); ylabel('% neurons p<0.05');
    xticks(1:nRats); xticklabels(ratNames);
    title([metrics{m} ' : sig-neuron fraction']);
    yline(5,'--k','\alpha=0.05','LabelHorizontalAlignment','left');
end

fprintf('\nRat-level p-values (rows: rats, cols: metrics):\n');
disp(pRat);

fprintf('Fraction sig. neurons (%%, rows: rats, cols: metrics):\n');
disp(fracSig*100);
end
