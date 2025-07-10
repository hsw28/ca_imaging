function temp(method)
% -------------------------------------------------------------------------
% FR–vs-speed correlations in the CS window (0–2 s) across three days
% (day–2, day–1, day–0) for each rat.
%
% Multiple‐comparison control:
%   method = 'fdr'  (default) → Benjamini–Hochberg, q = 0.05
%          = 'bonf'           → Bonferroni,          α = 0.05
%
% (1) Per‐neuron:
%     • Pearson r → pPear
%     • Shuffle test (nShuff perms) → pShuf
%     • MC‐correction on pPear & pShuf
%
% (2) Population‐level:
%     • pPopAll/Inc: compare |mean(r)| to null of mean(r) from all shuffles
%
% Two pools: ALL neurons vs INCLUDED (≥ minSpk spikes in trace window)
% -------------------------------------------------------------------------
if nargin<1, method='fdr'; end
method = lower(method);
assert(ismember(method,{'fdr','bonf'}),'method must be ''fdr'' or ''bonf''');

% params
qFDR    = 0.05;
alphaFW = 0.05;
minSpk  = 5;
win     = [0 2];
binSize = 1/7.5;
nShuff  = 5;
ratNames= {'rat0222','rat0307','rat0313','rat0314','rat0816'};

fprintf('\n=====  plotFRvsSpeedSummary   (%s correction)  =====\n',upper(method));

nRats    = numel(ratNames);
all      = makeBlankStruct(nRats);
inc      = makeBlankStruct(nRats);
pPopAll  = nan(nRats,1);
pPopInc  = nan(nRats,1);

for r = 1:nRats
    fprintf('Processing %s...\n', ratNames{r});
    rat      = evalin('base',ratNames{r});
    dates    = autoDateList(rat);
    idx      = find(strcmp(dates,rat.An));
    days     = dates(idx-2:idx);

    % per‐rat collectors
    r_all    = [];  p_all    = [];  pShufA = [];  shMat_all = [];
    r_inc    = [];  p_inc    = [];  pShufI = [];  shMat_inc = [];

    for d = 1:3
        day     = days{d};
        spkMat  = rat.Ca_peaks.(['CA_peaks_' day]);
        posRaw  = rat.pos.(['pos_' day]);
        csTimes = rat.CS_times.(['CS_' day]);

        % speed trace
        V       = ca_velocity(posRaw');
        speed   = interp1(V(2,:),V(1,:),posRaw(:,1),'linear','extrap');
        speed(~isfinite(speed)) = 0;
        ts      = posRaw(:,1);

        nTr    = numel(csTimes);
        nBins  = round(diff(win)/binSize);

        % precompute
        binEdges  = nan(nTr,nBins+1);
        speedBins = nan(nTr,nBins);
        for t = 1:nTr
            e = linspace(csTimes(t)+win(1),csTimes(t)+win(2),nBins+1);
            binEdges(t,:) = e;
            for b = 1:nBins
                ix = ts>=e(b) & ts<e(b+1);
                speedBins(t,b) = mean(speed(ix));
            end
        end

        % shuffle template
        shSpeed = nan(nTr,nBins,nShuff);
        for s = 1:nShuff
          for t = 1:nTr
            shSpeed(t,:,s) = speedBins(t,randperm(nBins));
          end
        end

        % per‐neuron
        for ni = 1:size(spkMat,1)
            spk = spkMat(ni,:);
            spk = spk(~isnan(spk));

            rT            = nan(nTr,1);
            shuffleMeans  = zeros(1,nShuff);   % <-- reset per neuron
            nSpkTrc       = 0;

            for t = 1:nTr
                cnt      = histcounts(spk,binEdges(t,:));
                nSpkTrc  = nSpkTrc + sum(cnt);
                ok       = ~isnan(cnt) & ~isnan(speedBins(t,:));
                if nnz(ok)>=3
                    rT(t) = corr(cnt(ok)',speedBins(t,ok)','type','Pearson');
                    for s = 1:nShuff
                        shuffleMeans(s) = shuffleMeans(s) + ...
                            corr(cnt(ok)',shSpeed(t,ok,s)','type','Pearson');
                    end
                end
            end

            rT = rT(~isnan(rT));
            if numel(rT)<3, continue; end

            % per‐cell stats
            rMu             = mean(rT);
            [~,pPear]       = ttest(rT);
            shuffleMeans    = shuffleMeans ./ numel(rT);  % <<-- CORRECTED
            pSh             = mean(abs(shuffleMeans) >= abs(rMu));

            % store ALL
            r_all(end+1)       = rMu;
            p_all(end+1)       = pPear;
            pShufA(end+1)      = pSh;
            shMat_all(end+1,:) = shuffleMeans;            % <<-- CORRECTED

            % store INCLUDED
            if nSpkTrc>=minSpk
                r_inc(end+1)       = rMu;
                p_inc(end+1)       = pPear;
                pShufI(end+1)      = pSh;
                shMat_inc(end+1,:) = shuffleMeans;          % <<-- CORRECTED
            end
        end
    end

    % population‐level p
    if ~isempty(r_all)
        nullAll     = mean(shMat_all,1);
        pPopAll(r)  = mean(abs(nullAll) >= abs(mean(r_all)));
    end
    if ~isempty(r_inc)
        nullInc     = mean(shMat_inc,1);
        pPopInc(r)  = mean(abs(nullInc) >= abs(mean(r_inc)));
    end

    % MC‐correction
    [sigPA,sigSA] = sigTest(p_all, pShufA, method, qFDR, alphaFW);
    [sigPI,sigSI] = sigTest(p_inc, pShufI, method, qFDR, alphaFW);

    all = updateStruct(all, r, r_all, mean(shMat_all,2), sigPA, sigSA);
    inc = updateStruct(inc, r, r_inc, mean(shMat_inc,2), sigPI, sigSI);

    all.nIncl(r)=numel(r_all); all.nTotal(r)=numel(r_all); all.pAvsS(r)=pPopAll(r);
    inc.nIncl(r)=numel(r_inc); inc.nTotal(r)=numel(r_inc); inc.pAvsS(r)=pPopInc(r);
end

% console + figures
printConsole(inc, ratNames, method,'INCLUDED neuron-days');
printConsole(all, ratNames, method,'ALL neuron-days');
makeFigure(inc,'INCLUDED neuron-days',method,ratNames);
makeFigure(all ,'ALL neuron-days',     method,ratNames);
end


%% ——— HELPERS —————————————————————————————————————————————

function S = makeBlankStruct(n)
    S.meanCorr           = nan(n,1);
    S.semCorr            = nan(n,1);
    S.meanShuff          = nan(n,1);
    S.semShuff           = nan(n,1);
    S.sigPercent.Pearson = nan(n,1);
    S.sigPercent.Shuffle = nan(n,1);
    S.allCorrByRat       = cell(n,1);
    S.nIncl              = nan(n,1);
    S.nTotal             = nan(n,1);
    S.pAvsS              = nan(n,1);
end

function [sigP,sigS] = sigTest(pPear,pShuf,method,qFDR,alphaFW)
    switch method
      case 'bonf'
        thrP = alphaFW/numel(pPear);
        thrS = alphaFW/numel(pShuf);
      otherwise
        thrP = fdr_bh(pPear,qFDR);
        thrS = fdr_bh(pShuf,qFDR);
    end
    sigP = pPear < thrP;
    sigS = pShuf < thrS;
end

function thr = fdr_bh(p,q)
    p = sort(p(:)); m = numel(p);
    if m==0, thr=0; return; end
    k = find(p <= (1:m)'/m*q,1,'last');
    if isempty(k), thr=0; else thr=p(k); end
end

function S = updateStruct(S,r,C,Csh,sigP,sigS)
S.meanCorr(r)            = mean(C ,'omitnan');
S.semCorr (r)            = std (C ,'omitnan')/sqrt(max(1,numel(C)));
S.meanShuff(r)           = mean(Csh,'omitnan');
S.semShuff (r)           = std (Csh,'omitnan')/sqrt(max(1,numel(Csh)));
S.sigPercent.Pearson(r)  = 100*mean(sigP);
S.sigPercent.Shuffle(r)  = 100*mean(sigS);
S.allCorrByRat{r}        = C;
end

function printConsole(S,ratNames,method,label)
fprintf('\n############################################################\n');
fprintf('###  FR-vs-Speed SUMMARY PER RAT  (%s)   <%s>\n',...
        upper(method),label);
fprintf('############################################################\n');

hdr = ['\n-------------------  %-15s-------------------\n' ...
       '%-8s %-11s %-11s %-11s %-11s %-11s %-11s\n'];
printfn = @(lbl)fprintf(hdr,lbl,'Rat','nIncl/Total','p(A-v-S)', ...
                              'mean r','SD r','%Sig(Pear)','%Sig(shuf)');
printfn(label);

for k = 1:numel(ratNames)
    nI = S.nIncl(k);  nT = S.nTotal(k);
    r  = S.allCorrByRat{k};
    fprintf('%-8s %4d/%-6d %11.3g %11.3f %11.3f %11.1f %11.1f\n',...
        ratNames{k}, nI, nT, S.pAvsS(k), ...
        mean(r,'omitnan'), std(r,'omitnan'), ...
        S.sigPercent.Pearson(k), S.sigPercent.Shuffle(k));
end
fprintf('############################################################\n\n');
end

function makeFigure(S,figTitle,method,ratNames)
nR = numel(ratNames);
figure('Color','w','Position',[100 300 1600 420],'Name',figTitle);
sgtitle(sprintf('%s  (%s correction)',figTitle,upper(method)));

subplot(1,4,1); hold on;
S.meanShuff
bh = bar([S.meanCorr S.meanShuff],'grouped');
x1=bh(1).XEndPoints;  x2=bh(2).XEndPoints;
errorbar(x1,S.meanCorr ,S.semCorr ,'k.','CapSize',8);
errorbar(x2,S.meanShuff,S.semShuff,'k.','CapSize',8);
xticks(1:nR); xticklabels(ratNames); ylabel('Mean FR–spd r');
legend({'Actual','Shuffle'},'location','northwest'); title('Mean ± SEM r');

subplot(1,4,2);
boxData = S.allCorrByRat(~cellfun(@isempty,S.allCorrByRat));
gdat=[]; glab=[];
for i = 1:numel(boxData)
    gdat = [gdat; boxData{i}(:)];
    glab = [glab; i*ones(numel(boxData{i}),1)];
end
boxplot(gdat,glab,'Labels',ratNames); ylabel('Pearson r'); title('Neuron-wise');

subplot(1,4,3);
bar(S.sigPercent.Pearson); ylim([0 100]);
xticks(1:nR); xticklabels(ratNames);
ylabel('% sig (Pearson)'); title('Sig. neurons');

subplot(1,4,4);
bar(S.sigPercent.Shuffle); ylim([0 100]);
xticks(1:nR); xticklabels(ratNames);
ylabel('% sig (shuffle)'); title('Shuffle sig.');
end
