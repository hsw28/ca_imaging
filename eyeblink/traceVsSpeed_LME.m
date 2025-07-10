function traceVsSpeed_LME(method)
% -------------------------------------------------------------------------
% Mixed‐effects + residual‐jump plots:
%   spikes ~ SpeedZ + Trace + (1|Trial)
%   • 30 rows/trial (15 baseline + 15 trace bins)
%   • fits Poisson GLME, log link
%   • extracts Raw residuals, builds ΔR per neuron
%   • per‐rat FDR('fdr') or Bonferroni('bonf') on pTrace
%   • plots: %sig bar, ΔR box & swarm
% -------------------------------------------------------------------------
if nargin<1, method='fdr'; end
method = lower(method);
alphaFW = 0.05;  qFDR = 0.05;

ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
%ratNames = {'rat0307'};

preWin   = [-2 0];    trcWin = [0 2];
binSize  = 1/7.5;     nBins  = round(diff(trcWin)/binSize);
minSpk   = 5;

nR      = numel(ratNames);
statsAll = repmat(emptyStats,nR,1);
statsInc = statsAll;
deltaR_all = cell(nR,1);
deltaR_inc = cell(nR,1);

fprintf('\n=== LME + residual ΔR  (%s) ===\n',upper(method));

for r = 1:nR
    rat   = evalin('base',ratNames{r});
    dL    = autoDateList(rat);
    idx   = find(strcmp(dL,rat.An));
    days  = dL(idx-2:idx);

    pAll = [];  nAll = [];  DRall = [];
    pInc = [];  nInc = [];  DRinc = [];

    for d = 1:3
        day    = days{d};
        spkM   = rat.Ca_peaks.(['CA_peaks_' day]);
        ts     = rat.pos.(['pos_' day])(:,1);
        xy     = rat.pos.(['pos_' day])(:,2:3);
        csOn   = rat.CS_times.(['CS_' day]);
        nTr    = numel(csOn);

        % compute speed trace
        v    = ca_velocity([ts.'; xy(:,1).'; xy(:,2).']);
        spd  = interp1(v(2,:),v(1,:),ts,'linear','extrap');
        spd(~isfinite(spd)) = 0;

        % precompute edges
        baseE = arrayfun(@(t)linspace(csOn(t)+preWin(1), csOn(t)+preWin(2), nBins+1), ...
                         1:nTr,'uni',0);
        trcE  = arrayfun(@(t)linspace(csOn(t)+trcWin(1),csOn(t)+trcWin(2),nBins+1),...
                         1:nTr,'uni',0);

        for ni = 1:size(spkM,1)
            spk = spkM(ni,:); spk = spk(~isnan(spk));
            rowsPerTrial = 2*nBins;
            Y       = nan(nTr*rowsPerTrial,1);
            SpeedZ  = nan(size(Y));
            Trace   = nan(size(Y));
            Trial   = nan(size(Y));
            nSpkTrc = 0;
            row = 0;

            for t = 1:nTr
                cntB = histcounts(spk,baseE{t});
                cntT = histcounts(spk,trcE{t});
                spdB = arrayfun(@(b)mean(spd(ts>=baseE{t}(b)&ts<baseE{t}(b+1))),1:nBins);
                spdT = arrayfun(@(b)mean(spd(ts>=trcE{t}(b)&ts<trcE{t}(b+1))),1:nBins);
                nSpkTrc = nSpkTrc + sum(cntT);
                zB = zscore(spdB);
                zT = zscore(spdT);

                for b = 1:nBins
                    row = row+1;
                    Y(row)      = cntB(b);
                    SpeedZ(row) = zB(b);
                    Trace(row)  = 0;
                    Trial(row)  = t;
                end
                for b = 1:nBins
                    row = row+1;
                    Y(row)      = cntT(b);
                    SpeedZ(row) = zT(b);
                    Trace(row)  = 1;
                    Trial(row)  = t;
                end
            end

            if all(isnan(Y)), continue; end

            Ttbl = table(Y,SpeedZ,Trace,categorical(Trial),...
                         'VariableNames',{'Spk','SpeedZ','Trace','Trial'});

            % fit generalized LME
            try
                mdl = fitglme(Ttbl, ...
                    'Spk ~ SpeedZ + Trace + (1|Trial)', ...
                    'Distribution','Poisson','Link','log');
            catch
                continue
            end

            % extract p for Trace effect
            ct = mdl.Coefficients;
            pTrace = ct.pValue(strcmp(ct.Name,'Trace'));

            % raw residuals (Y - fitted μ)
            res = residuals(mdl,'ResidualType','Raw');

            % compute ΔR per trial
            dR = nan(nTr,1);
            for t = 1:nTr
                idxs = (t-1)*rowsPerTrial + (1:rowsPerTrial);
                baseRes = res(idxs(1:nBins));
                trcRes  = res(idxs(nBins+1:end));
                dR(t)   = mean(trcRes) - mean(baseRes);
            end
            meanDR = mean(dR,'omitnan');

            % collect ALL
            pAll   = [pAll;   pTrace];
            nAll   = [nAll;   nSpkTrc];
            DRall  = [DRall;  meanDR];
            % INCLUDED
            if nSpkTrc >= minSpk
                pInc  = [pInc;  pTrace];
                nInc  = [nInc;  nSpkTrc];
                DRinc = [DRinc; meanDR];
            end
        end
    end

    deltaR_all{r} = DRall;
    deltaR_inc{r} = DRinc;
    statsAll(r)   = finishStatsLME(pAll,method,alphaFW,qFDR,nAll);
    statsInc(r)   = finishStatsLME(pInc,method,alphaFW,qFDR,nInc);
end

% print & plot
summarizeLME(statsAll,ratNames,'ALL neurons (LME)');
summarizeLME(statsInc,ratNames,'INCLUDED neurons (LME)');
plotResiduals(deltaR_all, ratNames, 'ΔR (ALL neurons, LME)');
plotResiduals(deltaR_inc, ratNames, 'ΔR (INCLUDED neurons, LME)');

% print ΔR residuals
fprintf('\n--- ΔR residuals (spikes/bin) ---\n');
for r = 1:nR
    drA = deltaR_all{r};
    drI = deltaR_inc{r};
    fprintf('%s ALL: mean=%.3f  SD=%.3f  N=%d\n', ...
      ratNames{r}, mean(drA,'omitnan'), std(drA,'omitnan'), numel(drA));
    fprintf('%s INC: mean=%.3f  SD=%.3f  N=%d\n', ...
      ratNames{r}, mean(drI,'omitnan'), std(drI,'omitnan'), numel(drI));
end
end

%% helpers

function [sig,adj] = mcCorrect(p,method,alph,q)
    if isempty(p), sig=false(size(p)); adj=p; return; end
    switch method
      case 'bonf'
        thr = alph/numel(p); sig = p<thr; adj = p;
      otherwise
        adj = mafdr(p,'BHFDR',true); sig = adj<q;
    end
end

function S = finishStatsLME(p,method,alph,q,spk)
    [sig,adj] = mcCorrect(p,method,alph,q);
    S = struct('sigTrace',sig, ...
               'pctSigTrace',100*mean(sig), ...
               'nNeurons',numel(p), ...
               'spikes',spk, ...
               'adjP',adj);
end

function summarizeLME(st,ratN,header)
    fprintf('\n--- %s ---\n',header);
    for k=1:numel(ratN)
        s = st(k);
        fprintf('%s : %4d / %4d  (%.1f%%)\n',...
          ratN{k}, sum(s.sigTrace), s.nNeurons, s.pctSigTrace);
    end
end

function plotResiduals(deltaRcell, ratN, ttl)
    nR = numel(ratN);
    figure('Name',ttl,'Color','w','Position',[200 300 600 400]);
    sgtitle(ttl,'FontWeight','bold'); hold on
    allDR=[]; grp=[];
    for k=1:nR
      allDR = [allDR; deltaRcell{k}];
      grp   = [grp;   k*ones(numel(deltaRcell{k}),1)];
    end
    if ~isempty(allDR)
      boxplot(allDR,grp,'Labels',ratN,'Symbol','k.');
      swarmchart(grp,allDR,8,'filled','MarkerFaceAlpha',0.3);
      yline(0,'k--');
      ylabel('mean Δ residual'); xlabel('Rat');
    end
end

function S = emptyStats
    S = struct('sigTrace',[],'pctSigTrace',NaN,'nNeurons',0,'spikes',[],'adjP',[]);
end
