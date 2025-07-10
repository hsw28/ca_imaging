function traceVsSpeed_RESIDUAL(method,nPerm)
% -------------------------------------------------------------------------
%   Speed-corrected residual jump  (trace – baseline)  with *matching*
%   15-bin baseline and 15-bin trace windows.
%
%   – For every trial:  30 rows
%         rows 1-15  : baseline  (-2 → 0 s, 133-ms bins)
%         rows 16-30 : trace     (0 → 2 s, 133-ms bins)
%   – Ordinary least-squares:  spk ~ speed   (row-wise)
%   – Residual per row      :  r = spk – ŷ
%   – ΔRtrial = mean(rTrace) – mean(rBase)
%   – Per neuron (≥3 trials):
%         • two-tailed paired-t  on ΔRtrial
%         • sign-flip permutation (nPerm flips)
%   – Per rat multiple-comparison:
%         method = 'fdr'  (default, BH q = 0.05)
%                = 'bonf' (Bonferroni α = 0.05)
%   – Two pools:  ALL neurons  vs  INCLUDED (≥ minSpk spikes in **trace**)
%
%   Figures (per pool):
%         1. % significant (t vs flip) grouped bars
%         2. Box-plots of mean ΔR per rat
%         3. Swarm of all ΔR pooled
%         4. Histogram of adjusted p / q values
%
%   Dependencies:  rat structs,  autoDateList,  ca_velocity,  mafdr
% -------------------------------------------------------------------------

if nargin < 1, method = 'fdr';   end, method = lower(method);
if nargin < 2, nPerm  = 500;     end                       % flips / neuron
alphaFW = 0.05;  qFDR = 0.05;                             % thresholds

% ---------------- USER SETTINGS -----------------------------------------
ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};

preWin   = [-2 0];    trcWin   = [0 2];
binSize  = 1/7.5;     nBins    = round(diff(trcWin)/binSize); % 15
minSpk   = 5;                                           % inclusion rule
% ------------------------------------------------------------------------

nR = numel(ratNames);
statsAll = repmat(emptyStats,nR,1);
statsInc = repmat(emptyStats,nR,1);
allDR = cell(nR,1);   incDR = cell(nR,1);                % for plots

fprintf('\n=====  Residual trace-vs-baseline  (%s, %d flips)  =====\n',...
        upper(method),nPerm);

% ================= RAT LOOP =============================================
for r = 1:nR
  r
    rat = evalin('base',ratNames{r});
    dL  = autoDateList(rat);
    days= dL(find(strcmp(dL,rat.An))-2 : find(strcmp(dL,rat.An)));

    pT = [];      pPerm = [];      nSpk_all = [];  DR_all = [];
    pT_i = [];    pPerm_i = [];    nSpk_inc = [];  DR_inc = [];

    % -------------- DAY LOOP -------------------------------------------
    for d = 1:3
      d
        day   = days{d};
        spkM  = rat.Ca_peaks.(['CA_peaks_' day]);
        ts    = rat.pos.(['pos_' day])(:,1);
        xy    = rat.pos.(['pos_' day])(:,2:3);
        csT   = rat.CS_times.(['CS_' day]);

        % speed trace
        v     = ca_velocity([ts.';xy(:,1).';xy(:,2).']);
        speed = interp1(v(2,:),v(1,:),ts,'linear','extrap'); speed(~isfinite(speed))=0;

        nTr   = numel(csT);

        % Precompute edges ------------------------------------------------
        baseE = arrayfun(@(t) linspace(csT(t)+preWin(1),csT(t)+preWin(2),nBins+1), ...
                         1:nTr,'uni',0);  baseE = vertcat(baseE{:});  % nTr×(nBins+1)
        trcE  = arrayfun(@(t) linspace(csT(t)+trcWin(1),csT(t)+trcWin(2),nBins+1), ...
                         1:nTr,'uni',0);  trcE  = vertcat(trcE{:});

        % ---------------- NEURON LOOP ----------------------------------
        for ni = 1:size(spkM,1)
            spk = spkM(ni,:); spk = spk(~isnan(spk));

            dR = nan(nTr,1);   nSpkTrc = 0;


            for t = 1:nTr %% trial loop
                % baseline rows
                cntB = histcounts(spk,baseE(t,:));
                spdB = arrayfun(@(b) ...
                       mean(speed(ts>=baseE(t,b)&ts<baseE(t,b+1))),1:nBins);

                % trace rows
                cntT = histcounts(spk,trcE(t,:));
                spdT = arrayfun(@(b) ...
                       mean(speed(ts>=trcE(t,b)&ts<trcE(t,b+1))),1:nBins);

                Y = [cntB cntT].';                         % 30 × 1
                X = [ones(30,1) , [spdB spdT].'];          % 30 × 2

                B = X\Y;             resid = Y - X*B;      % residuals

                dR(t) = max(resid(nBins+1:end)) - max(resid(1:nBins));
                nSpkTrc = nSpkTrc + sum(cntT);
            end

            if nnz(~isnan(dR)) < 3, continue, end
            meanDR = mean(dR,'omitnan');

            % parametric paired-t (two-tailed)
            [~,p] = ttest(dR);

            % permutation p (sign-flip)
            flips = 2*(rand(nPerm,numel(dR))>0.5) - 1;     % ±1 mask
            permMeans = mean(flips .* dR',2,'omitnan');
            pSF = mean(abs(permMeans) >= abs(meanDR));

            % ---- store ALL --------------------------------------------
            pT      = [pT ; p];
            pPerm   = [pPerm ; pSF];
            nSpk_all= [nSpk_all ; nSpkTrc];
            DR_all  = [DR_all  ; meanDR];

            % ---- INCLUDED ---------------------------------------------
            if nSpkTrc >= minSpk
                pT_i    = [pT_i ; p];
                pPerm_i = [pPerm_i ; pSF];
                nSpk_inc= [nSpk_inc ; nSpkTrc];
                DR_inc  = [DR_inc  ; meanDR];
            end
        end
    end % day

    allDR{r} = DR_all;  incDR{r} = DR_inc;

    statsAll(r) = packStats(pT ,pPerm ,method,alphaFW,qFDR,nSpk_all);
    statsInc(r) = packStats(pT_i,pPerm_i,method,alphaFW,qFDR,nSpk_inc);
end % rat

% --------------- Figures & Console --------------------------------------

summarize(statsAll,ratNames,'ALL neurons');
summarize(statsInc,ratNames,'INCLUDED neurons');
makeFigures(statsAll,ratNames,'ALL neurons',    allDR);
makeFigures(statsInc,ratNames,'INCLUDED neurons',incDR);
end
% ========================================================================
%                               HELPERS
% ========================================================================
function [sig,adj] = mcCorrect(p,method,a,q)
if isempty(p), sig=false(size(p)); adj=p; return, end
switch method
    case 'bonf',  thr = a/numel(p); sig = p<thr; adj=p;
    otherwise,   adj = mafdr(p,'BHFDR',true); sig = adj<q;
end
end
% ------------------------------------------------------------------------
function S = packStats(pT,pSF,method,a,q,spk)
[sigT ,adjT ] = mcCorrect(pT ,method,a,q);
[sigSF,adjSF] = mcCorrect(pSF,method,a,q);

S.sigT   = sigT;    S.sigSF = sigSF;
S.pctT   = 100*mean(sigT);   S.pctSF = 100*mean(sigSF);
S.nNeurons = numel(sigT);
S.spikes   = spk;
S.adjT   = adjT;    S.adjSF = adjSF;
end
% ------------------------------------------------------------------------
function S = emptyStats
S = struct('sigT',[],'sigSF',[], ...
           'pctT',NaN,'pctSF',NaN, ...
           'nNeurons',0,'spikes',[], ...
           'adjT',[],'adjSF',[]);
end
% ------------------------------------------------------------------------
function makeFigures(st,ratN,ttl,deltaRcell)
nR = numel(ratN);
pctT  = arrayfun(@(s)s.pctT ,st);
pctSF = arrayfun(@(s)s.pctSF,st);

figure('color','w','Name',ttl,'Position',[100 300 1400 520]);
sgtitle(sprintf('%s  (baseline 15 bins)',ttl),'FontWeight','bold');

% 1. % significant grouped bar
subplot(3,1,1); bar([pctT(:) pctSF(:)],'grouped');
ylim([0 100]); ylabel('% significant'); grid on
xticks(1:nR); xticklabels(ratN);
legend({'t-test','sign-flip'},'location','northwest'); title('%Sig');

% 2. Box-plot ΔR per rat
subplot(3,1,2); hold on
g=[];lab=[];
for k=1:nR
  g=[g;deltaRcell{k}(:).*7.5];
  lab=[lab;k*ones(numel(deltaRcell{k}),1)];
end
if ~isempty(g), boxplot(g,lab,'Labels',ratN,'Symbol','k.'); end
yline(0,'k--'); ylabel('mean Δ residual (hz)'); grid on
title('Per-neuron ΔR');

% 3. Swarm pooled ΔR
subplot(3,1,3); hold on
pool = vertcat(deltaRcell{:}).*7.5;
if ~isempty(pool), swarmchart(ones(size(pool)),pool,'filled','MarkerFaceAlpha',0.4); end
yline(0,'k--'); xlim([0.5 1.5]); set(gca,'XTick',[]);
ylabel('mean Δ residual (hz)'); box on; title('All neurons pooled');

end
% ------------------------------------------------------------------------
function summarize(st,ratN,hdr)
fprintf('\n--- %s ---\n',hdr);
for k=1:numel(ratN)
    fprintf('%s : %4d / %4d  (T %4.1f%% | SF %4.1f%%)\n',ratN{k}, ...
        sum(st(k).sigT), st(k).nNeurons, st(k).pctT, st(k).pctSF);
end
pT  = arrayfun(@(s)s.pctT ,st);
pSF = arrayfun(@(s)s.pctSF,st);
fprintf('Grand mean  T = %.1f ± %.1f %%   |   SF = %.1f ± %.1f %%\n\n',...
        mean(pT,'omitnan'), std(pT,'omitnan'), ...
        mean(pSF,'omitnan'),std(pSF,'omitnan'));
end
