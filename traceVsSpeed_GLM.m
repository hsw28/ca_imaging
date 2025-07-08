function traceVsSpeed_GLM
% -------------------------------------------------------------------------
% GLM analysis:   Does firing in the CS-trace window persist after
% controlling for speed?
%
%  • Speed is ALWAYS computed via your supplied helper `ca_velocity.m`
%    (which internally calls `smoothpos`).
%  • Neurons with < minSpk spikes IN THE CS WINDOW are excluded from the
%    “included” analysis (but still counted for the “all neurons” stats).
%  • Per-rat Benjamini–Hochberg FDR (q=.05) on Trace coefficients.
%  • Two figures + console summaries:   ALL cells vs INCLUDED cells.
%
% -------------------------------------------------------------------------

% ---------------- USER PARAMETERS ----------------------------------------
ratNames  = {'rat0222','rat0307','rat0313','rat0314','rat0816'};

win       = [0 2];           % CS-trace window (s)
binSize   = 1/7.5;           % bin size (s)
nBins     = round(diff(win)/binSize);
minSpk    = 5;               % ≥ spikes in CS window to be “included”
alphaFW   = 0.05;            % desired FDR level (per rat)
% -------------------------------------------------------------------------

nRats   = numel(ratNames);
fmtHdr  = '\n=== %s =================================\n';

statsAll = repmat(emptyStats,nRats,1);
statsInc = repmat(emptyStats,nRats,1);

for r = 1:nRats
    fprintf(fmtHdr,ratNames{r});
    rat     = evalin('base',ratNames{r});      % struct already in workspace
    dates   = autoDateList(rat);
    idx     = find(strcmp(dates,rat.An));
    days    = dates(idx-2:idx);                % last 3 training days

    % Collectors for this rat ---------------------------------------------
    betaT_all = [];  betaS_all = [];  pT_all = [];   nSpk_all = [];
    betaT_inc = [];  betaS_inc = [];  pT_inc = [];   nSpk_inc = [];

    % ===================== DAY LOOP ======================================
    for d = 1:3
        dayStr = days{d};

        spikeMat = rat.Ca_peaks.(['CA_peaks_' dayStr]);
        ts       = rat.pos.(['pos_' dayStr])(:,1);
        posXY    = rat.pos.(['pos_' dayStr])(:,2:3);
        csTimes  = rat.CS_times.(['CS_' dayStr]);

        % ---------- velocity via helper (returns 2×N: [vel ; time]) -------
        velDat   = ca_velocity([ ts.' ; posXY(:,1).' ; posXY(:,2).' ]);
        speed    = interp1(velDat(2,:), velDat(1,:), ts,'linear','extrap');
        speed(~isfinite(speed)) = 0;

        % ---------- speed per bin per trial -------------------------------
        nTrials  = numel(csTimes);
        speedBins = nan(nTrials,nBins);
        binEdgs   = nan(nTrials,nBins+1);

        for t = 1:nTrials
            edg                = linspace(csTimes(t)+win(1), csTimes(t)+win(2), nBins+1);
            binEdgs(t,:)       = edg;
            for b = 1:nBins
                inBin                  = ts>=edg(b) & ts<edg(b+1);
                speedBins(t,b)         = mean(speed(inBin));
            end
        end

        % =================== NEURON LOOP =================================
        for ni = 1:size(spikeMat,1)
            spk = spikeMat(ni,:);  spk = spk(~isnan(spk));

            betaS_tr = nan(nTrials,1);
            betaT_tr = nan(nTrials,1);
            pT_tr    = nan(nTrials,1);

            % loop trials --------------------------------------------------
            for t = 1:nTrials
                counts = histcounts(spk, binEdgs(t,:));          % 1×nBins
                good   = ~isnan(counts) & ~isnan(speedBins(t,:));

                if nnz(good) >= 3
                    X = [ zscore(speedBins(t,good)).'  ones(nnz(good),1) ]; % [speedZ TraceFlag]
                    Y = counts(good).';

                    try
                      [B,~,stats] = glmfit(X,Y,'poisson','constant','off',...
                                              'options',statset('glmfit','MaxIter',300));

                      if stats.iterations >= 200
                          % mark as non-converged & skip collecting
                          continue     % simply ignore this neuron
                      end
                                              betaS_tr(t) = B(1);
                        betaT_tr(t) = B(2);

                        % Wald p-value for Trace coef
                        [~,pTmp]   = glmval(B,X,'log',struct('Constant','off'));
                        pT_tr(t)   = pTmp;
                    catch
                        % ill-conditioned → leave NaN
                    end
                end
            end % trial

            % ---------- aggregate across trials --------------------------
            betaS_all = [betaS_all ; mean(betaS_tr,'omitnan')];
            betaT_all = [betaT_all ; mean(betaT_tr,'omitnan')];
            pT_all    = [pT_all    ; mean(pT_tr ,'omitnan')];

            % spike count in CS window
            nSpkWin = 0;
            for t = 1:nTrials
                nSpkWin = nSpkWin + sum( spk >= binEdgs(t,1) & spk < binEdgs(t,end) );
            end
            nSpk_all = [nSpk_all ; nSpkWin];

            % ------ inclusion list ---------------------------------------
            if nSpkWin >= minSpk
                betaS_inc = [betaS_inc ; mean(betaS_tr,'omitnan')];
                betaT_inc = [betaT_inc ; mean(betaT_tr,'omitnan')];
                pT_inc    = [pT_inc    ; mean(pT_tr ,'omitnan')];
                nSpk_inc  = [nSpk_inc ; nSpkWin];
            end
        end % neuron
    end % day

    % ---------------- FDR (per rat) --------------------------------------
    statsAll(r) = finishStats(betaT_all,betaS_all,pT_all,alphaFW,nSpk_all);
    statsInc(r) = finishStats(betaT_inc,betaS_inc,pT_inc,alphaFW,nSpk_inc);
end % rat

% ---------------- Figures & summaries ------------------------------------
makeFigures(statsAll,ratNames,'ALL neurons');
makeFigures(statsInc,ratNames,sprintf('INCLUDED (≥%d spikes)',minSpk));

summarize(statsAll,ratNames,'ALL neurons');
summarize(statsInc,ratNames,sprintf('INCLUDED (≥%d spikes)',minSpk));
end
% ===================== Helper sub-functions ==============================
function S = emptyStats
S.betaTrace   = [];  S.betaSpeed = [];
S.qTrace      = [];  S.sigTrace  = [];
S.pctSigTrace = NaN;
S.nNeurons    = 0;
S.spikes      = [];
end
% -------------------------------------------------------------------------
function S = finishStats(bT,bS,pT,alphaFW,spkCnt)
if isempty(pT)
    S = emptyStats;  return
end
qvals    = mafdr(pT,'BHFDR',true);
sigMask  = qvals < alphaFW;

S.betaTrace   = bT;
S.betaSpeed   = bS;
S.qTrace      = qvals;
S.sigTrace    = sigMask;
S.pctSigTrace = 100*mean(sigMask);
S.nNeurons    = numel(pT);
S.spikes      = spkCnt;
end
% -------------------------------------------------------------------------
function makeFigures(stats,ratNames,figTitle)
pct = arrayfun(@(s)s.pctSigTrace,stats);

figure('color','w','name',figTitle,'pos',[200 300 600 380]);
bar(pct); ylim([0 100]);
xticks(1:numel(ratNames)); xticklabels(ratNames);
ylabel('% Trace coeff. sig (q<0.05)');
title(figTitle,'interpreter','none');
end
% -------------------------------------------------------------------------
function summarize(stats,ratNames,header)
fprintf('\n--- %s ---\n',header);
for k = 1:numel(ratNames)
    fprintf('%s : %4d / %4d  (%.1f%%)\n',ratNames{k}, ...
        sum(stats(k).sigTrace), stats(k).nNeurons, stats(k).pctSigTrace);
end
allPct = arrayfun(@(s)s.pctSigTrace,stats);
fprintf('Grand mean = %.1f ± %.1f %%\n',mean(allPct,'omitnan'),std(allPct,'omitnan'));
end
