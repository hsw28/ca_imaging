function traceVsSpeed_GLM
% TRACEVSSPEED_GLM --------------------------------------------------------
% Quantifies **trace-period firing** while controlling for running speed.
%
% ───────── WHAT THE SCRIPT DOES ──────────────────────────────────────────
% * For each rat (last 3 training days) and every neuron:
%       1.  Bin the 2 s *trace* window (CS→US) into 15 × 133 ms bins.
%       2.  Collapse the preceding 2 s *baseline* window (–2→0 s) into
%           **one row** whose spike count is the *sum* of those 15 bins and
%           whose speed is the *mean* speed in that window.
%           → 16 design-matrix rows per trial: 1 baseline  + 15 trace.
%       3.  Build predictors per row
%              – z-scored speed
%              – centred Trace dummy  (baseline = −0.47, trace = +0.53)
%              – intercept (explicit column so we can switch it off in glmfit)
%       4.  Fit a **Poisson GLM** with log link
%              log λ = β₀ + β_speed · speedZ + β_trace · Tracẽ
%       5.  Correct standard errors and p-values with
%              quasi-Poisson dispersion φ = deviance / d.f.e.
%              (β̂ stays unbiased; SE ← SE·√φ ; p from z = β̂/SE)
%       6.  Average β_trace and its p across trials for that neuron.
% * Multiple-comparison: per-rat Benjamini–Hochberg FDR (q = 0.05).
% * Two summary figures:  all neurons vs neurons with ≥ minSpk spikes in
%   the **trace** window.
%
% ───────── WHY “ONE BASELINE ROW” WORKS ──────────────────────────────────
%   Collapsing the 15 baseline bins into a single aggregate row
%   ① retains the *trace-vs-baseline* contrast (β_trace is still
%      log(rate_trace / rate_baseline) after speed correction),
%   ② greatly reduces sparsity (fewer all-zero rows),
%   ③ stabilises IRLS and all but eliminates “ill-conditioned weights”
%      warnings without losing interpretability.
%
% Dependencies
%   ·  rat structs in base workspace with fields
%                  Ca_peaks.*, pos.*, CS_times.*, An
%   ·  helper functions:  autoDateList, ca_velocity
%   ·  MATLAB Statistics & ML Toolbox  (glmfit, mafdr)
%
% ------------------------------------------------------------------------

% ---------------- USER PARAMETERS ----------------------------------------
ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};

preWin  = [-2 0];          % baseline (s rel. CS)  → collapsed to 1 row
trcWin  = [ 0 2];          % trace window (s rel. CS)
binSize = 1/7.5;           % Ca frame 133 ms
nBins   = round(diff(trcWin)/binSize);   % 15 trace bins

minSpk  = 5;               % ≥ spikes in TRACE window to be “included”
alphaFW = 0.05;            % per-rat FDR threshold
% -------------------------------------------------------------------------

% ---------------- HOUSE-KEEPING -----------------------------------------
warning('off','stats:glmfit:IterationLimit');   % quiet convergence chatter
warning('off',lastwarnID('stats:glmfit'));      % quiet ill-conditioned msg

opts = statset('MaxIter',200);                 % IRLS cap

nRats = numel(ratNames);
statsAll = repmat(emptyStats,nRats,1);
statsInc = repmat(emptyStats,nRats,1);

fmtHdr = '\n=== %s =================================\n';

% ===================== PER-RAT LOOP ======================================
for r = 1:nRats
    fprintf(fmtHdr,ratNames{r});
    rat  = evalin('base',ratNames{r});
    dStr = autoDateList(rat);
    days = dStr(find(strcmp(dStr,rat.An))-2 : find(strcmp(dStr,rat.An)));

    % per-rat collectors
    betaT_all = [];  betaS_all = [];  pT_all = [];   nSpk_all = [];
    betaT_inc = [];  betaS_inc = [];  pT_inc = [];   nSpk_inc = [];

    % ---------------- PER-DAY LOOP --------------------------------------
    for d = 1:3
        day = days{d};

        spikeMat = rat.Ca_peaks.(['CA_peaks_' day]);   % cells × times
        ts       = rat.pos.(['pos_' day])(:,1);        % time vector
        posXY    = rat.pos.(['pos_' day])(:,2:3);      % x,y
        csTimes  = rat.CS_times.(['CS_' day]);         % CS onsets

        % ---- speed trace ------------------------------------------------
        velDat = ca_velocity([ ts.' ; posXY(:,1).' ; posXY(:,2).' ]);
        speed  = interp1(velDat(2,:),velDat(1,:),ts,'linear','extrap');
        speed(~isfinite(speed)) = 0;

        nTrials = numel(csTimes);

        % ---- pre-compute edges & speed bins -----------------------------
        preEdges = nan(nTrials,2);                 % collapsed baseline
        trcEdges = nan(nTrials,nBins+1);
        speedTrc = nan(nTrials,nBins);
        speedPre = nan(nTrials,1);                 % mean speed baseline

        for t = 1:nTrials
            % baseline single interval
            preEdges(t,:) = [csTimes(t)+preWin(1), csTimes(t)+preWin(2)];
            idxB = ts>=preEdges(t,1) & ts<preEdges(t,2);
            speedPre(t) = mean(speed(idxB));

            % trace 15 bins
            trcEdges(t,:) = linspace(csTimes(t)+trcWin(1), csTimes(t)+trcWin(2), nBins+1);
            for b = 1:nBins
                idxT = ts>=trcEdges(t,b) & ts<trcEdges(t,b+1);
                speedTrc(t,b) = mean(speed(idxT));
            end
        end

        % ---------------- PER-NEURON LOOP -------------------------------
        for ni = 1:size(spikeMat,1)
            spk = spikeMat(ni,:); spk = spk(~isnan(spk));

            betaS_tr = nan(nTrials,1);   % β_speed per trial
            betaT_tr = nan(nTrials,1);   % β_trace per trial
            pT_tr    = nan(nTrials,1);   % p-val per trial

            for t = 1:nTrials
                % ----- counts -------------------------------------------
                cntPre = sum(spk >= preEdges(t,1) & spk < preEdges(t,2)); % scalar
                cntTrc = histcounts(spk,trcEdges(t,:));                   % 1×15

                Y      = [cntPre; cntTrc(:)];                            % 16×1
                spdAll = [speedPre(t); speedTrc(t,:).'];                 % 16×1

                Trace  = [0; ones(nBins,1)];                             % baseline=0
                Trace  = Trace - mean(Trace);                            % centre (−.47,+.53)

                ok = ~isnan(Y) & ~isnan(spdAll);
                if nnz(ok) < 4,  continue, end   % need enough rows

                X = [ zscore(spdAll(ok)) , Trace(ok) , ones(nnz(ok),1) ];

                % Poisson GLM (intercept in X, so constant off)
                [B,dev,st] = glmfit(X(:,1:2),Y(ok),'poisson','link','log',...
                                     'constant','off','options',opts);

                % quasi-Poisson (φ) correction
                phi  = max(dev/st.dfe,1);
                seQ  = st.se*sqrt(phi);
                zQ   = B./seQ;
                pQ   = 2*(1-normcdf(abs(zQ)));

                betaS_tr(t) = B(1);
                betaT_tr(t) = B(2);
                pT_tr(t)    = pQ(2);
            end % trial loop

            % ---- aggregate across trials --------------------------------
            betaS_all = [betaS_all; mean(betaS_tr,'omitnan')];
            betaT_all = [betaT_all; mean(betaT_tr,'omitnan')];
            pT_all    = [pT_all   ; mean(pT_tr ,'omitnan')];

            % spike inclusion count in TRACE window
            nSpkWin = 0;
            for t = 1:nTrials
                nSpkWin = nSpkWin + sum(spk >= trcEdges(t,1) & spk < trcEdges(t,end));
            end
            nSpk_all = [nSpk_all; nSpkWin];

            % ---- inclusion subset ---------------------------------------
            if nSpkWin >= minSpk
                betaS_inc = [betaS_inc; mean(betaS_tr,'omitnan')];
                betaT_inc = [betaT_inc; mean(betaT_tr,'omitnan')];
                pT_inc    = [pT_inc   ; mean(pT_tr ,'omitnan')];
                nSpk_inc  = [nSpk_inc ; nSpkWin];
            end
        end % neuron loop
    end % day loop

    % ---------------- FDR per rat ----------------------------------------
    statsAll(r) = finishStats(betaT_all,betaS_all,pT_all,alphaFW,nSpk_all);
    statsInc(r) = finishStats(betaT_inc,betaS_inc,pT_inc,alphaFW,nSpk_inc);
end % rat loop

% ---------------- PLOTS & CONSOLE ---------------------------------------
makeFigures(statsAll,ratNames,'ALL neurons');
makeFigures(statsInc,ratNames,sprintf('INCLUDED (≥%d spikes in Trace)',minSpk));

summarize(statsAll,ratNames,'ALL neurons');
summarize(statsInc,ratNames,sprintf('INCLUDED (≥%d spikes in Trace)',minSpk));
end
% ========================================================================
%                               HELPERS
% ========================================================================
function id = lastwarnID(startsWithStr)
   % helper to fetch unknown glmfit warning id for this MATLAB
   lastwarn(''); warning('query'); x=lastwarn; [~,id]=lastwarn;
   if contains(id,startsWithStr), return, else id=''; end
end
function S = emptyStats
S = struct('betaTrace',[],'betaSpeed',[],'qTrace',[],'sigTrace',[], ...
           'pctSigTrace',NaN,'nNeurons',0,'spikes',[]);
end
function S = finishStats(bT,bS,pT,alphaFW,spkCnt)
if isempty(pT), S = emptyStats; return; end
q      = mafdr(pT,'BHFDR',true);
sig    = q < alphaFW;
S      = struct('betaTrace',bT,'betaSpeed',bS,'qTrace',q, ...
                'sigTrace',sig,'pctSigTrace',100*mean(sig), ...
                'nNeurons',numel(pT),'spikes',spkCnt);
end
function makeFigures(stats,ratNames,ttl)
pct = arrayfun(@(s)s.pctSigTrace,stats);
figure('color','w','name',ttl,'pos',[200 300 600 380]);
bar(pct); ylim([0 100]); xticks(1:numel(ratNames)); xticklabels(ratNames);
ylabel('% Trace coeff. sig (q<0.05)'); title(ttl,'interpreter','none');
end
function summarize(stats,ratNames,header)
fprintf('\n--- %s ---\n',header);
for k = 1:numel(ratNames)
    fprintf('%s : %4d / %4d  (%.1f%%)\n', ratNames{k}, ...
        sum(stats(k).sigTrace), stats(k).nNeurons, stats(k).pctSigTrace);
end
g = arrayfun(@(s)s.pctSigTrace,stats);
fprintf('Grand mean = %.1f ± %.1f %%\n',mean(g,'omitnan'),std(g,'omitnan'));
end
