function traceVsSpeed_RESIDUAL(method,nPerm)
%% NOTEE: doesnt work great bc:
%your baseline and trace windows don’t have identical speed distributions AND
%your bins have very few spikes (so residuals are dominated by noise)

% traceVsSpeed_RESIDUAL  Speed-corrected ΔR per cell, with FDR/Bonferroni
%
% USAGE:
%    traceVsSpeed_RESIDUAL          % default: FDR, 500 perms
%    traceVsSpeed_RESIDUAL('bonf',1000)
%
% DEPENDENCIES: autoDateList, ca_velocity, mafdr
% -------------------------------------------------------------------------

if nargin<1, method = 'fdr'; end
if nargin<2, nPerm  = 5;   end
method = lower(method);

alphaFW = 0.05;    % family-wise Bonferroni
qFDR    = 0.05;    % FDR threshold

ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
preWin   = [-2 0];    trcWin = [ 0 2];
binSize  = 1/7.5;     nBins  = round(diff(trcWin)/binSize);
minSpk   = 5;         % inclusion threshold

nR = numel(ratNames);
statsAll = repmat(emptyStats,nR,1);
statsInc = repmat(emptyStats,nR,1);
allDR    = cell(nR,1);
incDR    = cell(nR,1);

fprintf('\n=== traceVsSpeed_RESIDUAL  (%s, %d perms) ===\n',...
        upper(method),nPerm);

for r = 1:nR
  rat = evalin('base',ratNames{r});
  dL  = autoDateList(rat);
  days= dL(find(strcmp(dL,rat.An))-2 : find(strcmp(dL,rat.An)));

  DRall = [];   spkAll = [];
  DRinc = [];   spkInc = [];

  %— LOOP DAYS ----------------------------------------------------------
  for d=1:3
    day   = days{d};
    spkM  = rat.Ca_peaks.(['CA_peaks_' day]);
    pos   = rat.pos.     (['pos_'      day]);
    csT   = rat.CS_times.(['CS_'       day]);

    ts    = pos(:,1);
    xy    = pos(:,2:3);

    % SPEED
    v     = ca_velocity([ts.'; xy(:,1).'; xy(:,2).']);
    speed = interp1(v(2,:),v(1,:),ts,'linear','extrap');
    speed(~isfinite(speed)) = 0;

    % PRUNE CS TIMES with incomplete windows
    valid = csT+preWin(1) >= ts(1) & csT+trcWin(2) <= ts(end);
    csT   = csT(valid);
    nTr   = numel(csT);
    if nTr==0, continue; end

    % PRECOMPUTE BIN‐EDGES
    baseE = arrayfun(@(t) linspace(csT(t)+preWin(1), csT(t)+preWin(2), nBins+1), ...
                     (1:nTr)','uni',false);
    trcE  = arrayfun(@(t) linspace(csT(t)+trcWin(1),csT(t)+trcWin(2),nBins+1), ...
                     (1:nTr)','uni',false);

    %— LOOP NEURONS -------------------------------------------------------
    for ci = 1:size(spkM,1)
      spk = spkM(ci,:);
      spk = spk(~isnan(spk));

      dR_trials = zeros(nTr,1)*nan;
      trcSpks   = zeros(nTr,1);

      %— LOOP TRIALS ------------------------------------------------------
      for t = 1:nTr
        edgesB = baseE{t};
        edgesT = trcE{t};
        % skip invalid edges
        if numel(edgesB)<2 || any(~isfinite(edgesB)) || any(diff(edgesB)<=0) ...
         || numel(edgesT)<2 || any(~isfinite(edgesT)) || any(diff(edgesT)<=0)
          continue
        end

        cntB = histcounts(spk, edgesB);
        cntT = histcounts(spk, edgesT);
        spdB = arrayfun(@(b) mean(speed(ts>=edgesB(b)   & ts<edgesB(b+1))), 1:nBins);
        spdT = arrayfun(@(b) mean(speed(ts>=edgesT(b)   & ts<edgesT(b+1))), 1:nBins);

        Y = [cntB cntT].';
        X = [ones(2*nBins,1) [spdB spdT].'];

        beta   = X\Y;
        resid  = Y - X*beta;
        dR_trials(t) = mean(resid(nBins+1:end)) - mean(resid(1:nBins));
        trcSpks(t)   = sum(cntT);

        [~, p] = ttest(dR_trials);      % two-tailed test of mean(dR_trials) ≠ 0


      end

      validIdx = ~isnan(dR_trials);
      if nnz(validIdx)<3
        continue
      end

      dRmean = mean(dR_trials(validIdx));
      totalSp= sum(trcSpks(validIdx));

      DRall = [DRall; dRmean];
      spkAll= [spkAll; totalSp];
      if totalSp>=minSpk
        DRinc = [DRinc; dRmean];
        spkInc= [spkInc; totalSp];
      end
    end
  end % days

  allDR{r} = DRall;
  incDR{r} = DRinc;

  %— SIGN-FLIP P-VALUES -----------------------------------------------
  pSF_all = arrayfun(@(x) mean(abs((rand(nPerm,1)>0.5)*2-1 .* x)>=abs(x)), DRall);
  pSF_inc = arrayfun(@(x) mean(abs((rand(nPerm,1)>0.5)*2-1 .* x)>=abs(x)), DRinc);

  statsAll(r) = packStats(pSF_all,method,alphaFW,qFDR,spkAll);
  statsInc(r) = packStats(pSF_inc,method,alphaFW,qFDR,spkInc);
end

%— SUMMARIZE & PLOT ------------------------------------------------------
summarize(statsAll,ratNames,'ALL neurons');
summarize(statsInc,ratNames,'INCLUDED neurons');
makeFigures(statsAll,ratNames,'ALL neurons',    allDR);
makeFigures(statsInc,ratNames,'INCLUDED neurons',incDR);
end


%% HELPERS ================================================================
function [sig,adj] = mcCorrect(p,method,a,q)
  if isempty(p)
    sig = false(size(p));
    adj = p;
    return
  end
  switch method
    case 'bonf'
      thr    = a/numel(p);
      sig    = p<thr;
      adj    = p;
    otherwise  % FDR
      adj    = mafdr(p,'BHFDR',true);
      sig    = adj<q;
  end
end

function S = packStats(pSF,method,a,q,spk)
  [sigSF,adjSF]  = mcCorrect(pSF,method,a,q);
  S.sigSF        = sigSF;
  S.pctSF        = 100*mean(sigSF);
  S.nNeurons     = numel(sigSF);
  S.spikes       = spk;
  S.adjSF        = adjSF;
end

function S = emptyStats
  S = struct('sigSF',[],'pctSF',NaN,'nNeurons',0,'spikes',[],'adjSF',[]);
end

function makeFigures(st,ratN,titleStr,DRcell)
  nR = numel(ratN);
  pct= arrayfun(@(x)x.pctSF,st);

  figure('Name',titleStr,'Color','w','Position',[200 300 1200 380]);
  sgtitle(titleStr,'FontWeight','bold');

  % 1: % sig
  subplot(1,3,1)
  bar(pct)
  ylim([0 100]), ylabel('% sig'), xticks(1:nR), xticklabels(ratN)
  title('% significant'), grid on

  % 2: box of ΔR
  subplot(1,3,2), hold on
  allV=[];grp=[];
  for k=1:nR
    allV=[allV; DRcell{k}];
    grp  =[grp;   k*ones(numel(DRcell{k}),1)];
  end
  if ~isempty(allV)
    boxplot(allV,grp,'Labels',ratN,'Symbol','k.')
  end
  yline(0,'k--'), ylabel('ΔR (spk/bin)'), title('per neuron ΔR'), grid on

  % 3: swarm
  subplot(1,3,3), hold on
  pool = vertcat(DRcell{:});
  if ~isempty(pool)
    swarmchart(ones(size(pool)),pool,'filled','MarkerFaceAlpha',0.4)
  end
  yline(0,'k--'), set(gca,'XTick',[])
  ylabel('ΔR (spk/bin)'), title('all neurons'), box on
end

function summarize(st,ratN,header)
  fprintf('\n--- %s ---\n',header);
  for i=1:numel(st)
    fprintf('%s: %3d / %3d   (%%SF = %4.1f%%)\n',ratN{i}, ...
            sum(st(i).sigSF),st(i).nNeurons,st(i).pctSF);
  end
  m=mean([st.pctSF]); sd=std([st.pctSF]);
  fprintf('Mean%%SF = %.1f ± %.1f%%\n\n',m,sd);
end
