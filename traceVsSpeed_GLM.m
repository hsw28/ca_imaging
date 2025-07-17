function traceVsSpeed_GLM
%note: DATA TOO SPARSE FOR THIS
%
%
%Quantifies **trace‐period firing** while controlling for running speed.
%
% ───────── WHAT THE SCRIPT DOES ──────────────────────────────────────────
% * For each rat (last 3 training days) and every neuron:
%       1.  Bin the 2 s *trace* window (CS→US) into 15 × 133 ms bins.
%       2.  Collapse the preceding 2 s *baseline* window (–2→0 s) into
%           **one row** whose spike count is the *sum* of those 15 bins and
%           whose speed is the *mean* speed in that window.
%           → 16 design‐matrix rows per trial: 1 baseline  + 15 trace.
%       3.  Build predictors per row
%              – z‐scored speed
%              – centred Trace dummy  (baseline = –15/16, trace = +1/16)
%              – intercept
%       4.  Fit a **Poisson GLM** with log link, then do quasi‐Poisson
%           correction on stderr/p‐values.
%       5.  Average β_trace and its p across trials for that neuron.
% * Multiple‐comparison: per‐rat Benjamini–Hochberg FDR (q = 0.05).
% * Two summary figures:  all neurons vs neurons with ≥ minSpk spikes in
%   the **trace** window.
%
% Dependencies
%   ·  rat structs in base workspace with fields
%                  Ca_peaks.*, pos.*, CS_times.*, An
%   ·  helper functions:  autoDateList, ca_velocity
%   ·  MATLAB Statistics & ML Toolbox  (glmfit, mafdr)
% ------------------------------------------------------------------------

%% USER PARAMETERS
ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};

preWin   = [-2  0];     % baseline window rel. CS
trcWin   = [ 0  2];     % trace  window rel. CS
binSize  = 1/7.5;       % 133 ms bins
nBins    = round(diff(trcWin)/binSize);  % 15

minSpk   = 5;           % min spikes in trace to “include”
alphaFW  = 0.05;        % FDR q‐level

%% QUIET WARNINGS
warning('off','stats:glmfit:IterationLimit');
warning('off','stats:glmfit:IllConditioned');

%% LOOP OVER RATS
nRats    = numel(ratNames);
statsAll = repmat(emptyStats, nRats,1);
statsInc = repmat(emptyStats, nRats,1);

for r = 1:nRats
  fprintf('\n=== Rat %s ===\n', ratNames{r});
  rat = evalin('base', ratNames{r});
  dates = autoDateList(rat);
  idx   = find(strcmp(dates, rat.An));
  days  = dates(idx-2:idx);

  % collectors
  betaT_all = []; pT_all = []; nSpk_all = [];
  betaT_inc = []; pT_inc = []; nSpk_inc = [];

  for d = 1:3
    day     = days{d};
    spkMat  = rat.Ca_peaks.(['CA_peaks_' day]);  % cells×times
    posData = rat.pos.      (['pos_'      day]);  % [t×[time x y]]
    csTimes = rat.CS_times. (['CS_'       day]);

    ts      = posData(:,1);
    xy      = posData(:,2:3);

    %% PRE‐COMPUTE SPEED
    vel    = ca_velocity([ts.'; xy(:,1).'; xy(:,2).']);
    speed  = interp1(vel(2,:), vel(1,:), ts, 'linear','extrap');
    speed(~isfinite(speed)) = 0;

    nTrials = numel(csTimes);
    % build edge arrays
    preEdges = csTimes(:) + preWin;
    trcEdges = cell(nTrials,1);
    speedPre = nan(nTrials,1);
    speedTrc = nan(nTrials,nBins);
    for t=1:nTrials
      be = preEdges(t,:);
      idxB = ts>=be(1)&ts<be(2);
      speedPre(t)=mean(speed(idxB));
      te = linspace(csTimes(t)+trcWin(1),csTimes(t)+trcWin(2),nBins+1);
      trcEdges{t}=te;
      for b=1:nBins
        idxT = ts>=te(b)&ts<te(b+1);
        speedTrc(t,b)=mean(speed(idxT));
      end
    end

    %% ONE GLM PER CELL
    [nCells,~]=size(spkMat);
    for ci=1:nCells
      spk = spkMat(ci,:);
      spk = spk(~isnan(spk));

      Nrows = nTrials*(1+nBins);
      Y     = zeros(Nrows,1);
      spd   = zeros(Nrows,1);
      Trc   = zeros(Nrows,1);

      row=0; traceSpkCount=0;
      for t=1:nTrials
        % baseline row
        be = preEdges(t,:);
        row=row+1;
        Y(row)=sum(spk>=be(1)&spk<be(2));
        spd(row)=speedPre(t);
        Trc(row)=-nBins/16;  % centre: sum(trace)=+15/16

        % trace bins
        te = trcEdges{t};
        for b=1:nBins
          row=row+1;
          Y(row)=sum(spk>=te(b)&spk<te(b+1));
          spd(row)=speedTrc(t,b);
          Trc(row)=+1/16;
          traceSpkCount = traceSpkCount + Y(row);
        end
      end

      % valid rows
      ok = ~isnan(Y)&~isnan(spd);
      if nnz(ok)<10, continue, end

      % design
      X = [ zscore(spd(ok)), Trc(ok), ones(nnz(ok),1) ];
      % glmfit
      [B, dev, ST] = glmfit(X(:,1:2), Y(ok), 'poisson','link','log','constant','off','options',statset('MaxIter',200));
      phi = max(dev/ST.dfe,1);
      seQ = ST.se*sqrt(phi);
      zQ  = B./seQ;
      pQ  = 2*(1-normcdf(abs(zQ)));

      % collect
      betaT_all(end+1,1)=B(2);
      pT_all   (end+1,1)=pQ(2);
      nSpk_all (end+1,1)=traceSpkCount;

      if traceSpkCount>=minSpk
        betaT_inc(end+1,1)=B(2);
        pT_inc   (end+1,1)=pQ(2);
        nSpk_inc (end+1,1)=traceSpkCount;
      end
    end
  end

  % FDR per rat
  statsAll(r)=finishStats(betaT_all, pT_all, alphaFW, nSpk_all);
  statsInc(r)=finishStats(betaT_inc, pT_inc, alphaFW, nSpk_inc);
end

%% PLOT & PRINT
makeFigures(statsAll, ratNames, 'ALL cells');
makeFigures(statsInc, ratNames, 'INCLUDED cells');

summarize(statsAll, ratNames, 'ALL cells');
summarize(statsInc, ratNames, 'INCLUDED cells');
end


%%==========================================================================%%
function S=emptyStats
  S=struct('betaTrace',[],'qTrace',[],'sigTrace',[],'pctSigTrace',nan,'nCells',0,'spikes',[]);
end
function S=finishStats(bT,pT,alpha,nSpk)
  if isempty(pT), S=emptyStats; return; end
  q=mafdr(pT,'BHFDR',true);
  h=q<alpha;
  S=struct('betaTrace',bT,'qTrace',q,'sigTrace',h, ...
           'pctSigTrace',100*mean(h),'nCells',numel(pT),'spikes',nSpk);
end
function makeFigures(stats, names, ttl)
  pct=arrayfun(@(s)s.pctSigTrace,stats);
  figure('Name',ttl,'Color','w');
  bar(pct); ylim([0 100]);
  xticks(1:numel(names)); xticklabels(names);
  ylabel('% sig β\_trace'); title(ttl);
end
function summarize(stats,names,ttl)
  fprintf('\n--- %s ---\n',ttl);
  for k=1:numel(names)
    s=stats(k);
    fprintf('%s : %4d / %4d (%.1f%%)\n', names{k}, sum(s.sigTrace), s.nCells, s.pctSigTrace);
  end
  g=arrayfun(@(s)s.pctSigTrace,stats);
  fprintf('Grand mean: %.1f%% ± %.1f\n', mean(g), std(g));
end
