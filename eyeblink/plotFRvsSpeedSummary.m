function plotFRvsSpeedSummary(method)
% -------------------------------------------------------------------------
% FR–vs-speed correlations across trials in the CS window (0–2 s) across three days
% (day–2, day–1, day–0) for each rat.
%
% Correlation is computed across trials (one Pearson r per neuron).
%
% Multiple‐comparison control:
%   method = 'fdr'  (default) → Benjamini–Hochberg, q = 0.05
%          = 'bonf'          → Bonferroni,        α = 0.05
%
% Population‐level test: compare mean(r) to shuffle‐based null.
%
% Needs in the base workspace:
%   • rat structs   (fields: Ca_peaks, pos, CS_times, An)
%   • helper fns    ca_velocity.m   autoDateList.m
% -------------------------------------------------------------------------
if nargin<1, method='fdr'; end
method = lower(method);
assert(ismember(method,{'fdr','bonf'}),'method must be ''fdr'' or ''bonf''');

% — user parameters —
qFDR     = 0.05;      % FDR level
alphaFW  = 0.05;      % family-wise for Bonferroni
minSpk   = 5;         % min spikes in each CS trial to include
win      = [0 2];     % CS window (s)
nShuff   = 5;       % # of shuffle replicates (permute trial speeds)
ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};

fprintf('\n===== plotFRvsSpeedSummary (%s correction) =====\n',upper(method));
nRats = numel(ratNames);

% prepare output structs
inc = makeBlankStruct(nRats);
all = makeBlankStruct(nRats);

for r = 1:nRats
  fprintf('Rat %s:\n',ratNames{r});
  rat   = evalin('base',ratNames{r});
  dates = autoDateList(rat);
  idx   = find(strcmp(dates,rat.An));
  days  = dates(idx-2:idx);

  % per‐rat collectors
  rAll    = [];  rInc    = [];
  shAll   = [];  shInc   = [];
  pAll    = [];  pInc    = [];

  for d=1:3
    day     = days{d};
    spkMat  = rat.Ca_peaks.(['CA_peaks_' day]);   % neurons × event times
    posRaw  = rat.pos.(['pos_' day]);
    csTimes = rat.CS_times.(['CS_' day]);

    % compute speed at each timestamp
    V     = ca_velocity(posRaw');
    speed = interp1(V(2,:),V(1,:),posRaw(:,1),'linear','extrap');
    speed(~isfinite(speed)) = 0;
    ts    = posRaw(:,1);

    % per‐trial FR & mean speed
    nTr    = numel(csTimes);
    fr     = nan(nTr,1);
    spd    = nan(nTr,1);
    for t=1:nTr
      t1 = csTimes(t)+win(1);
      t2 = csTimes(t)+win(2);
      % firing rate in CS window
      spks    = spkMat(:,:) ;  % we'll index per neuron below
      % we'll recompute per-neuron
    end

    % loop over neurons
    for ni = 1:size(spkMat,1)
      spks = spkMat(ni,:);
      spks = spks(~isnan(spks));

      % count total spikes in all CS trials
      nSpkWin = 0;
      for t=1:numel(csTimes)
        t1 = csTimes(t)+win(1);
        t2 = csTimes(t)+win(2);
        nSpkWin = nSpkWin + sum(spks>=t1 & spks<=t2);
      end
      passed = (nSpkWin>=minSpk);

      % build trial‐wise FR & mean speed
      for t=1:nTr
        t1 = csTimes(t)+win(1);
        t2 = csTimes(t)+win(2);
        fr(t)  = sum(spks>=t1 & spks<=t2)/diff(win);
        idxSp  = ts>=t1 & ts<=t2;
        spd(t) = mean(speed(idxSp));
      end

      % correlate FR vs. speed across trials
      if numel(fr)>=5 && numel(spd)>=5
        rval = corr(fr,spd,'type','Pearson');
        % get p‐value
        [~,pval] = corr(fr,spd,'type','Pearson');
      else
        rval = NaN;
        pval = NaN;
      end

      % shuffle null: permute trial speeds
      sh = nan(nShuff,1);
      for s=1:nShuff
        spdSh = spd(randperm(nTr));
        if numel(fr)>=5
          sh(s) = corr(fr,spdSh,'type','Pearson');
        end
      end

      % store ALL neurons
      rAll(end+1)     = rval;
      pAll(end+1)     = pval;
      shAll(end+1,:)  = sh';

      % store INCLUDED neurons
      if passed
        rInc(end+1)     = rval;
        pInc(end+1)     = pval;
        shInc(end+1,:)  = sh';
      end
    end
  end

  % multiple‐comparison per‐neuron
  [sigAllP,sigAllS] = sigTest(pAll, mean(abs(shAll-mean(shAll,2)),2)<abs(rAll), method,qFDR,alphaFW);
  [sigIncP,sigIncS] = sigTest(pInc, mean(abs(shInc-mean(shInc,2)),2)<abs(rInc), method,qFDR,alphaFW);

  % population‐level
  popNullAll = mean(shAll,2);
  obsAll     = mean(rAll,'omitnan');
  pPopAll(r) = mean(abs(popNullAll) >= abs(obsAll));

  popNullInc = mean(shInc,2);
  obsInc     = mean(rInc,'omitnan');
  pPopInc(r) = mean(abs(popNullInc) >= abs(obsInc));

  % write into structs


  all = updateStruct(all,r, rAll, mean(popNullAll), sigAllP, sigAllS);
  all.pAvsS(r)=pPopAll(r);
  inc = updateStruct(inc,r, rInc, mean(popNullInc), sigIncP, sigIncS);
  inc.pAvsS(r)=pPopInc(r);

  all.nIncl(r)=numel(rAll);    all.nTotal(r)=numel(rAll);
  inc.nIncl(r)=numel(rInc);    inc.nTotal(r)=numel(rAll);
end

% final tables & figures
printConsole(inc, ratNames, method, 'INCLUDED neuron-days');
printConsole(all, ratNames, method, 'ALL neuron-days');
makeFigure(inc,'INCLUDED neuron-days',method,ratNames);
makeFigure(all ,'ALL neuron-days',     method,ratNames);
end


%% HELPERS BELOW (unchanged from original)
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
    p = sort(p(:));
    m = numel(p);
    if m == 0
        thr = 0;
        return;
    end
    k = find(p <= (1:m)'/m*q, 1, 'last');
    if isempty(k)
        thr = 0;
    else
        thr = p(k);
    end
end

function S = updateStruct(S,r,C,Csh,sigP,sigS)
  S.meanCorr(r)            = mean(C,'omitnan');
  S.semCorr (r)            = std(C,'omitnan')/sqrt(max(1,numel(C)));
  S.meanShuff(r)           = Csh;
  S.semShuff (r)           = 0;  % popNull single value
  S.sigPercent.Pearson(r)  = 100*nanmean(sigP(:));
  S.sigPercent.Shuffle(r)  = 100*nanmean(sigS(:));
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
      C = S.allCorrByRat{k};
      fprintf('%-8s %4d/%-6d %11.3g %11.3f %11.3f %11.1f %11.1f\n',...
          ratNames{k}, S.nIncl(k), S.nTotal(k), S.pAvsS(k), ...
          mean(C,'omitnan'), std(C,'omitnan'), ...
          S.sigPercent.Pearson(k), S.sigPercent.Shuffle(k));
  end
  fprintf('############################################################\n\n');

end

function makeFigure(S,figTitle,method,ratNames)
  nR = numel(ratNames);
  figure('Color','w','Position',[100 300 1600 420],'Name',figTitle);
  sgtitle(sprintf('%s  (%s correction)',figTitle,upper(method)));

  subplot(1,4,1); hold on;
  bh = bar([S.meanCorr S.meanShuff],'grouped');
  errorbar(bh(1).XEndPoints,S.meanCorr,S.semCorr,'k.','CapSize',8);
  errorbar(bh(2).XEndPoints,S.meanShuff,S.semShuff,'k.','CapSize',8);
  xticks(1:nR); xticklabels(ratNames); ylabel('Mean FR–speed r');
  legend({'Actual','Shuffle'},'location','northwest'); title('Mean ± SEM r');
  fprintf('REMEMBER FOR ALL CELLS THE FIRST GRAPH WILL BE WRONG BUT ALSO IRRELEVANT BC YOU CANT AVERAGE FOR NON INCLUDED CELLS ANYWAY')

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
  ylabel('%Sig (Pearson)'); title('Sig. neurons');

  subplot(1,4,4);
  bar(S.sigPercent.Shuffle); ylim([0 100]);
  xticks(1:nR); xticklabels(ratNames);
  ylabel('%Sig (shuffle)'); title('Shuffle sig.');
end
