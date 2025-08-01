function plotTaskVsNonTask(ratNames)
% plotTaskVsNonTask   Pre- vs Task-epoch rates (±MI/bps adjustment)
%
%   plotTaskVsNonTask({'rat0222','rat0307',...})
%
% Rows: one per rat, plus pooled row
% Cols:
%   1) normPre vs normTask
%   2) ratePre vs rateTask
%   3) residPre_MI vs residTask_MI
%   4) residPre_bps vs residTask_bps

  % epoch definitions (relative to CS)
  epochPre  = [-5, 0];
  epochTask = [ 0, 2];

  nRats = numel(ratNames);
  figure('Color','w','Position',[100 100 1200 300*(nRats+1)]);

  % pooled accumulators
  ALL.normPre       = [];
  ALL.normTask      = [];
  ALL.ratePre       = [];
  ALL.rateTask      = [];
  ALL.residPre_MI   = [];
  ALL.residTask_MI  = [];
  ALL.residPre_bps  = [];
  ALL.residTask_bps = [];

  for r = 1:nRats
    rat   = evalin('base', ratNames{r});
    days  = autoDateList(rat);
    iAn   = find(strcmp(days, rat.An),1);
    days  = days(iAn-2:iAn);  % last 3 days

    % --- build per-cell lists across days ----------------------------
    spikesAll = {};    % each entry = vector of spike times for one cell
    MI_spat   = [];    % spatial MI per cell
    MI_task   = [];    % task MI per cell
    B_spat    = [];    % spatial bits/s per cell
    B_task    = [];    % task bits/s per cell
    csAll     = [];    % all CS times

    sessionDur = 0;
    for d = 1:3
      day = days{d};
      % get data
      S   = rat.Ca_peaks.(sprintf('CA_peaks_%s',day));       % n_d×T
      M   = rat.ratemask.(sprintf('ratemask_%s',day))==1;    % 1×n_d
      cs  = rat.CS_times.(sprintf('CS_%s',day))(:);          % nTrials×1
      miS = rat.MI_noCSUS.(sprintf('MI_%s',day));            % n_d×1
      miT = rat.MI_CSUS15.(sprintf('MI_%s',day));            % n_d×1
      Bt  = rat.bitsper.(sprintf('MI_%s',day));   Bt = Bt(2,:)';  % n_d×1
      Bt_task = rat.bitsperCSUS.(sprintf('MI_%s',day));
      Bt_task = Bt_task(2,:)';                             % n_d×1

      % record session duration
      posT = rat.pos.(sprintf('pos_%s',day))(:,1);
      sessionDur = sessionDur + (posT(end)-posT(1));

      %get velocities
      pos = rat.pos.(sprintf('pos_%s',day));
      vel = ca_velocity(pos);


      good = find(M);
      csAll = [csAll; cs];

      for iCell = 1:numel(good)
        c = good(iCell);
        st = S(c,:); st = st(~isnan(st)&st>0);
        spikesAll{end+1,1} = st;               %#ok<AGROW>
        MI_spat(end+1,1)   = miS(c);           %#ok<AGROW>
        MI_task(end+1,1)   = miT(c);           %#ok<AGROW>
        B_spat(end+1,1)    = Bt(c);            %#ok<AGROW>
        B_task(end+1,1)    = Bt_task(c);       %#ok<AGROW>
      end
    end

    N = numel(spikesAll);

    % --- compute Pre/Task rates & normalize --------------------------
    ratePre  = nan(N,1);
    rateTask = nan(N,1);
    meanRate = nan(N,1);
    for i = 1:N
      st = spikesAll{i};
      meanRate(i) = numel(st) / sessionDur;
      cntP = 0; cntT = 0;
      for tt = 1:numel(csAll)
        cntP = cntP + sum(st>=csAll(tt)+epochPre(1)  & st<csAll(tt)+epochPre(2));
        cntT = cntT + sum(st>=csAll(tt)+epochTask(1)& st<csAll(tt)+epochTask(2));
      end
      ratePre(i)  = cntP  / (numel(csAll)*diff(epochPre));
      rateTask(i) = cntT  / (numel(csAll)*diff(epochTask));
    end
    normPre  = ratePre  ./ meanRate;
    normTask = rateTask ./ meanRate;




    % regress raw rates on spatial MI
    mask = ~isnan(MI_spat) & ~isnan(ratePre);
    Xmi  = [ones(sum(mask),1), MI_spat(mask)];
    b_pre = regress(ratePre(mask),  Xmi);
    residPre_MI = ratePre(mask) - Xmi*b_pre;

    b_task = regress(rateTask(mask), Xmi);
    residTask_MI = rateTask(mask) - Xmi*b_task;

    % --- residualize Pre & Task against spatial bits/s -------
    maskB = ~isnan(B_spat) & ~isnan(ratePre);
    Xb   = [ones(sum(maskB),1), B_spat(maskB)];
    b_preB = regress(ratePre(maskB),  Xb);
    residPre_bps = ratePre(maskB) - Xb*b_preB;

    b_taskB = regress(rateTask(maskB), Xb);
    residTask_bps = rateTask(maskB) - Xb*b_taskB;

    % --- store for pooled row ----------------------------------------
    ALL.normPre       = [ALL.normPre;       normPre];
    ALL.normTask      = [ALL.normTask;      normTask];
    ALL.ratePre       = [ALL.ratePre;       ratePre];
    ALL.rateTask      = [ALL.rateTask;      rateTask];
    ALL.residPre_MI   = [ALL.residPre_MI;   residPre_MI];
    ALL.residTask_MI  = [ALL.residTask_MI;  residTask_MI];
    ALL.residPre_bps  = [ALL.residPre_bps;  residPre_bps];
    ALL.residTask_bps = [ALL.residTask_bps; residTask_bps];

    % --- plot this rat’s 4 columns -----------------------------------
    for col = 1:4
      ax = subplot(nRats+1,4,(r-1)*4+col); hold(ax,'on'), axis(ax,'square');
      switch col
        case 1  % raw normalized
          x = normPre;   y = normTask;   ttl = 'norm Pre vs Task';
        case 2  % raw absolute
          x = ratePre;   y = rateTask;   ttl = 'rate Pre vs Task (Hz)';
        case 3  % MI‐adjusted residuals
          x = residPre_MI;   y = residTask_MI;   ttl = 'MI‐adj residuals';
        case 4  % bits/s‐adjusted residuals
          x = residPre_bps;  y = residTask_bps;  ttl = 'bps‐adj residuals';
      end
      scatter(ax,x,y,15,'filled');
      m = max([x; y]); plot(ax,[0 m],[0 m],'k--');
      [rV,pV] = corr(x,y,'Rows','pairwise');
      text(ax,0.05*m,0.9*m,sprintf('r=%.2f, p=%.3f',rV,pV),'FontSize',10);
      xlabel(ax,'Pre'); ylabel(ax,'Task');
      title(ax,ttl);
    end
  end

  % --- pooled “all cells” row ----------------------------------------
  row = nRats+1;
  for col = 1:4
    ax = subplot(nRats+1,4,(row-1)*4+col); hold(ax,'on'), axis(ax,'square');
    switch col
      case 1
        x = ALL.normPre;   y = ALL.normTask;  ttl = 'All: norm';
      case 2
        x = ALL.ratePre;   y = ALL.rateTask;  ttl = 'All: Hz';
      case 3
        x = ALL.residPre_MI;  y = ALL.residTask_MI; ttl = 'All: MI‐resid';
      case 4
        x = ALL.residPre_bps; y = ALL.residTask_bps; ttl = 'All: bps‐resid';
    end
    scatter(ax,x,y,15,'filled','r');
    m = max([x; y]); plot(ax,[0 m],[0 m],'k--');
    [rV,pV] = corr(x,y,'Rows','pairwise');
    text(ax,0.05*m,0.9*m,sprintf('r=%.2f, p=%.3f',rV,pV),'FontSize',10);
    xlabel(ax,'Pre'); ylabel(ax,'Task');
    title(ax,ttl);
  end
end
