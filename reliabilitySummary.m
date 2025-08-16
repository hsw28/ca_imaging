function reliabilitySummary(ratNames)
% reliabilitySummary  Compare spike‐count variability in place field vs trace
%
% Splits each cell’s spikes into in-trial vs out-of-trial based on a 2 s window after each CS.
% gets smoothed rate maps for those two conditions (2.5 cm spatial bins).
% For non task time, we threshold each map at mean + 1·std and combines them into a binary mask of “significant” spatial bins for each cell.
% then for each cell:
  % for each spike, you
  % find its (x,y) bin,
  % check whether that bin sits inside the cell’s place-field mask,
  % record a 1 (in-field) or 0 (out-of-field).
  %countsTr: for each CS onset, you count how many spikes fall in the 0–2 s task window.
  %then compute variability (Fano) and reliability (split-half) for task period and non task period
  % graphs variability (Fano factor) and  Reliability (split‐half correlation) for each cell in task period vs place field




%   reliabilitySummary({'rat0222','rat0307',...})

nRats = numel(ratNames);
figure('Color','w','Position',[100 100 1200 300*(nRats+1)]);

allFanoPF   = [];
allFanoTr   = [];
allSplitPF  = [];
allSplitTr  = [];

for r = 1:nRats
  rat   = evalin('base', ratNames{r});
  dates = autoDateList(rat);
  iAn   = find(strcmp(dates,rat.An),1);
  days  = dates(iAn-2:iAn);

  fanoPF = [];
  fanoTr = [];
  spPF   = [];
  spTr   = [];

  for d = 1:3
    D       = days{d};
    S       = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));     % nCells×T
    good    = find(rat.ratemask.(sprintf('ratemask_%s',D))==1);
    pos     = rat.pos.(sprintf('pos_%s',D));               % Tpos×3
    csTimes = rat.CS_times.(sprintf('CS_%s',D))(:);
    t_pos   = pos(:,1);
    xy      = pos(:,2:3);


%    S = S(good,:);



    % --- build PF mask grid and get edges ---
    [maskGrid, xEdges, yEdges] = computePFmaskGrid2p5(S, pos, csTimes, 1);
    [nY,nX,~] = size(maskGrid);

    for c = good(:)'
      st = S(c,:);
      st = st(~isnan(st)&st>0);
      if numel(st)<3, continue, end

      % bin spike positions
      idx_spk = interp1(t_pos, (1:numel(t_pos))', st, 'nearest','extrap');
      bx      = discretize(xy(idx_spk,1), xEdges);
      by      = discretize(xy(idx_spk,2), yEdges);
      thisMask = maskGrid(:,:,c);

      % identify PF visits
      validSpk = bx>=1 & bx<=nX & by>=1 & by<=nY;
      vPF = thisMask(sub2ind([nY,nX], by(validSpk), bx(validSpk)));
      % collapse runs: treat each spike as one count
      % can refine to count per visit if desired
      countsPF = vPF;

      % trace counts
      countsTr = arrayfun(@(t0) sum(st>=t0 & st< t0+2), csTimes);

           % only accumulate when we have >1 spike in both PF and Trace
           if numel(countsPF)>1 && numel(countsTr)>1
             fanoPF(end+1) = var(countsPF)/mean(countsPF);
             fanoTr(end+1) = var(countsTr)/mean(countsTr);
             spPF(end+1)   = splitHalfCorr(countsPF);
             spTr(end+1)   = splitHalfCorr(countsTr);
           end
    end
  end

  % aggregate
  allFanoPF  = [allFanoPF;  fanoPF(:)];
  allFanoTr  = [allFanoTr;  fanoTr(:)];
  allSplitPF = [allSplitPF; spPF(:)];
  allSplitTr = [allSplitTr; spTr(:)];


  % plotting per rat
  ax1 = subplot(nRats+1,2, (r-1)*2+1);
  hold(ax1,'on');
  %axis(ax1,'square');
  scatter(ax1,fanoPF, fanoTr, 12, 'filled');
  axis([0 max(fanoTr), 0 max(fanoTr)])
  plot(ax1, xlim, xlim,'k--');
  xlabel(ax1,'Fano Place Field'); ylabel(ax1,'Fano Trace'); title(ax1, ratNames{r});

  ax2 = subplot(nRats+1,2, r*2);
  hold(ax2,'on');
  %axis(ax2,'square');
  scatter(ax2,spPF, spTr, 12,'filled');
  plot(ax2, xlim, xlim,'k--');
  xlabel(ax2,'Split-half Place Field'); ylabel(ax2,'Split-half Trace');
end

% pooled row
ax1 = subplot(nRats+1,2, nRats*2+1);
hold(ax1,'on'); %axis(ax1,'square');
scatter(ax1,allFanoPF,allFanoTr,12,'r','filled');
  axis([0 max(allFanoTr), 0 max(allFanoTr)])
plot(ax1, xlim, xlim,'k--');
xlabel(ax1,'Fano Place Field'); ylabel(ax1,'Fano Trace'); title(ax1,'All pooled');
ax2 = subplot(nRats+1,2, nRats*2+2);
hold(ax2,'on'); %axis(ax2,'square');
scatter(ax2,allSplitPF,allSplitTr,12,'r','filled'); plot(ax2, xlim, xlim,'k--');
xlabel(ax2,'Split-half Place Field'); ylabel(ax2,'Split-half Trace'); title(ax2,'All pooled');
end

function [maskGrid,xEdges,yEdges] = computePFmaskGrid2p5(spikeMat,pos,cs_times,N)
% reuse CA_normalizePosData grid

nCells = size(spikeMat,1);
t_pos = pos(:,1);
isInTrial = false(size(t_pos)); for t=cs_times(:)', isInTrial = isInTrial | (t_pos>=t & t_pos< t+2); end

makeNewGrid = 0;
for c=1:nCells
  st = spikeMat(c,:); st = st(~isnan(st)&st>0);
  if numel(st)<3, continue; end
  idx_spk = interp1(t_pos,(1:numel(t_pos))', st, 'nearest','extrap');
  inTr = isInTrial(idx_spk);
  if size(st(inTr),2)==0 || size(st(~inTr),2) ==0
    continue
  end
  rateT = CA_normalizePosData(st(inTr), pos, 2.5, 1);
  rateN = CA_normalizePosData(st(~inTr), pos, 2.5, 1);
  if makeNewGrid == 0
    [xEdges,yEdges] = size(rateT);
    maskGrid = false(xEdges,yEdges,nCells);
  end
  makeNewGrid = 1;
  rt = rateT;
  rn=rateN;
  maskGrid(:,:,c) = (rt>nanmean(rt(:))+N*nanstd(rt(:))) | (rn>nanmean(rn(:))+N*nanstd(rn(:)));
end

end

function r = splitHalfCorr(v)
  n   = numel(v);
  idx = randperm(n);
  h   = floor(n/2);
  idx = idx(1:2*h);                  % only use the first 2*h spikes
  r   = corr( v(idx(1:h)), ...
              v(idx(h+1:2*h)), ...
              'Rows','pairwise');
end
