function R = epochModulation_venn(ratNames, varargin)
% epochModulation_venn
% Detect per-cell firing-rate modulation in CS / Trace / US / Post epochs
% compared with FULL non-trial time (optional speed restriction),
% optionally limited to place cells, and visualize overlaps (UpSet or Venn).
%
% Usage:
%   R = epochModulation_venn({'rat0222','rat0314'}, ...
%       'Alpha',0.05, ...
%       'RestrictSpeed',false, ...
%       'PlotMode','venn', ...
%       'VennEpochs',{'CS','Trace','US','Post'}, ...
%       'VennPercent',true, ...
%       'PlaceOnly',true);

% ---------------- parameters ----------------
p = inputParser;
p.addParameter('Alpha',0.05);
p.addParameter('RestrictSpeed',false);
p.addParameter('SpeedMin',0.0);
p.addParameter('SpeedMax',4);
p.addParameter('SpeedSource','auto');
p.addParameter('PlotMode','venn');
p.addParameter('VennEpochs',{'CS','Trace','US','Post'});
p.addParameter('CollapsePreIntoOne', true, @islogical);  % CS∪Trace∪US vs Post
p.addParameter('VennPercent',true,@islogical);
p.addParameter('PercentBase','kept',@(s)any(strcmpi(s,{'kept','union'})));
p.addParameter('Verbose',true);
% place-cell options
p.addParameter('PlaceOnly',true,@islogical);
p.addParameter('PlaceSpec',struct('Struct','MI_wCSUS_shuff','Prefix','MI_', ...
    'Col',3,'Thresh',0.95,'Comparator','>'));
p.parse(varargin{:});
opt = p.Results;

labels = {'CS','Trace','US','Post'};
epochs = [0.00 0.25; 0.25 0.75; 0.75 0.85; 0.85 2.00];
trialWin = [0 6.0];

nR = numel(ratNames);
ALL_sets = [];

for r = 1:nR
    rat = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates,rat.An),1);
    days  = dates(max(1,idx-2):idx);

    for d = 1:numel(days)
        dayStr = days{d};
        spk = rat.Ca_peaks.(sprintf('CA_peaks_%s',dayStr));
        csTimes = rat.CS_times.(sprintf('CS_%s',dayStr));
        mask = rat.ratemask.(sprintf('ratemask_%s',dayStr))==1;

        % optional: place cells only
        if opt.PlaceOnly
            [pcMask, okPC] = fetchPlaceMask(rat, dayStr, opt.PlaceSpec);
            if okPC
                mask = mask & pcMask;
            else
                warning('[%s %s] PlaceOnly requested but MI table missing.', ratNames{r}, dayStr);
            end
        end

        if isempty(csTimes), continue; end
        sessMin = min(csTimes)-60; sessMax = max(csTimes)+60;
        excl = [csTimes(:)+trialWin(1), csTimes(:)+trialWin(2)];
        pool = complementIntervals([sessMin sessMax], excl);

        poolDur = sum(max(0,pool(:,2)-pool(:,1)));
        if poolDur<=0, continue; end

        nEp=size(epochs,1);
        epAbs=cell(nEp,1); epDur=zeros(nEp,1);
        for e=1:nEp
            epAbs{e}=[csTimes(:)+epochs(e,1),csTimes(:)+epochs(e,2)];
            epDur(e)=diff(epochs(e,:));
        end
        nTrials=numel(csTimes); nCells=size(spk,1);
        pvals=nan(nCells,nEp); dirsgn=zeros(nCells,nEp,'int8'); keep=false(nCells,1);

        for c=1:nCells
            if ~mask(c), continue; end
            st=spk(c,:); st=st(~isnan(st)&st>0);
            if isempty(st), continue; end
            keep(c)=true;
            poolCounts=intervalCountsPerRow(st,pool);
            r_pool=sum(poolCounts)/poolDur;
            for e=1:nEp
                counts=intervalCountsPerRow(st,epAbs{e});
                obsK=sum(counts); expT=epDur(e)*nTrials;
                lambda=r_pool*expT;
                pvals(c,e)=poisson_two_sided(obsK,lambda);
                if pvals(c,e)<=opt.Alpha
                    dirsgn(c,e)=int8(sign(double(obsK)-lambda));
                end
            end
        end

        sigs=pvals<=opt.Alpha;
        ALL_sets=[ALL_sets; sigs(keep,:)];
        % --- Optional collapse: (CS ∪ Trace ∪ US) vs Post ---
          plotSets   = ALL_sets;
          plotLabels = {'CS','Trace','US','Post'};
          if opt.CollapsePreIntoOne
              if size(ALL_sets,2) < 4
                  warning('CollapsePreIntoOne requested, but ALL_sets has <4 columns. Skipping collapse.');
              else
                  plotSets   = [ any(ALL_sets(:,1:3), 2), ALL_sets(:,4) ];  % [N x 2]
                  plotLabels = {'CS/Trace/US', 'Post'};
              end
          end
        if opt.Verbose
            fprintf('[%s %s] %d cells kept, sig CS=%d Trace=%d US=%d Post=%d\n',...
                ratNames{r},dayStr,nnz(keep),nnz(sigs(:,1)),nnz(sigs(:,2)),nnz(sigs(:,3)),nnz(sigs(:,4)));
        end
    end
end

% ---------- plotting ----------
switch lower(opt.PlotMode)
    case 'upset'
        figure('Color','w');
        upsetPlotGeneric(plotSets, plotLabels);   % <— NEW generic upset (handles 2–4)
      case 'venn'
      % Build plotSets / plotLabels first (handles CollapsePreIntoOne)
      plotSets   = ALL_sets;
      plotLabels = {'CS','Trace','US','Post'};
      if isfield(opt,'CollapsePreIntoOne') && opt.CollapsePreIntoOne
          if size(ALL_sets,2) < 4
              warning('CollapsePreIntoOne requested, but ALL_sets has <4 columns. Skipping collapse.');
          else
              plotSets   = [ any(ALL_sets(:,1:3), 2), ALL_sets(:,4) ];  % [N x 2]
              plotLabels = {'CS/Trace/US','Post'};
          end
      end

      % If user specified VennEpochs and we didn't collapse, honor them.
      if ~(isfield(opt,'CollapsePreIntoOne') && opt.CollapsePreIntoOne)
          [epMask, labsK] = selectEpochSubset(plotLabels, opt.VennEpochs, 4);
          setsK  = plotSets(:, epMask);
      else
          setsK  = plotSets;     % already collapsed to 2 sets
          labsK  = plotLabels;
      end

      kSel = size(setsK,2);
      if kSel < 2
          warning('Venn requires ≥2 sets; got %d. Skipping plot.', kSel);
          return
      end

      % denominator for % labels
      switch lower(opt.PercentBase)
          case 'union', denom = max(1, sum(any(setsK,2)));
          otherwise,     denom = size(setsK,1);  % 'kept'
      end
      asPct = opt.VennPercent;

      % Build region labels and call your venn() function
      if kSel == 2
          C   = computeExclusiveCounts(setsK);
          lbl = { fmtCount(C.A_only,denom,asPct), ...
                  fmtCount(C.B_only,denom,asPct), ...
                  fmtCount(C.AB,     denom,asPct) };
          col = parula(5);
          col = col([2,5],:);
          venn(2, 'sets', labsK, 'labels', lbl, 'colors', col, ...
                  'alpha', 0.5, 'edgeC', [1 1 1], 'edgeW', 3);

      elseif kSel == 3
          C   = computeExclusiveCounts(setsK);
          lbl = { fmtCount(C.A_only,denom,asPct), fmtCount(C.B_only,denom,asPct), fmtCount(C.C_only,denom,asPct), ...
                  fmtCount(C.AB,denom,asPct),     fmtCount(C.AC,denom,asPct),     fmtCount(C.BC,denom,asPct), ...
                  fmtCount(C.ABC,denom,asPct) };
          col = parula(3);
          venn(3, 'sets', labsK, 'labels', lbl, 'colors', col, ...
                  'alpha', 0.5, 'edgeC', [1 1 1], 'edgeW', 3);

      elseif kSel == 4
          C   = computeExclusiveCounts4(setsK);
          lbl = { fmtCount(C.A_only,denom,asPct), fmtCount(C.B_only,denom,asPct), ...
                  fmtCount(C.C_only,denom,asPct), fmtCount(C.D_only,denom,asPct), ...
                  fmtCount(C.AB,denom,asPct),     fmtCount(C.AC,denom,asPct),     fmtCount(C.AD,denom,asPct), ...
                  fmtCount(C.BC,denom,asPct),     fmtCount(C.BD,denom,asPct),     fmtCount(C.CD,denom,asPct), ...
                  fmtCount(C.ABC,denom,asPct),    fmtCount(C.ABD,denom,asPct),    fmtCount(C.ACD,denom,asPct), ...
                  fmtCount(C.BCD,denom,asPct),    fmtCount(C.ABCD,denom,asPct) };
          col = parula(4);
          venn(4, 'sets', labsK, 'labels', lbl, 'colors', col, ...
                  'alpha', 0.5, 'edgeC', [1 1 1], 'edgeW', 3);
      end
end
R.opts=opt;
end
% ---------------- helpers ----------------
function upsetPlotGeneric(sets, labels)
% Works for 2–4 sets
K = size(sets,2);
if K < 2 || K > 4
    title(sprintf('UpSet: supports 2–4 sets (got %d)', K)); return
end
bits = packBitsN(sets);       % packs to 2^K-1 combos (1..(2^K-1))
u = (1:(2^K - 1))';           % exclude 0 (no-set)
cnt = arrayfun(@(b) sum(bits==b), u);
[cnt, ord] = sort(cnt, 'descend');
u = u(ord); u = u(cnt>0); cnt = cnt(cnt>0);
M = numel(u);

subplot(2,1,1); bar(1:M, cnt, 0.85); ylabel('# cells'); set(gca,'XTick',[]);
title(sprintf('UpSet (%d sets)', K)); box off

subplot(2,1,2); cla; hold on
for i=1:M
    bvec = bitget(u(i), 1:K);  % little-endian
    for s=1:K
        if bvec(s)
            plot(i, K-s+1, 'ko', 'MarkerFaceColor','k');
        else
            plot(i, K-s+1, 'o',  'Color',[0.8 0.8 0.8]);
        end
    end
    idx = find(bvec);
    if ~isempty(idx)
        y = K - idx + 1;
        plot([i i], [min(y) max(y)], '-', 'Color','k', 'LineWidth',1);
    end
end
set(gca,'YTick',1:K,'YTickLabel',fliplr(labels));
xlim([0.5 M+0.5]); ylim([0.5 K+0.5]); xlabel('Intersection'); box off
end

function bits = packBitsN(sets)
% sets: [N x K], K<=8 (here we use up to 4)
K = size(sets,2);
s = uint8(sets>0);
bits = zeros(size(sets,1),1,'uint8');
for k=1:K
    bits = bitor(bits, bitshift(s(:,k), k-1));
end
end

function counts=intervalCountsPerRow(st,win)
if isempty(st)||isempty(win),counts=zeros(size(win,1),1);return,end
N=size(win,1);counts=zeros(N,1);
for i=1:N,counts(i)=sum(st>=win(i,1)&st<win(i,2));end
end

function p=poisson_two_sided(k,lambda)
if lambda<=0, p=double(k>0); return, end
Fc=poisscdf(k,lambda);Fg=1-poisscdf(k-1,lambda);
p=2*min(Fc,Fg);p=min(1,p);
end

function [pcMask,ok]=fetchPlaceMask(rat,dayStr,spec)
ok=false;pcMask=[];
fld=[spec.Prefix dayStr];
if isfield(rat,spec.Struct)&&isfield(rat.(spec.Struct),fld)
    MI=rat.(spec.Struct).(fld);
    ok=true;
else,return,end
if size(MI,2)<spec.Col,ok=false;return,end
v=MI(:,spec.Col);
switch spec.Comparator
    case '>',pcMask=v>spec.Thresh;
    case '<',pcMask=v<spec.Thresh;
end
pcMask=pcMask(:)>0;
end

function s=fmtCount(n,denom,asPercent)
if asPercent, s=sprintf('%.1f%%',100*n/max(1,denom));
else, s=sprintf('%d',n);
end
end

function C=computeExclusiveCounts(sets)
K=size(sets,2);
if K==2
    A=sets(:,1);B=sets(:,2);AB=A&B;
    C=struct('A_only',sum(A&~B),'B_only',sum(B&~A),'AB',sum(AB));
elseif K==3
    A=sets(:,1);B=sets(:,2);Cc=sets(:,3);
    ABC=A&B&Cc;
    AB=(A&B)&~Cc;AC=(A&Cc)&~B;BC=(B&Cc)&~A;
    A1=A&~B&~Cc;B1=B&~A&~Cc;C1=Cc&~A&~B;
    C=struct('A_only',sum(A1),'B_only',sum(B1),'C_only',sum(C1), ...
             'AB',sum(AB),'AC',sum(AC),'BC',sum(BC),'ABC',sum(ABC));
end
end

function C=computeExclusiveCounts4(sets)
A=sets(:,1)>0;B=sets(:,2)>0;C1=sets(:,3)>0;D=sets(:,4)>0;
C.A_only=sum(A&~B&~C1&~D);
C.B_only=sum(B&~A&~C1&~D);
C.C_only=sum(C1&~A&~B&~D);
C.D_only=sum(D&~A&~B&~C1);
C.AB=sum(A&B&~C1&~D);
C.AC=sum(A&C1&~B&~D);
C.AD=sum(A&D&~B&~C1);
C.BC=sum(B&C1&~A&~D);
C.BD=sum(B&D&~A&~C1);
C.CD=sum(C1&D&~A&~B);
C.ABC=sum(A&B&C1&~D);
C.ABD=sum(A&B&D&~C1);
C.ACD=sum(A&C1&D&~B);
C.BCD=sum(B&C1&D&~A);
C.ABCD=sum(A&B&C1&D);
end

function C=complementIntervals(full,excl)
a=full(1);b=full(2);E=excl;E=E(~any(isnan(E),2),:);
E=sortrows(E,1);E(:,1)=max(E(:,1),a);E(:,2)=min(E(:,2),b);
E=E(E(:,1)<E(:,2),:);C=[];t=a;
for i=1:size(E,1)
    if E(i,1)>t,C=[C;t E(i,1)];end
    t=max(t,E(i,2));
end
if t<b,C=[C;t b];end
end

function upsetPlot4(sets,labels)
bits=packBits4(sets);
u=(1:15)';cnt=arrayfun(@(b)sum(bits==b),u);
[cnt,ord]=sort(cnt,'descend');u=u(ord);u=u(cnt>0);cnt=cnt(cnt>0);
K=numel(u);
subplot(2,1,1);bar(1:K,cnt,0.8);ylabel('#cells');xlim([0.5 K+0.5]);set(gca,'XTick',[]);
subplot(2,1,2);hold on
for i=1:K
    b=dec2bin(u(i),4)=='1';
    for s=1:4
        if b(5-s),plot(i,s,'ko','MarkerFaceColor','k');
        else,plot(i,s,'o','Color',[0.8 0.8 0.8]);end
    end
end
set(gca,'YTick',1:4,'YTickLabel',fliplr(labels));xlim([0.5 K+0.5]);ylim([0.5 4.5]);
xlabel('Intersection');box off
end

function bits=packBits4(sets)
s=uint8(sets);
bits=bitor(bitor(bitor(bitshift(s(:,1),0),bitshift(s(:,2),1)),bitshift(s(:,3),2)),bitshift(s(:,4),3));
end

function [mask,labsOut]=selectEpochSubset(allLabs,subset,maxK)
labs=cellstr(allLabs);want=cellstr(subset);mask=false(1,numel(labs));
for i=1:numel(want)
    j=find(strcmpi(labs,want{i}),1);
    if ~isempty(j),mask(j)=true;end
end
if nnz(mask)>maxK,mask(find(mask,maxK,'first'))=true;end
labsOut=labs(mask);
end
