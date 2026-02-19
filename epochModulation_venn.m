function R = epochModulation_venn(ratNames, varargin)
% epochModulation_venn (Option 2)
%
% Bars: fraction of KEPT cells per epoch (rat-equalized optional), plus Place bar.
% Epoch overlaps/Venn: computed ONLY within full-window task-modulated cells, using effect size (zMat).
% Task (any epoch) for overlap diagrams: full-window task-modulated (from rat.mod.mod_<dayStr>).
%
% !!!!!!!! Required upstream: plotProportionModulated(...,'SaveMod',true) to populate rat.mod.mod_<dayStr>.
%
% External deps (not included here):
%   - autoDateList(rat)
%   - venn(...) plotting function (same one you currently call)

% ---------------- parameters ----------------
p = inputParser;

% Full-window task gate
p.addParameter('TaskAlpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);

% Epoch membership rule within task-modulated cells (effect-size thresholds)
p.addParameter('RelThresh',0.50,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('AbsZMin',1.0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);

% Plot options
p.addParameter('PlotMode','venn',@(s)ischar(s)||isstring(s));
p.addParameter('CollapsePreIntoOne',true,@islogical);
p.addParameter('VennEpochs',{'CS','Trace','US','Post'});
p.addParameter('VennPercent',true,@islogical);
p.addParameter('PercentBase','kept',@(s)any(strcmpi(s,{'kept','union'})));
p.addParameter('EqualWeightRats',true,@islogical);
p.addParameter('EqualizeBars',false,@islogical);
p.addParameter('Verbose',true,@islogical);
p.addParameter('Tail','right',@(s) any(strcmpi(s,{'two','right','left'})));

% Place-cell options
p.addParameter('PlaceOnly',false,@islogical);
p.addParameter('PlaceSpec',struct('Struct','MI_noCSUS15_shuff','Prefix','MI_', ...
    'Col',3,'Thresh',0.95,'Comparator','>'));

p.parse(varargin{:});
opt = p.Results;

Tail = lower(p.Results.Tail);

% ---------------- fixed epoch definitions ----------------
labels = {'CS','Trace','US','Post'};
epochs = [0.00 0.25; 0.25 0.75; 0.75 0.85; 0.85 2.00];   % 4 epochs in s post-CS
trialWin = [0 6.0]; % exclusion around CS for baseline pool

K = numel(labels);
epochDur = diff(epochs,1,2);

nR = numel(ratNames);

% ---------------- accumulators ----------------
% Task-modulated cells only (for epoch overlaps / venn membership)
ALL_sets_task   = [];    % [Ntask x 4] epoch memberships
ALL_dir_task    = [];    % [Ntask x 4] sign of z
ALL_z_task      = [];    % [Ntask x 4] zMat
ALL_place_task  = [];    % [Ntask x 1] place membership (NaN if missing)
ALL_rat_task    = [];    % [Ntask x 1] rat id (1..nR)
ALL_dom_task    = [];    % [Ntask x 1] dominant epoch (1..4)

% Kept-cell bookkeeping (for denominators and Place bar, and Task-vs-Place)
ALL_rat_kept      = [];  % [Nkept x 1]
ALL_place_kept    = [];  % [Nkept x 1]
ALL_taskflag_kept = [];  % [Nkept x 1] full-window task-modulated (tail-aware)

% NEW: kept-aligned pre/post membership among task-modulated cells
ALL_pre_kept      = [];  % [Nkept x 1] in any of CS/Trace/US epochs (within task-modulated)
ALL_post_kept     = [];  % [Nkept x 1] in Post epoch (within task-modulated)
% ---------------- main loop ----------------
for r = 1:nR
    rat = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates,rat.An),1);
    if isempty(idx), idx = numel(dates); end
    days  = dates(max(1,idx-2):idx);

    for d = 1:numel(days)
        dayStr = days{d};

        spk     = rat.Ca_peaks.(sprintf('CA_peaks_%s',dayStr));
        csTimes = rat.CS_times.(sprintf('CS_%s',dayStr));
        if isempty(csTimes), continue; end

        mask0 = rat.ratemask.(sprintf('ratemask_%s',dayStr))==1;

        % Place membership (track even if PlaceOnly=false)
        placeThisDay = nan(size(mask0));
        [pcMask, okPC] = fetchPlaceMask(rat, dayStr, opt.PlaceSpec);
        if okPC
            placeThisDay = pcMask(:);
        else
            if opt.PlaceOnly
                warning('[%s %s] PlaceOnly requested but MI missing.', ratNames{r}, dayStr);
            end
        end

        % Apply PlaceOnly restriction if requested
        mask = mask0;
        if opt.PlaceOnly && okPC
            mask = mask & pcMask;
        end

        % Load full-window task p-values saved by plotProportionModulated (tail-aware)
        taskP = nan(size(mask0));
        fldMod = sprintf('mod_%s', dayStr);

        if isfield(rat,'mod') && isfield(rat.mod, fldMod)
            modObj = rat.mod.(fldMod);

            % New format: struct with p_two / p_right / p_left
            if isstruct(modObj) && isfield(modObj,'p_two')
                switch Tail
                    case 'two'
                        taskP = modObj.p_two;
                    case 'right'
                        taskP = modObj.p_right;
                    case 'left'
                        taskP = modObj.p_left;
                end

            % Backward compatible: numeric vector (assume two-sided)
            else
                taskP = modObj;
            end
        else
            warning('[%s %s] Missing rat.mod.%s (run plotProportionModulated with SaveMod=true).', ...
                ratNames{r}, dayStr, fldMod);
        end

        taskSig0 = isfinite(taskP) & (taskP <= opt.TaskAlpha);

        taskSig = taskSig0;
        if opt.PlaceOnly && okPC
            taskSig = taskSig & pcMask;
        end

        % Build baseline pool = complement of [CS + trialWin] windows
        sessMin = min(csTimes)-60;
        sessMax = max(csTimes)+60;
        excl = [csTimes(:)+trialWin(1), csTimes(:)+trialWin(2)];
        pool = complementIntervals([sessMin sessMax], excl);
        poolDur = sum(max(0,pool(:,2)-pool(:,1)));
        if poolDur<=0, continue; end

        % Precompute absolute epoch windows
        nTrials = numel(csTimes);
        epAbs = cell(K,1);
        for e = 1:K
            epAbs{e} = [csTimes(:)+epochs(e,1), csTimes(:)+epochs(e,2)];
        end

        nCells = size(spk,1);
        zMat   = nan(nCells,K);
        dirsgn = zeros(nCells,K,'int8');
        keep   = false(nCells,1);
        preKeep  = false(nCells,1);
        postKeep = false(nCells,1);

        % Compute zMat for all kept cells
        for c = 1:nCells
            if ~mask(c), continue; end
            st = spk(c,:);
            st = st(~isnan(st) & st>0);
            if isempty(st), continue; end
            keep(c) = true;

            poolCounts = intervalCountsPerRow(st, pool);
            r_pool     = sum(poolCounts)/poolDur;

            for e = 1:K
                counts = intervalCountsPerRow(st, epAbs{e});
                obsK   = sum(counts);
                expT   = epochDur(e)*nTrials;
                lambda = r_pool*expT;

                if lambda > 0
                    zMat(c,e) = (obsK - lambda) ./ sqrt(lambda);
                    dirsgn(c,e) = int8(sign(zMat(c,e)));
                else
                    zMat(c,e) = NaN;
                    dirsgn(c,e) = int8(0);
                end
            end
        end

        % Task-modulated subset for epoch membership
      taskKeep = keep & taskSig;
      idxTask  = find(taskKeep);

      if opt.Verbose
          fprintf('[%s %s] kept=%d | task(fullwin)=%d\n', ratNames{r}, dayStr, nnz(keep), numel(idxTask));
      end

      % Default: pre/post are all-false unless we populate below
      % (you already initialized preKeep/postKeep = false(nCells,1) earlier)
      if ~isempty(idxTask)

          zTask   = zMat(idxTask,:);
          dirTask = dirsgn(idxTask,:);

          setsTask = false(numel(idxTask), K);
          domEpoch = zeros(numel(idxTask),1);

          for ii = 1:numel(idxTask)
              z = zTask(ii,:);
              if all(~isfinite(z)), continue; end
              [zmax, dom] = max(abs(z));
              domEpoch(ii) = dom;
              if ~isfinite(zmax) || zmax<=0, continue; end

              thr = max(opt.AbsZMin, opt.RelThresh*zmax);
              setsTask(ii,:) = abs(z) >= thr;
              setsTask(ii,dom) = true;
          end

          % pre = any of CS/Trace/US (cols 1:3), post = Post (col 4)
          preTask  = any(setsTask(:,1:3), 2);
          postTask = setsTask(:,4) > 0;

          preKeep(idxTask)  = preTask;
          postKeep(idxTask) = postTask;

          % Accumulate TASK-only arrays
          ALL_sets_task  = [ALL_sets_task;  setsTask];                 %#ok<AGROW>
          ALL_dir_task   = [ALL_dir_task;   dirTask];                  %#ok<AGROW>
          ALL_z_task     = [ALL_z_task;     zTask];                    %#ok<AGROW>
          ALL_place_task = [ALL_place_task; placeThisDay(idxTask)];    %#ok<AGROW>
          ALL_rat_task   = [ALL_rat_task;   r*ones(numel(idxTask),1)]; %#ok<AGROW>
          ALL_dom_task   = [ALL_dom_task;   domEpoch];                 %#ok<AGROW>

          % Save dominant epoch per cell for this day (optional)
          epochVec = nan(nCells,1);
          epochVec(keep & ~taskSig) = 0;
          epochVec(idxTask) = domEpoch;

          if ~isfield(rat,'epoch') || ~isstruct(rat.epoch)
              rat.epoch = struct();
          end
          rat.epoch.(sprintf('epoch_%s', dayStr)) = epochVec;
          assignin('base', ratNames{r}, rat);
      end

      % NOW append kept bookkeeping (after preKeep/postKeep are correct)
      keptIdx = find(keep);
      if ~isempty(keptIdx)
          ALL_rat_kept      = [ALL_rat_kept;      r*ones(numel(keptIdx),1)]; %#ok<AGROW>
          ALL_place_kept    = [ALL_place_kept;    placeThisDay(keptIdx)];    %#ok<AGROW>
          ALL_taskflag_kept = [ALL_taskflag_kept; taskSig(keptIdx)];         %#ok<AGROW>
          ALL_pre_kept      = [ALL_pre_kept;      preKeep(keptIdx)];         %#ok<AGROW>
          ALL_post_kept     = [ALL_post_kept;     postKeep(keptIdx)];        %#ok<AGROW>
      end

      % If idxTask was empty, we already appended kept with pre/post all-false and move on
      if isempty(idxTask)
          continue
      end
    end
end

% ---------------- guardrails ----------------
if isempty(ALL_rat_kept)
    warning('No kept cells found.');
    R.opts = opt;
    return
end

% ---------------- build bar data (Option 2) ----------------
% Numerator: epoch membership among TASK cells
% Denominator: KEPT cells
if opt.EqualizeBars
    stackCounts  = ratEqualEpochStackCounts_option2(ALL_sets_task, ALL_rat_task, ALL_rat_kept, K);
    placeFrac    = ratEqualPlaceFrac(ALL_place_kept, ALL_rat_kept);
else
    stackCounts  = pooledEpochStackCounts_option2(ALL_sets_task, ALL_rat_kept, K);
    placeFrac    = pooledPlaceFrac(ALL_place_kept);
end

% Append Place as its own bar row, drawn as "Only this epoch"
stackCounts_with_place = [stackCounts; zeros(1,K)];
stackCounts_with_place(end,1) = placeFrac;

% Time-normalize (Place gets ~0 since duration=Inf)
epochDur_with_place = [epochDur; Inf];
Ynorm = bsxfun(@rdivide, stackCounts_with_place, epochDur_with_place);

printNumEpochsSummary(ALL_sets_task);

% ---------------- plotting ----------------
switch lower(string(opt.PlotMode))
    case "upset"
        figure('Color','w');
        upsetPlotGeneric(ALL_sets_task, labels);

    otherwise % "venn"
        % ---- stacked bars: fraction of kept ----
        figure('Color','w');
        subplot(1,2,1);
        bh1 = bar(1:(K+1), stackCounts_with_place, 'stacked');
        ylabel(ternary(opt.EqualizeBars,'Fraction of kept (rat-equalized)','Fraction of kept'));
        set(gca,'XTick',1:(K+1),'XTickLabel',[labels {'Place'}]);
        box off;

        legStr = cell(1,K);
        for k = 1:K
            if k == 1
                legStr{k} = 'Only this epoch';
            else
                legStr{k} = sprintf('+%d other epoch%s', k-1, ternary(k-1==1,'','s'));
            end
        end
        legend(bh1, legStr, 'Location','bestoutside');
        title('Epoch membership within task-modulated cells (normalized to kept)');

        subplot(1,2,2);
        bh2 = bar(1:(K+1), Ynorm, 'stacked');
        ylabel(ternary(opt.EqualizeBars,'(Fraction of kept)/s (rat-equalized)','(Fraction of kept)/s'));
        set(gca,'XTick',1:(K+1),'XTickLabel',[labels {'Place'}]);
        box off;
        legend(bh2, legStr, 'Location','bestoutside');
        title('Epoch membership density (normalized by epoch duration)');

        % ---- Venn (task cells only) ----
        if ~isempty(ALL_sets_task)
            plotSets   = ALL_sets_task;
            plotLabels = {'CS','Trace','US','Post'};

            if opt.CollapsePreIntoOne
                plotSets   = [ any(plotSets(:,1:3),2), plotSets(:,4) ];
                plotLabels = {'CS/Trace/US','Post'};
            end

            if ~opt.CollapsePreIntoOne
                [epMask, labsK] = selectEpochSubset(plotLabels, opt.VennEpochs);
                setsK = plotSets(:,epMask);
            else
                setsK = plotSets;
                labsK = plotLabels;
            end

            kSel = size(setsK,2);
            if kSel >= 2 && kSel <= 4
                if opt.EqualWeightRats
                    if kSel == 2
                        F = ratEqualExclusiveFrac2(setsK, ALL_rat_task, opt.PercentBase);
                        lbl = {sprintf('%.1f%%',100*F.A_only),sprintf('%.1f%%',100*F.B_only),sprintf('%.1f%%',100*F.AB)};
                    elseif kSel == 3
                        F = ratEqualExclusiveFrac3(setsK, ALL_rat_task, opt.PercentBase);
                        lbl = {sprintf('%.1f%%',100*F.A_only),sprintf('%.1f%%',100*F.B_only),sprintf('%.1f%%',100*F.C_only), ...
                               sprintf('%.1f%%',100*F.AB),sprintf('%.1f%%',100*F.AC),sprintf('%.1f%%',100*F.BC), ...
                               sprintf('%.1f%%',100*F.ABC)};
                    else
                        F = ratEqualExclusiveFrac4(setsK, ALL_rat_task, opt.PercentBase);
                        lbl = {sprintf('%.1f%%',100*F.A_only),sprintf('%.1f%%',100*F.B_only),sprintf('%.1f%%',100*F.C_only),sprintf('%.1f%%',100*F.D_only), ...
                               sprintf('%.1f%%',100*F.AB),sprintf('%.1f%%',100*F.AC),sprintf('%.1f%%',100*F.AD), ...
                               sprintf('%.1f%%',100*F.BC),sprintf('%.1f%%',100*F.BD),sprintf('%.1f%%',100*F.CD), ...
                               sprintf('%.1f%%',100*F.ABC),sprintf('%.1f%%',100*F.ABD),sprintf('%.1f%%',100*F.ACD),sprintf('%.1f%%',100*F.BCD),sprintf('%.1f%%',100*F.ABCD)};
                    end
                else
                    if strcmpi(opt.PercentBase,'union')
                        denom = max(1, sum(any(setsK,2)));
                    else
                        denom = size(setsK,1);
                    end
                    asPct = opt.VennPercent;

                    if kSel == 2
                        C = computeExclusiveCounts(setsK);
                        lbl = {fmtCount(C.A_only,denom,asPct), fmtCount(C.B_only,denom,asPct), fmtCount(C.AB,denom,asPct)};
                    elseif kSel == 3
                        C = computeExclusiveCounts(setsK);
                        lbl = {fmtCount(C.A_only,denom,asPct), fmtCount(C.B_only,denom,asPct), fmtCount(C.C_only,denom,asPct), ...
                               fmtCount(C.AB,denom,asPct), fmtCount(C.AC,denom,asPct), fmtCount(C.BC,denom,asPct), fmtCount(C.ABC,denom,asPct)};
                    else
                        C = computeExclusiveCounts4(setsK);
                        lbl = {fmtCount(C.A_only,denom,asPct), fmtCount(C.B_only,denom,asPct), fmtCount(C.C_only,denom,asPct), fmtCount(C.D_only,denom,asPct), ...
                               fmtCount(C.AB,denom,asPct), fmtCount(C.AC,denom,asPct), fmtCount(C.AD,denom,asPct), ...
                               fmtCount(C.BC,denom,asPct), fmtCount(C.BD,denom,asPct), fmtCount(C.CD,denom,asPct), ...
                               fmtCount(C.ABC,denom,asPct), fmtCount(C.ABD,denom,asPct), fmtCount(C.ACD,denom,asPct), fmtCount(C.BCD,denom,asPct), fmtCount(C.ABCD,denom,asPct)};
                    end
                end

                figure('Color','w'); hold on;
                if kSel == 2
                    col = parula(5); col = col([2,5],:);
                    venn(2, 'sets', labsK, 'labels', lbl, 'colors', col, 'alpha', 0.5, 'edgeC', [1 1 1], 'edgeW', 3);
                elseif kSel == 3
                    col = parula(3);
                    venn(3, 'sets', labsK, 'labels', lbl, 'colors', col, 'alpha', 0.5, 'edgeC', [1 1 1], 'edgeW', 3);
                else
                    col = parula(4);
                    venn(4, 'sets', labsK, 'labels', lbl, 'colors', col, 'alpha', 0.5, 'edgeC', [1 1 1], 'edgeW', 3);
                end
            end
        end

        % ---- Proportional Task-vs-Place (kept denominator; Task = full-window) ----
        taskMask = ALL_taskflag_kept;
        size(ALL_taskflag_kept)

        placeMask = false(size(ALL_place_kept));
        ok = isfinite(ALL_place_kept);
        placeMask(ok) = ALL_place_kept(ok) > 0;

        valid = true(size(taskMask)); % keep all kept cells
        A = taskMask(valid);
        B = placeMask(valid);
        ratIDs = ALL_rat_kept(valid);

        % ---- Proportional Pre/Post/Place (kept denominator) ----
        % Pre/Post are defined among task-modulated cells via epoch membership,
        % then stored aligned to kept cells.
        preMask  = ALL_pre_kept(:)  > 0;
        postMask = ALL_post_kept(:) > 0;

        placeMask = false(size(ALL_place_kept));
        ok = isfinite(ALL_place_kept);
        placeMask(ok) = ALL_place_kept(ok) > 0;

        ratIDs = ALL_rat_kept(:);

        figure('Color','w');
        drawProportionalVenn3(preMask, postMask, placeMask, ratIDs, ...
            {'Tone/Trace/Shock','Post','Place'}, opt.PercentBase, opt.EqualWeightRats);

            taskMask = ALL_taskflag_kept(:) > 0;

figure('Color','w');
drawProportionalVenn2(taskMask, placeMask, ratIDs, {'Task (full window)','Place'}, ...
    opt.PercentBase, opt.EqualWeightRats);
end

% ---------------- output struct ----------------
R.opts                 = opt;
R.labels               = labels;
R.epochs               = epochs;
R.epochDur             = epochDur;
R.ALL_sets_task        = ALL_sets_task;
R.ALL_dir_task         = ALL_dir_task;
R.ALL_z_task           = ALL_z_task;
R.ALL_place_task       = ALL_place_task;
R.ALL_rat_task         = ALL_rat_task;
R.ALL_dom_task         = ALL_dom_task;
R.ALL_rat_kept         = ALL_rat_kept;
R.ALL_place_kept       = ALL_place_kept;
R.ALL_taskflag_kept    = ALL_taskflag_kept;
R.stackCounts_with_place = stackCounts_with_place;
R.Ynorm                = Ynorm;
R.ALL_pre_kept  = ALL_pre_kept;
R.ALL_post_kept = ALL_post_kept;

end

% =====================================================================
% Helpers (epochModulation_venn)
% =====================================================================

function y = ternary(cond,a,b)
if cond, y = a; else, y = b; end
end

function [mask, labs] = selectEpochSubset(allLabels, wanted)
mask = ismember(allLabels, wanted);
labs = allLabels(mask);
if nnz(mask) < 2
    mask = true(size(allLabels));
    labs = allLabels;
end
end

function counts = intervalCountsPerRow(st, win)
if isempty(st) || isempty(win)
    counts = zeros(size(win,1),1);
    return
end
N = size(win,1);
counts = zeros(N,1);
for i = 1:N
    counts(i) = sum(st>=win(i,1) & st<win(i,2));
end
end

function C = complementIntervals(full, excl)
a = full(1); b = full(2);
E = excl;
E = E(~any(isnan(E),2),:);
E = sortrows(E,1);
E(:,1) = max(E(:,1),a);
E(:,2) = min(E(:,2),b);
E = E(E(:,1)<E(:,2),:);

C = [];
t = a;
for i = 1:size(E,1)
    if E(i,1) > t
        C = [C; t E(i,1)]; %#ok<AGROW>
    end
    t = max(t, E(i,2));
end
if t < b
    C = [C; t b]; %#ok<AGROW>
end
end

function [pcMask, ok] = fetchPlaceMask(rat, dayStr, spec)
ok = false; pcMask = [];
fld = [spec.Prefix dayStr];

if isfield(rat, spec.Struct) && isfield(rat.(spec.Struct), fld)
    MI = rat.(spec.Struct).(fld);
    ok = true;
else
    return
end

if size(MI,2) < spec.Col
    ok = false;
    return
end

v = MI(:,spec.Col);
switch spec.Comparator
    case '>'
        pcMask = v > spec.Thresh;
    case '<'
        pcMask = v < spec.Thresh;
    otherwise
        pcMask = v > spec.Thresh;
end
pcMask = pcMask(:) > 0;
end

function s = fmtCount(n, denom, asPercent)
if asPercent
    s = sprintf('%.1f%%', 100*n/max(1,denom));
else
    s = sprintf('%d', n);
end
end

function C = computeExclusiveCounts(sets)
K = size(sets,2);
if K==2
    A = sets(:,1); B = sets(:,2); AB = A & B;
    C = struct('A_only',sum(A&~B), 'B_only',sum(B&~A), 'AB',sum(AB));
elseif K==3
    A = sets(:,1); B = sets(:,2); Cc = sets(:,3);
    ABC = A & B & Cc;
    AB  = (A & B)  & ~Cc;
    AC  = (A & Cc) & ~B;
    BC  = (B & Cc) & ~A;
    A1  = A  & ~B  & ~Cc;
    B1  = B  & ~A  & ~Cc;
    C1  = Cc & ~A  & ~B;
    C = struct('A_only',sum(A1),'B_only',sum(B1),'C_only',sum(C1), ...
               'AB',sum(AB),'AC',sum(AC),'BC',sum(BC),'ABC',sum(ABC));
else
    error('computeExclusiveCounts supports 2 or 3 sets.');
end
end

function C = computeExclusiveCounts4(sets)
A = sets(:,1)>0; B = sets(:,2)>0; C1 = sets(:,3)>0; D = sets(:,4)>0;

C.A_only = sum(A & ~B & ~C1 & ~D);
C.B_only = sum(B & ~A & ~C1 & ~D);
C.C_only = sum(C1 & ~A & ~B & ~D);
C.D_only = sum(D & ~A & ~B & ~C1);

C.AB = sum(A & B & ~C1 & ~D);
C.AC = sum(A & C1 & ~B & ~D);
C.AD = sum(A & D & ~B & ~C1);
C.BC = sum(B & C1 & ~A & ~D);
C.BD = sum(B & D & ~A & ~C1);
C.CD = sum(C1 & D & ~A & ~B);

C.ABC  = sum(A & B & C1 & ~D);
C.ABD  = sum(A & B & D  & ~C1);
C.ACD  = sum(A & C1 & D & ~B);
C.BCD  = sum(B & C1 & D & ~A);
C.ABCD = sum(A & B & C1 & D);
end

% ==========================================================
% Option-2 bar helpers:
% Numerator from TASK epoch-membership sets; denominator = KEPT cells
% ==========================================================
function S = pooledEpochStackCounts_option2(setsTask, ratKept, K)
% setsTask: [Ntask x K] logical epoch membership among task-modulated cells
% ratKept : [Nkept x 1] (only used for denom length)
S = zeros(K,K);
denom = max(1, numel(ratKept));
for e = 1:K
    idx = find(setsTask(:,e));
    if isempty(idx), continue; end
    othersCount = sum(setsTask(idx,:),2) - 1; % 0..K-1
    for k = 0:(K-1)
        S(e,k+1) = sum(othersCount == k) / denom;
    end
end
end

function S = ratEqualEpochStackCounts_option2(setsTask, ratTask, ratKept, K)
% ratTask: [Ntask x 1] rat ID for each row of setsTask
% ratKept: [Nkept x 1] rat ID for each kept cell
rats = unique(ratKept(:)); rats = rats(isfinite(rats));
nR = numel(rats);

Srat = zeros(nR, K, K); % [rat x epoch x overlapBin] as fractions of kept

for i = 1:nR
    rr = rats(i);
    denom = max(1, sum(ratKept==rr));   % kept cells in that rat
    idxR_task = (ratTask==rr);
    setsR = setsTask(idxR_task,:);

    for e = 1:K
        idx = find(setsR(:,e));
        if isempty(idx), continue; end
        othersCount = sum(setsR(idx,:),2) - 1;
        for k = 0:(K-1)
            Srat(i,e,k+1) = sum(othersCount == k) / denom;
        end
    end
end

S = squeeze(mean(Srat,1,'omitnan')); % [K x K]
end

function f = pooledPlaceFrac(placeKept)
ok = isfinite(placeKept);
if ~any(ok), f = 0; return; end
f = mean(placeKept(ok) > 0);
end

function f = ratEqualPlaceFrac(placeKept, ratKept)
rats = unique(ratKept(:)); rats = rats(isfinite(rats));
vals = nan(numel(rats),1);
for i = 1:numel(rats)
    idx = (ratKept==rats(i));
    pk  = placeKept(idx);
    ok  = isfinite(pk);
    if any(ok)
        vals(i) = mean(pk(ok) > 0);
    else
        vals(i) = NaN;
    end
end
f = mean(vals,'omitnan');
end

% ==========================================================
% Rat-equalized exclusive fractions for Venn labels
% ==========================================================
function F = ratEqualExclusiveFrac2(sets2, ratIDs, percentBase)
A = sets2(:,1)>0; B = sets2(:,2)>0;
rats = unique(ratIDs(:)); rats = rats(isfinite(rats));
vals = nan(numel(rats),3); % A_only, B_only, AB
for i = 1:numel(rats)
    idx = (ratIDs==rats(i));
    Ai = A(idx); Bi = B(idx);
    if strcmpi(percentBase,'union')
        denom = max(1, sum(Ai|Bi));
    else
        denom = max(1, numel(Ai));
    end
    vals(i,1) = sum(Ai & ~Bi)/denom;
    vals(i,2) = sum(Bi & ~Ai)/denom;
    vals(i,3) = sum(Ai &  Bi)/denom;
end
F.A_only = mean(vals(:,1),'omitnan');
F.B_only = mean(vals(:,2),'omitnan');
F.AB     = mean(vals(:,3),'omitnan');
end

function F = ratEqualExclusiveFrac3(sets3, ratIDs, percentBase)
A = sets3(:,1)>0; B = sets3(:,2)>0; C = sets3(:,3)>0;
rats = unique(ratIDs(:)); rats = rats(isfinite(rats));
vals = nan(numel(rats),7);
for i = 1:numel(rats)
    idx = (ratIDs==rats(i));
    Ai=A(idx); Bi=B(idx); Ci=C(idx);
    if strcmpi(percentBase,'union')
        denom = max(1, sum(Ai|Bi|Ci));
    else
        denom = max(1, numel(Ai));
    end
    vals(i,1) = sum(Ai & ~Bi & ~Ci)/denom;
    vals(i,2) = sum(Bi & ~Ai & ~Ci)/denom;
    vals(i,3) = sum(Ci & ~Ai & ~Bi)/denom;
    vals(i,4) = sum(Ai &  Bi & ~Ci)/denom;
    vals(i,5) = sum(Ai &  Ci & ~Bi)/denom;
    vals(i,6) = sum(Bi &  Ci & ~Ai)/denom;
    vals(i,7) = sum(Ai &  Bi &  Ci)/denom;
end
F.A_only = mean(vals(:,1),'omitnan');
F.B_only = mean(vals(:,2),'omitnan');
F.C_only = mean(vals(:,3),'omitnan');
F.AB     = mean(vals(:,4),'omitnan');
F.AC     = mean(vals(:,5),'omitnan');
F.BC     = mean(vals(:,6),'omitnan');
F.ABC    = mean(vals(:,7),'omitnan');
end

function F = ratEqualExclusiveFrac4(sets4, ratIDs, percentBase)
A=sets4(:,1)>0; B=sets4(:,2)>0; C=sets4(:,3)>0; D=sets4(:,4)>0;
rats = unique(ratIDs(:)); rats = rats(isfinite(rats));
vals = nan(numel(rats),15);

for i = 1:numel(rats)
    idx = (ratIDs==rats(i));
    Ai=A(idx); Bi=B(idx); Ci=C(idx); Di=D(idx);
    if strcmpi(percentBase,'union')
        denom = max(1, sum(Ai|Bi|Ci|Di));
    else
        denom = max(1, numel(Ai));
    end

    vals(i, 1)=sum(Ai & ~Bi & ~Ci & ~Di)/denom;
    vals(i, 2)=sum(Bi & ~Ai & ~Ci & ~Di)/denom;
    vals(i, 3)=sum(Ci & ~Ai & ~Bi & ~Di)/denom;
    vals(i, 4)=sum(Di & ~Ai & ~Bi & ~Ci)/denom;

    vals(i, 5)=sum(Ai & Bi & ~Ci & ~Di)/denom;
    vals(i, 6)=sum(Ai & Ci & ~Bi & ~Di)/denom;
    vals(i, 7)=sum(Ai & Di & ~Bi & ~Ci)/denom;
    vals(i, 8)=sum(Bi & Ci & ~Ai & ~Di)/denom;
    vals(i, 9)=sum(Bi & Di & ~Ai & ~Ci)/denom;
    vals(i,10)=sum(Ci & Di & ~Ai & ~Bi)/denom;

    vals(i,11)=sum(Ai & Bi & Ci & ~Di)/denom;
    vals(i,12)=sum(Ai & Bi & Di & ~Ci)/denom;
    vals(i,13)=sum(Ai & Ci & Di & ~Bi)/denom;
    vals(i,14)=sum(Bi & Ci & Di & ~Ai)/denom;
    vals(i,15)=sum(Ai & Bi & Ci & Di)/denom;
end

F.A_only = mean(vals(:, 1),'omitnan');
F.B_only = mean(vals(:, 2),'omitnan');
F.C_only = mean(vals(:, 3),'omitnan');
F.D_only = mean(vals(:, 4),'omitnan');
F.AB     = mean(vals(:, 5),'omitnan');
F.AC     = mean(vals(:, 6),'omitnan');
F.AD     = mean(vals(:, 7),'omitnan');
F.BC     = mean(vals(:, 8),'omitnan');
F.BD     = mean(vals(:, 9),'omitnan');
F.CD     = mean(vals(:,10),'omitnan');
F.ABC    = mean(vals(:,11),'omitnan');
F.ABD    = mean(vals(:,12),'omitnan');
F.ACD    = mean(vals(:,13),'omitnan');
F.BCD    = mean(vals(:,14),'omitnan');
F.ABCD   = mean(vals(:,15),'omitnan');
end

% ==========================================================
% Proportional Task-vs-Place (2-circle)
% ==========================================================
function drawProportionalVenn2(A, B, ratIDs, names, percentBase, equalWeight)
A = A(:)>0; B = B(:)>0;

if equalWeight
    F = ratEqualExclusiveFrac2([A B], ratIDs, percentBase);
    fA  = F.A_only + F.AB;
    fB  = F.B_only + F.AB;
    fAB = F.AB;
else
    fA  = mean(A);
    fB  = mean(B);
    fAB = mean(A & B);
end

labA  = sprintf('%.1f%%',100*fA);
labB  = sprintf('%.1f%%',100*fB);
labAB = sprintf('%.1f%%',100*fAB);

rA = sqrt(max(1e-12,fA) / pi);
rB = sqrt(max(1e-12,fB) / pi);
d  = solveCircleDistanceForOverlap(rA, rB, fAB);

cla; hold on; axis equal; axis off;
cA = [0, 0];
cB = [d, 0];
drawCircle(cA, rA, 0.35);
drawCircle(cB, rB, 0.35);

title('Task vs Place (area-proportional)');
text(cA(1)-0.6*rA, cA(2), sprintf('%s\n%s', names{1}, labA), ...
    'HorizontalAlignment','center','FontSize',11);
text(cB(1)+0.6*rB, cB(2), sprintf('%s\n%s', names{2}, labB), ...
    'HorizontalAlignment','center','FontSize',11);

mid = (cA + cB)/2;
text(mid(1), mid(2), labAB, 'HorizontalAlignment','center', ...
    'FontSize',12,'FontWeight','bold');

pad = 0.25*max(rA,rB);
xlim([min(cA(1)-rA,cB(1)-rB)-pad, max(cA(1)+rA,cB(1)+rB)+pad]);
ylim([-max(rA,rB)-pad, max(rA,rB)+pad]);
end

function d = solveCircleDistanceForOverlap(r1, r2, targetArea)
maxOverlap = pi*min(r1,r2)^2;
targetArea = min(max(targetArea,0), maxOverlap);

if targetArea == 0
    d = r1 + r2 + 0.2*max([r1,r2]);
    return;
end
if abs(targetArea - maxOverlap) < 1e-12
    d = abs(r1 - r2);
    return;
end

lo = abs(r1 - r2) + 1e-12;
hi = (r1 + r2) - 1e-12;

f = @(x) circleOverlapArea(r1,r2,x) - targetArea;

if f(lo) < 0, d = lo; return; end
if f(hi) > 0, d = hi; return; end

try
    d = fzero(f, [lo hi]);
catch
    for it = 1:60
        mid = 0.5*(lo+hi);
        if f(mid) > 0, lo = mid; else, hi = mid; end
    end
    d = 0.5*(lo+hi);
end
end

function A = circleOverlapArea(r1,r2,d)
if d >= r1 + r2
    A = 0; return;
end
if d <= abs(r1 - r2)
    A = pi*min(r1,r2)^2; return;
end
alpha = 2*acos( (d^2 + r1^2 - r2^2) / (2*d*r1) );
beta  = 2*acos( (d^2 + r2^2 - r1^2) / (2*d*r2) );
A = 0.5*r1^2*(alpha - sin(alpha)) + 0.5*r2^2*(beta - sin(beta));
end

function drawCircle(center, radius, faceAlpha)
t = linspace(0,2*pi,400);
x = center(1) + radius*cos(t);
y = center(2) + radius*sin(t);
patch(x,y,1,'FaceAlpha',faceAlpha,'EdgeColor',[1 1 1],'LineWidth',2);
end

function drawProportionalVenn3(A, B, C, ratIDs, names, percentBase, equalWeight)
A=A(:)>0; B=B(:)>0; C=C(:)>0;

if equalWeight
    F = ratEqualExclusiveFrac3([A B C], ratIDs, percentBase);
    fA   = F.A_only + F.AB + F.AC + F.ABC;
    fB   = F.B_only + F.AB + F.BC + F.ABC;
    fC   = F.C_only + F.AC + F.BC + F.ABC;
    fAB  = F.AB + F.ABC;
    fAC  = F.AC + F.ABC;
    fBC  = F.BC + F.ABC;
    fABC = F.ABC;
else
    fA   = mean(A); fB = mean(B); fC = mean(C);
    fAB  = mean(A&B); fAC = mean(A&C); fBC = mean(B&C);
    fABC = mean(A&B&C);
end

labA   = sprintf('%.1f%%',100*fA);
labB   = sprintf('%.1f%%',100*fB);
labC   = sprintf('%.1f%%',100*fC);
labAB  = sprintf('%.1f%%',100*fAB);
labAC  = sprintf('%.1f%%',100*fAC);
labBC  = sprintf('%.1f%%',100*fBC);
labABC = sprintf('%.1f%%',100*fABC);

rA = sqrt(max(1e-12,fA)/pi);
rB = sqrt(max(1e-12,fB)/pi);
rC = sqrt(max(1e-12,fC)/pi);

dAB = solveCircleDistanceForOverlap(rA, rB, fAB);
dAC = solveCircleDistanceForOverlap(rA, rC, fAC);
dBC = solveCircleDistanceForOverlap(rB, rC, fBC);

cA = [0,0];
cB = [dAB,0];
[cC, ok] = thirdPointFromDistances(cA, cB, dAC, dBC);
if ~ok
    cC = [dAB/2, max([rA,rB,rC])*0.9 + 0.25];
end

cla; hold on; axis equal; axis off;
drawCircle(cA, rA, 0.30);
drawCircle(cB, rB, 0.30);
drawCircle(cC, rC, 0.30);

title('Pre/Post/Place (pairwise area-proportional)');
text(cA(1)-0.7*rA, cA(2)-0.1*rA, sprintf('%s\n%s', names{1}, labA), ...
    'HorizontalAlignment','center','FontSize',11);
text(cB(1)+0.7*rB, cB(2)-0.1*rB, sprintf('%s\n%s', names{2}, labB), ...
    'HorizontalAlignment','center','FontSize',11);
text(cC(1), cC(2)+0.85*rC, sprintf('%s\n%s', names{3}, labC), ...
    'HorizontalAlignment','center','FontSize',11);

mAB = (cA+cB)/2;
mAC = (cA+cC)/2;
mBC = (cB+cC)/2;
ctr = (cA+cB+cC)/3;

text(mAB(1), mAB(2), labAB,  'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
text(mAC(1), mAC(2), labAC,  'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
text(mBC(1), mBC(2), labBC,  'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
text(ctr(1), ctr(2), labABC, 'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');

pad = 0.35*max([rA,rB,rC]);
xmin = min([cA(1)-rA, cB(1)-rB, cC(1)-rC]) - pad;
xmax = max([cA(1)+rA, cB(1)+rB, cC(1)+rC]) + pad;
ymin = min([cA(2)-rA, cB(2)-rB, cC(2)-rC]) - pad;
ymax = max([cA(2)+rA, cB(2)+rB, cC(2)+rC]) + pad;
xlim([xmin xmax]); ylim([ymin ymax]);
end


function [cC, ok] = thirdPointFromDistances(cA, cB, dAC, dBC)
ok = true;
AB = norm(cB - cA);
if AB < 1e-12
    ok = false; cC = [NaN NaN]; return;
end
if (dAC + dBC < AB) || (abs(dAC - dBC) > AB)
    ok = false; cC = [NaN NaN]; return;
end
x = (dAC^2 - dBC^2 + AB^2) / (2*AB);
y2 = max(0, dAC^2 - x^2);
y = sqrt(y2);
u = (cB - cA) / AB;
v = [-u(2), u(1)];
cC = cA + x*u + y*v;
end

function printEpochBarStats(labels, stackCounts_with_place, Ynorm)

K = numel(labels);
overlapLabels = cell(1,K);
for k = 1:K
    if k==1
        overlapLabels{k} = 'Only this epoch';
    else
        overlapLabels{k} = sprintf('+%d other epoch%s', k-1, ternary(k-1==1,'','s'));
    end
end

fprintf('\n============================================\n');
fprintf('STACKED BAR STATS (POOLED ACROSS ALL CELLS)\n');
fprintf('============================================\n');

for e = 1:(K+1)   % includes Place as last row
    if e <= K
        name = labels{e};
    else
        name = 'Place';
    end

    totalFrac = sum(stackCounts_with_place(e,:));
    totalNorm = sum(Ynorm(e,:));

    fprintf('\n%s:\n', name);
    fprintf('  Total fraction of kept: %.4f (%.2f%%)\n', ...
        totalFrac, 100*totalFrac);

    if isfinite(totalNorm)
        fprintf('  Total density (fraction/s): %.4f\n', totalNorm);
    end

    for k = 1:K
        if stackCounts_with_place(e,k) > 0
            fprintf('    %-18s: %.4f (%.2f%%)\n', ...
                overlapLabels{k}, ...
                stackCounts_with_place(e,k), ...
                100*stackCounts_with_place(e,k));
        end
    end
end

fprintf('\n============================================\n\n');
end

function printNumEpochsSummary(ALL_sets_task)

if isempty(ALL_sets_task)
    fprintf('\nNum-epoch summary: no task-modulated cells.\n');
    return
end

nEpochs = sum(ALL_sets_task>0, 2);   % 1..4 (should never be 0 for rows that exist)
N = numel(nEpochs);

counts = zeros(1,4);
for k = 1:4
    counts(k) = sum(nEpochs == k);
end
pct = 100 * counts / max(1,N);

fprintf('\n============================================\n');
fprintf('TASK CELLS: # epochs each cell belongs to (POOLED)\n');
fprintf('Denominator = task-modulated cells (N=%d)\n', N);
fprintf('============================================\n');
fprintf('1 epoch : %6d  (%.2f%%)\n', counts(1), pct(1));
fprintf('2 epochs: %6d  (%.2f%%)\n', counts(2), pct(2));
fprintf('3 epochs: %6d  (%.2f%%)\n', counts(3), pct(3));
fprintf('4 epochs: %6d  (%.2f%%)\n', counts(4), pct(4));
fprintf('============================================\n\n');

end
