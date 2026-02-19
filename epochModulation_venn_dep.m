function R = epochModulation_venn_(ratNames, varargin)
% epochModulation_venn
% Detect per-cell firing-rate modulation in CS / Trace / US / Post epochs
% vs FULL non-trial time, visualize overlaps (UpSet or Venn), and save
% dominant epoch label per cell (rat.epoch.epoch_<dayStr>).
%
% ADDITIONS:
%   - Track place membership for every kept cell (even when PlaceOnly=false)
%   - Two extra proportional-circle diagrams:
%       (1) Task (any epoch) vs Place (2 circles)
%       (2) Pre(CS∪Trace∪US) vs Post vs Place (3 circles; pairwise-proportional)
%   - EqualWeightRats: Venn labels are rat-equalized (each rat contributes equally)
%   - EqualizeBars: stacked bar panels are rat-equalized (each rat contributes equally)

% ---------------- parameters ----------------
p = inputParser;
p.addParameter('Alpha',0.05);
p.addParameter('RestrictSpeed',false);
p.addParameter('SpeedMin',0.0);
p.addParameter('SpeedMax',4);
p.addParameter('SpeedSource','auto');
p.addParameter('PlotMode','venn');
%p.addParameter('PlotMode','bar');
p.addParameter('VennEpochs',{'CS','Trace','US','Post'});
p.addParameter('CollapsePreIntoOne', true, @islogical);  % CS∪Trace∪US vs Post
p.addParameter('VennPercent',true,@islogical);
p.addParameter('PercentBase','kept',@(s)any(strcmpi(s,{'kept','union'})));
p.addParameter('EqualWeightRats', true, @islogical);     % NEW: rat-equalized Venn numbers
p.addParameter('EqualizeBars', false, @islogical);       % NEW: rat-equalized bar panels
p.addParameter('Verbose',true);

% place-cell options
p.addParameter('PlaceOnly',false,@islogical);
p.addParameter('PlaceSpec',struct('Struct','MI_noCSUS15_shuff','Prefix','MI_', ...
    'Col',3,'Thresh',0.95,'Comparator','>'));
p.parse(varargin{:});
opt = p.Results;



labels = {'CS','Trace','US','Post'};
epochs = [0.00 0.25; 0.25 0.75; 0.75 0.85; 0.85 2.00];   % 4 epochs
trialWin = [0 6.0];

nR = numel(ratNames);
ALL_sets   = [];    % [Ncells_total x 4] logical: significant vs non-trial
ALL_dir    = [];    % [Ncells_total x 4] int8 direction (up/down)
ALL_pvals  = [];    % [Ncells_total x 4] Poisson p-values
ALL_lambda = [];    % [Ncells_total x 4] expected counts lambda
ALL_zlike  = [];    % [Ncells_total x 4] (obs - lambda)/sqrt(lambda)
ALL_place  = [];    % [Ncells_total x 1] logical place membership (NaN if MI missing)
ALL_rat    = [];    % [Ncells_total x 1] rat index for each kept cell

for r = 1:nR
    rat = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates,rat.An),1);
    days  = dates(max(1,idx-2):idx);



    for d = 1:numel(days)
        dayStr = days{d};
        spk    = rat.Ca_peaks.(sprintf('CA_peaks_%s',dayStr));
        csTimes= rat.CS_times.(sprintf('CS_%s',dayStr));
        mask   = rat.ratemask.(sprintf('ratemask_%s',dayStr))==1;

        % --- place membership tracking (even when PlaceOnly=false) ---
        placeThisDay = nan(size(mask));  % NaN means MI missing
        [pcMask, okPC] = fetchPlaceMask(rat, dayStr, opt.PlaceSpec);
        if okPC
            placeThisDay = pcMask(:);
        else
            if opt.PlaceOnly
                warning('[%s %s] PlaceOnly requested but MI table missing.', ...
                        ratNames{r}, dayStr);
            end
        end

        % optional: restrict analysis to place cells only
        if opt.PlaceOnly && okPC
            mask = mask & pcMask;
        end

        if isempty(csTimes), continue; end
        sessMin = min(csTimes)-60;
        sessMax = max(csTimes)+60;
        excl = [csTimes(:)+trialWin(1), csTimes(:)+trialWin(2)];
        pool = complementIntervals([sessMin sessMax], excl);

        poolDur = sum(max(0,pool(:,2)-pool(:,1)));
        if poolDur<=0, continue; end

        nEp   = size(epochs,1);
        epAbs = cell(nEp,1);
        epDur = zeros(nEp,1);
        for e = 1:nEp
            epAbs{e} = [csTimes(:)+epochs(e,1), csTimes(:)+epochs(e,2)];
            epDur(e) = diff(epochs(e,:));
        end

        nTrials = numel(csTimes);
        nCells  = size(spk,1);
        pvals     = nan(nCells,nEp);
        dirsgn    = zeros(nCells,nEp,'int8');
        lambdaMat = nan(nCells,nEp);
        zMat      = nan(nCells,nEp);
        keep      = false(nCells,1);

        for c = 1:nCells
            if ~mask(c), continue; end
            st = spk(c,:);
            st = st(~isnan(st) & st>0);
            if isempty(st), continue; end
            keep(c) = true;

            poolCounts = intervalCountsPerRow(st, pool);
            r_pool     = sum(poolCounts)/poolDur;

            for e = 1:nEp
                counts = intervalCountsPerRow(st, epAbs{e});
                obsK   = sum(counts);
                expT   = epDur(e)*nTrials;
                lambda = r_pool*expT;

                pvals(c,e)     = poisson_two_sided(obsK, lambda);
                lambdaMat(c,e) = lambda;

                if lambda > 0
                    zMat(c,e) = (obsK - lambda) ./ sqrt(lambda);
                else
                    zMat(c,e) = NaN;
                end

                if pvals(c,e) <= opt.Alpha
                    dirsgn(c,e) = int8(sign(double(obsK) - lambda));
                end
            end

            % After filling pvals (nCells x nEp)
            alpha = opt.Alpha;
            sigs = false(size(pvals));  % corrected significance calls

            for c = 1:nCells
                pv = pvals(c,:);
                ok = isfinite(pv);
                if nnz(ok) < 1, continue; end

                % Holm-Bonferroni across the epochs for this cell
                [pSort, ord] = sort(pv(ok), 'ascend');
                m = numel(pSort);

                pass = false(1,m);
                for k = 1:m
                    if pSort(k) <= alpha / (m - k + 1)
                        pass(k) = true;
                    else
                        break; % Holm: once you fail, all larger p fail
                    end
                end

                % Map back to epoch indices
                idxOk = find(ok);
                sigIdxLocal = ord(pass);
                sigEpochs = idxOk(sigIdxLocal);
                sigs(c, sigEpochs) = true;
            end
        end

        %sigs = pvals <= opt.Alpha;
        % Fraction of kept cells significant in ANY epoch
        anyEpoch = any(sigs(keep,:), 2);
        fracAny  = mean(anyEpoch);

        if opt.Verbose
            fprintf('  Any-epoch modulated (this session): %.2f%%\n', 100*fracAny);
        end

        ALL_sets   = [ALL_sets;   sigs(keep,:)];          %#ok<AGROW>
        ALL_dir    = [ALL_dir;    dirsgn(keep,:)];        %#ok<AGROW>
        ALL_pvals  = [ALL_pvals;  pvals(keep,:)];         %#ok<AGROW>
        ALL_lambda = [ALL_lambda; lambdaMat(keep,:)];     %#ok<AGROW>
        ALL_zlike  = [ALL_zlike;  zMat(keep,:)];          %#ok<AGROW>
        ALL_place  = [ALL_place;  placeThisDay(keep)];    %#ok<AGROW>
        ALL_rat    = [ALL_rat;    r * ones(nnz(keep),1)]; %#ok<AGROW>

        if opt.Verbose
            fprintf('[%s %s] %d cells kept, sig CS=%d Trace=%d US=%d Post=%d\n',...
                ratNames{r}, dayStr, nnz(keep), ...
                nnz(sigs(:,1)), nnz(sigs(:,2)), nnz(sigs(:,3)), nnz(sigs(:,4)));
        end

        % ===== Save dominant epoch label per cell for this day =====
        % Codes: 1 = CS, 2 = Trace, 3 = US, 4 = Post
        %        0 = analyzed but not significantly modulated
        %      NaN = not analyzed / excluded cell
        epochVec = nan(nCells,1);
        for c = 1:nCells
            if ~keep(c), continue; end
            pv = pvals(c,:);
            sigMask = pv <= opt.Alpha & isfinite(pv);

            if ~any(sigMask)
                epochVec(c) = 0;
            else
                zc = zMat(c,:);
                zc(~sigMask) = NaN;
                [~, bestIdx] = max(abs(zc)); % 1..4
                epochVec(c) = bestIdx;
            end
        end

        if ~isfield(rat,'epoch') || ~isstruct(rat.epoch)
            rat.epoch = struct();
        end
        fld = sprintf('epoch_%s', dayStr);
        rat.epoch.(fld) = epochVec;

        assignin('base', ratNames{r}, rat);
    end
end

anyEpoch_all = any(ALL_sets, 2);
fracAny_all  = mean(anyEpoch_all);

fprintf('\nOverall any-epoch fraction: %.2f%%\n', 100*fracAny_all);

if isempty(ALL_sets)
    warning('No cells kept for epochModulation_venn.');
    R.opts = opt;
    return;
end

epochDur = diff(epochs,1,2);   % [4 x 1]
K = size(ALL_sets,2);
N = size(ALL_sets,1);

% ---------- epoch-wise significance diagnostics ----------
fprintf('\n=== Epoch-wise significance diagnostics ===\n');
for e = 1:K
    nSig   = nnz(ALL_sets(:,e));
    nUp    = nnz(ALL_dir(:,e) > 0);
    nDown  = nnz(ALL_dir(:,e) < 0);

    lam_e    = ALL_lambda(:,e);
    testMask = isfinite(lam_e) & lam_e > 0;
    nTested  = nnz(testMask);
    expNull  = opt.Alpha * nTested;

    fracSig  = nSig / max(1,N);
    ratioExp = nSig / max(1,expNull);

    fprintf('Epoch %s:\n', labels{e});
    fprintf('  Sig cells: %d (%.2f%% of %d; up=%d, down=%d)\n', ...
        nSig, 100*fracSig, N, nUp, nDown);
    fprintf('  Tested cells (lambda>0): %d, expected under null ≈ %.1f\n', ...
        nTested, expNull);
    fprintf('  Observed / expected (null) ≈ %.2f\n', ratioExp);
end
fprintf('=========================================\n\n');

% ---------- bar panel data (pooled vs rat-equalized) ----------
if opt.EqualizeBars
    stackCounts  = ratEqualEpochStackCounts(ALL_sets, ALL_rat); % [K x K], fractions per rat then mean
    prepostStack = ratEqualPrePostStack(ALL_sets, ALL_rat);     % [2 x 2], fractions per rat then mean
else
    stackCounts  = pooledEpochStackCounts(ALL_sets);            % [K x K], raw counts
    prepostStack = pooledPrePostStack(ALL_sets);                % [2 x 2], raw counts
end

% Time-normalized stack (same transform for either)
Ynorm = bsxfun(@rdivide, stackCounts, epochDur);

% ==========================================================
%   Plotting
% ==========================================================
switch lower(opt.PlotMode)
    case 'upset'
        figure('Color','w');
        upsetPlotGeneric(ALL_sets, labels);

    case 'venn'
        % FIGURE 1: stacked bars
        figure('Color','w');

        subplot(3,1,1);
        bh1 = bar(1:K, stackCounts, 'stacked');
        ylabel(ternary(opt.EqualizeBars,'Fraction of cells (rat-equalized)',' # cells'));
        set(gca,'XTick',1:K,'XTickLabel',labels);
        box off;
        legStr = cell(1,K);
        for k = 1:K
            if k == 1
                legStr{k} = 'Only this epoch';
            else
                legStr{k} = sprintf('+%d other epoch%s', ...
                    k-1, ternary(k-1==1,'','s'));
            end
        end
        legend(bh1, legStr, 'Location','bestoutside');
        title(ternary(opt.EqualizeBars,'Epoch modulation overlap (rat-equalized)','Epoch modulation overlap (counts)'));

        subplot(3,1,2);
        bh2 = bar(1:2, prepostStack, 'stacked');
        ylabel(ternary(opt.EqualizeBars,'Fraction of cells (rat-equalized)',' # cells'));
        set(gca,'XTick',1:2,'XTickLabel',{'CS/Trace/US','Post'});
        box off;
        legend(bh2, {'Only this bin','Both bins'}, 'Location','bestoutside');
        title(ternary(opt.EqualizeBars,'CS/Trace/US vs Post overlap (rat-equalized)','CS/Trace/US vs Post overlap'));

        subplot(3,1,3);
        bh3 = bar(1:K, Ynorm, 'stacked');
        ylabel(ternary(opt.EqualizeBars,'(Fraction per s) rat-equalized','Cells / s of epoch window'));
        set(gca,'XTick',1:K,'XTickLabel',labels);
        box off;
        legend(bh3, legStr, 'Location','bestoutside');
        title(ternary(opt.EqualizeBars,'Epoch modulation overlap / s (rat-equalized)','Epoch modulation overlap (time-normalized)'));

        % FIGURE 2: Venn (labels can be rat-equalized)
        plotSets   = ALL_sets;
        plotLabels = {'CS','Trace','US','Post'};

        if opt.CollapsePreIntoOne
            if size(ALL_sets,2) < 4
                warning('CollapsePreIntoOne requested, but ALL_sets has <4 columns. Skipping collapse.');
            else
                plotSets   = [ any(ALL_sets(:,1:3), 2), ALL_sets(:,4) ];  % [N x 2]
                plotLabels = {'CS/Trace/US','Post'};
            end
        end

        if ~opt.CollapsePreIntoOne
            [epMask, labsK] = selectEpochSubset(plotLabels, opt.VennEpochs);
            setsK  = plotSets(:, epMask);
        else
            setsK  = plotSets;
            labsK  = plotLabels;
        end

        kSel = size(setsK,2);
        if kSel < 2
            warning('Venn requires ≥2 sets; got %d. Skipping Venn plot.', kSel);
            R.opts = opt;
            return
        end

        % --- venn label strings ---
        if opt.EqualWeightRats
            if kSel == 2
                F = ratEqualExclusiveFrac2(setsK, ALL_rat, opt.PercentBase);
                lbl = { sprintf('%.1f%%',100*F.A_only), ...
                        sprintf('%.1f%%',100*F.B_only), ...
                        sprintf('%.1f%%',100*F.AB) };
            elseif kSel == 3
                F = ratEqualExclusiveFrac3(setsK, ALL_rat, opt.PercentBase);
                lbl = { sprintf('%.1f%%',100*F.A_only), sprintf('%.1f%%',100*F.B_only), sprintf('%.1f%%',100*F.C_only), ...
                        sprintf('%.1f%%',100*F.AB),     sprintf('%.1f%%',100*F.AC),     sprintf('%.1f%%',100*F.BC), ...
                        sprintf('%.1f%%',100*F.ABC) };
            else
                F = ratEqualExclusiveFrac4(setsK, ALL_rat, opt.PercentBase);
                F
                lbl = { sprintf('%.1f%%',100*F.A_only), sprintf('%.1f%%',100*F.B_only), ...
                        sprintf('%.1f%%',100*F.C_only), sprintf('%.1f%%',100*F.D_only), ...
                        sprintf('%.1f%%',100*F.AB),     sprintf('%.1f%%',100*F.AC),     sprintf('%.1f%%',100*F.AD), ...
                        sprintf('%.1f%%',100*F.BC),     sprintf('%.1f%%',100*F.BD),     sprintf('%.1f%%',100*F.CD), ...
                        sprintf('%.1f%%',100*F.ABC),    sprintf('%.1f%%',100*F.ABD),    sprintf('%.1f%%',100*F.ACD), ...
                        sprintf('%.1f%%',100*F.BCD),    sprintf('%.1f%%',100*F.ABCD) };


            end
        else
            switch lower(opt.PercentBase)
                case 'union'
                    denom = max(1, sum(any(setsK,2)));
                otherwise
                    denom = size(setsK,1);
            end
            asPct = opt.VennPercent;

            if kSel == 2
                C = computeExclusiveCounts(setsK);
                lbl = { fmtCount(C.A_only,denom,asPct), ...
                        fmtCount(C.B_only,denom,asPct), ...
                        fmtCount(C.AB,     denom,asPct) };
            elseif kSel == 3
                C = computeExclusiveCounts(setsK);
                lbl = { fmtCount(C.A_only,denom,asPct), fmtCount(C.B_only,denom,asPct), fmtCount(C.C_only,denom,asPct), ...
                        fmtCount(C.AB,denom,asPct),     fmtCount(C.AC,denom,asPct),     fmtCount(C.BC,denom,asPct), ...
                        fmtCount(C.ABC,denom,asPct) };
            else
                C = computeExclusiveCounts4(setsK);
                lbl = { fmtCount(C.A_only,denom,asPct), fmtCount(C.B_only,denom,asPct), ...
                        fmtCount(C.C_only,denom,asPct), fmtCount(C.D_only,denom,asPct), ...
                        fmtCount(C.AB,denom,asPct),     fmtCount(C.AC,denom,asPct),     fmtCount(C.AD,denom,asPct), ...
                        fmtCount(C.BC,denom,asPct),     fmtCount(C.BD,denom,asPct),     fmtCount(C.CD,denom,asPct), ...
                        fmtCount(C.ABC,denom,asPct),    fmtCount(C.ABD,denom,asPct),    fmtCount(C.ACD,denom,asPct), ...
                        fmtCount(C.BCD,denom,asPct),    fmtCount(C.ABCD,denom,asPct) };
            end
        end

        figure('Color','w'); hold on;
        if kSel == 2
            col = parula(5); col = col([2,5],:);
            venn(2, 'sets', labsK, 'labels', lbl, 'colors', col, ...
                 'alpha', 0.5, 'edgeC', [1 1 1], 'edgeW', 3);
        elseif kSel == 3
            col = parula(3);
            venn(3, 'sets', labsK, 'labels', lbl, 'colors', col, ...
                 'alpha', 0.5, 'edgeC', [1 1 1], 'edgeW', 3);
        elseif kSel == 4
            col = parula(4);
            venn(4, 'sets', labsK, 'labels', lbl, 'colors', col, ...
                 'alpha', 0.5, 'edgeC', [1 1 1], 'edgeW', 3);
        end

        % ==========================================================
        %   EXTRA proportional-size diagrams (rat-equalized geometry optional)
        % ==========================================================
        placeOK = isfinite(ALL_place);
        if any(placeOK)
            placeMask = false(size(ALL_place));
            placeMask(placeOK) = ALL_place(placeOK) > 0;

            % (A) Task vs Place (2 circles)
            taskMask = any(ALL_sets, 2);
            valid2 = placeOK;
            A = taskMask(valid2);
            B = placeMask(valid2);
            ratIDs2 = ALL_rat(valid2);

            figure('Color','w');
            drawProportionalVenn2(A, B, ratIDs2, {'Task (any epoch)','Place'}, ...
                opt.PercentBase, opt.EqualWeightRats);

            % (B) Pre/Post/Place (3 circles)
            preAny  = any(ALL_sets(:,1:3), 2);
            postAny = ALL_sets(:,4);

            valid3 = placeOK;
            A3 = preAny(valid3);
            B3 = postAny(valid3);
            C3 = placeMask(valid3);
            ratIDs3 = ALL_rat(valid3);

            figure('Color','w');
            drawProportionalVenn3(A3, B3, C3, ratIDs3, {'CS/Trace/US','Post','Place'}, ...
                opt.PercentBase, opt.EqualWeightRats);
        else
            warning('No valid place labels (MI missing) so skipping proportional Place overlap plots.');
        end
end

% Expose everything
R.opts         = opt;
R.labels       = labels;
R.epochDur     = epochDur;
R.ALL_sets     = ALL_sets;
R.ALL_dir      = ALL_dir;
R.ALL_pvals    = ALL_pvals;
R.ALL_lambda   = ALL_lambda;
R.ALL_zlike    = ALL_zlike;
R.ALL_place    = ALL_place;
R.ALL_rat      = ALL_rat;
R.stackCounts  = stackCounts;
R.prepostStack = prepostStack;
R.Ynorm        = Ynorm;
end

% =====================================================================
% Helpers
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

function upsetPlotGeneric(sets, labels)
K = size(sets,2);
if K < 2 || K > 4
    title(sprintf('UpSet: supports 2–4 sets (got %d)', K)); return
end
bits = packBitsN(sets);
u = (1:(2^K - 1))';
cnt = arrayfun(@(b) sum(bits==b), u);
[cnt, ord] = sort(cnt, 'descend');
u = u(ord); u = u(cnt>0); cnt = cnt(cnt>0);
M = numel(u);

subplot(2,1,1);
bar(1:M, cnt, 0.85);
ylabel('# cells');
set(gca,'XTick',[]);
title(sprintf('UpSet (%d sets)', K));
box off

subplot(2,1,2);
cla; hold on
for i=1:M
    bvec = bitget(u(i), 1:K);
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
for i=1:N
    counts(i)=sum(st>=win(i,1)&st<win(i,2));
end
end

function p=poisson_two_sided(k,lambda)
if lambda<=0, p=double(k>0); return, end
Fc=poisscdf(k,lambda);
Fg=1-poisscdf(k-1,lambda);
p=2*min(Fc,Fg);
p=min(1,p);
end

function [pcMask,ok]=fetchPlaceMask(rat,dayStr,spec)
ok=false;pcMask=[];
fld=[spec.Prefix dayStr];
if isfield(rat,spec.Struct) && isfield(rat.(spec.Struct),fld)
    MI=rat.(spec.Struct).(fld);
    ok=true;
else
    return
end
if size(MI,2) < spec.Col
    ok=false;
    return
end
v=MI(:,spec.Col);
switch spec.Comparator
    case '>', pcMask = v > spec.Thresh;
    case '<', pcMask = v < spec.Thresh;
    otherwise, pcMask = v > spec.Thresh;
end
pcMask = pcMask(:) > 0;
end

function s=fmtCount(n,denom,asPercent)
if asPercent
    s=sprintf('%.1f%%',100*n/max(1,denom));
else
    s=sprintf('%d',n);
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
a=full(1);b=full(2);
E=excl;
E=E(~any(isnan(E),2),:);
E=sortrows(E,1);
E(:,1)=max(E(:,1),a);
E(:,2)=min(E(:,2),b);
E=E(E(:,1)<E(:,2),:);
C=[];
t=a;
for i=1:size(E,1)
    if E(i,1)>t
        C=[C;t E(i,1)]; %#ok<AGROW>
    end
    t=max(t,E(i,2));
end
if t<b
    C=[C;t b]; %#ok<AGROW>
end
end

% ==========================================================
% Bar helpers: pooled vs rat-equalized
% ==========================================================
function S = pooledEpochStackCounts(ALL_sets)
K = size(ALL_sets,2);
S = zeros(K,K);
for e = 1:K
    idx = find(ALL_sets(:,e));
    if isempty(idx), continue; end
    othersCount = sum(ALL_sets(idx,:),2) - 1;
    for k = 0:(K-1)
        S(e,k+1) = sum(othersCount == k);
    end
end
end

function S = ratEqualEpochStackCounts(ALL_sets, ratIDs)
K = size(ALL_sets,2);
rats = unique(ratIDs(:)); rats = rats(isfinite(rats));
nR = numel(rats);

Srat = zeros(nR, K, K);  % [rat x epoch x overlapBin] as fractions within rat
for i = 1:nR
    idxR = ratIDs == rats(i);
    setsR = ALL_sets(idxR,:);
    denom = max(1, size(setsR,1));
    for e = 1:K
        idx = find(setsR(:,e));
        if isempty(idx), continue; end
        othersCount = sum(setsR(idx,:),2) - 1;
        for k = 0:(K-1)
            Srat(i,e,k+1) = sum(othersCount == k) / denom;
        end
    end
end
S = squeeze(mean(Srat,1,'omitnan')); % [K x K] mean fraction across rats
end

function P = pooledPrePostStack(ALL_sets)
preMask  = any(ALL_sets(:,1:3), 2);
postMask = ALL_sets(:,4);
P = [ sum(preMask  & ~postMask), sum(preMask  &  postMask); ...
      sum(postMask & ~preMask), sum(postMask &  preMask) ];
end

function P = ratEqualPrePostStack(ALL_sets, ratIDs)
rats = unique(ratIDs(:)); rats = rats(isfinite(rats));
nR = numel(rats);

vals = zeros(nR,2,2);  % [rat x (pre/post) x (only/both)] fractions within rat
for i = 1:nR
    idxR = ratIDs == rats(i);
    setsR = ALL_sets(idxR,:);
    denom = max(1, size(setsR,1));

    preMask  = any(setsR(:,1:3), 2);
    postMask = setsR(:,4);

    vals(i,1,1) = sum(preMask  & ~postMask) / denom;
    vals(i,1,2) = sum(preMask  &  postMask) / denom;
    vals(i,2,1) = sum(postMask & ~preMask) / denom;
    vals(i,2,2) = sum(postMask &  preMask) / denom;
end
P = squeeze(mean(vals,1,'omitnan')); % [2 x 2]
end

% ==========================================================
% Rat-equalized exclusive fractions (Venn labels + proportional geometry)
% ==========================================================
function F = ratEqualExclusiveFrac2(sets2, ratIDs, percentBase)
A = sets2(:,1)>0; B = sets2(:,2)>0;
rats = unique(ratIDs(:)); rats = rats(isfinite(rats));
nR = numel(rats);
vals = nan(nR,3); % [A_only, B_only, AB]
for i = 1:nR
    idx = ratIDs==rats(i);
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
nR = numel(rats);
vals = nan(nR,7); % Aonly,Bonly,Conly,AB,AC,BC,ABC
for i = 1:nR
    idx = ratIDs==rats(i);
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
nR = numel(rats);
vals = nan(nR,15);
for i = 1:nR
    idx = ratIDs==rats(i);
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
% Proportional circle diagrams (rat-equalized geometry optional)
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

if f(lo) < 0
    d = lo; return;
end
if f(hi) > 0
    d = hi; return;
end

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

function drawCircle(center, radius, faceAlpha)
t = linspace(0,2*pi,400);
x = center(1) + radius*cos(t);
y = center(2) + radius*sin(t);
patch(x,y,1,'FaceAlpha',faceAlpha,'EdgeColor',[1 1 1],'LineWidth',2);
end
