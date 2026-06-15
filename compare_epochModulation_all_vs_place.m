function S = compare_epochModulation_all_vs_place(R_all, R_pc)
% Compare epoch proportions between:
%   R_all = epochModulation_venn(...,'PlaceOnly',false)
%   R_pc  = epochModulation_venn(...,'PlaceOnly',true)
%
% Paired t-test: per-rat epoch fractions
% KS test: pooled cell-level epoch-membership distributions

labels = R_all.labels;
K = numel(labels);

% ---------- per-rat fractions (better for t-test) ----------
Fall = perRatTotalEpochFractions(R_all.ALL_sets_task, R_all.ALL_rat_task, R_all.ALL_rat_kept, K);
Fpc  = perRatTotalEpochFractions(R_pc.ALL_sets_task,  R_pc.ALL_rat_task,  R_pc.ALL_rat_kept,  K);

fprintf('\n============================================\n');
fprintf('ALL CELLS vs PLACE-ONLY\n');
fprintf('Per-rat epoch fractions (paired t-test)\n');
fprintf('============================================\n');
fprintf('%-10s %-10s %-10s %-10s %-10s\n', 'Epoch','meanAll','meanPlace','tstat','p_t');

S.ttest = struct();
for e = 1:K
    xa = Fall(:,e);
    xb = Fpc(:,e);

    good = isfinite(xa) & isfinite(xb);
    xa = xa(good);
    xb = xb(good);

    if numel(xa) >= 2
        [~,p,~,stats] = ttest(xa, xb);
        tstat = stats.tstat;
    else
        p = NaN; tstat = NaN;
    end

    fprintf('%-10s %-10.4f %-10.4f %-10.4f %-10.3g\n', ...
        labels{e}, mean(xa,'omitnan'), mean(xb,'omitnan'), tstat, p);

    S.ttest.(labels{e}).all   = xa;
    S.ttest.(labels{e}).place = xb;
    S.ttest.(labels{e}).p     = p;
    S.ttest.(labels{e}).tstat = tstat;
end

% ---------- pooled cell-level distributions (for KS) ----------
% Here I compare the binary epoch-membership distributions across cells:
% is the cell in this epoch-set or not?
fprintf('\n============================================\n');
fprintf('ALL CELLS vs PLACE-ONLY\n');
fprintf('Pooled cell-level membership (KS test)\n');
fprintf('============================================\n');
fprintf('%-10s %-10s %-10s\n', 'Epoch','KS stat','p_KS');

S.ks = struct();
for e = 1:K
    xa = double(R_all.ALL_sets_task(:,e) > 0);
    xb = double(R_pc.ALL_sets_task(:,e)  > 0);

    if numel(xa) >= 2 && numel(xb) >= 2
        [~,pks,ksstat] = kstest2(xa, xb);
    else
        pks = NaN; ksstat = NaN;
    end

    fprintf('%-10s %-10.4f %-10.3g\n', labels{e}, ksstat, pks);

    S.ks.(labels{e}).p      = pks;
    S.ks.(labels{e}).ksstat = ksstat;
end

% ---------- optional omnibus summaries ----------
% Compare distribution of "number of epochs per task cell"
na = sum(R_all.ALL_sets_task > 0, 2);
nb = sum(R_pc.ALL_sets_task  > 0, 2);

if numel(na) >= 2 && numel(nb) >= 2
    [~,pksN,ksN] = kstest2(na, nb);
else
    pksN = NaN; ksN = NaN;
end

fprintf('\nNum-epochs-per-cell KS: stat = %.4f, p = %.3g\n', ksN, pksN);

S.numEpochs.ksstat = ksN;
S.numEpochs.p      = pksN;
end


function F = perRatTotalEpochFractions(setsTask, ratTask, ratKept, K)
% Per-rat total fraction of kept cells belonging to each epoch
rats = unique(ratKept(:));
rats = rats(isfinite(rats));

F = nan(numel(rats), K);

for i = 1:numel(rats)
    rr = rats(i);

    denom = sum(ratKept == rr);   % kept cells in this rat
    if denom <= 0
        continue;
    end

    setsR = setsTask(ratTask == rr, :);
    if isempty(setsR)
        F(i,:) = 0;
        continue;
    end

    for e = 1:K
        F(i,e) = sum(setsR(:,e) > 0) / denom;
    end
end
end
