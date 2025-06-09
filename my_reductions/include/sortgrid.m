function f = sortgrid(allAUROC)


sz = size(allAUROC);
diffs = nan(sz);

pos = 0;
neg = 0;
% Loop through all combinations
for i = 1:sz(1)
    for j = 1:sz(2)
        for k = 1:sz(3)
            for l = 1:sz(4)
                for m = 1:sz(5)
                    auc = allAUROC{i,j,k,l,m};
                    if isempty(auc) || numel(auc) < 15, continue; end
                    early = mean(auc(1:6), 'omitnan');
                    late  = mean(auc(7:14), 'omitnan');

                    diffs(i,j,k,l,m) = late - early;
                    diffs(i,j,k,l,m) = auc(14);
                end
            end
        end
    end
end
pos
neg

% Find top N largest increases
N = 10;
diff_vals = diffs(:);
[sortedVals, sortedIdx] = maxk(diff_vals, N);
[topI, topJ, topK, topL, topM] = ind2sub(sz, sortedIdx);

figure
hold on
% Display
fprintf('Top %d increases in AUROC (late - early):\n', N);
for idx = 10:10
    fprintf('(%d,%d,%d,%d,%d): Δ=%.3f\n', topI(idx), topJ(idx), topK(idx), topL(idx), topM(idx), sortedVals(idx));
    vals = cell2mat(allAUROC(topI(idx), topJ(idx), topK(idx), topL(idx), topM(idx)));
    plot(vals, 'LineWidth', 1.5)
end

f = [topI, topJ, topK, topL, topM];
