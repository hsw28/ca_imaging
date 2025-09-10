function [R,f] = plotSpeedDependenceSummary(ratNames, varargin)
% runAndPlotSpeedDependence  Wrapper that calls speedDependence_taskVsRun
% and immediately produces summary plots.
%
% Usage:
%   [R,f] = plotSpeedDependenceSummary({'rat0222','rat0307','rat0313','rat0314','rat0816'});
%
% Inputs:
%   ratNames : cellstr of rat variable names in base workspace
%   varargin : name/value options passed to speedDependence_taskVsRun
%
% Outputs:
%   R : results struct from speedDependence_taskVsRun
%   f : handle to summary figure

% ---------- run analysis ----------
R = speedDependence_taskVsRun(ratNames, varargin{:});

% ---------- plot summary ----------
f = figure('Color','w','Position',[100 100 1200 800]);
pooledColors = [0.2 0.2 0.2];
nR = numel(ratNames);

rt = R.pooled.r_task(:);
rr = R.pooled.r_run(:);
dr = R.pooled.delta(:);

ratIdx = [];
if isfield(R,'r_task') && ~isempty(R.r_task{1})
    ratIdx = cell2mat(arrayfun(@(i) i*ones(numel(R.r_task{i}),1), (1:nR)', 'uni', 0));
end

% --- A: scatter r_task vs r_run
subplot(2,3,1); hold on;
ok = isfinite(rt) & isfinite(rr);
if ~isempty(ratIdx)
    cmap = lines(nR);
    for i = 1:nR
        ii = ok & (ratIdx==i);
        scatter(rr(ii), rt(ii), 10, cmap(i,:), 'filled', ...
            'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2);
    end
else
    scatter(rr(ok), rt(ok), 10, pooledColors, 'filled', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.2);
end
plot([-1 1],[-1 1],'k--');
xlim([-1 1]); ylim([-1 1]); axis square;
xlabel('r_{run} (speed–rate)'); ylabel('r_{task} (speed–rate)');
title('Task vs Run speed–dependence (pooled)');

% --- B: Δr histogram
subplot(2,3,2); hold on;
edges = -0.6:0.03:0.6;
histogram(dr(isfinite(dr)), edges, 'FaceAlpha',0.7,'EdgeColor','none');
plot([0 0], ylim, 'k--');
xlabel('\Delta r = r_{task} - r_{run}');
ylabel('# cells');
title('\Delta r distribution');

% --- C: Fraction significant bars
subplot(2,3,3); hold on;
vals = nan(nR+1,3);
for i = 1:nR
    vals(i,1) = R.fracSig_task(i);
    vals(i,2) = R.fracSig_run(i);
    vals(i,3) = R.fracSig_delta(i);
end
vals(end,1) = R.pooled.fracSig_task;
vals(end,2) = R.pooled.fracSig_run;
vals(end,3) = R.pooled.fracSig_delta;
bar(vals,'grouped');
legend({'Task','Run','\Delta'},'Location','northwest'); legend boxoff;
set(gca,'XTick',1:nR+1,'XTickLabel',[ratNames(:); {'pooled'}], ...
    'XTickLabelRotation',30);
ylabel('fraction significant (FDR q=0.05)');
title('FDR‐significant proportions');

% --- D: Speed‐matched rate ratios
subplot(2,3,[4 6]); hold on;
boxes = cell(1,nR+1);
for i = 1:nR
    if i<=numel(R.rateRatio_match) && ~isempty(R.rateRatio_match{i})
        boxes{i} = R.rateRatio_match{i}(:);
    else
        boxes{i} = NaN;
    end
end
pooled = cell2mat(cellfun(@(v) v(:), R.rateRatio_match(:), 'uni', false));
boxes{end} = pooled;
positions = 1:(nR+1);
for i = 1:numel(boxes)
    yi = boxes{i}(isfinite(boxes{i}));
    if isempty(yi), continue; end
    boxchart(i*ones(size(yi)), yi, 'BoxFaceAlpha',0.4, 'MarkerStyle','.');
end
plot(xlim, [1 1], 'k--');
set(gca,'XTick',positions,'XTickLabel',[ratNames(:); {'pooled'}], ...
    'XTickLabelRotation',30);
ylabel('Geo-mean rate ratio: trial / matched run');
title('Trial vs speed-matched non-task rate');
end
