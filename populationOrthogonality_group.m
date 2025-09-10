function G = populationOrthogonality_group(ratNames, varargin)
% populationOrthogonality_group
% Run populationOrthogonality for each rat, then build one group summary figure.
%
% Example:
%   G = populationOrthogonality_group( ...
%         {'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%         'NSplits',15,'Mode','zscore','Similarity','pearson', ...
%         'DoStats',true,'NBoot',500,'NPerm',500, ...
%         'GroupSaveFig',true,'GroupFigName','all_rats_group.png');

% ---------------- Parse group options; pass through the rest ----------------
ratNames = cellstr(ratNames);
p = inputParser;  p.KeepUnmatched = true;
addParameter(p,'GroupSaveFig',false,@islogical);
addParameter(p,'GroupFigName','',@(s) ischar(s) || isstring(s));
parse(p,varargin{:});
grpSave = p.Results.GroupSaveFig;
grpName = string(p.Results.GroupFigName);

% Pass-through args for populationOrthogonality
un = p.Unmatched;
unNames  = fieldnames(un);
unValues = struct2cell(un);
passArgs = reshape([unNames.'; unValues.'], 1, []);  % name/value alternating

% ---------------- Per-rat analysis (produces your per-rat figs) ----------------
nR = numel(ratNames);
perRat = cell(nR,1);
for k = 1:nR
    perRat{k} = populationOrthogonality(ratNames{k}, passArgs{:});
end

% ---------------- Collect per-rat pieces ----------------
nSplits = perRat{1}.params.NSplits;
for k = 2:nR
    assert(perRat{k}.params.NSplits == nSplits, 'All rats must use the same NSplits.');
end

simStr = lower(perRat{1}.params.Similarity);
if strcmp(simStr,'cosine'), clim = [0 1]; else, clim = [-1 1]; end

simMats = nan(nSplits,nSplits,nR);
lagMat  = nan(nR, nSplits);
within  = nan(nR,1);
between = nan(nR,1);
deltaWB = nan(nR,1);
mantelR = nan(nR,1);
rsaSame = nan(nR,1); rsaSameCI = nan(nR,2);
rsaProx = nan(nR,1); rsaProxCI = nan(nR,2);

for k = 1:nR
    simMats(:,:,k) = perRat{k}.pooled.simMat;

    if isfield(perRat{k},'stats') && isfield(perRat{k}.stats,'lag') && isfield(perRat{k}.stats.lag,'mean')
        lagMat(k,:) = perRat{k}.stats.lag.mean;
    end

    if isfield(perRat{k}.stats,'epoch')
        within(k)  = perRat{k}.stats.epoch.within_mean;
        between(k) = perRat{k}.stats.epoch.between_mean;
        deltaWB(k) = perRat{k}.stats.epoch.diff;
    end

    if isfield(perRat{k}.stats,'mantel') && isfield(perRat{k}.stats.mantel,'r')
        mantelR(k) = perRat{k}.stats.mantel.r;
    end

    if isfield(perRat{k}.stats,'rsa')
        rS = perRat{k}.stats.rsa;

        if isfield(rS,'beta_same'), rsaSame(k) = rS.beta_same; end
        if isfield(rS,'beta_prox'), rsaProx(k) = rS.beta_prox; end

        % Accept 1×2 or 2×2 CI formats
        if isfield(rS,'beta_same_CI')
            ci = rS.beta_same_CI;
            if isequal(size(ci),[2 2])
                rsaSameCI(k,:) = ci(:,1).';     % take column 1 (same)
            else
                rsaSameCI(k,:) = ci(:).';       % 1×2 or 2×1 -> row
            end
        end
        if isfield(rS,'beta_prox_CI')
            ci = rS.beta_prox_CI;
            if isequal(size(ci),[2 2])
                rsaProxCI(k,:) = ci(:,2).';     % take column 2 (prox)
            else
                rsaProxCI(k,:) = ci(:).';
            end
        end
    end
end

simMean = mean(simMats,3,'omitnan');

% Across-rat within > between (rat is the unit)
pWB = NaN;

    [~,pWB] = ttest(within, between);


% ---------------- Group figure ----------------
nColsHeat = min(5, nR+1);      % heatmaps per row (incl. mean)
nRowsHeat = ceil((nR+1)/nColsHeat);
nCols     = nColsHeat;
nRows     = nRowsHeat + 2;     % +2 rows for lag + within/between

fig = figure('Color','w','Position',[80 80 320*nCols 260*nRows]);
tl  = tiledlayout(nRows, nCols, 'Padding','compact','TileSpacing','compact');

% (1) Heatmaps per rat
%for k = 1:nR
%    nexttile; imagesc(simMats(:,:,k), clim); axis square
%    colormap(parula); colorbar;
%    title(ratNames{k},'Interpreter','none');
%    xticks(1:nSplits); yticks(1:nSplits); xlabel('Split'); ylabel('Split');
%end

% (2) Across-rat mean heatmap
nexttile; imagesc(simMean, clim); axis square
colormap(parula); colorbar; title('Across-rat mean');
xticks(1:nSplits); yticks(1:nSplits); xlabel('Split'); ylabel('Split');

% (3) Lag overlay (span whole row)
nexttile([1 2]); hold on
cols = lines(nR);
x = 0:(nSplits-1);
hLag = gobjects(1,nR);
for k = 1:nR
    if all(isfinite(lagMat(k,:)))
        hLag(k) = plot(x, lagMat(k,:), '-', 'Color', cols(k,:), 'LineWidth', 1.2);
    end
end
% across-rat mean ± SEM
mLag  = mean(lagMat, 1, 'omitnan');
sLag  = std(lagMat, 0, 1, 'omitnan') ./ max(1, sqrt(sum(isfinite(lagMat),1)));
fill([x fliplr(x)], [mLag-sLag fliplr(mLag+sLag)], [0 0 0], ...
     'FaceAlpha',0.08, 'EdgeColor','none');
plot(x, mLag, 'k-', 'LineWidth', 2);
xlabel('Lag (bins)'); ylabel(sprintf('Mean %s similarity', perRat{1}.params.Similarity));
title('Lag curves (per rat & mean ± SEM)'); box off
if any(isgraphics(hLag))
    legend(hLag(isgraphics(hLag)), ratNames(isgraphics(hLag)), 'Location','eastoutside');
end

% (4) Across-rat within vs between + paired dots (span ~half row)
nexttile([1 ceil(nCols/2)]); hold on
muW = mean(within,'omitnan');  muB = mean(between,'omitnan');
seW = std(within,'omitnan')/sqrt(max(1,sum(isfinite(within))));
seB = std(between,'omitnan')/sqrt(max(1,sum(isfinite(between))));
bar([1 2], [muW muB], 0.6, 'FaceColor',[.72 .76 .92]);
errorbar([1 2], [muW muB], [seW seB], 'k.', 'LineWidth',1.2);
set(gca,'XTick',[1 2],'XTickLabel',{'Within','Between'});
ylabel(sprintf('Mean %s similarity', perRat{1}.params.Similarity));
title(sprintf('Across-rat: Within > Between?  p = %.3g', pWB)); box off
% per-rat dots + paired lines
scatter(ones(nR,1)*0.85, within, 28, cols, 'filled', 'MarkerFaceAlpha',0.7);
scatter(ones(nR,1)*2.15, between, 28, cols, 'filled', 'MarkerFaceAlpha',0.7);
for k=1:nR
    if isfinite(within(k)) && isfinite(between(k))
        plot([0.85 2.15],[within(k) between(k)],'-','Color',cols(k,:));
    end
end

% (5) Text panel with Mantel & RSA summary (span remaining columns)
axTxt = nexttile([1 max(1, nCols - ceil(nCols/2))]);
axis(axTxt,'off'); xlim([0 1]); ylim([0 1]);
txt = sprintf(['Mantel r (across rats): mean=%.2f, SD=%.2f\n' ...
               'RSA  \\beta_{prox}:      mean=%.3f, SD=%.3f\n' ...
               'RSA  \\beta_{same}:      mean=%.3f, SD=%.3f'], ...
               mean(mantelR,'omitnan'), std(mantelR,'omitnan'), ...
               mean(rsaProx,'omitnan'), std(rsaProx,'omitnan'), ...
               mean(rsaSame,'omitnan'), std(rsaSame,'omitnan'));
text(0,0.9, txt, 'FontSize',11, 'Interpreter','tex');

title(tl,'Population orthogonality: group summary');

% Save if requested
if grpSave
    if strlength(grpName)==0, grpName = "population_orthogonality_group.png"; end
    exportgraphics(fig, grpName, 'Resolution', 300);
end

% ---------------- Pack outputs ----------------
G = struct;
G.perRat   = perRat;
G.simMats  = simMats;
G.simMean  = simMean;
G.lagMat   = lagMat;
G.within   = within;
G.between  = between;
G.delta    = deltaWB;
G.mantelR  = mantelR;
G.rsaSame  = rsaSame;  G.rsaSameCI = rsaSameCI;
G.rsaProx  = rsaProx;  G.rsaProxCI = rsaProxCI;
G.fig      = fig;
end
