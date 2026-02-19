function compareVels_allRats(ratNames,winSecs, type)
% compareVels_allRats  Call compareVels_perRat for several animals and plot together
%
%   compareVels_allRats({'rat0222','rat0307',...},[0 2])
%
%   * ratNames – cell array of variable names (strings/chars) that exist in base workspace
%   * winSecs  – 1×2 vector, window relative to CS onset (default [0 2])
%
%   Produces a grand box-plot with all “in–trial” and “out-of-trial” velocity
%   samples for every rat (colour-coded), plus prints mean ± SD for each group.

if nargin<1, error('Pass a cell array of rat names'); end
if nargin<2, winSecs = [0 2]; end
if nargin<3
  type=1 %compares to all non trial times
  type = 2; %compres to prev 2 sec
end

figure('Color','w');

nRats      = numel(ratNames);
velInCell  = cell(nRats,1);
velOutCell = cell(nRats,1);
meanIn     = nan(nRats,1);
meanOut    = nan(nRats,1);
pvals      = nan(nRats,1);

% ---------- gather data ---------------------------------------------------
for r = 1:nRats
    if type ==1
      [mIn,mOut,p,vin,vout] = compareVels_perRat(ratNames{r},winSecs);
      fprintf('type1')
    else
      [mIn,mOut,p,vin,vout] = compareVels_perRat2(ratNames{r},winSecs);
            fprintf('type2')
    end
    velInCell{r}  = vin;
    velOutCell{r} = vout;
    meanIn(r)     = mIn;
    meanOut(r)    = mOut;
    pvals(r)      = p;
end

% ---------- build long vectors for boxplot --------------------------------
allVals   = [];
allGroups = strings(0);

for r = 1:nRats
    thisRat = ratNames{r};
    allVals   = [allVals ; velOutCell{r} ; velInCell{r}]; %#ok<AGROW>
    allGroups = [allGroups ;
                 repmat(thisRat + " Out" , numel(velOutCell{r}),1) ;
                 repmat(thisRat + " In"  , numel(velInCell{r}) ,1)]; %#ok<AGROW>
end

% ---------- grand BOX PLOT (original) -----------------------------
figure('Color','w','Position',[200 250 700 480]);
boxplot(allVals,allGroups,'PlotStyle','traditional','Symbol','k.');
ylabel('Velocity (cm/s)');
title('Velocity inside vs. outside CS window – all rats (boxplot)');

% ---------- grand VIOLIN PLOT (new; consistent capping) ------------------
figure('Color','w','Position',[200 250 1000 560]); hold on
set(gca,'Clipping','on');
ylabel('Velocity (cm/s)');
title('Velocity inside vs. outside CS window – all rats (violin)');

% Colors: same for all animals
COL_OUT = [0.85 0.33 0.10];
COL_IN  = [0.10 0.45 0.90];

% Layout
sep=1.6; offset=0.25; halfW=0.35; nbins=128; alphaFill=0.35;
xRat = (1:nRats)*sep;

% ----------------- choose consistent cap across figures ------------------
Cap.mode   = 'pct';      % 'pct' or 'abs'
Cap.value  = 99;         % percentile if 'pct' (e.g., 95/97/99) or numeric if 'abs'
Cap.source = 'in';       % 'in' | 'out' | 'all'

% Pool for cap
poolIn  = cell2mat(cellfun(@(v) v(:), velInCell,  'UniformOutput',false));
poolOut = cell2mat(cellfun(@(v) v(:), velOutCell, 'UniformOutput',false));
switch lower(Cap.source)
    case 'in',  pool = poolIn;
    case 'out', pool = poolOut;
    otherwise,  pool = [poolIn; poolOut];
end
pool = pool(isfinite(pool));

% Compute y-range once and use it for all violins in this figure
if strcmpi(Cap.mode,'pct')
    yMax = prctile(pool, Cap.value);
else
    yMax = Cap.value;
end
yMin = max(0, min(pool));     % adjust if negative velocities possible
ylim([yMin yMax]);
% -------------------------------------------------------------------------

for r = 1:nRats
    gOut = [ratNames{r} ' Out'];  yOut = allVals(strcmp(allGroups,gOut));
    gIn  = [ratNames{r} ' In' ];  yIn  = allVals(strcmp(allGroups,gIn));

    xOut = xRat(r) - offset;
    xIn  = xRat(r) + offset;

    drawViolin(yOut, xOut, halfW, COL_OUT, nbins, alphaFill, [yMin yMax]);
    drawSummary(yOut, xOut, [yMin yMax]);

    drawViolin(yIn,  xIn,  halfW, COL_IN,  nbins, alphaFill, [yMin yMax]);
    drawSummary(yIn,  xIn,  [yMin yMax]);
end

xlim([xRat(1)-sep*0.7, xRat(end)+sep*0.7]);
set(gca,'XTick',xRat,'XTickLabel',ratNames,'TickDir','out','Box','off');
legend({'mean','median'},'Location','northeastoutside');


for r = 1:nRats
    fprintf('%-8s  %7.2f  %7.2f  %g\n', ratNames{r}, meanOut(r), meanIn(r), pvals(r));
end
end

% ================= helper: draw one violin =================
function drawViolin(y, xCenter, halfWidth, colorRGB, nbins, alphaFill, yRange)
    y = y(isfinite(y));
    if numel(y)<2
        plot([xCenter xCenter], [min(y) max(y)], 'Color', colorRGB, 'LineWidth', 1.2, 'Clipping','on');
        return
    end
    % limit density to display range (percentile/absolute cap)
    lo = max(min(y), yRange(1));
    hi = min(max(y), yRange(2));
    if hi <= lo, hi = lo + eps; end
    yGrid = linspace(lo, hi, nbins);
    [f, yi] = ksdensity(y, yGrid, 'Function','pdf');   % density from full data
    if max(f) > 0
        f = f / max(f) * halfWidth;
    end
    X = [xCenter - f, fliplr(xCenter + f)];
    Y = [yi, fliplr(yi)];
    patch('XData',X,'YData',Y,'FaceColor',colorRGB,'FaceAlpha',alphaFill,...
          'EdgeColor',colorRGB*0.6,'LineWidth',1.0,'Clipping','on');
end

% ============== helper: draw mean & median markers ==============

function drawSummary(y, xCenter, yRange)
    y = y(isfinite(y));
    mMed  = median(y);
    mMean = mean(y);
    % clamp markers to visible range so they don't draw outside the axes
    mMedP  = min(max(mMed,  yRange(1)), yRange(2));
    mMeanP = min(max(mMean, yRange(1)), yRange(2));
    plot([xCenter-0.08 xCenter+0.08], [mMeanP mMeanP], 'k-', 'LineWidth',1.2);          % mean (short black line)
    plot(xCenter, mMedP, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4);                   % median (black dot)
end



% ============== optional helper: jittered scatter ==============
function jitterScatter(y, xCenter, w, c)
    if isempty(y), return; end
    xj = xCenter + (rand(size(y))-0.5)*2*w;
    plot(xj, y, '.', 'Color', c, 'MarkerSize', 6);
end
