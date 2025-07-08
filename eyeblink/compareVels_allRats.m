function compareVels_allRats(ratNames,winSecs)
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

nRats      = numel(ratNames);
velInCell  = cell(nRats,1);
velOutCell = cell(nRats,1);
meanIn     = nan(nRats,1);
meanOut    = nan(nRats,1);
pvals      = nan(nRats,1);

% ---------- gather data ---------------------------------------------------
for r = 1:nRats
    [mIn,mOut,p,vin,vout] = compareVels_perRat(ratNames{r},winSecs);
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

% ---------- grand plot ----------------------------------------------------
figure('Color','w','Position',[200 250 700 480]);
boxplot(allVals,allGroups,'PlotStyle','traditional','Symbol','k.');
ylabel('Velocity (cm/s)');
title('Velocity inside vs. outside CS window – all rats');

% print a tiny table in command window
fprintf('\n%-8s  meanOut  meanIn   p-value\n', 'Rat');
for r = 1:nRats
    fprintf('%-8s  %7.2f  %7.2f  %g\n', ratNames{r}, meanOut(r), meanIn(r), pvals(r));
end
end
