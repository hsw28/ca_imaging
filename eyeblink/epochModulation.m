function epochModulation(ratNames, k)
% epochModulation Build 4-epoch modulation indices and cluster them
%%%%% NOT THE MAIN WAY TO FIND OUT SIG
%
%   epochModulationFingerprints({'rat0222', ...}, 4)
%
% INPUTS
%   ratNames : cell array of rat variable names (in base workspace)
%   k        : k-means cluster count (default 4)
%
%% ex: epochModulation({'rat0222','rat0307','rat0313','rat0314','rat0816'})

if nargin < 2, k = 4; end

% ---- epoch boundaries relative to CS onset ----------------------------
edges = [-2 0     0.25 0.75 0.85 2.00];   % 5 epochs → 4 modulation scores
nEpoch = 5;

allIdx   = [];          % [nCellsTotal × 4] modulation matrix
cellRat  = [];          % rat ID per cell (for prevalence plot)
labels   = ["CS" "Trace" "US" "Post"];

for r = 1:numel(ratNames)
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);

    for d = 1:3
        dayStr   = days{d};
        spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',dayStr));
        csTimes  = rat.CS_times.(sprintf('CS_%s',dayStr));
        maskCell = rat.ratemask.(sprintf('ratemask_%s',dayStr)) == 1;

        nCells   = size(spk,1);

        % --- NEW: per-day epoch vector, default NaN for cells NOT analyzed -------
          epochVec = nan(nCells,1);   % NaN = not analyzed

          for c = 1:nCells
              if ~maskCell(c), continue, end

              st = spk(c,:);
              st = st(~isnan(st) & st>0);
              if isempty(st), continue, end

              % ----- compute firing rates per epoch -----
              fr = zeros(1,numel(edges)-1);
              for e = 1:numel(edges)-1
                  cnt = 0;
                  for t = 1:numel(csTimes)
                      t0  = csTimes(t) + edges(e);
                      t1  = csTimes(t) + edges(e+1);
                      cnt = cnt + sum(st >= t0 & st < t1);
                  end
                  fr(e) = cnt / ((edges(e+1)-edges(e)) * numel(csTimes)); % Hz
              end

              baseFR = fr(1);
              if baseFR == 0
                  continue   % remains NaN
              end

              % ----- Δ-modulation -----
              modIdx = (fr(2:end) - baseFR) ./ (fr(2:end) + baseFR);   % 1×4

              allIdx  = [allIdx ; modIdx];
              cellRat = [cellRat;  r];

              % =======================================================
              % NEW LOGIC FOR EPOCH LABEL
              % =======================================================
              thr = 0.05;                       % threshold for modulation
              [mx, eDom] = max(abs(modIdx));    % best epoch
              if mx < thr
                  epochVec(c) = 0;              % analyzed but NOT modulated
              else
                  epochVec(c) = eDom;           % 1=CS,2=Trace,3=US,4=Post
              end
          end


    end

    % --- NEW: write updated rat struct back to base workspace --------------
    assignin('base', ratNames{r}, rat);
end

% ---------- k-means clustering & deterministic ordering -----------------
rng default                                     % reproducible
[clusterIdx,C] = kmeans(allIdx,k,'Replicates',10);   % C : k × 4

% =======================================================================
% Build row order so heat-map shows: CS-dominant block first, then Trace,
% then US, then Post.  Inside each block rows go from strong suppression
% (Δ ≈ –1, blue) to strong excitation (Δ ≈ +1, yellow).
% -----------------------------------------------------------------------

% 1) dominant epoch for each cell (1 = CS, 2 = Trace, 3 = US, 4 = Post)
[~, domEpoch] = max(abs(allIdx), [], 2);

ord = [];
for e = 1:4                                   % 1=CS 2=Trace 3=US 4=Post
    rows = find(domEpoch == e);               % cells whose dominant epoch = e
    if isempty(rows), continue, end

    % signed sort: +1 (yellow) → 0 (green) → –1 (blue)
    [~, sub] = sort(allIdx(rows, e), 'descend');

    ord = [ord; rows(sub)];                   %#ok<AGROW>
end



% ------ heat-map --------------------------------------------------------
figure('Color','w');
subplot(1,4,1)
imagesc(allIdx(ord ,:)); colormap(parula); caxis([-1 1])
xticks(1:4); xticklabels({'CS','Trace','US','Post'});
ylabel('Cells (CS → Trace → US → Post)'); title('Modulation fingerprints');
colorbar;


% ---------- (panel-2) centroid Δ modulation with error bars ------------
subplot(1,4,2); cla; hold on
cols = lines(k);                      % one colour per cluster
x    = 1:4;                           % CS  Trace  US  Post

for j = 1:k
    rows = (clusterIdx == j);         % cells belonging to this cluster

    mu  = C(j,:);                     % 1 × 4  (centroid means already)
    sem = std(allIdx(rows,:),0,1,'omitnan') ./ sqrt(sum(rows));

    % mean line with circular markers
    plot(x, mu, '-o', ...
         'Color', cols(j,:), ...
         'MarkerFaceColor', cols(j,:), ...
         'LineWidth', 1.8);

    % vertical error bars (± SEM)
  %  errorbar(x, mu, sem, ...
  %           'Color', cols(j,:), ...
  %           'LineStyle','none', ...
  %           'LineWidth',1.2);
end


yline(0,'k');                         % baseline at zero
xlim([0.7 4.3]); ylim([-1 1]);
xticks(x); xticklabels({'CS','Trace','US','Post'});
ylabel('Δ modulation (mean ± SEM)');
title('Cluster centroids');
box off
legend("C1","C2","C3","C4",'Location','eastoutside');



% (c) Cluster prevalence per rat
subplot(1,4,3);
nRats = numel(ratNames);
prev  = zeros(nRats, k);

for r = 1:nRats
    for j = 1:k
        prev(r,j) = mean( clusterIdx(cellRat==r) == j );
    end
end

% 1) across-rat mean (each rat = 1 datum)
prevMean = mean(prev,1,'omitnan');        % 1 × k

% 2) append as extra row
prevPlot = [prev; prevMean];              % (nRats+1) × k

% 3) labels
xLabs = [ratNames, {'All'}];

% 4) draw
subplot(1,4,3); cla
bh = bar(prevPlot,'stacked');             % stacked bars
xlim([0.5  (size(prevPlot,1)+0.5)])
xticks(1:numel(xLabs)); xticklabels(xLabs)
ylabel('Fraction of cells');
title('Cluster prevalence');
box off
legend("C1","C2","C3","C4",'Location','EastOutside');

sgtitle('4-epoch modulation and cell-type clustering');


% -------- build +/– counts per cluster & epoch ------------------------
posFrac = zeros(k,4); negFrac = zeros(k,4);
for j = 1:k
    rows = (clusterIdx == j);
    posFrac(j,:) = mean(allIdx(rows,:) >  0, 1, 'omitnan');
    negFrac(j,:) = mean(allIdx(rows,:) <  0, 1, 'omitnan');
end

% -------- plot --------------------------------------------------------

% ---------------- panel-4 :  + vs – within each epoch -------------------
subplot(1,4,4); cla; hold on

epochCols = lines(4);                    % colours for CS,Trace,US,Post
epochLbls = {'CS','Trace','US','Post'};

% dominant epoch of *each cell* (1‒4)
[~, domEpoch] = max(abs(allIdx), [], 2); % allIdx is nCells × 4

posFrac = zeros(1,4);
negFrac = zeros(1,4);

for e = 1:4                               % loop over epochs
    rows =  domEpoch == e;                % cells whose *own* epoch is e
    v    =  allIdx(rows, e);              % their Δ-modulation in epoch e
    v    =  v(~isnan(v));                 % drop NaNs

    if ~isempty(v)
        posFrac(e) = mean(v > 0);
        negFrac(e) = mean(v < 0);
    end
end

x = 1:4;                                  % CS … Post on x-axis

% positive bars
bhPos        = bar(x,  posFrac, 0.6, ...
                   'FaceColor','flat', 'EdgeColor','none');
bhPos.CData  = epochCols;                 % one RGB triplet per bar

% negative bars (transparent & below zero)
bhNeg        = bar(x, -negFrac, 0.6, ...
                   'FaceColor','flat', 'EdgeColor','none', ...
                   'FaceAlpha',0.30);
bhNeg.CData  = epochCols;                 % same colours

% axis cosmetics
yline(0,'k');
xticks(x); xticklabels(epochLbls);
ylim([-1 1]); ylabel('Fraction of cells (↑ pos, ↓ neg)');

title('+ vs – in each cell’s *own* epoch');
box off

% legend swatches: opaque = positive, translucent = negative
patch(NaN,NaN,epochCols(1,:),'EdgeColor','none');      % dummy for legend
patch(NaN,NaN,epochCols(1,:),'EdgeColor','none','FaceAlpha',0.30);
legend({'Positive','Negative'},'Location','eastoutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRINT OUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =========================================================
%  (a)  Centroid table  (panel-2)
% =========================================================
fprintf('\n-- Cluster-centroid Δ modulation (mean ± SEM) --\n');
hdr = '%7s  %8s  %8s  %8s  %8s\n';
fprintf(hdr, 'Clust', 'CS', 'Trace', 'US', 'Post');

for j = 1:k
    rows = (clusterIdx == j);
    mu   = mean(allIdx(rows,:), 1, 'omitnan');
    sem  = std(allIdx(rows,:), 0, 1, 'omitnan') ./ sqrt(sum(rows));
    fprintf(hdr, sprintf('C%d',j), ...
            sprintf('%.3f±%.3f',mu(1),sem(1)), ...
            sprintf('%.3f±%.3f',mu(2),sem(2)), ...
            sprintf('%.3f±%.3f',mu(3),sem(3)), ...
            sprintf('%.3f±%.3f',mu(4),sem(4)));
end

% =========================================================
%  (b)  Cluster-prevalence Chi² across rats  (panel-3)
% =========================================================
% -------- percentages of + / – cells in each cluster -------------------
labels = {'CS+','CS-','Trace+','Trace-','US+','US-','Post+','Post-'};
fprintf('\n%% of cells in each sign/epoch category:\n');
hdr = '%6s';
for j = 1:numel(labels), hdr = [hdr '  %7s']; end
fprintf([hdr '\n'], 'Cluster', labels{:});

for c = 1:k
    rows = (clusterIdx == c);
    tab  = [ mean(allIdx(rows,1)>0)  mean(allIdx(rows,1)<0) ...
             mean(allIdx(rows,2)>0)  mean(allIdx(rows,2)<0) ...
             mean(allIdx(rows,3)>0)  mean(allIdx(rows,3)<0) ...
             mean(allIdx(rows,4)>0)  mean(allIdx(rows,4)<0) ] * 100;
    fmt  = '%6s';
    for j = 1:numel(tab), fmt = [fmt '  %6.1f']; end
    fprintf([fmt '\n'], ['C' num2str(c)], tab);
end


fprintf('\nCluster   Epoch   %%Up   %%Down\n');
for j = 1:k
    fprintf('   C%-2d     %-5s  %5.1f  %5.1f\n', ...
        j, labels{domEpoch(j)}, posFrac(j)*100, negFrac(j)*100);
end
