function SpikesPerBin(ratNames)
% For each rat and each of its last 3 analysed days,
% build a 15×nCells spike-count matrix (2 s after CS, 0.133 s bins)
% and store it in rat.spikesperbin.cpb_<date>.
% ex: SpikesPerBin({'rat0222','rat0307','rat0313','rat0314','rat0816'})

%win          = [0 2];         % 2-s CS-trace window
%nBins        = 15;
%binEdges     = linspace(win(1), win(2), nBins+1);   % 0 : 0.133 : 2.0
%binWidth     = diff(binEdges(1:2));

winFull   = [-1 2];                 % 3-second window
binWidth  = (1/7.5);               % keep the 0–2 s resolution
binEdges  = winFull(1):binWidth:winFull(2);   % edges –1 : 0.133 : 2
nBins     = numel(binEdges)-1;      % 23 bins (7 pre-CS, 15 post-CS, 1 overlap)

for r = 1:numel(ratNames)
    rat   = evalin('base', ratNames{r});          % pull from workspace
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);                     % day-2, day-1, day-0

    for d = 1:3
        dayStr   = days{d};
        spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',dayStr));   % n×T
        csTimes  = rat.CS_times.(sprintf('CS_%s',dayStr));         % trials
        maskCell = rat.ratemask.(sprintf('ratemask_%s',dayStr));   % 1×n

        % --- sizes ----------------------------------------------------------
        ratemask = (rat.ratemask.(sprintf('ratemask_%s',dayStr)) == 1);  % 1×nTotal

        % --- sizes ----------------------------------------------------------
        nTotal   = size(spk,1);                       % total cells for that day
        inclMask = ratemask == 1;                     % logical row-vector

        nBins    = nBins;
        spikesperbin = zeros(nBins, nTotal);          % ←  start at 0, not NaN
        inclIdx  = find(inclMask);                    % indices of cells we keep

        % ---------- accumulate spikes ---------------------------------------
        for t = 1:numel(csTimes)
            t0 = csTimes(t);                          % CS onset
            for jj = 1:numel(inclIdx)
                c   = inclIdx(jj);
                st  = spk(c,:);  st = st(~isnan(st) & st>0);
                rel = st - t0;
                rel = rel(rel>=0 & rel<2);            % 0–2 s after CS
                if isempty(rel),  continue, end
                cnt = histcounts(rel, binEdges);      % 15-bin histogram
                spikesperbin(:,c) = spikesperbin(:,c) + cnt(:);
              end
            end

            % ---------- convert counts → Hz  ------------------------------------
            spikesperbin = spikesperbin / (binWidth * numel(csTimes));

        % ---------- blank out the excluded columns --------------------------
        spikesperbin(:, ~inclMask) = NaN;




        % store back into the rat struct
        if ~isfield(rat,'spikesperbin'), rat.spikesperbin = struct(); end
        rat.spikesperbin.(sprintf('cpb_%s',dayStr)) = spikesperbin';
    end

    % push updated struct back to base workspace
    assignin('base', ratNames{r}, rat);
    fprintf('Stored spikesperbin for %s\n', ratNames{r});

end  % ← closes the outer loop over rats

% ---------- figure: mean ± SEM with pre-CS ---------------------------
nRats   = numel(ratNames);
tAxis   = binEdges(1:end-1) + binWidth/2;   % bin centres (-1 … 1.933)

figure('Color','w','Position',[200 200 1200 300]);

grandMat = [];                             % collect for all-rats panel

for rr = 1:nRats
    rat   = evalin('base', ratNames{rr});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);

    M = [];                                % concat three days
    for d = 1:3
        M = [M; rat.spikesperbin.(sprintf('cpb_%s',days{d}))]; %#ok<AGROW>
    end

    mu  = nanmean(M,1);                    % 1×23
    sem = nanstd(M,0,1) ./ sqrt(sum(~isnan(M),1));

    subplot(1,nRats+1,rr); hold on;

    % grey shading of the trace interval (0 – 2 s)
    fill([0 2 2 0], [0 0 max(mu+sem)*1.1 max(mu+sem)*1.1], ...
         [.9 .9 .9], 'EdgeColor','none');

    bar(tAxis, mu, 1, 'FaceColor',[.5 .5 .8], 'EdgeColor','none');
    errorbar(tAxis, mu, sem, '.k', 'CapSize',4);

    line([0 0], ylim, 'Color','r','LineStyle','--');  % CS marker
    line([.75 .75], ylim, 'Color','blue','LineStyle','--');  % US marker

    xlabel('Time from CS (s)'); ylabel('Spikes / bin');
    title(ratNames{rr}); box off;

    grandMat = [grandMat; mu];              %#ok<AGROW>
end

% ----------- combined panel ------------------------------------------
grandMu  = nanmean(grandMat,1);
grandSem = nanstd (grandMat,0,1) ./ sqrt(sum(~isnan(grandMat),1));

subplot(1,nRats+1,nRats+1); hold on;
fill([0 2 2 0], [0 0 max(grandMu+grandSem)*1.1 max(grandMu+grandSem)*1.1], ...
     [.9 .9 .9], 'EdgeColor','none');
bar(tAxis, grandMu, 1, 'FaceColor',[.8 .5 .5], 'EdgeColor','none');
errorbar(tAxis, grandMu, grandSem, '.k', 'CapSize',4);
line([0 0], ylim, 'Color','r','LineStyle','--');
line([.75 .75], ylim, 'Color','blue','LineStyle','--');  % US marker

xlabel('Time from CS (s)'); ylabel('Spikes / bin');
title('All rats'); box off;

sgtitle('Mean spikes per 0.133-s bin  (-1 → +2 s, CS at 0)');

end
