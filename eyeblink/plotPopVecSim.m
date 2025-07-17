function plotPopVecSim()
% plotPopVecSim  Population‐vector similarity: Trial vs Random Windows
%
%   plotPopVecSim()
%
% Computes and compares two similarity distributions across rats and days:
%   • TT: trial‐window vs. trial‐window (conditioning periods)
%   • RR: random non‐trace windows vs. random non‐trace windows
%
% Hypothesis: TT similarities are higher, reflecting population responsiveness.
%
% Expects in base workspace:
%   ratXXX.Ca_peaks.CA_peaks_<date>   = {nCells×1} cell array of spike times
%   ratXXX.pos.pos_<date>             = [nTime×3] [time, x, y]
%   ratXXX.CS_times.CS_<date>         = CS onsets
%   autoDateList(ratXXX)              = function returning dates

    %% USER PARAMETERS
    ratNames  = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
    win       = [0 2];    % CS window relative to CS onset (s)
    minSpikes = 5;        % min total spikes in all CS windows per cell

    TT_all = [];
    RR_all = [];

    for r = 1:numel(ratNames)
        rat   = evalin('base', ratNames{r});
        dates = autoDateList(rat);
        idx   = find(strcmp(dates, rat.An),1);
        days  = dates(idx-2:idx);  % day−2, day−1, day−0

        for d = 1:3
            day    = days{d};
            spk    = rat.Ca_peaks.(['CA_peaks_' day]);
            posMat = rat.pos.     (['pos_'      day]);
            ts     = posMat(:,1);
            cs     = rat.CS_times.(['CS_' day]);

            % pass minSpikes into helper
            [FRt, FRr] = popVecSim(spk, ts, cs, win, minSpikes);

            % z-score normalize across cells
            muT = mean(FRt,1);   sT = std(FRt,0,1);  sT(sT==0)=1;
            FRt = (FRt - muT) ./ sT;
            muR = mean(FRr,1);   sR = std(FRr,0,1);  sR(sR==0)=1;
            FRr = (FRr - muR) ./ sR;

            % compute cosine similarity matrices
            nTrials = size(FRt,1);
            normT   = vecnorm(FRt,2,2);
            Ctt     = (FRt * FRt') ./ (normT * normT');
            normR   = vecnorm(FRr,2,2);
            Crr     = (FRr * FRr') ./ (normR * normR');

            % extract upper-triangle values
            idxUT   = triu(true(nTrials),1);
            TT_all  = [TT_all; Ctt(idxUT)];
            RR_all  = [RR_all; Crr(idxUT)];
        end
    end

    % Plot distributions
    figure('Color','w','Position',[200 200 800 400]);
    boxplot([TT_all; RR_all], ...
            [ones(size(TT_all)); 2*ones(size(RR_all))], ...
            'Labels', {'TT','RR'});
    ylabel('Cosine similarity');
    title('Population‐vector similarity: Trial vs Random');

    % Console summary
    fprintf('\n=== Similarity summary (mean ± SD) ===\n');
    fprintf('TT: %.3f ± %.3f\n', mean(TT_all,'omitnan'), std(TT_all,'omitnan'));
    fprintf('RR: %.3f ± %.3f\n', mean(RR_all,'omitnan'), std(RR_all,'omitnan'));
end

%%---------------------------------------------------------------------------%%
function [FRt,FRr] = popVecSim(spikeCell, ts, csTimes, win, minSpikes)
% popVecSim  Build FR matrices for CS vs. random windows
%
%   [FRt,FRr] = popVecSim(spikeCell, ts, csTimes, win, minSpikes)
%   - spikeCell: {nCells×1} cell array of spike times, or
%                numeric [nCells×nS] spike-time matrix (zeros/NaNs padded)
%   - ts:        [nBins×1] timestamps (s)
%   - csTimes:   [nTrials×1] CS onset times (s)
%   - win:       [t0 t1] CS window relative to CS onset (s)
%   - minSpikes: filter threshold on total spikes per cell in CS windows

    % default minSpikes
    if nargin<5, minSpikes = 5; end

    % 1) Convert numeric matrix to cell array of spike times
    if ~iscell(spikeCell)
        [nCells, ~] = size(spikeCell);
        tmp = cell(nCells,1);
        for c = 1:nCells
            st = spikeCell(c,:);
            st = st(~isnan(st) & st>0);  % drop padding
            tmp{c} = st(:);
        end
        spikeCell = tmp;
    end

    % 2) parameters
    nTrials   = numel(csTimes);
    nCells    = numel(spikeCell);
    windowLen = diff(win);
    dt        = median(diff(ts));

    % 3) build mask for all CS windows on ts
    maskCS = false(size(ts));
    for t = 1:nTrials
        maskCS = maskCS | (ts>=csTimes(t)+win(1) & ts<csTimes(t)+win(2));
    end
    outs = ts(~maskCS);  % non‐CS timepoints

    % 4) FR in each CS window
    FRt = nan(nTrials, nCells);
    for t = 1:nTrials
        t0 = csTimes(t) + win(1);
        t1 = t0 + windowLen;
        for c = 1:nCells
            st = spikeCell{c};
            FRt(t,c) = sum(st>=t0 & st<t1) / windowLen;
        end
    end

    % 5) FR in random non‐CS windows
    bins     = round(windowLen / dt);
    maxStart = numel(outs) - bins + 1;
    starts   = randi(maxStart, nTrials, 1);

    FRr = nan(nTrials, nCells);
    for t = 1:nTrials
        t0 = outs(starts(t));
        t1 = t0 + windowLen;
        for c = 1:nCells
            st = spikeCell{c};
            FRr(t,c) = sum(st>=t0 & st<t1) / windowLen;
        end
    end

    % 6) filter cells by total spikes in CS windows
    totalSpikes = sum(FRt,1) * windowLen;  % back to counts
    keepCells   = totalSpikes >= minSpikes;
    FRt = FRt(:, keepCells);
    FRr = FRr(:, keepCells);
end
