function slidingWindowLatency(ratNames, varargin)
% slidingWindowLatency  Per-rat and global latency summaries (12 panels)
%
%   slidingWindowLatency({'rat0222','rat0307',…})
%
% Optional name–value pairs:
%   'Win'   [t0 t1] analysis window  (default [-5 2]  s)
%   'Bin'   scalar  sliding window   (default 0.020 s  = 20 ms)
%   'Alpha' scalar  per-bin α        (default 0.01)

% ------------ parameters ------------------------------------------------
p = inputParser;
addParameter(p,'Win',[-5 2],@(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'Bin',0.133,@isscalar);
addParameter(p,'Alpha',0.01,@isscalar);
parse(p,varargin{:});
win   = p.Results.Win;
binW  = p.Results.Bin;
alpha = p.Results.Alpha;

edges = win(1):binW:win(2);                    % sliding-window edges
nBins = numel(edges)-1;

% -------- containers ----------------------------------------------------
% store latencies per rat so we can plot later
latPerRat   = cell(numel(ratNames)+1,1);       % +1 row for “All rats”
bandPerRat  = cell(numel(ratNames)+1,1);       % 1 = CS, 2 = Trace

fprintf('Scanning %d rats …\n', numel(ratNames));

for r = 1:numel(ratNames)
    rat   = evalin('base', ratNames{r});
    days  = autoDateList(rat);
    idx   = find(strcmp(days,rat.An),1);
    days  = days(idx-2:idx);                   % last 3 analysed days

    latencies = [];  whichBand = [];           % reset per rat

    for d = 1:3
        spk  = rat.Ca_peaks.(sprintf('CA_peaks_%s',days{d}));
        cs   = rat.CS_times.(sprintf('CS_%s',days{d}));
        keep = rat.ratemask.(sprintf('ratemask_%s',days{d})) == 1;
        if isempty(cs), continue, end

        for c = find(keep)'
            st = spk(c,:); st = st(~isnan(st)&st>0);
            if isempty(st), continue, end

            % ---------- FR matrix  (trials × bins) -----------------------
            FR = nan(numel(cs), nBins);
            for t = 1:numel(cs)
                for b = 1:nBins
                    FR(t,b) = sum(st >= cs(t)+edges(b) & st < cs(t)+edges(b+1)) / binW;
                end
            end

            % ---------- 5-s baseline (-5–0 s) ---------------------------
            baseIdx = edges >= -5 & edges < 0;
            baseFR  = mean(FR(:,baseIdx), 2, 'omitnan');     % per-trial baseline

            % ---------- paired t-test baseline vs each later bin ---------
            sig = false(1,nBins);
            for b = 1:nBins
                [~,p] = ttest(FR(:,b), baseFR, 'Alpha', alpha);
                sig(b) = p < alpha;
            end

            first = find(sig & edges(1:end-1) >= 0, 1, 'first');
            if isempty(first), continue, end

            latencies(end+1) = edges(first)*1000;            %#ok<AGROW>
            if edges(first) < 0.25
                whichBand(end+1) = 1;                        % CS
            elseif edges(first) < 0.75
                whichBand(end+1) = 2;                        % Trace
            else
                whichBand(end+1) = 0;                        % ignore
            end
        end
    end

    latPerRat{r}  = latencies;
    bandPerRat{r} = whichBand;
end

% ---------------- combine all rats into last “global” slot -------------
latPerRat{end}  = horzcat(latPerRat{:});
bandPerRat{end} = horzcat(bandPerRat{:});
ratLabels       = [ratNames, {'All'}];

% --------------------- PLOTS (6×2 grid) --------------------------------
nRows = numel(ratLabels);  % 6 (5 rats + global)
figure('Color','w','Position',[100 50 900 120*nRows]);

for rr = 1:nRows
    lat  = latPerRat{rr};
    band = bandPerRat{rr};

    % ---- Histogram panel (left) ----------------------------------------
    ax1 = subplot(nRows,2,2*rr-1); hold(ax1,'on')
    csLat = lat(band==1);
    trLat = lat(band==2);

    histogram(ax1, csLat, 30, 'BinLimits',[0 750], ...
              'FaceColor',[0.3 0.5 0.9], 'EdgeColor','none');
    histogram(ax1, trLat, 30, 'BinLimits',[0 750], ...
              'FaceColor',[0.9 0.6 0.2], 'EdgeColor','none','FaceAlpha',0.7);
    xline(ax1,250,'--k'); xline(ax1,750,'--k');
    xlabel(ax1,'Latency (ms)'); ylabel(ax1,'# Neurons');
    title(ax1, sprintf('%s: CS & Trace latencies', ratLabels{rr}));

    % ---- Cumulative fraction panel (right) -----------------------------
    ax2 = subplot(nRows,2,2*rr); hold(ax2,'on')
    edgesCDF = 0:binW:win(2);
    cdfVals  = arrayfun(@(x) mean(lat<=x*1000), edgesCDF);
    plot(ax2, edgesCDF*1000, cdfVals,'k-','LineWidth',1.8);
    xline(ax2,250,'--'); xline(ax2,750,'--');
    ylim(ax2,[0 1]); xlim(ax2,[0 win(2)*1000]);
    xlabel(ax2,'Time from CS (ms)'); ylabel(ax2,'Cumulative fraction');
    title(ax2, sprintf('%s: recruitment curve', ratLabels{rr}));
end

figure
rr=6;
lat  = latPerRat{rr};
band = bandPerRat{rr};

% ---- seperate graph of all ----------------------------------------
ax1 = subplot(1,2,1); hold(ax1,'on')
csLat = lat(band==1);
trLat = lat(band==2);

histogram(ax1, [csLat, trLat], 30, 'BinLimits',[0 750], ...
          'FaceColor',[0.3 0.5 0.9], 'EdgeColor','none', 'Normalization','probability');
%histogram(ax1, trLat, 30, 'BinLimits',[0 750], ...
%          'FaceColor',[0.9 0.6 0.2], 'EdgeColor','none','FaceAlpha',0.7);
xline(ax1,250,'--k'); xline(ax1,750,'--k');
axis([0 2000, 0 1])
xlabel(ax1,'Latency (ms)'); ylabel(ax1,'Fraction of Neurons');
title(ax1, sprintf('%s: CS & Trace latencies', ratLabels{rr}));

% ---- Cumulative fraction panel (right) -----------------------------
ax2 = subplot(1,2,2); hold(ax2,'on')
edgesCDF = 0:binW:win(2);
cdfVals  = arrayfun(@(x) mean(lat<=x*1000), edgesCDF);
plot(ax2, edgesCDF*1000, cdfVals,'k-','LineWidth',1.8);
xline(ax2,250,'--'); xline(ax2,750,'--');
ylim(ax2,[0 1]); xlim(ax2,[0 win(2)*1000]);
xlabel(ax2,'Time from CS (ms)'); ylabel(ax2,'Cumulative fraction');
title(ax2, sprintf('%s: recruitment curve', ratLabels{rr}));



sgtitle(sprintf('Sliding-window latency (bin = %.0f ms, α = %.3f)', ...
                 binW*1e3, alpha));

end
