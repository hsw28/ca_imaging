function plotCellComposite(animal, neuronIDs, dateStr, varargin)
% plotCellComposite
% Rows per cell:
%   1) CS-aligned spike raster (dots) with epoch bars
%   2) Spatial rate maps computed via *legacy logic*:
%        - Non-task: CA_normalizePosData(highspeedspikes, ...) then imgaussfilt(.75)
%        - Task    : CA_normalizePosData(CSUSspikes,      ...) then imgaussfilt(.75)
%   3) Tone responses: PSTH and mean±SEM calcium
%
% plotCellComposite(rat0816, [2], '2022_11_04', 'CLimNon', [0.25 0.55], 'CLimTask', [0.0 0.4]);
% Columns: 2 per cell → [Non-task (high-speed) | Task (CSUS)]

% ---------- options ----------
p = inputParser;
p.addParameter('Win',        [-1 2]);     % legacy: applies to both if WinPlot/WinTask not set
p.addParameter('WinPlot',    [-1 2]);    % for raster/PSTH/calcium
p.addParameter('WinTask',    [0 2]);     % defines CSUS windows (task)
p.addParameter('VelThresh',   4);        % LEGACY high-speed velocity threshold (cm/s)
p.addParameter('MapBin',      2.5);      % LEGACY 'dim' in CA_normalizePosData
p.addParameter('SmoothSigma', 0.75);     % imgaussfilt sigma for maps
p.addParameter('BinSize',     0.05);     % PSTH bin (s)
p.addParameter('HeatmapCLim', []);       % [low high] for both maps
p.addParameter('CLimNon',     []);       % [low high] just for Non-task
p.addParameter('CLimTask',    []);       % [low high] just for Task
p.addParameter('TaskFixedMax',0.4);      % default Task map upper bound (legacy)
p.addParameter('RasterStyle','line');   % 'line' | 'rect' | 'dot'
p.addParameter('TickHeight', 0.8);      % trial-height of a tick/rect (in trials)
p.addParameter('RectWidth',  0.03);     % seconds
p.parse(varargin{:});

WinLegacy    = p.Results.Win;
WinPlot      = p.Results.WinPlot;
WinTask      = p.Results.WinTask;

if ~ismember('Win', p.UsingDefaults)
    if ismember('WinPlot', p.UsingDefaults), WinPlot = WinLegacy; end
    if ismember('WinTask', p.UsingDefaults), WinTask = WinLegacy; end
end


VelThresh    = p.Results.VelThresh;
MapBin       = p.Results.MapBin;
SmoothSigma  = p.Results.SmoothSigma;
BinSize      = p.Results.BinSize;

HeatmapCLim  = p.Results.HeatmapCLim;
CLimNon      = p.Results.CLimNon;
CLimTask     = p.Results.CLimTask;
TaskFixedMax = p.Results.TaskFixedMax;

RasterStyle = lower(p.Results.RasterStyle);
TickH       = p.Results.TickHeight;
RectW       = p.Results.RectWidth;

% If user passed Win but didn’t pass the new ones, use Win for both
if ismember('WinPlot', p.UsingDefaults), WinPlot = WinLegacy; end
if ismember('WinTask', p.UsingDefaults), WinTask = WinLegacy; end

% ---------- unpack ----------
pos      = animal.pos.(['pos_' dateStr]);           % [t, x, y]
tPos     = pos(:,1);
CS_times = animal.CS_times.(['CS_'  dateStr])(:);   % [nTrials x 1]
CA_peaks = animal.Ca_peaks.(['CA_peaks_' dateStr]); % [nNeurons x times-as-row]
ca_ts    = animal.Ca_ts.(['CA_time_' dateStr]);

have_traces = isfield(animal,'Ca_traces') && isfield(animal.Ca_traces, ['CA_traces_' dateStr]);
if have_traces
    CA_traces = animal.Ca_traces.(['CA_traces_' dateStr]); % [nNeurons x time]
else
    CA_traces = [];
end

% CA time in seconds, frame-sampled
if size(ca_ts,2) > 1
    ca_ts = ca_ts(:,2) ./ 1000;
    ca_ts = ca_ts(2:2:end);
end
Fs = 1 / max(eps, median(diff(ca_ts)));

% Epochs + colors (match Fig 1b/1c)
epochNames = {'CS','Trace','US','Post'};
epochWin   = {[0 0.25], [0.25 0.75], [0.75 0.85], [0.85 min(2.0, WinTask(2))]};
epochCols  = {[0 0 1], [0 0.7 0], [1 0 0], [0 0 0]};  % blue, green, red, black

% ---------- figure layout ----------
nCells = numel(neuronIDs);
nRows  = 3 * nCells;
nCols  = 2 * nCells;

figure('Color','w','Position',[80 80 1400 900]);

for c = 1:nCells
    nid = neuronIDs(c);
    spk_t = CA_peaks(nid,:); spk_t = spk_t(~isnan(spk_t));

    % ---------- LEGACY spike sets for maps (high-speed vs CSUS) ----------
    [CSUSspikes, highspeedspikes] = buildLegacySpikeSets(spk_t, pos, CS_times, WinTask, VelThresh);

    % ---------- Row 1: CS-aligned raster (WinPlot) ----------
    rowRaster = (c-1)*3 + 1;
    colL = (c-1)*2 + 1;  colR = (c-1)*2 + 2;
    axL = subplot(nRows,nCols,sub2ind([nCols nRows], colL, rowRaster));
    axR = subplot(nRows,nCols,sub2ind([nCols nRows], colR, rowRaster));
    posL = get(axL,'Position'); posR = get(axR,'Position');
    delete(axL); delete(axR);
    axRaster = axes('Position',[posL(1) posL(2) (posR(1)+posR(3)-posL(1)) posL(4)]); hold(axRaster,'on');

    % align spikes to CS using WinPlot  ✅ (build BEFORE raster switch)
    nTrials = numel(CS_times);
    alignedSpikes = cell(nTrials,1);
    for tr = 1:nTrials
        cs = CS_times(tr);
        rel = spk_t(spk_t >= cs + WinPlot(1) & spk_t <= cs + WinPlot(2)) - cs;
        alignedSpikes{tr} = rel(:);
    end

    switch RasterStyle
        case 'line'   % vertical tick per spike (classic raster)
            for tr = 1:numel(CS_times)
                r = alignedSpikes{tr};
                if isempty(r), continue; end
                for k = 1:numel(r)
                    plot(axRaster, [r(k) r(k)], [tr - TickH/2, tr + TickH/2], ...
                        'k-', 'LineWidth', 1.5);
                end
            end

        case 'rect'   % filled block per spike
            for tr = 1:numel(CS_times)
                r = alignedSpikes{tr};
                if isempty(r), continue; end
                yB = tr - TickH/2;
                for k = 1:numel(r)
                    rectangle(axRaster, 'Position', [r(k)-RectW/2, yB, RectW, TickH], ...
                        'FaceColor', 'k', 'EdgeColor', 'none');
                end
            end

        otherwise     % 'dot' (fallback)
            for tr = 1:numel(CS_times)
                r = alignedSpikes{tr};
                if isempty(r), continue; end
                plot(axRaster, r, tr*ones(size(r)), 'k.', 'MarkerSize', 8);
            end
    end

    yTop = numel(CS_times) + 1;
    for e = 1:numel(epochWin)
        xr = epochWin{e};
        patch(axRaster,[xr(1) xr(2) xr(2) xr(1)], [yTop-0.5 yTop-0.5 yTop+0.5 yTop+0.5], ...
              epochCols{e}, 'EdgeColor','none', 'FaceAlpha',0.25);
        text(axRaster, mean(xr), yTop+0.8, epochNames{e}, 'Color', epochCols{e}, ...
             'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
    end

    for tr = 1:numel(CS_times)
        r = alignedSpikes{tr};
        if isempty(r), continue; end
        plot(axRaster, r, tr*ones(size(r)), 'k.', 'MarkerSize',8);
    end
    xline(axRaster, 0, '--', 'Color',[1 0 0], 'LineWidth',1);
    xlim(axRaster, WinPlot);
    ylim(axRaster, [0.5 numel(CS_times)+1.5]);
    ylabel(axRaster, 'Trial');
    title(axRaster, sprintf('Neuron %d — CS-aligned raster', nid));
    set(axRaster,'Box','off');

    % ---------- Row 2: Rate maps via LEGACY CA_normalizePosData ----------
    rowMaps = (c-1)*3 + 2;

    rate_non = computeRateMapLegacy(highspeedspikes, pos, MapBin, SmoothSigma);
    rate_tsk = computeRateMapLegacy(CSUSspikes,     pos, MapBin, SmoothSigma);

    % CLim selection (legacy defaults + overrides)
    [loNon, hiNon] = legacyCLim(rate_non, 'non',  [], [], [], []);
    [loTsk, hiTsk] = legacyCLim(rate_tsk, 'task', [], [], [], TaskFixedMax);

    if ~isempty(HeatmapCLim)
        loNon = HeatmapCLim(1); hiNon = HeatmapCLim(2);
        loTsk = loNon;          hiTsk = hiNon;
    end
    if ~isempty(CLimNon),  loNon = CLimNon(1);  hiNon = CLimNon(2); end
    if ~isempty(CLimTask), loTsk = CLimTask(1); hiTsk = CLimTask(2); end

    ax1 = subplot(nRows,nCols,sub2ind([nCols nRows], colL, rowMaps));
    imagesc(rate_non, [loNon hiNon]); axis off; colormap(parula);
    title(sprintf('Neuron %d — Non-task (high-speed)', nid)); colorbar;

    ax2 = subplot(nRows,nCols,sub2ind([nCols nRows], colR, rowMaps));
    imagesc(rate_tsk, [loTsk hiTsk]); axis off; colormap(parula);
    title('Task (CSUS)'); colorbar;

    % ---------- Row 3: Tone responses (PSTH | mean±SEM calcium) ----------
    rowTone = (c-1)*3 + 3;

    % PSTH
    axP = subplot(nRows,nCols,sub2ind([nCols nRows], colL, rowTone)); hold(axP,'on');
    edges = WinPlot(1):BinSize:WinPlot(2);
    allRel = cell2mat(alignedSpikes);
    counts = histcounts(allRel, edges);
    psth = counts / max(1,(numel(CS_times)*BinSize));
    ctrs = edges(1:end-1) + BinSize/2;

    trW = epochWin{2};  yl = [0 max(eps, max(psth)*1.1)];
    patch(axP,[trW(1) trW(2) trW(2) trW(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          [0 0.7 0], 'FaceAlpha',0.08, 'EdgeColor','none');
    bar(axP, ctrs, psth, 1, 'FaceColor',[0.25 0.25 0.25],'EdgeColor','none');
    xline(axP, 0, '--', 'Color',[1 0 0], 'LineWidth',1);
    xlim(axP, WinPlot); ylim(axP, yl);
    xlabel(axP,'Time from CS (s)'); ylabel(axP,'Spikes/s'); title(axP,'PSTH'); set(axP,'Box','off');

    % Calcium mean±SEM
    axC = subplot(nRows,nCols,sub2ind([nCols nRows], colR, rowTone)); hold(axC,'on');
    if have_traces && ~isempty(CA_traces)
        alignN = max(2, round((WinPlot(2)-WinPlot(1))*Fs));
        caAligned = nan(numel(CS_times), alignN);
        trc = CA_traces(nid,:);  % 1 x time
        for i = 1:numel(CS_times)
            cs = CS_times(i);
            mask = ca_ts >= cs+WinPlot(1) & ca_ts <= cs+WinPlot(2);
            seg  = trc(mask);
            L = min(alignN, numel(seg));
            if L >= 2, caAligned(i,1:L) = seg(1:L); end
        end
        tAligned = linspace(WinPlot(1), WinPlot(2), alignN);
        m = nanmean(caAligned,1);
        s = nanstd( caAligned,[],1) ./ max(1,sum(~isnan(caAligned(:,1))));
        yl2 = [min(m - s), max(m + s)]; if ~all(isfinite(yl2)) || yl2(2)<=yl2(1), yl2 = [0 1]; end

        patch(axC,[trW(1) trW(2) trW(2) trW(1)], [yl2(1) yl2(1) yl2(2) yl2(2)], ...
              [0 0.7 0], 'FaceAlpha',0.08, 'EdgeColor','none');

        fill(axC, [tAligned fliplr(tAligned)], [m+s fliplr(m-s)], [0.6 0.1 0.8], ...
             'FaceAlpha',0.25, 'EdgeColor','none');
        plot(axC, tAligned, m, 'Color',[0.6 0.1 0.8], 'LineWidth',1.5);

        xline(axC, 0, '--', 'Color',[1 0 0], 'LineWidth',1);
        xlim(axC, WinPlot); ylim(axC, yl2);
        xlabel(axC,'Time from CS (s)'); ylabel(axC,'Calcium (a.u.)');
        title(axC,'Mean \pm SEM calcium'); set(axC,'Box','off');
    else
        axis(axC,'off'); text(axC,0.5,0.5,'(No CA\_traces for this session)','HorizontalAlignment','center');
    end
end

%sgtitle(sprintf('%s — %s', getFieldOr(animal,'name','animal'), strrep(dateStr,'_','-')), 'FontWeight','bold');
%set(gcf,'Renderer','painters');
%outfile = sprintf('Fig_1cStyle_composite_%s.svg', dateStr);
%print(gcf, outfile, '-dsvg', '-r600');
%fprintf('Saved %s\n', outfile);
end

% ---------- LEGACY HELPERS ----------
function [CSUSspikes, highspeedspikes] = buildLegacySpikeSets(peaks_time, pos, cs_times, winTask, velthreshold)
% STRICT legacy replication of your plot_place_vs_CSUS logic.

t_pos = pos(:,1);

% ---- CSUS spikes (absolute times), plus CSUS position indices ----
CSUSspikes = [];
csusPosIdx = [];
for i = 1:numel(cs_times)
    cs = cs_times(i);
    relSpikes = peaks_time(peaks_time >= cs + winTask(1) & peaks_time <= cs + winTask(2));
    CSUSspikes = [CSUSspikes, relSpikes]; %#ok<AGROW>

    csusIdx = find(t_pos >= cs + winTask(1) & t_pos <= cs + winTask(2));
    csusPosIdx = [csusPosIdx, csusIdx']; %#ok<AGROW>
end

% ---- velocity threshold (NO index offset; use exactly goodtime = pos(goodvel,1)) ----
vel = ca_velocity(pos);              % assume vel(1,:) is speed
goodvel = find(vel(1,:) >= velthreshold);
goodtime = pos(goodvel, 1);          % legacy line (no +1)
% goodpos  = pos(goodvel,:);         % not needed here

% (legacy code removes csusPosIdx from indices AFTER computing goodtime)
goodvel = setdiff(goodvel, csusPosIdx); %#ok<NASGU> % kept only to mirror your flow

% ---- high-speed spikes: nearest-time match within 1/15 s ----
tol = 1/15;  % ~0.067 s
highspeedspikes = [];
for ii = 1:numel(peaks_time)
    if isnan(peaks_time(ii)), continue; end
    [minValue_vel, closestIndex] = min(abs(peaks_time(ii) - goodtime)); %#ok<ASGLU>
    if ~isempty(minValue_vel) && minValue_vel <= tol
        highspeedspikes(end+1) = peaks_time(ii); %#ok<AGROW>
    end
end

% legacy exclusion of CSUS spikes from high-speed
highspeedspikes = setdiff(highspeedspikes, CSUSspikes);
end

function rate = computeRateMapLegacy(spikeTimes, pos, dim, smoothSigma)
% EXACT legacy mapping:
%   rate = CA_normalizePosData(spikeTimes, pos, dim, 1.000);
%   rate(isnan(rate)) = 0;
%   rate = imgaussfilt(rate, .75);

rate = CA_normalizePosData(spikeTimes, pos, dim, 1.000);
rate(isnan(rate)) = 0;

% force sigma to your legacy default if none provided
if nargin < 4 || isempty(smoothSigma)
    smoothSigma = 0.75;
end
if smoothSigma > 0
    rate = imgaussfilt(rate, smoothSigma);
end
end

function v = getFieldOr(s, f, d)
if isstruct(s) && isfield(s,f), v = s.(f); else, v = d; end
end

function [lo, hi] = legacyCLim(rate, mode, heatmapLOW, heatmapHIGH, useExplicit, taskFixedMax)
vals = rate(~isnan(rate));
if isempty(vals), lo = 0; hi = 1; return, end

sortedDesc = sort(vals,'descend');
kTop = max(1, ceil(numel(sortedDesc)*0.01));
maxratefive = min(sortedDesc(1:kTop));

firstZero = find(sortedDesc==0, 1, 'first');
if isempty(firstZero), nz = sortedDesc; else, nz = sortedDesc(1:firstZero-1); end
if isempty(nz)
    lo_est = 0;
else
    idxLo = max(1, floor(numel(nz)*0.18));
    lo_est = max(nz(idxLo:end));
end

switch lower(mode)
    case 'non'   % high-speed panel
        lo = lo_est; hi = maxratefive;
        if ~isfinite(lo) || ~isfinite(hi) || hi<=lo, lo = 0; hi = max(1,maxratefive); end
    case 'task'  % CSUS panel
        lo = 0.0; if nargin<6 || isempty(taskFixedMax), taskFixedMax = 0.4; end
        hi = taskFixedMax;
    otherwise, lo = 0; hi = 1;
end
end
