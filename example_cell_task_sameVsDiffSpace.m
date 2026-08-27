function OUT = example_cells_task_sameVsDiffSpace(ratVar, dayLabel, varargin)
% EXAMPLE_CELLS_TASK_SAMEVSDIFFSPACE
% Multi-cell illustration (one row per cell) of SAME-space vs DIFF-space,
% time-locked within TraceWin. Uses temporal bin-matched (b↔b) trial correlations.
%
% Usage:
%   OUT = example_cells_task_sameVsDiffSpace('rat0816','2022_11_03');
%   OUT = example_cells_task_sameVsDiffSpace('rat0816','2022_11_03','CellIdx',[12 39 51]);
%
% Key options (name/value):
%   'CellIdx'        []        % vector of cell indices to plot; if empty, auto-pick
%   'TopN'           10        % how many cells to show (if auto-picking)
%   'TraceWin'       [0 2]     % seconds relative to CS
%   'TimeBins'       12        % temporal bins across TraceWin
%   'GridRC'         [2 2]     % spatial grid (rows, cols)
%   'UseSpeedMask'   true
%   'VelThresh'      4         % cm/s
%   'MinFrames'      2         % min frames in intersection (temporal-bin ∩ spatial-bin)
%   'MinBinsCorr'    4         % min overlapping temporal bins for a trial–trial corr
%   'MinTaskSpikes'  20        % MIN spikes for the cell across ALL task frames (filter)
%   'CellNorm'       'none'    % 'none'|'demean'|'zscore' for frame vectors S
%   'Colormap'       'parula'
%   'ShowPSTH'       true      % add same-place/different-place PSTH lines
%   'PSTHWindow'     [0 2]     % seconds relative to CS
%   'PSTHBin'        0.05      % PSTH bin width (s)
%   'Verbose'        true
%
% Auto-pick logic:
%   Rank valid cells by |mean(SAME_r) - mean(DIFF_r)| (smaller is better),
%   so SAME and DIFF are ~equal (per your request), then keep TopN.

% -------------------- options --------------------
p = inputParser;
addParameter(p,'CellIdx',[]);
addParameter(p,'TopN',10);
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'TimeBins',15);
addParameter(p,'GridRC',[2 2]);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'VelThresh',4);
addParameter(p,'MinFrames',0); %dont think this does what i think, cant be >2
addParameter(p,'MinBinsCorr',10);
addParameter(p,'MinTaskSpikes',35);
addParameter(p,'CellNorm','none');
addParameter(p,'Colormap','parula');
addParameter(p,'ShowPSTH',true);
addParameter(p,'PSTHWindow',[0 2]);
addParameter(p,'PSTHBin',0.05);
addParameter(p,'PSTHSmoothBins',1);
addParameter(p,'Verbose',true);
parse(p,varargin{:});
o = p.Results; B = o.TimeBins;

assert(evalin('base',sprintf('exist(''%s'',''var'')',ratVar))==1, 'Rat variable %s not found.', ratVar);
RAT   = evalin('base', ratVar);
spk   = RAT.Ca_peaks.(sprintf('CA_peaks_%s', dayLabel));
posd  = RAT.pos.(sprintf('pos_%s', dayLabel));
csVec = RAT.CS_times.(sprintf('CS_%s', dayLabel));

% ---- position & frame timebase ----
[t, x, y] = coerce_pos(posd);
if max(t,[],'omitnan')>1e4, t = t/1000; end
if numel(t)<5, error('Too few position frames.'); end

% ---- spikes → frame histogram per cell ----
cellSpk = to_cell_spikes(spk);
Nc      = numel(cellSpk);
S       = spikes_to_matrix(cellSpk, t);
S       = normalize_cells(S, o.CellNorm);

% ---- optional speed mask ----
if o.UseSpeedMask
    v = speed_cm_per_s(t,x,y);
    speedOK = (v >= o.VelThresh);
else
    speedOK = true(size(t));
end

% ---- grid & spatial bin per frame ----
GridR = o.GridRC(1); GridC = o.GridRC(2);
xedges = linspace(min(x,[],'omitnan'), max(x,[],'omitnan'), GridC+1);
yedges = linspace(min(y,[],'omitnan'), max(y,[],'omitnan'), GridR+1);
bx = discretize(x, xedges);
by = discretize(y, yedges);
K  = GridR*GridC;
bSpace = nan(size(t));
m = isfinite(bx) & isfinite(by);
bSpace(m) = sub2ind([GridR,GridC], by(m), bx(m));

% ---- temporal-bin frame indices for each trial ----
cs = double(csVec(:)); if isempty(cs), error('No CS times found.'); end
tw0 = o.TraceWin(1); tw1 = o.TraceWin(2);
nT  = numel(cs);

framesBy = repmat(struct('idx',[]), nT, B);
taskFramesMask = false(size(t));  % for MinTaskSpikes

for tr = 1:nT
    idx_all = find(t >= cs(tr)+tw0 & t < cs(tr)+tw1);
    if isempty(idx_all), continue; end
    taskFramesMask(idx_all) = true;
    e = round(linspace(1, numel(idx_all)+1, B+1));
    e(1)=1; e(end)=numel(idx_all)+1;
    for b = 1:B
        if e(b) >= e(b+1), continue; end
        seg = idx_all(e(b):e(b+1)-1);
        seg = seg(speedOK(seg));
        if numel(seg) >= o.MinFrames
            framesBy(tr,b).idx = seg(:);
        end
    end
end

% =================== compute SAME/DIFF per cell ======================
cell_same_mu = nan(Nc,1);
cell_diff_mu = nan(Nc,1);
cell_ok      = false(Nc,1);
cell_spikes_task = sum(S,2);  % total spikes per cell across ALL frames
cell_spikes_task = cell_spikes_task(:)'; % row

% limit to task frames only for MinTaskSpikes
task_spikes = sum(S(:, taskFramesMask), 2);
task_spikes = task_spikes(:);

for c = 1:Nc
    if task_spikes(c) < o.MinTaskSpikes
        continue;  % fails spike-count gate during task period
    end
    [same_r, diff_r] = corr_same_vs_diff_for_cell(c, S, bSpace, framesBy, B, K, o.MinBinsCorr);
    if ~isempty(same_r) || ~isempty(diff_r)
        cell_same_mu(c) = mean(same_r,'omitnan');
        cell_diff_mu(c) = mean(diff_r,'omitnan');
        cell_ok(c) = true;
    end
end

% ===== choose cells to plot =====
if isempty(o.CellIdx)
    idxValid = find(cell_ok);
    idxValid1 = find(cell_same_mu>0);
    idxValid2 = find(cell_diff_mu>0);
    idxValid3 = intersect(idxValid1,idxValid2);
    idxValid = intersect(idxValid, idxValid3);
    if isempty(idxValid), error('No cells passed MinTaskSpikes and correlation gates.'); end
    % rank by closeness of SAME vs DIFF (small |Δr| => more equal)
    delta = abs(cell_same_mu(idxValid) + cell_diff_mu(idxValid));
    %[~, ord] = sort(delta, 'ascend');
    [~, ord] = sort(delta, 'descend');
    cells = idxValid(ord);
    cells = cells(1:min(o.TopN, numel(cells)));
else
    cells = o.CellIdx(:)';
    cells = cells(cells>=1 & cells<=Nc);
    if isempty(cells), error('CellIdx outside valid range 1..%d.', Nc); end
    % still enforce MinTaskSpikes gate on requested cells
    keep = task_spikes(cells) >= o.MinTaskSpikes;
    if ~any(keep), error('Requested CellIdx all failed MinTaskSpikes gate.'); end
    cells = cells(keep);
end

% ===== Verbose printing of which cells we’ll plot (row order) =====
if o.Verbose
    fprintf('[%s | %s] plotting %d cells (row order): %s\n', ...
        ratVar, dayLabel, numel(cells), mat2str(cells));
    for ii = 1:numel(cells)
        fprintf('  Row %2d -> cell %d (task spikes = %d)\n', ii, cells(ii), task_spikes(cells(ii)));
    end
end

% ============================ plot =============================
R = numel(cells); C = 3 + (o.ShowPSTH ~= 0); % SAME heat | DIFF heat | bars | PSTH
figure('Color','w','Position',[100 80 1300 max(300, 140*R)]);
tiledlayout(R, C, 'TileSpacing','compact','Padding','compact');
cm = o.Colormap; if ischar(cm), cm = feval(cm, 256); end
frameDt = median(diff(t),'omitnan');
if ~isfinite(frameDt) || frameDt<=0, frameDt = max(eps, mean(diff(t),'omitnan')); end

% compute a per-figure color range so rows are comparable
% (based on all values used in SAME/DIFF averages below)
allVals = [];
for iRow = 1:R
    cIdx = cells(iRow);
    [meanSame, meanDiff] = averaged_same_diff_profiles(cIdx, S, bSpace, framesBy, B, K);
    allVals = [allVals; meanSame(:)./frameDt; meanDiff(:)./frameDt]; %#ok<AGROW>
end
clim = [min(allVals(isfinite(allVals)),[],'omitnan'), max(allVals(isfinite(allVals)),[],'omitnan')];
if ~all(isfinite(clim)) || clim(2)<=clim(1), clim = [0 1]; end

for iRow = 1:R
    cIdx = cells(iRow);

    % (1) SAME trial-averaged temporal profile heat (1×B)
    [meanSame, meanDiff] = averaged_same_diff_profiles(cIdx, S, bSpace, framesBy, B, K);
    meanSameRate = meanSame ./ frameDt;
    meanDiffRate = meanDiff ./ frameDt;
    nexttile; imagesc(1:B, 1, meanSameRate); set(gca,'YTick',[]); axis tight
    colormap(cm); cb = colorbar; ylabel(cb,'Events/s'); caxis(clim);
    title(sprintf('Cell %d  SAME space (avg trials)', cIdx));
    xlabel('Temporal bin');

    % (2) DIFF trial-averaged temporal profile heat (1×B)
    nexttile; imagesc(1:B, 1, meanDiffRate); set(gca,'YTick',[]); axis tight
    colormap(cm); cb = colorbar; ylabel(cb,'Events/s'); caxis(clim);
    title('DIFF space (avg other bins)');
    xlabel('Temporal bin');

    % (3) bar + jitter: SAME vs DIFF mean r for this cell
    [same_r, diff_r] = corr_same_vs_diff_for_cell(cIdx, S, bSpace, framesBy, B, K, o.MinBinsCorr);
    same_mu = mean(same_r,'omitnan'); diff_mu = mean(diff_r,'omitnan');

    nexttile; hold on; box on; yline(0,'k:');
    bar(1, same_mu, 0.65, 'FaceColor',[0.30 0.60 1.00], 'EdgeColor','k');
    bar(2, diff_mu, 0.65, 'FaceColor',[0.85 0.40 0.20], 'EdgeColor','k');
    % jittered points
    if ~isempty(same_r)
        scatter(1 + 0.10*(rand(size(same_r))-0.5), same_r, 16, 'k', 'filled', 'MarkerFaceAlpha',0.35);
    end
    if ~isempty(diff_r)
        scatter(2 + 0.10*(rand(size(diff_r))-0.5), diff_r, 16, 'k', 'filled', 'MarkerFaceAlpha',0.35);
    end
    xlim([0.5 2.5]); xticks([1 2]); xticklabels({'SAME','DIFF'});
    ylabel('Trial–trial temporal corr (r)');
    title(sprintf('Cell %d  mean r: %.2f vs %.2f  |  |Δ|=%.2f', cIdx, same_mu, diff_mu, abs(same_mu-diff_mu)));

    if o.ShowPSTH
        kRef = reference_space_bin(bSpace, framesBy, B, K);
        [tp, samePSTH, diffPSTH, nSameFrames, nDiffFrames] = same_diff_place_psth( ...
            cellSpk{cIdx}, cs, t, bSpace, kRef, o.PSTHWindow, o.PSTHBin, o.PSTHSmoothBins);
        nexttile; hold on
        plot(tp, samePSTH, 'Color',[0.30 0.60 1.00], 'LineWidth',1.8);
        plot(tp, diffPSTH, 'Color',[0.85 0.40 0.20], 'LineWidth',1.8);
        xline(0,'k--','LineWidth',1.0);
        xline(o.TraceWin(2),'--','Color',[0.65 0.15 0.10],'LineWidth',1.0);
        yline(0,':','Color',[0.65 0.65 0.65]);
        yl = [0 max([samePSTH(:); diffPSTH(:)],[],'omitnan')*1.15];
        if ~all(isfinite(yl)) || yl(2)<=yl(1), yl = [0 1]; end
        ylim(yl); xlim(o.PSTHWindow);
        xlabel('Time from CS (s)'); ylabel('Events/s');
        title(sprintf('PSTH same vs diff place (bin %d)', kRef));
        legend({sprintf('same place (n=%d frames)',nSameFrames), ...
                sprintf('different place (n=%d frames)',nDiffFrames)}, ...
               'Location','northoutside','Orientation','horizontal');
        set(gca,'Box','off');
    end
end

sgtitle(sprintf('%s | %s | Trace %.1f–%.1fs | Grid %dx%d | speed\\geq%.1f=%d | MinTaskSpikes=%d', ...
    ratVar, strrep(dayLabel,'_','\_'), o.TraceWin(1), o.TraceWin(2), ...
    GridR, GridC, o.VelThresh, o.UseSpeedMask, o.MinTaskSpikes));

% --------------- pack outputs ---------------
OUT = struct();
OUT.cells = cells(:);
OUT.meta  = struct('rat',ratVar,'day',dayLabel,'TraceWin',o.TraceWin,'TimeBins',B, ...
                   'GridRC',o.GridRC,'UseSpeedMask',o.UseSpeedMask,'VelThresh',o.VelThresh, ...
                   'MinFrames',o.MinFrames,'MinBinsCorr',o.MinBinsCorr,'MinTaskSpikes',o.MinTaskSpikes);

end

% ================= helper functions =================

function [same_r, diff_r] = corr_same_vs_diff_for_cell(c, S, bSpace, framesBy, B, K, minBins)
% Compute SAME and DIFF trial–trial correlations for cell c (temporal b matched).
same_r = []; diff_r = [];
nT = size(framesBy,1);
V  = nan(nT, K, B);
for tr = 1:nT
    for b = 1:B
        idx = framesBy(tr,b).idx;
        if isempty(idx), continue; end
        bs  = bSpace(idx);
        for k = 1:K
            hit = (bs == k);
            if any(hit)
                V(tr,k,b) = mean(S(c, idx(hit)), 2, 'omitnan');
            end
        end
    end
end
% SAME: (t1,b,k) vs (t2,b,k)
for k = 1:K
    trs = find(any(isfinite(squeeze(V(:,k,:))), 2));
    if numel(trs) < 2, continue; end
    for i1 = 1:numel(trs)-1
        v1 = squeeze(V(trs(i1),k,:));
        for i2 = i1+1:numel(trs)
            v2 = squeeze(V(trs(i2),k,:));
            m  = isfinite(v1) & isfinite(v2);
            if nnz(m) >= minBins
                r = corr(v1(m), v2(m), 'type','Pearson');
                if isfinite(r), same_r(end+1,1) = r; end %#ok<AGROW>
            end
        end
    end
end
% DIFF: (t1,b,k1) vs (t2,b,k2), k1≠k2
for k1 = 1:K
    tr1 = find(any(isfinite(squeeze(V(:,k1,:))), 2));
    if isempty(tr1), continue; end
    for k2 = 1:K
        if k2==k1, continue; end
        tr2 = find(any(isfinite(squeeze(V(:,k2,:))), 2));
        if isempty(tr2), continue; end
        for i1 = 1:numel(tr1)
            v1 = squeeze(V(tr1(i1),k1,:));
            for i2 = 1:numel(tr2)
                if tr2(i2)==tr1(i1), continue; end
                v2 = squeeze(V(tr2(i2),k2,:));
                m  = isfinite(v1) & isfinite(v2);
                if nnz(m) >= minBins
                    r = corr(v1(m), v2(m), 'type','Pearson');
                    if isfinite(r), diff_r(end+1,1) = r; end %#ok<AGROW>
                end
            end
        end
    end
end
end

function [meanSame, meanDiff] = averaged_same_diff_profiles(cIdx, S, bSpace, framesBy, B, K)
% Build trial-averaged temporal profiles (length B) for SAME and DIFF for this cell.
% We pick kRef as the spatial bin with the most task frames.
kRef = reference_space_bin(bSpace, framesBy, B, K);
nT = size(framesBy,1);

meanSame = nan(1,B);
meanDiff = nan(1,B);
for b = 1:B
    fSame = [];
    fDiff = [];
    for tr = 1:nT
        idx = framesBy(tr,b).idx;
        if isempty(idx), continue; end
        bs = bSpace(idx);
        fSame = [fSame; idx(bs==kRef)]; %#ok<AGROW>
        fDiff = [fDiff; idx(bs~=kRef & isfinite(bs))]; %#ok<AGROW>
    end
    if ~isempty(fSame), meanSame(b) = mean(S(cIdx, fSame), 2, 'omitnan'); end
    if ~isempty(fDiff), meanDiff(b) = mean(S(cIdx, fDiff), 2, 'omitnan'); end
end
end

function kRef = reference_space_bin(bSpace, framesBy, B, K)
% Pick the spatial bin with the most task frames; this is the SAME place.
nT = size(framesBy,1);
countK = zeros(K,1);
for k = 1:K
    for b = 1:B
        for tr = 1:nT
            idx = framesBy(tr,b).idx;
            if isempty(idx), continue; end
            countK(k) = countK(k) + nnz(bSpace(idx) == k);
        end
    end
end
[~, kRef] = max(countK);
if ~isfinite(kRef) || countK(kRef)==0, kRef = 1; end
end

function [tp, sameRate, diffRate, nSameFrames, nDiffFrames] = same_diff_place_psth( ...
    events, cs, t, bSpace, kRef, psthWin, bin, smoothBins)
% Occupancy-normalized event-rate PSTHs split by SAME place versus DIFF place.
edges = psthWin(1):bin:psthWin(2);
tp = edges(1:end-1) + bin/2;
dt = median(diff(t),'omitnan');
if ~isfinite(dt) || dt<=0, dt = max(eps, mean(diff(t),'omitnan')); end

sameEvt = zeros(1,numel(tp));
diffEvt = zeros(1,numel(tp));
sameOcc = zeros(1,numel(tp));
diffOcc = zeros(1,numel(tp));
nSameFrames = 0;
nDiffFrames = 0;

events = double(events(:));
events = events(isfinite(events) & events>0);
for tr = 1:numel(cs)
    frameIdx = find(t >= cs(tr)+edges(1) & t < cs(tr)+edges(end));
    if ~isempty(frameIdx)
        relFrame = t(frameIdx) - cs(tr);
        [~,~,frameBin] = histcounts(relFrame, edges);
        sameFrame = bSpace(frameIdx) == kRef;
        diffFrame = isfinite(bSpace(frameIdx)) & bSpace(frameIdx) ~= kRef;
        for b = 1:numel(tp)
            sameOcc(b) = sameOcc(b) + nnz(frameBin == b & sameFrame) * dt;
            diffOcc(b) = diffOcc(b) + nnz(frameBin == b & diffFrame) * dt;
        end
        nSameFrames = nSameFrames + nnz(sameFrame);
        nDiffFrames = nDiffFrames + nnz(diffFrame);
    end

    relEvt = events(events >= cs(tr)+edges(1) & events < cs(tr)+edges(end)) - cs(tr);
    if isempty(relEvt), continue; end
    evtFrame = nearest_frame_index(t, relEvt + cs(tr));
    [~,~,evtBin] = histcounts(relEvt, edges);
    sameEvent = bSpace(evtFrame) == kRef;
    diffEvent = isfinite(bSpace(evtFrame)) & bSpace(evtFrame) ~= kRef;
    for b = 1:numel(tp)
        sameEvt(b) = sameEvt(b) + nnz(evtBin == b & sameEvent);
        diffEvt(b) = diffEvt(b) + nnz(evtBin == b & diffEvent);
    end
end

sameRate = sameEvt ./ max(sameOcc, eps);
diffRate = diffEvt ./ max(diffOcc, eps);
sameRate(sameOcc <= 0) = NaN;
diffRate(diffOcc <= 0) = NaN;

if smoothBins > 0
    w = ones(1, 2*round(smoothBins)+1);
    w = w ./ sum(w);
    sameRate = smooth_with_nan(sameRate, w);
    diffRate = smooth_with_nan(diffRate, w);
end
end

function y = smooth_with_nan(x, w)
good = isfinite(x);
x0 = x;
x0(~good) = 0;
den = conv(double(good), w, 'same');
y = conv(x0, w, 'same') ./ max(den, eps);
y(den <= 0) = NaN;
end

function idx = nearest_frame_index(tFrames, tEvents)
idx = round(interp1(tFrames, 1:numel(tFrames), tEvents, 'linear', 'extrap'));
idx = max(1, min(numel(tFrames), idx));
end

function [t,x,y] = coerce_pos(pos)
if istable(pos)
    vn = lower(string(pos.Properties.VariableNames));
    t = pos{:, find(ismember(vn, {'t','time','ts'}),1)};
    x = pos{:, find(ismember(vn, {'x','xpos','x_cm','xcm','xsmooth','posx','x_smooth'}),1)};
    y = pos{:, find(ismember(vn, {'y','ypos','y_cm','ycm','ysmooth','posy','y_smooth'}),1)};
elseif isnumeric(pos) && size(pos,2)>=3
    t = pos(:,1); x = pos(:,2); y = pos(:,3);
else
    f = lower(string(fieldnames(pos)));
    t = pos.(f{find(ismember(f,{'t','time','ts'}),1)});
    x = pos.(f{find(ismember(f,{'x','xpos','x_cm','xcm','xsmooth','posx','x_smooth'}),1)});
    y = pos.(f{find(ismember(f,{'y','ypos','y_cm','ycm','ysmooth','posy','y_smooth'}),1)});
end
t = double(t(:)); x = double(x(:)); y = double(y(:));
end

function v = speed_cm_per_s(t, x, y)
dt = diff(t); dt(end+1,1) = median(dt(dt>0),'omitnan');
dx = [diff(x); 0]; dy = [diff(y); 0];
v = hypot(dx,dy) ./ max(dt, eps);
v(~isfinite(v)) = 0;
end

function cellSpk = to_cell_spikes(spk)
if iscell(spk), cellSpk = spk; return; end
[nc, ~] = size(spk);
cellSpk = cell(nc,1);
for c = 1:nc
    st = double(spk(c,:));
    st = st(isfinite(st) & st > 0);
    cellSpk{c} = st(:);
end
end

function S = spikes_to_matrix(spikeCell, t)
t = double(t(:));
if numel(t)<2, S=zeros(numel(spikeCell), numel(t), 'single'); return; end
dt = median(diff(t),'omitnan'); if ~isfinite(dt)||dt<=0, dt = max(eps, mean(diff(t),'omitnan')); end
edges = [t - dt/2; t(end) + dt/2];
nc = numel(spikeCell); S = zeros(nc, numel(t), 'single');
for c=1:nc
    st = spikeCell{c}; if isempty(st), continue; end
    S(c,:) = histcounts(st, edges);
end
end

function S = normalize_cells(S, mode)
switch lower(mode)
    case 'none'
    case 'demean'
        mu = mean(S,2,'omitnan'); S = S - mu;
    case 'zscore'
        mu = mean(S,2,'omitnan'); sd = std(S,0,2,'omitnan'); sd(sd==0|~isfinite(sd)) = 1; S = (S - mu) ./ sd;
    otherwise
        error('CellNorm must be ''none'',''demean'',''zscore''.');
end
end
