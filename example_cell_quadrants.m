function OUT = example_cell_quadrants(ratVar, dayLabel, varargin)
% EXAMPLE_CELL_QUADRANTS
% Pick Top N cells that look "place-y" WITHOUT task (one quadrant dominates)
% but become spread out WITH task (uniform), and plot them.
%
% Also adds two panels (right side):
%  - Per-quadrant CTRL vs WITH-task PV correlations (population metric)
%  - Overall CTRL vs WITH-task correlation
% If a single CellIdx is provided, the final panel shows that cell's
% WITH-vs-WITHOUT rate-map correlation (so it reflects the example cell).
%
% Example:
%   OUT = example_cell_quadrants('rat0816','2022_11_04', ...
%       'TopN',6,'ScoreMode','peakspread','MinOccSec',0,'MinEvents',0);

% ---------------- options ----------------
p = inputParser;
addParameter(p,'TopN',6);
addParameter(p,'CellIdx',[]);                % [] => auto-pick TopN; scalar => show that cell only
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'BufferPost',0);
addParameter(p,'VelThresh',4);
addParameter(p,'UseSpeedMask',true);
addParameter(p,'GridRC',[2 2]);              % quadrants by default
addParameter(p,'MinOccSec',0.10);
addParameter(p,'MinEvents',5);
addParameter(p,'Colormap','parula');
addParameter(p,'ScoreMode','peakspread');    % 'peakspread' | 'cosine' | 'corr' | 'l2'
addParameter(p,'SavePNG','');                % '' => no save
% peakspread knobs
addParameter(p,'PeakSpreadMinOccFrac',0.3);
addParameter(p,'PeakSpreadPseudo',1e-8);
addParameter(p,'PeakSpreadDropMode','diff'); % 'ratio' | 'diff'
addParameter(p,'PeakSpreadAlpha',[1 1 1]);   % [conc, even, drop] weights
% correlation panels
addParameter(p,'ShowCorrPanels',true);
addParameter(p,'CorrTimeBins',15);           % for PV-with corr (population)
parse(p,varargin{:});
o = p.Results; TopN = o.TopN;

assert(evalin('base',sprintf('exist(''%s'',''var'')',ratVar))==1, 'Rat variable %s not found.', ratVar);
RAT   = evalin('base', ratVar);
spk   = RAT.Ca_peaks.(sprintf('CA_peaks_%s', dayLabel));
posd  = RAT.pos.(sprintf('pos_%s', dayLabel));
csVec = RAT.CS_times.(sprintf('CS_%s', dayLabel));

% ---- position & condition masks ----
[tPos, xPos, yPos] = coerce_pos(posd);
dt = median(diff(tPos),'omitnan'); if ~isfinite(dt)||dt<=0, dt = max(eps, mean(diff(tPos),'omitnan')); end
isTaskCore  = false(size(tPos));
isTaskExcl  = false(size(tPos));
for j = 1:numel(csVec)
    w0 = csVec(j) + o.TraceWin(1);
    w1 = csVec(j) + o.TraceWin(2);
    isTaskCore = isTaskCore | (tPos >= w0 & tPos < w1);
    isTaskExcl = isTaskExcl | (tPos >= w0 & tPos < (w1 + max(0,o.BufferPost)));
end
isWithout = ~isTaskExcl;

if o.UseSpeedMask
    v = speed_cm_per_s(tPos, xPos, yPos);
    speedOK = (v >= o.VelThresh);
else
    speedOK = true(size(tPos));
end
maskWith    = isTaskCore & speedOK;
maskWithout = isWithout  & speedOK;

% ---- grid & occupancy (sec) ----
[edges, K] = make_grid_edges(xPos, yPos, o.GridRC);
occ_with    = occupancy2D(tPos, xPos, yPos, maskWith,    edges, dt);
occ_without = occupancy2D(tPos, xPos, yPos, maskWithout, edges, dt);

% ---- spikes → per-cell event times ----
cellSpk = to_cell_spikes(spk);
nC      = numel(cellSpk);
xAt = @(tt) interp1_nonan(tPos, xPos, tt);
yAt = @(tt) interp1_nonan(tPos, yPos, tt);

evt_with    = zeros(K, nC);
evt_without = zeros(K, nC);

% (only used if we later want per-quadrant time-resolved stuff)
evt_with_times    = cell(K, nC);
evt_without_times = cell(K, nC);

for c = 1:nC
    ts = double(cellSpk{c});
    if isempty(ts), continue; end
    ts = ts(ts>=tPos(1) & ts<=tPos(end));
    if isempty(ts), continue; end

    % Build trial masks on spike times
    inCore = false(size(ts));
    inWout = true(size(ts));
    for j = 1:numel(csVec)
        w0 = csVec(j) + o.TraceWin(1);
        w1 = csVec(j) + o.TraceWin(2);
        inCore = inCore | (ts >= w0 & ts < w1);
        inWout = inWout & ~(ts >= w0 & ts < (w1 + max(0,o.BufferPost)));
    end
    if o.UseSpeedMask
        vEvt  = speed_cm_per_s(tPos, xPos, yPos, ts);
        inCore = inCore & (vEvt >= o.VelThresh);
        inWout = inWout & (vEvt >= o.VelThresh);
    end

    % --- WITH (correct mapping of submask back to ts) ---
    if any(inCore)
        idxCore = find(inCore);
        bx = discretize(xAt(ts(inCore)), edges.x);
        by = discretize(yAt(ts(inCore)), edges.y);
        k  = sub2ind([numel(edges.y)-1, numel(edges.x)-1], by, bx);
        k  = k(isfinite(k));
        evt_with(:,c) = accum_k(k, K);

        for kk = 1:K
            m_sub = isfinite(bx) & isfinite(by) & ...
                sub2ind([numel(edges.y)-1, numel(edges.x)-1], by, bx) == kk;
            if any(m_sub)
                evt_with_times{kk,c} = ts( idxCore(m_sub) );
            end
        end
    end

    % --- WITHOUT (correct mapping of submask back to ts) ---
    if any(inWout)
        idxWout = find(inWout);
        bx = discretize(xAt(ts(inWout)), edges.x);
        by = discretize(yAt(ts(inWout)), edges.y);
        k  = sub2ind([numel(edges.y)-1, numel(edges.x)-1], by, bx);
        k  = k(isfinite(k));
        evt_without(:,c) = accum_k(k, K);

        for kk = 1:K
            m_sub = isfinite(bx) & isfinite(by) & ...
                sub2ind([numel(edges.y)-1, numel(edges.x)-1], by, bx) == kk;
            if any(m_sub)
                evt_without_times{kk,c} = ts( idxWout(m_sub) );
            end
        end
    end
end

% ---- convert to rates (Hz) with per-bin normalization ----
occW = max(occ_with(:),    eps);  % K×1
occU = max(occ_without(:), eps);  % K×1
% implicit expansion (R2016b+) or use bsxfun
rate_with_all    = evt_with   ./ occW;
rate_without_all = evt_without ./ occU;

% ---- quality gates ----
minOcc = max(o.MinOccSec, eps);
validOcc = (occ_with >= minOcc) & (occ_without >= minOcc);

validCells = false(nC,1);
for c=1:nC
    if sum(evt_with(:,c)) >= o.MinEvents && sum(evt_without(:,c)) >= o.MinEvents ...
            && any(validOcc)
        validCells(c) = true;
    end
end

% ---- scoring ----
scores = -inf(nC,1);
for c = 1:nC
    if ~validCells(c), continue; end
    a = rate_with_all(:,c);     % WITH
    b = rate_without_all(:,c);  % WITHOUT
    a = a(:); b = b(:);

    switch lower(o.ScoreMode)
        case 'cosine'
            aa=a; bb=b; aa(~validOcc)=0; bb(~validOcc)=0;
            an = aa./max(norm(aa),eps); bn = bb./max(norm(bb),eps);
            scores(c) = 1 - sum(an.*bn);

        case 'corr'
            aa=a; bb=b; aa(~validOcc)=0; bb(~validOcc)=0;
            if all(aa==0)||all(bb==0)
                scores(c) = -inf;
            else
                r = corr(aa, bb, 'type','Pearson','rows','complete');
                scores(c) = 1 - r;
            end

        case 'l2'
            aa=a; bb=b; aa(~validOcc)=0; bb(~validOcc)=0;
            scores(c) = norm(aa - bb, 2);

        case 'peakspread'
            % 1) Concentration WITHOUT task: max probability mass in a single bin
            bb = b; bb(~validOcc) = 0;
            sumB = sum(bb); if sumB<=0, scores(c) = -inf; continue; end
            pB = bb ./ sumB; conc = max(pB);

            % 2) Evenness WITH task: normalized entropy (0..1)
            aa = a; aa(~validOcc) = 0;
            sumA = sum(aa); if sumA<=0, scores(c) = -inf; continue; end
            pA = aa ./ sumA;
            H  = -nansum(pA .* log(pA + o.PeakSpreadPseudo));
            Hmax = log(max(1, nnz(validOcc)));
            even = min(1, H / max(Hmax, eps));

            % 3) Drop in dominant bin when task is present
            [~, kStar] = max(pB);
            rWout = bb(kStar); rWith = aa(kStar);
            switch lower(o.PeakSpreadDropMode)
                case 'ratio'
                    drop = max(0, (rWout+eps)/(rWith+eps) - 1);
                    dropN = min(1, drop/2);   % ~1 at ~3x drop
                case 'diff'
                    drop = max(0, rWout - rWith);
                    dropN = drop / (prctile(b(b>0), 90)+eps);
                    dropN = max(0, min(1, dropN));
                otherwise
                    error('PeakSpreadDropMode must be ''ratio'' or ''diff''.');
            end

            w = o.PeakSpreadAlpha(:).'; w = w / max(sum(w),eps);
            scores(c) = (conc^w(1)) * (even^w(2)) * (dropN^w(3));

        otherwise
            error('Unknown ScoreMode.');
    end
end

% ---- choose cells to plot ----
if ~isempty(o.CellIdx)
    % force a specific cell (or vector)
    cellIdx = o.CellIdx(:);
    scoreTop = scores(cellIdx);
    TopN = numel(cellIdx);
else
    [scoreSorted, idxSorted] = sort(scores, 'descend');
    keep = isfinite(scoreSorted) & scoreSorted > -inf;
    idxSorted = idxSorted(keep); scoreSorted = scoreSorted(keep);
    if isempty(idxSorted)
        warning('No cells passed quality/scoring gates. Try lowering MinOccSec/MinEvents.');
        OUT = struct('cellIdx',[],'scores',[]); return
    end
    TopN   = min(TopN, numel(idxSorted));
    cellIdx  = idxSorted(1:TopN);
    scoreTop = scoreSorted(1:TopN);
end

% ---- print summary of chosen examples ----
fprintf('\nTop %d example cell(s) (ranked by %s):\n', numel(cellIdx), o.ScoreMode);
T = table((1:numel(cellIdx))', cellIdx(:), scoreTop(:), ...
    'VariableNames', {'Rank','CellIndex','Score'});
disp(T);

% ---- prepare plotting data ----
R_with    = reshape(rate_with_all(:,cellIdx),    o.GridRC(1), o.GridRC(2), TopN);
R_without = reshape(rate_without_all(:,cellIdx), o.GridRC(1), o.GridRC(2), TopN);
mx = prctile([R_with(:); R_without(:)], 99);

% ---- optional PV-based CTRL vs WITH correlations (population) ----
% Only show these when we aren't focusing on a single cell
showPopCorr = o.ShowCorrPanels && (numel(cellIdx) > 1 || isempty(o.CellIdx));
if showPopCorr
    [r_with_k, r_without_k, r_with_overall, r_without_overall] = ...
        local_compute_C_quadrant_corr(RAT, dayLabel, ...
            'TraceWin',o.TraceWin,'VelThresh',o.VelThresh,'UseSpeedMask',o.UseSpeedMask, ...
            'GridRC',o.GridRC,'TimeBins',o.CorrTimeBins,'BufferPost',o.BufferPost);
end

% ---- layout: WITHOUT | WITH | rates | (pop) per-quad corr | (pop or cell) overall corr
nCols = 3 + (o.ShowCorrPanels ~= 0) + (o.ShowCorrPanels ~= 0);
if ~showPopCorr
    % For single-cell mode: only add the final "overall correlation" panel for the cell
    nCols = 4;
end

figure('Color','w','Position',[60 50 1200 max(500,240*TopN)]);
tiledlayout(TopN, nCols, 'Padding','compact','TileSpacing','compact');

for i = 1:TopN
    % WITHOUT
    nexttile; imagesc(R_without(:,:,i), [0 mx]); axis image off
    set(gca,'YDir','normal'); colormap(o.Colormap); colorbar
    title(sprintf('WITHOUT task | Cell %d', cellIdx(i)));

    % WITH
    nexttile; imagesc(R_with(:,:,i), [0 mx]); axis image off
    set(gca,'YDir','normal'); colormap(o.Colormap); colorbar
    title(sprintf('WITH task (CS..%.1fs) | Cell %d', o.TraceWin(2), cellIdx(i)));

    % Rates bars
    nexttile; hold on
    rw = R_with(:,:,i); ru = R_without(:,:,i);
    rw = rw(:); ru = ru(:); nb = numel(rw); bw = 0.38; x = 1:nb;
    bar(x-bw/2, ru, bw, 'FaceColor',[0.6 0.6 0.6], 'EdgeColor','k');
    bar(x+bw/2, rw, bw, 'FaceColor',[0.25 0.45 0.95], 'EdgeColor','k');
    ylabel('Rate (Hz)'); xlabel('Quadrant'); box on; xlim([0.5 nb+0.5]); yline(0,'k-');
    title(sprintf('Score = %.3f (%s)', scoreTop(i), o.ScoreMode));

    % Per-quad PV corr (population)
    if showPopCorr
        nexttile; hold on
        bw2 = 0.38; x2 = 1:numel(r_with_k);
        bar(x2-bw2/2, r_without_k, bw2, 'FaceColor',[0.65 0.65 0.65], 'EdgeColor','k');
        bar(x2+bw2/2, r_with_k,    bw2, 'FaceColor',[0.10 0.60 0.95], 'EdgeColor','k');
        ylim([-0.05 1]); yline(0,'k:');
        xticks(x2); xlabel('Quadrant'); ylabel('PV corr (r)');
        legend({'CTRL split-half','WITH task'},'Location','northoutside','Orientation','horizontal');
        title('Population PV corr by quadrant');
    end

    % Overall corr panel
    nexttile; hold on
    if showPopCorr
        % overall population bars
        bar(1, r_with_overall,    0.6, 'FaceColor',[0.10 0.60 0.95], 'EdgeColor','k');
        bar(2, r_without_overall, 0.6, 'FaceColor',[0.65 0.65 0.65], 'EdgeColor','k');
        xlim([0.5 2.5]); xticks([1 2]); xticklabels({'WITH','CTRL'});
        ylabel('PV corr (r)'); yline(0,'k:'); title('Population overall PV corr');
    else
      % ---- NEW: single-cell within-condition reliabilities (two bars) ----
      % Build even/odd frame masks inside each condition
      evenodd = false(size(tPos));
      evenodd(2:2:end) = true;      % even-indexed frames true

      % Base masks for conditions
      mask_with_all = maskWith;         % task frames (speed-filtered)
      mask_wout_all = maskWithout;      % control frames (speed-filtered)

      % Split each into even/odd sets
      mask_with_even = mask_with_all &  evenodd;
      mask_with_odd  = mask_with_all & ~evenodd;
      mask_wout_even = mask_wout_all &  evenodd;
      mask_wout_odd  = mask_wout_all & ~evenodd;

      % Build two rate maps per condition (even vs odd) for THIS cell
      cidx = cellIdx(i);                         % current cell index
      ts_cell = double(cellSpk{cidx});           % spike times for this cell

      % WITH task
      rm_with_even = rate_map_for_mask(ts_cell, mask_with_even, tPos, xPos, yPos, edges);
      rm_with_odd  = rate_map_for_mask(ts_cell, mask_with_odd,  tPos, xPos, yPos, edges);
      % WITHOUT task
      rm_wout_even = rate_map_for_mask(ts_cell, mask_wout_even, tPos, xPos, yPos, edges);
      rm_wout_odd  = rate_map_for_mask(ts_cell, mask_wout_odd,  tPos, xPos, yPos, edges);

      % Correlations (within condition)
      r_with_reli  = safe_corr(rm_with_even, rm_with_odd,  2, 1e-12);
      r_wout_reli  = safe_corr(rm_wout_even, rm_wout_odd,  2, 1e-12);
      if ~isfinite(r_with_reli), r_with_reli = 0; end
      if ~isfinite(r_wout_reli), r_wout_reli = 0; end

      % Plot two bars: WITH and WITHOUT within-condition reliability
      bar(1, r_with_reli, 0.6, 'FaceColor',[0.25 0.45 0.95], 'EdgeColor','k'); hold on
      bar(2, r_wout_reli, 0.6, 'FaceColor',[0.60 0.60 0.60], 'EdgeColor','k');
      xlim([0.5 2.5]); xticks([1 2]); xticklabels({'WITH','WITHOUT'});
      ylabel('Within-condition rate-map corr (r)'); yline(0,'k:');
      title(sprintf('Cell %d reliability (even vs odd)', cidx));
    end

    % Optional per-cell PNG (one figure per cell)
    if ~isempty(o.SavePNG)
        if ~exist(o.SavePNG,'dir'), mkdir(o.SavePNG); end
        print(gcf, fullfile(o.SavePNG, sprintf('%s_%s_cell%04d.png', ratVar, dayLabel, cellIdx(i))), '-dpng','-r200');
    end
end

sgtitle(sprintf('%s  |  %s  |  Grid %dx%d  |  speed\\geq%.1f cm/s=%d', ...
    ratVar, strrep(dayLabel,'_','\_'), o.GridRC(1), o.GridRC(2), o.VelThresh, o.UseSpeedMask));

% ---- pack outputs ----
OUT = struct();
OUT.cellIdx      = cellIdx(:);
OUT.scores       = scoreTop(:);
OUT.rate_with    = rate_with_all(:,cellIdx);
OUT.rate_without = rate_without_all(:,cellIdx);
OUT.occ_with     = occ_with(:);
OUT.occ_without  = occ_without(:);
OUT.meta = struct('rat',ratVar,'day',dayLabel,'GridRC',o.GridRC,'TraceWin',o.TraceWin, ...
                  'VelThresh',o.VelThresh,'UseSpeedMask',o.UseSpeedMask,'ScoreMode',o.ScoreMode);

end

% ======================= helpers ========================
function [t,x,y] = coerce_pos(pos)
if istable(pos)
    vn = lower(string(pos.Properties.VariableNames));
    t = pos{:, find(ismember(vn, {'t','time','ts'}),1)};
    x = pos{:, find(ismember(vn, {'x','xpos','x_cm','xcm','xsmooth','posx','x_smooth'}),1)};
    y = pos{:, find(ismember(vn, {'y','ypos','y_cm','ycm','ysmooth','posy','y_smooth'}),1)};
elseif isstruct(pos)
    f = lower(string(fieldnames(pos)));
    t = pos.(f{find(ismember(f,{'t','time','ts'}),1)});
    x = pos.(f{find(ismember(f,{'x','xpos','x_cm','xcm','xsmooth','posx','x_smooth'}),1)});
    y = pos.(f{find(ismember(f,{'y','ypos','y_cm','ycm','ysmooth','posy','y_smooth'}),1)});
elseif isnumeric(pos) && size(pos,2)>=3
    t = pos(:,1); x = pos(:,2); y = pos(:,3);
else
    error('Unsupported pos format.');
end
t = double(t(:)); x = double(x(:)); y = double(y(:));
if max(t,[],'omitnan')>1e4, t = t/1000; end
end

function v = speed_cm_per_s(t, x, y, tQuery)
dt = diff(t); dt(end+1,1) = median(dt(dt>0),'omitnan');
dx = [diff(x); 0]; dy = [diff(y); 0];
spd = hypot(dx,dy) ./ max(dt, eps);
if nargin<4, v = spd; else, v = interp1_nonan(t, spd, tQuery); end
end

function val = interp1_nonan(t, y, tq)
[tU, ia] = unique(t,'stable'); yU = y(ia);
val = interp1(tU, yU, tq, 'linear','extrap');
end

function [edges, K] = make_grid_edges(x, y, rc)
edges.x = linspace(min(x,[],'omitnan'), max(x,[],'omitnan'), rc(2)+1);
edges.y = linspace(min(y,[],'omitnan'), max(y,[],'omitnan'), rc(1)+1);
K = rc(1)*rc(2);
end

function occ = occupancy2D(t, x, y, mask, edges, dt)
xM = x(mask); yM = y(mask);
bx = discretize(xM, edges.x); by = discretize(yM, edges.y);
K  = (numel(edges.x)-1) * (numel(edges.y)-1);
k  = sub2ind([numel(edges.y)-1, numel(edges.x)-1], by, bx);
k  = k(isfinite(k));
occ = accumarray(k, dt, [K,1], @sum, 0);
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

function v = accum_k(k, K)
v = accumarray(k, 1, [K,1], @sum, 0);
end

function r = safe_corr(a,b,minN,minStd)
if nargin<3, minN=3; end
if nargin<4, minStd=1e-12; end
a=a(:); b=b(:);
m = isfinite(a)&isfinite(b);
a=a(m); b=b(m);
if numel(a)<minN, r=NaN; return; end
sa=std(a); sb=std(b);
if ~isfinite(sa)||~isfinite(sb)||sa<minStd||sb<minStd, r=0; return; end
r = corr(a,b,'type','Pearson');
end

% ---------- Population PV correlations by quadrant + overall ----------
function [r_with_k, r_without_k, r_with_overall, r_without_overall] = ...
    local_compute_C_quadrant_corr(RAT, dayLabel, varargin)

q = inputParser;
addParameter(q,'TraceWin',[0 2]);
addParameter(q,'BufferPost',0);
addParameter(q,'GridRC',[2 2]);
addParameter(q,'TimeBins',15);
addParameter(q,'VelThresh',4);
addParameter(q,'UseSpeedMask',true);
parse(q,varargin{:});
o = q.Results;

% Pull standardized containers similar to run_task_to_space_interference
posd = RAT.pos.(sprintf('pos_%s', dayLabel));
spk  = RAT.Ca_peaks.(sprintf('CA_peaks_%s', dayLabel));
ts   = RAT.Ca_ts.(sprintf('Ca_ts_%s', dayLabel));
cs   = RAT.CS_times.(sprintf('CS_%s', dayLabel));

% t-aligned spikes matrix
t = double(ts(:));
[x, y] = interp_pos_local(posd, t);
S = spikes_to_matrix_local(spk, t);  % [Nc x T]
S = normalize_cells_local(S, 'demean');

% speed mask (optional)
if o.UseSpeedMask
    v = speed_cm_per_s_local(posd);                % at pos timestamps
    v_i = interp1(double(posd.t(:)), double(v(:)), t, 'linear','extrap');
    speed_ok = (v_i >= o.VelThresh);
else
    speed_ok = true(size(t));
end

% grid edges for this day
edges = build_grid_edges_single_day_local(posd, o.GridRC);
GridR = o.GridRC(1); GridC = o.GridRC(2); K = GridR*GridC;

% ---------- CTRL PV per bin (non-task frames) ----------
is_taskbuf = false(size(t));
for j = 1:numel(cs)
    is_taskbuf = is_taskbuf | (t >= cs(j)+o.TraceWin(1) & t <= cs(j)+o.TraceWin(2)+o.BufferPost);
end
use_ctrl = ~is_taskbuf & speed_ok;

% Precompute control-frame indices once
idx_ctrl = find(use_ctrl);

% Bin assignments for control frames
[~, b_ctrl] = pos2bin_local(x(idx_ctrl), y(idx_ctrl), edges);

PV_ctrl = nan(size(S,1), K);
for k = 1:K
    idxk = find(b_ctrl == k);             % indices into idx_ctrl
    if numel(idxk) >= 2
        cols = idx_ctrl(idxk);            % actual columns in S
        PV_ctrl(:,k) = mean(S(:, cols), 2, 'omitnan');
    end
end

% split-half within CTRL for each quadrant
r_without_k = nan(K,1);
for k = 1:K
    idxk = find(b_ctrl == k);             % into idx_ctrl
    if numel(idxk) < 4, continue; end
    cols_even = idx_ctrl( idxk(2:2:end) );
    cols_odd  = idx_ctrl( idxk(1:2:end) );
    v1 = mean(S(:, cols_even), 2, 'omitnan');
    v2 = mean(S(:, cols_odd ), 2, 'omitnan');
    r_without_k(k) = safe_corr(v1, v2, 3, 1e-12);
end
r_without_overall = nanmean(r_without_k);

% ---------- WITH-task PV per (bin,time-slice) then average ----------
r_with_k = nan(K,1);
B = o.TimeBins;
r_k_accum = nan(K, B);

for tr = 1:numel(cs)
    w0 = cs(tr) + o.TraceWin(1);
    w1 = cs(tr) + o.TraceWin(2);
    idx_all = find(t >= w0 & t < w1);
    if isempty(idx_all), continue; end
    idx_all = idx_all(speed_ok(idx_all));
    if isempty(idx_all), continue; end

    edges_idx = round(linspace(1, numel(idx_all)+1, B+1));
    edges_idx(1) = 1; edges_idx(end) = numel(idx_all)+1;

    for bti = 1:(numel(edges_idx)-1)
        if edges_idx(bti) >= edges_idx(bti+1), continue; end
        seg = idx_all(edges_idx(bti):edges_idx(bti+1)-1);

        [~, b_with] = pos2bin_local(x(seg), y(seg), edges);
        for k = 1:K
            hit = (b_with == k);
            if nnz(hit) >= 2 && all(isfinite(PV_ctrl(:,k)))
                pv_task = mean(S(:, seg(hit)), 2, 'omitnan');
                r = safe_corr(PV_ctrl(:,k), pv_task, 3, 1e-12);
                r_k_accum(k,bti) = nanmean([r_k_accum(k,bti), r]);
            end
        end
    end
end

% average across time-bins then summarize by quadrant + overall
for k = 1:K
    rv = r_k_accum(k,:);
    rv = rv(isfinite(rv));
    if ~isempty(rv)
        r_with_k(k) = tanh(mean(atanh(max(min(rv,0.99999),-0.99999))));
    end
end
r_with_overall = nanmean(r_with_k);
end

% ---- minimal locals mirroring your pipeline ----
function [x_i, y_i] = interp_pos_local(posd, t)
tt = double(posd.t(:)); xx = double(posd.x(:)); yy = double(posd.y(:));
[ttu, ia] = unique(tt, 'stable'); xxu = xx(ia); yyu = yy(ia);
x_i = interp1(ttu, xxu, t, 'linear','extrap');
y_i = interp1(ttu, yyu, t, 'linear','extrap');
end

function v = speed_cm_per_s_local(posd)
t  = double(posd.t(:)); x  = double(posd.x(:)); y  = double(posd.y(:));
n  = min([numel(t), numel(x), numel(y)]);
dt = diff(t(1:n)); dt(end+1,1) = median(dt(dt>0),'omitnan');
dx = [diff(x(1:n)); 0]; dy = [diff(y(1:n)); 0];
v = hypot(dx,dy) ./ max(dt, eps);
end

function [rc_idx, k] = pos2bin_local(x, y, edges)
cx = discretize(x, edges.x); cy = discretize(y, edges.y);
GridR = numel(edges.y)-1; GridC = numel(edges.x)-1;
bad = isnan(cx) | isnan(cy) | cx<1 | cy<1 | cx>GridC | cy>GridR;
cx(bad) = NaN; cy(bad) = NaN;
rc_idx = [cy, cx]; k = nan(size(x));
m = isfinite(cx) & isfinite(cy);
if any(m), k(m) = sub2ind([GridR, GridC], cy(m), cx(m)); end
end

function edges = build_grid_edges_single_day_local(posd, GridRC)
allx = posd.x(:); ally = posd.y(:);
allx = allx(isfinite(allx));  ally = ally(isfinite(ally));
if isempty(allx) || isempty(ally)
    edges.x = linspace(0, 1, GridRC(2)+1);
    edges.y = linspace(0, 1, GridRC(1)+1);
else
    edges.x = linspace(min(allx), max(allx), GridRC(2)+1);
    edges.y = linspace(min(ally), max(ally), GridRC(1)+1);
end
end

function S = spikes_to_matrix_local(daySpikes, t)
% Accepts cell-of-vectors or numeric [nc x ?] with spike times.
daySpikes = unwrap_spike_container_local(daySpikes);
if iscell(daySpikes), Nc=numel(daySpikes);
elseif isnumeric(daySpikes), Nc=size(daySpikes,1);
else, error('Unsupported spikes container');
end
t=double(t(:));
if numel(t)<2, S=zeros(Nc, numel(t), 'single'); return; end
dt=median(diff(t)); if ~isfinite(dt)||dt<=0, dt=max(eps, mean(diff(t),'omitnan')); end
edges=[t - dt/2; t(end)+dt/2];
S=zeros(Nc, numel(t), 'single');
for c=1:Nc
    st = extract_cell_spikes_local(daySpikes, c);
    if isempty(st), continue; end
    S(c,:) = histcounts(st, edges);
end
end

function obj = unwrap_spike_container_local(S)
if isstruct(S)
    f=fieldnames(S); if isempty(f), obj=[]; return; end
    pick=[]; for j=1:numel(f)
        name=lower(f{j});
        if contains(name,'peak') || contains(name,'spike') || contains(name,'ca_peaks'), pick=j; break, end
    end
    if isempty(pick), pick=1; end
    obj = S.(f{pick});
else
    obj=S;
end
end

function st = extract_cell_spikes_local(container, c)
if iscell(container), st=container{c};
elseif isnumeric(container), if c>size(container,1), st=[]; return; end, st=container(c,:).';
else, st=[]; end
st=double(st(:)); st=st(isfinite(st) & st>0);
end

function S = normalize_cells_local(S, mode)
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

function idx = nearest_frame_index(t_frames, t_events)
% Map each event time to the nearest frame index (1..numel(t_frames))
idx = round(interp1(t_frames, 1:numel(t_frames), t_events, 'linear', 'extrap'));
idx = max(1, min(numel(t_frames), idx));
end

function rateK = rate_map_for_mask(ts_cell, mask_frames, tPos, xPos, yPos, edges)
% Build a single-cell rate map (Kx1) using only frames where mask_frames==true.
% Occupancy comes from masked frames; events restricted to times whose nearest
% frame is in the mask. Returned as Hz (events/sec).
dt = median(diff(tPos),'omitnan'); if ~isfinite(dt)||dt<=0, dt = max(eps, mean(diff(tPos),'omitnan')); end

% Occupancy per bin
xM = xPos(mask_frames); yM = yPos(mask_frames);
bx = discretize(xM, edges.x);
by = discretize(yM, edges.y);
GridR = numel(edges.y)-1; GridC = numel(edges.x)-1; K = GridR*GridC;
kOcc = sub2ind([GridR, GridC], by, bx);
kOcc = kOcc(isfinite(kOcc));
occ = accumarray(kOcc, dt, [K,1], @sum, 0);

% Event counts per bin (restrict events to masked frames)
if isempty(ts_cell), rateK = zeros(K,1); return; end
evt_idx = nearest_frame_index(tPos, ts_cell);
keep_evt = mask_frames(evt_idx);
te = ts_cell(keep_evt);

if isempty(te)
    rateK = zeros(K,1);
else
    % Bin events by position at event times
    xE = interp1(tPos, xPos, te, 'linear', 'extrap');
    yE = interp1(tPos, yPos, te, 'linear', 'extrap');
    bxE = discretize(xE, edges.x);
    byE = discretize(yE, edges.y);
    kEvt = sub2ind([GridR, GridC], byE, bxE);
    kEvt = kEvt(isfinite(kEvt));
    evt = accumarray(kEvt, 1, [K,1], @sum, 0);
    % Convert to Hz per bin (events/sec), with per-bin occupancy normalization
    occ_safe = occ; occ_safe(occ_safe<=0) = Inf;
    rateK = evt ./ occ_safe;
end
end
