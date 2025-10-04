function T = grid_trialVsSpeedMatched(grid, varargin)
% Grid over plotFRvsSpeedMatched (internally uses trialVsSpeedMatched).
% Early-quit rule: if ANY (rat,day) has < MinDayCells paired cells, the worker returns early and we skip the combo.

% ---------- options ----------
start_idx = 8;
if ~isempty(varargin)
    for k = 1:2:numel(varargin)
        if strcmpi(varargin{k}, 'start_idx')
            start_idx = varargin{k+1};
        end
    end
end

% ---------- defaults ----------
if nargin<1 || isempty(grid), grid = struct(); end
if ~isfield(grid,'CompareMode'),       grid.CompareMode = {'bins','seconds'}; end
if ~iscell(grid.CompareMode),          grid.CompareMode = cellstr(grid.CompareMode); end

if ~isfield(grid,'binSize'),           grid.binSize = [2 1 .66 .5 .25 1/7.5]; end   % for bins-mode
if ~isfield(grid,'AllowReuseControl'), grid.AllowReuseControl = [true]; end
if ~isfield(grid,'MinPairsFrac'),      grid.MinPairsFrac = [0.25 .5 .66]; end
if ~isfield(grid,'test'),              grid.test = {'ttest'}; end
if ~iscell(grid.test),                 grid.test = cellstr(grid.test); end
if ~isfield(grid,'alpha'),             grid.alpha = 0.05; end
if ~isfield(grid,'win'),               grid.win = {[0 2]}; end
if ~iscell(grid.win),                  grid.win = {grid.win}; end
if ~isfield(grid,'ratNames'),          grid.ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'}; end
if ~isfield(grid,'doPlots'),           grid.doPlots = false; end
if ~isfield(grid,'MinDayCells'),       grid.MinDayCells = 50; end  % EARLY-QUIT THRESHOLD

% vectorize
v_MODE  = grid.CompareMode(:).';
v_BS    = grid.binSize(:).';
v_REUSE = grid.AllowReuseControl(:).';
v_MPF   = grid.MinPairsFrac(:).';
v_TEST  = grid.test(:).';
v_A     = grid.alpha(:).';
v_WIN   = grid.win(:).';

% ---------- build jobs ----------
jobs = struct('CompareMode',{},'binSize',{},'AllowReuse',{}, ...
              'MinPairsFrac',{},'test',{},'alpha',{},'win',{});
for iM = 1:numel(v_MODE)
    mode = lower(string(v_MODE{iM}));
    switch mode
        case "bins"
            if isempty(v_BS), bsz = 1/7.5; else, bsz = v_BS; end
            for iB = 1:numel(bsz)
              for iR = 1:numel(v_REUSE)
                for iP = 1:numel(v_MPF)
                  for iT = 1:numel(v_TEST)
                    for iA = 1:numel(v_A)
                      for iW = 1:numel(v_WIN)
                        jobs(end+1) = struct( ... %#ok<AGROW>
                          'CompareMode','bins', ...
                          'binSize',     bsz(iB), ...
                          'AllowReuse',  v_REUSE(iR), ...
                          'MinPairsFrac',v_MPF(iP), ...
                          'test',        v_TEST{iT}, ...
                          'alpha',       v_A(iA), ...
                          'win',         v_WIN{iW});
                      end
                    end
                  end
                end
              end
            end

        case "seconds"
            for iR = 1:numel(v_REUSE)
              for iP = 1:numel(v_MPF)
                for iT = 1:numel(v_TEST)
                  for iA = 1:numel(v_A)
                    for iW = 1:numel(v_WIN)
                      jobs(end+1) = struct( ... %#ok<AGROW>
                        'CompareMode','seconds', ...
                        'binSize',     NaN, ...
                        'AllowReuse',  v_REUSE(iR), ...
                        'MinPairsFrac',v_MPF(iP), ...
                        'test',        v_TEST{iT}, ...
                        'alpha',       v_A(iA), ...
                        'win',         v_WIN{iW});
                    end
                  end
                end
              end
            end

        otherwise
            error('Unknown CompareMode: %s', mode);
    end
end

N = numel(jobs);
start_idx = max(1, min(N, round(ifempty(start_idx,1))));
if start_idx>1
    fprintf('Resuming at combo %d/%d\n', start_idx, N);
end

rows = cell(N, 15);

% ---------- run ----------
for k = start_idx:N
    P = jobs(k);

    % header
    if strcmpi(P.CompareMode,'bins')
        tag = sprintf('MODE=bins, binSize=%.4g s, reuse=%d, MinPairsFrac=%.3g, test=%s, alpha=%.3g, win=[%.3g %.3g]', ...
            P.binSize, P.AllowReuse, P.MinPairsFrac, P.test, P.alpha, P.win(1), P.win(2));
    else
        tag = sprintf('MODE=seconds, reuse=%d, MinPairsFrac=%.3g, test=%s, alpha=%.3g, win=[%.3g %.3g]', ...
            P.AllowReuse, P.MinPairsFrac, P.test, P.alpha, P.win(1), P.win(2));
    end
    fprintf('Combo %d/%d: %s\n', k, N, tag);

    % worker call (plotting off)
    try
        Res = plotFRvsSpeedMatched(P.CompareMode, ...
                'ratNames',          grid.ratNames, ...
                'win',               P.win, ...
                'binSize',           ifempty(P.binSize, 1/7.5), ...
                'alpha',             P.alpha, ...
                'test',              P.test, ...
                'MinPairsFrac',      P.MinPairsFrac, ...
                'AllowReuseControl', P.AllowReuse, ...
                'MinDayCells',       grid.MinDayCells, ...   % NEW
                'EarlyQuitOnLowN',   true, ...               % NEW
                'doPlots',           false);
    catch ME
        warning('  -> combo failed: %s', ME.message);
        continue;
    end

    % Early-quit?
    if isfield(Res,'earlyQuit') && Res.earlyQuit
        info = Res.lowN;
        fprintf('  -> early-quit: %s %s had n=%d (< %d). Skipping this combo.\n\n', ...
            info.rat, info.dayLabel, info.nCells, grid.MinDayCells);
        continue;
    end

    % Summaries across all (rat,day) cells
    [nR, nD] = size(Res.dayTestsByRat);
    dayNs    = nan(nR,nD);
    dayDelta = nan(nR,nD);
    dayP     = nan(nR,nD);
    for r=1:nR
        for d=1:nD
            dayNs(r,d)    = Res.dayTestsByRat(r,d).nCells;
            dayDelta(r,d) = Res.dayTestsByRat(r,d).deltaMeanHz;
            dayP(r,d)     = Res.dayTestsByRat(r,d).p;
        end
    end
    daySigFrac = mean(dayP(:) < P.alpha, 'omitnan');

    % console recap (optional)
    lbl = {'day-2','day-1','day 0'};
    for r=1:nR
        for d=1:nD
            if dayNs(r,d)>0 && isfinite(dayP(r,d))
                fprintf('    [%s %s] n=%d  Δ=%.4f Hz  p=%.3g\n', ...
                    grid.ratNames{r}, lbl{d}, dayNs(r,d), dayDelta(r,d), dayP(r,d));
            end
        end
    end
    if isfield(Res,'allRatsPooled') && isfinite(Res.allRatsPooled.p)
        fprintf('    [pooled ALL] n=%d  Δ=%.4f Hz  p=%.3g\n', ...
            Res.allRatsPooled.nCells, Res.allRatsPooled.deltaMeanHz, Res.allRatsPooled.p);
    end
    fprintf('    -> meanPctIncluded=%.2f%%, day p̃=%.3g, day sig frac=%.2f\n\n', ...
        Res.meanPctIncluded, median(dayP(:),'omitnan'), daySigFrac);

    % store row
    rows{k,1}  = char(P.CompareMode);
    rows{k,2}  = P.binSize;
    rows{k,3}  = P.AllowReuse;
    rows{k,4}  = P.MinPairsFrac;
    rows{k,5}  = P.test;
    rows{k,6}  = P.alpha;
    rows{k,7}  = sprintf('[%.3g %.3g]', P.win(1), P.win(2));
    rows{k,8}  = Res.meanPctIncluded;
    rows{k,9}  = daySigFrac;
    rows{k,10} = dayNs;                      % matrix nRats x 3
    rows{k,11} = dayDelta;                   % matrix nRats x 3
    rows{k,12} = dayP;                       % matrix nRats x 3
    % pooled summaries
    rows{k,13} = Res.allRatsPooled.nCells;
    rows{k,14} = Res.allRatsPooled.deltaMeanHz;
    rows{k,15} = Res.allRatsPooled.p;
end

% keep only filled rows
keep = ~cellfun(@isempty, rows(:,1));
rows = rows(keep,:);

T = cell2table(rows, 'VariableNames', { ...
    'CompareMode','binSize','AllowReuse','MinPairsFrac','test','alpha','win', ...
    'meanPctIncluded','daySigFrac','day_nCells','day_deltaHz','day_p', ...
    'pooled_n','pooled_deltaHz','pooled_p'});
end

% ---- helpers ----
function out = ifempty(x, alt)
    if isempty(x) || (isnumeric(x)&&any(~isfinite(x(:)))&&isscalar(alt)), out = alt; else, out = x; end
end
