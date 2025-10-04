function T = gridSig_speedBinMatched(ratNames, grid, varargin)
% GRID search wrapper for run_speedBinMatched.
% Supports sweeping *either* fixed SpeedBinWidth or SpeedBinCount.
% Prints per-combo day-level stats (paired t across cells) for each rat/day.
%
% Example:
%   grid.SpeedBinMode  = {'width','count'};
%   grid.SpeedBinWidth = [1 2 3];
%   grid.SpeedBinCount = [5 7 9];
%   grid.SpeedRange    = {'trial','all'};
%   grid.MinDurPerBin  = [0.5 1 2];
%   grid.MinBins       = [5 7 9];
%   grid.SpeedEdges    = {[]};   % explicit edges still supported (edges win)
%   grid.test          = {'ttest'};
%   grid.win           = {[0 2]};
%   grid.binSize       = [1/7.5 1 2];
%   T = gridSig_speedBinMatched([], grid);

if nargin<1 || isempty(ratNames)
    ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end
if nargin<2 || isempty(grid), grid = struct(); end

set(0,'DefaultFigureVisible','off');

% ----- optional args -----
start_idx = 1;
if ~isempty(varargin)
    for ii = 1:2:numel(varargin)
        if strcmpi(varargin{ii},'start_idx')
            start_idx = varargin{ii+1};
        end
    end
end

% ----- defaults -----
if ~isfield(grid,'SpeedBinMode'),   grid.SpeedBinMode   = {'width'}; end
if ~iscell(grid.SpeedBinMode),      grid.SpeedBinMode   = cellstr(grid.SpeedBinMode); end

if ~isfield(grid,'SpeedBinWidth'),  grid.SpeedBinWidth  = [1 2 3 4]; end
if ~isfield(grid,'SpeedBinCount'),  grid.SpeedBinCount  = [1 2 3 4 5]; end

if ~isfield(grid,'SpeedRange'),     grid.SpeedRange     = {'trial'}; end
if ~iscell(grid.SpeedRange),        grid.SpeedRange     = cellstr(grid.SpeedRange); end

if ~isfield(grid,'SpeedEdges'),     grid.SpeedEdges     = {[]}; end
if ~iscell(grid.SpeedEdges),        grid.SpeedEdges     = num2cell(grid.SpeedEdges); end

if ~isfield(grid,'MinDurPerBin'),   grid.MinDurPerBin   = [1/7.5 .25]; end
if ~isfield(grid,'MinBins'),        grid.MinBins        = [1 2 3 4 5]; end
if ~isfield(grid,'alpha'),          grid.alpha          = 0.05; end
if ~isfield(grid,'test'),           grid.test           = {'ttest'}; end
if ~iscell(grid.test),              grid.test           = cellstr(grid.test); end
if ~isfield(grid,'win'),            grid.win            = {[0 2]}; end
if ~iscell(grid.win),               grid.win            = {grid.win}; end
if ~isfield(grid,'binSize'),        grid.binSize        = [1]; end %%%NOT IN USE

% ----- build a job list (each item is a concrete parameter combo) -----
v_MODE  = grid.SpeedBinMode(:).';
v_SBW   = grid.SpeedBinWidth(:).';
v_SBC   = grid.SpeedBinCount(:).';
v_SR    = grid.SpeedRange(:).';
v_SE    = grid.SpeedEdges(:).';
v_MDPB  = grid.MinDurPerBin(:).';
v_MB    = grid.MinBins(:).';
v_A     = grid.alpha(:).';
v_TEST  = grid.test(:).';
v_WIN   = grid.win(:).';
v_BS    = grid.binSize(:).';

% Pre-init to avoid double->struct conversion errors
jobs = struct('Mode',{},'SpeedBinWidth',{},'SpeedBinCount',{},'SpeedEdges',{}, ...
              'SpeedRange',{},'MinDurPerBin',{},'MinBins',{},'alpha',{}, ...
              'test',{},'win',{},'binSize',{});

for m = 1:numel(v_MODE)
    mode = lower(string(v_MODE{m}));
    for iSE = 1:numel(v_SE)
      for iSR = 1:numel(v_SR)
        for iMD = 1:numel(v_MDPB)
          for iMB = 1:numel(v_MB)
            for iA  = 1:numel(v_A)
              for iT  = 1:numel(v_TEST)
                for iW  = 1:numel(v_WIN)
                  for iBS = 1:numel(v_BS)

                    switch mode
                      case "width"
                        widths = v_SBW(:).';
                        widths = widths(isfinite(widths) & widths>0);
                        if isempty(widths), continue; end
                        for iSBW = 1:numel(widths)
                          jobs(end+1) = struct( ... %#ok<AGROW>
                            'Mode',         'width', ...
                            'SpeedBinWidth', widths(iSBW), ...
                            'SpeedBinCount', [], ...
                            'SpeedEdges',    v_SE{iSE}, ...
                            'SpeedRange',    v_SR{iSR}, ...
                            'MinDurPerBin',  v_MDPB(iMD), ...
                            'MinBins',       v_MB(iMB), ...
                            'alpha',         v_A(iA), ...
                            'test',          v_TEST{iT}, ...
                            'win',           v_WIN{iW}, ...
                            'binSize',       v_BS(iBS) );
                        end

                      case "count"
                        counts = v_SBC(:).';
                        counts = counts(isfinite(counts) & counts>=1);
                        counts = round(counts);
                        if isempty(counts), continue; end
                        for iSBC = 1:numel(counts)
                          jobs(end+1) = struct( ... %#ok<AGROW>
                            'Mode',         'count', ...
                            'SpeedBinWidth', [], ...
                            'SpeedBinCount', counts(iSBC), ...
                            'SpeedEdges',    v_SE{iSE}, ...
                            'SpeedRange',    v_SR{iSR}, ...
                            'MinDurPerBin',  v_MDPB(iMD), ...
                            'MinBins',       v_MB(iMB), ...
                            'alpha',         v_A(iA), ...
                            'test',          v_TEST{iT}, ...
                            'win',           v_WIN{iW}, ...
                            'binSize',       v_BS(iBS) );
                        end

                      otherwise
                        error('Unknown SpeedBinMode: %s', mode);
                    end

                  end
                end
              end
            end
          end
        end
      end
    end
end

% ----- after building jobs -----
Njobs = numel(jobs);
if ~isscalar(start_idx) || ~isfinite(start_idx), start_idx = 1; end
start_idx = max(1, min(Njobs, round(start_idx)));
if start_idx > 1
    fprintf('Resuming at combo %d/%d\n', start_idx, Njobs);
end

rows = cell(Njobs, 14);

% ----- run the grid -----
for jobIdx = start_idx:Njobs
    P = jobs(jobIdx);

    % label bits
    if strcmpi(P.Mode,'width')
        modeLabel = sprintf('SBW=%.3g', P.SpeedBinWidth);
    else
        modeLabel = sprintf('SBC=%g',   P.SpeedBinCount);
    end
    seStr = shortStrSE(P.SpeedEdges);

    % pass-through args (empty means "use default inside")
    sbwArg = [];
    if ~isempty(P.SpeedBinWidth) && isfinite(P.SpeedBinWidth) && P.SpeedBinWidth > 0
        sbwArg = P.SpeedBinWidth;
    end
    sbcArg = [];
    if ~isempty(P.SpeedBinCount) && isfinite(P.SpeedBinCount) && P.SpeedBinCount >= 1
        sbcArg = round(P.SpeedBinCount);
    end

    % run
    Rk = run_speedBinMatched(ratNames, ...
        'SpeedEdges',    P.SpeedEdges, ...
        'SpeedBinWidth', sbwArg, ...
        'SpeedBinCount', sbcArg, ...
        'SpeedRange',    P.SpeedRange, ...
        'MinDurPerBin',  P.MinDurPerBin, ...
        'MinBins',       P.MinBins, ...
        'alpha',         P.alpha, ...
        'test',          P.test, ...
        'win',           P.win, ...
        'binSize',       P.binSize);

    % header for this combo
    fprintf('Combo %d/%d: %s, SE=%s, Range=%s, MinDur=%.3g, MinBins=%d, binSize=%.4g\n', ...
            jobIdx, Njobs, modeLabel, seStr, char(P.SpeedRange), ...
            P.MinDurPerBin, P.MinBins, P.binSize);

    % ---------- OPTION C: abort whole combo if ANY day has DL.delta < 0 ----------
    abort_combo = false;
    day_deltas = []; day_ps = []; day_ns = [];
    day_labels = {'day-2','day-1','day 0'};

    for rr = 1:numel(Rk)
        % guard early-exit shapes
        rname = '(unknown)';
        if isfield(Rk(rr),'rat') && ~isempty(Rk(rr).rat), rname = Rk(rr).rat; end
        if ~isfield(Rk(rr),'perDay') || isempty(Rk(rr).perDay) || ~iscell(Rk(rr).perDay)
            continue;
        end

        nDays = min(3, numel(Rk(rr).perDay));
        for d = 1:nDays
            S = Rk(rr).perDay{d};
            if isempty(S) || ~isfield(S,'dayLevel') || isempty(S.dayLevel), continue; end
            DL = S.dayLevel;
            if ~isfield(DL,'p') || ~isfinite(DL.p) || ~isfield(DL,'nCells') || DL.nCells < 2
                continue;
            end

            % check BEFORE printing/accumulating
            if isfield(DL,'delta') && isfinite(DL.delta) && DL.delta < 0
                fprintf('   [%s %s] Δ<0 — skipping this combo\n', rname, day_labels{d});
                abort_combo = true;
                break; % break day loop
            elseif DL.nCells<100
              fprintf('Combo %d/%d: skipping bc combo with not a good num of cells \n', jobIdx, Njobs);
            elseif isfield(DL,'delta') && isfinite(DL.delta) && DL.delta > 0
              fprintf('Combo %d/%d:   -> good news so far!\n', jobIdx, Njobs);
            end

            % normal print + accumulate
            fprintf('   [%s %s] day-level: n=%d  Δ=%.4f Hz  p=%.3g\n', ...
                    rname, day_labels{d}, DL.nCells, DL.delta, DL.p);
            day_deltas(end+1,1) = double(DL.delta); %#ok<AGROW>
            day_ps(end+1,1)     = double(DL.p);     %#ok<AGROW>
            day_ns(end+1,1)     = double(DL.nCells);%#ok<AGROW>
        end
        if abort_combo, break; end  % break rat loop
    end

    if abort_combo
        fprintf('Combo %d/%d:   -> skipped (Δ<0 encountered)\n\n', jobIdx, Njobs);
        continue;  % do NOT store a row; will be dropped by the filter below
    end

    % ---------- metrics & summary (only if we didn't abort) ----------
    % metric 1: average % sig across rats
    meanPct = NaN;
    if ~isempty(Rk) && all(isfield(Rk,{'pctSig'}))
        pctAll = vertcat(Rk.pctSig);
        if ~isempty(pctAll), meanPct = double(mean(pctAll,'omitnan')); end
    end

    % metric 2: “rats significant”
    rats_sig = 0; rats_tot = numel(Rk);
    for rr = 1:numel(Rk)
        if ~isfield(Rk(rr),'perDay') || ~iscell(Rk(rr).perDay), continue; end
        del = [];
        for d = 1:min(3,numel(Rk(rr).perDay))
            S = Rk(rr).perDay{d};
            if isempty(S) || ~isfield(S,'pVal'), continue; end
            tested = isfinite(S.pVal);
            if any(tested), del = [del; S.deltaRate(tested)]; end %#ok<AGROW>
        end
        if isempty(del), continue; end
        switch lower(P.test)
            case 'signrank', pRat = signrank(del,0);
            otherwise,       [~,pRat] = ttest(del,0);
        end
        if isfinite(pRat) && pRat < P.alpha, rats_sig = rats_sig + 1; end
    end
    ratsSigStr = sprintf('%d/%d', rats_sig, rats_tot);

    % day summaries (robust)
    if isempty(day_ps)
        fprintf('Combo %d/%d:   -> meanPct=%.1f%% | (no valid day-level stats)\n\n', ...
                jobIdx, Njobs, meanPct);
        dayDeltaMean = NaN; dayPMedian = NaN; daySigFrac = NaN; dayNtotal = 0;
    else
        dayDeltaMean = double(mean(day_deltas,'omitnan'));
        dayPMedian   = double(median(day_ps,'omitnan'));
        daySigFrac   = double(mean(day_ps < P.alpha,'omitnan'));
        dayNtotal    = double(nansum(day_ns));
        fprintf(['Combo %d/%d:   -> meanPct=%.1f%% | dayΔ_mean=%.4f Hz, ' ...
                 'day p̃ (median)=%.3g, days sig frac=%.2f, ΣnCells=%d\n\n'], ...
                 jobIdx, Njobs, meanPct, dayDeltaMean, dayPMedian, daySigFrac, round(dayNtotal));
    end

    % store row
    rows(jobIdx,:) = { ...
        char(P.Mode), ...
        ifempty(P.SpeedBinWidth, NaN), ...
        ifempty(P.SpeedBinCount, NaN), ...
        P.SpeedRange, ...
        seStr, ...
        P.MinDurPerBin, ...
        P.MinBins, ...
        P.binSize, ...
        P.alpha, ...
        P.test, ...
        sprintf('[%.3g %.3g]', P.win(1), P.win(2)), ...
        meanPct, ...
        ratsSigStr, ...
        struct('dayDeltaMean',dayDeltaMean,'dayPMedian',dayPMedian, ...
               'daySigFrac',daySigFrac,'totalCells',round(dayNtotal)) ...
    };
end

% keep only filled rows (skip skipped/empty)
filled = ~cellfun(@isempty, rows(:,1));
rows   = rows(filled,:);

T = cell2table(rows, 'VariableNames', { ...
    'Mode','SpeedBinWidth','SpeedBinCount','SpeedRange','SpeedEdges','MinDurPerBin','MinBins', ...
    'binSize','alpha','test','win','meanPctSig','ratsSignificant','dayLevelSummary'});
end

% ---- helpers ----
function out = ifempty(x, alt)
if isempty(x), out = alt; else, out = x; end
end

function s = shortStrSE(SE)
if isempty(SE), s='[]'; return; end
if numel(SE)<=5
    s=sprintf('[%s]',num2str(SE));
else
    s=sprintf('[%g..%g](n=%d)',SE(1),SE(end),numel(SE));
end
end
