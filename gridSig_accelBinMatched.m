function T = gridSig_accelBinMatched(ratNames, grid, varargin)
% GRID search wrapper for ACCELERATION bins: average % sig across all rats
%
% Example:
%   grid.AccelMode     = {'signed','magnitude'};   % sweep (optional)
%   grid.AccelBinWidth = [0.05 0.1 0.2];
%   grid.MinDurPerBin  = [0.5 1 2];
%   grid.MinBins       = [3 5 7];
%   grid.AccelEdges    = {[]};  % or {-2:0.2:2, 0:0.2:5} etc.
%   T = gridSig_accelBinMatched([], grid);
%
% Resume example:
%   T = gridSig_accelBinMatched([], grid, 'start_idx', 538);

if nargin<1 || isempty(ratNames)
    ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end
if nargin<2 || isempty(grid)
    grid = struct();
end

% ----- optional args -----
start_idx = 1;
if ~isempty(varargin)
    for k = 1:2:numel(varargin)
        if strcmpi(varargin{k},'start_idx')
            start_idx = varargin{k+1};
        end
    end
end

% ---- defaults (accel-flavored) ----
if ~isfield(grid,'AccelMode'),     grid.AccelMode     = {'signed'}; end  % 'signed' | 'magnitude'
if ~iscell(grid.AccelMode),        grid.AccelMode     = cellstr(grid.AccelMode); end
if ~isfield(grid,'AccelBinWidth'), grid.AccelBinWidth = [0.05 0.1 0.2 0.3 0.5]; end
if ~isfield(grid,'AccelEdges'),    grid.AccelEdges    = {[]}; end
if ~iscell(grid.AccelEdges),       grid.AccelEdges    = num2cell(grid.AccelEdges); end
if ~isfield(grid,'MinDurPerBin'),  grid.MinDurPerBin  = [0.5 1 2 3 4 5]; end
if ~isfield(grid,'MinBins'),       grid.MinBins       = [3 4 5 7 9]; end
if ~isfield(grid,'alpha'),         grid.alpha         = 0.05; end
if ~isfield(grid,'test'),          grid.test          = {'ttest'}; end
if ~iscell(grid.test),             grid.test          = cellstr(grid.test); end
if ~isfield(grid,'win'),           grid.win           = {[0 2]}; end
if ~iscell(grid.win),              grid.win           = {grid.win}; end
if ~isfield(grid,'binSize'),       grid.binSize       = [1/7.5 2/7.5 3/7.5 .5 4/7.5 5/7.5 6/7.5 1 11/7.5 2]; end

% make row arrays
v_MODE = grid.AccelMode(:).';
v_ABW  = grid.AccelBinWidth(:).';
v_AE   = grid.AccelEdges(:).';
v_MDPB = grid.MinDurPerBin(:).';
v_MB   = grid.MinBins(:).';
v_A    = grid.alpha(:).';
v_TEST = grid.test(:).';
v_WIN  = grid.win(:).';
v_BS   = grid.binSize(:).';

% cartesian product (note: MODE added)
[iMODE,iABW,iAE,iMD,iMB,iA,iTEST,iWIN,iBS] = ndgrid( ...
    1:numel(v_MODE), 1:numel(v_ABW), 1:numel(v_AE), 1:numel(v_MDPB), 1:numel(v_MB), ...
    1:numel(v_A), 1:numel(v_TEST), 1:numel(v_WIN), 1:numel(v_BS));
N = numel(iABW);

% sanitize start_idx
if ~isscalar(start_idx) || ~isfinite(start_idx), start_idx = 1; end
start_idx = max(1, min(N, round(start_idx)));
if start_idx>1
    fprintf('Resuming at combo %d/%d\n', start_idx, N);
end

rows = cell(N, 11);

for k=start_idx:N
    P.AccelMode     = v_MODE{iMODE(k)};
    P.AccelBinWidth = v_ABW(iABW(k));
    P.AccelEdges    = v_AE{iAE(k)};
    P.MinDurPerBin  = v_MDPB(iMD(k));
    P.MinBins       = v_MB(iMB(k));
    P.alpha         = v_A(iA(k));
    P.test          = v_TEST{iTEST(k)};
    P.win           = v_WIN{iWIN(k)};
    P.binSize       = v_BS(iBS(k));

    % Run without plots, quiet logging
    Rk = run_accelBinMatched(ratNames, ...
        'AccelMode',     P.AccelMode, ...
        'AccelBinWidth', P.AccelBinWidth, ...
        'AccelEdges',    P.AccelEdges, ...
        'MinDurPerBin',  P.MinDurPerBin, ...
        'MinBins',       P.MinBins, ...
        'alpha',         P.alpha, ...
        'test',          P.test, ...
        'win',           P.win, ...
        'binSize',       P.binSize);

    % --- average % sig over all rats
    pctAll  = vertcat(Rk.pctSig);
    meanPct = mean(pctAll,'omitnan');

    % --- “rats significant” (per-rat pooled across days)
    rats_sig = 0; rats_tot = numel(Rk);
    for rr = 1:numel(Rk)
        del = [];
        for d = 1:numel(Rk(rr).perDay)
            S = Rk(rr).perDay{d};
            if isempty(S), continue; end
            tested = isfinite(S.pVal);
            if any(tested)
                del = [del; S.deltaRate(tested)]; %#ok<AGROW>
            end
        end
        if isempty(del), continue; end
        switch lower(P.test)
            case 'signrank'
                pRat = signrank(del,0);
            otherwise
                [~,pRat] = ttest(del,0);
        end
        if isfinite(pRat) && pRat < P.alpha
            rats_sig = rats_sig + 1;
        end
    end

    rows(k, :) = { ...
        P.AccelMode, ...
        P.AccelBinWidth, ...
        shortStrSE(P.AccelEdges), ...
        P.MinDurPerBin, ...
        P.MinBins, ...
        P.binSize, ...
        P.alpha, ...
        P.test, ...
        P.win, ...
        meanPct, ...
        sprintf('%d/%d', rats_sig, rats_tot) ...
    };

    fprintf(['Combo %d/%d: MODE=%s, ABW=%.3g, AE=%s, MinDur=%.3g, MinBins=%d, ' ...
             'binSize=%.4g -> meanPct=%.1f%%\n\n'], ...
            k, N, P.AccelMode, P.AccelBinWidth, shortStrSE(P.AccelEdges), ...
            P.MinDurPerBin, P.MinBins, P.binSize, meanPct);
end

% keep only filled rows
filled = ~cellfun(@isempty, rows(:,1));
rows   = rows(filled,:);

T = cell2table(rows, 'VariableNames', { ...
    'AccelMode','AccelBinWidth','AccelEdges','MinDurPerBin','MinBins', ...
    'binSize','alpha','test','win','meanPctSig','ratsSignificant'});

end

function s = shortStrSE(SE)
if isempty(SE), s='[]'; return; end
if numel(SE)<=5
    s=sprintf('[%s]',num2str(SE));
else
    s=sprintf('[%g..%g](n=%d)',SE(1),SE(end),numel(SE));
end
end
