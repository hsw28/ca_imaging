function Rall = postShockPV_evenOdd(ratNames, varargin)
% postShockPV_evenOdd
%   Even/odd split PV-correlation to estimate how long post-US ensemble
%   patterns persist and potentially interfere with place coding.
%
% Usage (single rat):
%   R = postShockPV_evenOdd('rat0314', 'UseSpeedMask',true, ...);
%
% Usage (multiple rats):
%   Rall = postShockPV_evenOdd({'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%           'UseSpeedMask',true, 'DoPlot',true);
%
% When multiple rats are passed:
%   - Creates one figure with 2 columns per rat (PV corr heatmap; diag trace)
%   - Creates a second figure with pooled-across-rats results
%   - Returns Rall with per-rat results and pooled-across-rats fields.
%
% New option:
%   'Preprocess' : 'none' | 'demean' | 'zscore' | 'meanrate' (default 'none')
%                  Per-cell normalization applied to across all bins×trials
%                  X (cells × bins × trials) before forming PVs.

% ---------- args ----------
p = inputParser;
addParameter(p,'Days','An-2:An');
addParameter(p,'PostWin',[0 6],@(v)isnumeric(v)&&numel(v)==2&&v(2)>v(1));
%addParameter(p,'BinSec',(1/7.5)*3,@isscalar);
addParameter(p,'BinSec',(1/7.5),@isscalar);

addParameter(p,'MinEvents',3,@isscalar);
addParameter(p,'QFDR',0.05,@isscalar);
addParameter(p,'NShuff',1000,@(x)isscalar(x)&&x>=100);
addParameter(p,'DoPlot',true,@islogical);
addParameter(p,'Debug',false,@islogical);
addParameter(p,'UseSpeedMask',false,@islogical);
addParameter(p,'SpeedThresh',4,@isscalar);
addParameter(p,'SpeedMaskFrac',0.5,@(x)isscalar(x)&&x>=0&&x<=1);
addParameter(p,'Preprocess','demean',@(s) ischar(s)||isstring(s));

parse(p,varargin{:});
opts = p.Results;

% Normalize ratNames into a cellstr
if ischar(ratNames) || isstring(ratNames)
    ratNames = {char(ratNames)};
elseif ~iscell(ratNames)
    error('ratNames must be a char, string, or cell array of these.');
end

% Run per-rat and collect
perRat = cell(1,numel(ratNames));
pooled_Z_stack = [];   % collect Fisher-z mats across all rats (weighted sum later)
pooled_wts      = [];  % weights corresponding to each Z mat slice
tBins_master    = [];  % record tBins (must be the same across rats to pool C)
had_any_rat     = false;

% For per-rat subplot grid
if opts.DoPlot && numel(ratNames) > 1
    nR = numel(ratNames);
    figure('Color','w','Position',[100 100 1200 max(420, 240*nR)]);
end

for r = 1:numel(ratNames)
    rn = ratNames{r};
    Ri = run_one_rat(rn, opts);  % per-rat core analysis

    perRat{r} = Ri;
    if ~isempty(Ri)
        had_any_rat = true;
        % plot per-rat rows if multiple
        if opts.DoPlot && numel(ratNames) > 1
            % Row index
            row = r;
            % Heatmap
            subplot(numel(ratNames), 2, 2*(row-1)+1);
            imagesc(Ri.pooled.tBins, Ri.pooled.tBins, Ri.pooled.C, [-0.4 0.6]);
            axis image; colormap(gca, parula); colorbar; hold on;
            plot([Ri.pooled.tBins(1) Ri.pooled.tBins(end)], [Ri.pooled.tBins(1) Ri.pooled.tBins(end)], 'k--');
            title(sprintf('%s: PV corr (even vs odd)', rn), 'Interpreter','none');
            xlabel('Odd (s post-US)'); ylabel('Even (s post-US)');

            % Diagonal
            subplot(numel(ratNames), 2, 2*(row-1)+2);
            plot(Ri.pooled.tBins, Ri.pooled.diag_r, 'LineWidth', 1.5); hold on;
            yyaxis right; stairs(Ri.pooled.tBins, double(Ri.pooled.mask_fdr)*100, 'LineWidth', 1.2);
            ylabel('Sig (FDR)'); ylim([0 110]);
            yyaxis left; ylabel('r (diag)');
            xlabel('Time since US (s)');
            title(sprintf('%s: diag r (q=%.2g) | dur=%.2fs', rn, opts.QFDR, Ri.pooled.duration_sig_s), 'Interpreter','none');
            grid on;
        end

        % Accumulate pooled across-rats Fisher z
        if isempty(tBins_master)
            tBins_master = Ri.pooled.tBins;
        else
            % sanity: ensure same tBins
            if any(abs(Ri.pooled.tBins - tBins_master) > 1e-9)
                warning('tBins differ for %s; skipping from pooled-across-rats heatmap.', rn);
                Ri.skipFromGlobal = true; %#ok<*STRNU>
            end
        end

        if ~isfield(Ri,'skipFromGlobal') || ~Ri.skipFromGlobal
            pooled_Z_stack = cat(3, pooled_Z_stack, atanh(max(min(Ri.pooled.C,0.999999),-0.999999)));
            % weight by sqrt(total trials across that rat’s included days)
            w = 0;
            for dd = 1:numel(Ri.perDay)
                w = w + sqrt(Ri.perDay(dd).Ntrials);
            end
            pooled_wts = [pooled_wts; max(w, eps)];
        end
    end
end

if ~had_any_rat
    error('No valid rats analyzed. Check MinEvents, US detection, Days, or speed mask settings.');
end

% Build pooled-across-rats heatmap and diag/p-values if we have any Z
pooledAcross = struct();
if ~isempty(pooled_Z_stack)
    W = reshape(pooled_wts, 1,1,[]);
    Z = sum(pooled_Z_stack .* W, 3, 'omitnan') ./ sum(W,3,'omitnan');
    C_pooled = tanh(Z);
    diag_pooled_r = diag(C_pooled)';

    % Combine per-bin p across rats by Stouffer from per-rat p_diag
    % Gather per-rat p_diag on the (common) tBins
    pMat = [];
    ws   = [];
    for r = 1:numel(ratNames)
        Ri = perRat{r};
        if isempty(Ri) || (isfield(Ri,'skipFromGlobal') && Ri.skipFromGlobal), continue; end
        p = Ri.pooled.p_diag;
        if ~isempty(p)
            pMat = [pMat; p(:)']; %#ok<AGROW>
            % weight ~ sqrt(total trials across that rat’s days)
            w = 0;
            for dd = 1:numel(Ri.perDay), w = w + sqrt(Ri.perDay(dd).Ntrials); end
            ws = [ws; max(w, eps)]; %#ok<AGROW>
        end
    end
    if ~isempty(pMat)
        pMat(isnan(pMat)) = 1;
        Zs = -norminv(pMat/2);
        Wrep = repmat(ws,1,size(pMat,2));
        Zcomb = sum(Wrep .* Zs, 1) ./ sqrt(sum(Wrep.^2,1));
        p_pooled = 2*normcdf(-abs(Zcomb));
    else
        p_pooled = nan(1, numel(tBins_master));
    end

    % FDR and segments
    [~, mask_pool] = fdr_bh_local(replaceNaNwithOne(p_pooled), opts.QFDR);
    [seg_pool_s, dur_pool_s] = contiguousSegments(mask_pool, tBins_master, opts.BinSec);

    pooledAcross.C              = C_pooled;
    pooledAcross.diag_r         = diag_pooled_r;
    pooledAcross.p_diag         = p_pooled;
    pooledAcross.mask_fdr       = mask_pool;
    pooledAcross.segments_s     = seg_pool_s;
    pooledAcross.duration_sig_s = dur_pool_s;
    pooledAcross.tBins          = tBins_master;

    % Plot pooled-across-rats
    if opts.DoPlot
        figure('Color','w','Position',[100 100 1050 420]);
        subplot(1,2,1);
        imagesc(tBins_master, tBins_master, C_pooled, [-0.4 0.6]); axis image;
        colormap(gca, parula); colorbar; hold on;
        plot([tBins_master(1) tBins_master(end)], [tBins_master(1) tBins_master(end)], 'k--');
        title('POOLED across rats: PV corr (even vs odd)');
        xlabel('Odd (s post-US)'); ylabel('Even (s post-US)');

        subplot(1,2,2); hold on;
        plot(tBins_master, diag_pooled_r, 'LineWidth', 1.5);
        ylabel('r (diag)');
        xlabel('Time since US (s)');
        grid on;
    end
end

% Assemble return
Rall.perRat = perRat;
if ~isempty(fieldnames(pooledAcross))
    Rall.pooledAcrossRats = pooledAcross;
end
Rall.params = opts;

end % ==================== main ====================


% ==================== core per-rat ====================
function R = run_one_rat(ratName, opts)
% Returns a per-rat R struct in the same shape you already had for a single rat.

% ---------- resolve rat & days ----------
rat = evalin('base', ratName);
allDays = autoDateList(rat);

if ischar(opts.Days) || isstring(opts.Days)
    if strcmpi(opts.Days,'An-2:An')
        idx = find(strcmp(allDays, rat.An));
        assert(~isempty(idx) && idx>=3, 'Not enough days before An in %s.', ratName);
        daysToUse = allDays(idx-2:idx);
    else
        daysToUse = {char(opts.Days)};
    end
else
    daysToUse = opts.Days;   % assume cellstr
end

% Pull Ca_peaks for selected days
spike_struct = filterFieldsByDay(rat.Ca_peaks, daysToUse);
dayNames = fieldnames(spike_struct);

% Precompute bin edges & centers (relative to US)
tEdges = opts.PostWin(1):opts.BinSec:opts.PostWin(2);
Nbins  = numel(tEdges)-1;
tBins  = (tEdges(1:end-1)+tEdges(2:end))/2;

perDay = struct([]);
zMats  = [];    % Fisher-z per-day
diagPs = [];
wts    = [];

for d = 1:numel(dayNames)
    day = dayNames{d};
    useSpeedMask_day = opts.UseSpeedMask;

    CaMat = spike_struct.(day);
    % Convert Ca_peaks double (rows=cells, cols=event times) → cell of time vectors
    if iscell(CaMat)
        evByCell = CaMat;
    elseif isnumeric(CaMat)
        Ncells = size(CaMat,1);
        evByCell = cell(1,Ncells);
        for c = 1:Ncells
            v = CaMat(c,:);
            v = v(~isnan(v) & v>0);
            evByCell{c} = v(:)';   % seconds
        end
    else
        error('Unsupported Ca_peaks datatype for %s (got %s).', day, class(CaMat));
    end

    %%% NEW: apply ratemask==1 (BEFORE MinEvents) %%%
    Mfld = sprintf('ratemask_%s', day);
    if isfield(rat,'ratemask') && isstruct(rat.ratemask) && isfield(rat.ratemask, Mfld)
        rm = rat.ratemask.(Mfld) == 1;
        rm = rm(:);
        if numel(rm) == numel(evByCell)
            evByCell = evByCell(rm);
        end
    end
    %%% END NEW %%%

    % MinEvents per cell
    totEv = cellfun(@numel, evByCell);
    keepC = totEv >= opts.MinEvents;
    evByCell = evByCell(keepC);
    NcellsUsed = sum(keepC);
    if NcellsUsed < 3
        if opts.Debug, whySkipped(day, sprintf('MinEvents left %d cells', NcellsUsed)); end
        continue;
    end

    % US detection
    [usTimes, whyUS] = getUSTimes_guess(rat, day);
    if isempty(usTimes)
        if opts.Debug, whySkipped(day, strjoin(whyUS,' | ')); end
        continue;
    end
    Ntr = numel(usTimes);
    if Ntr < 6
        if opts.Debug, whySkipped(day, sprintf('Only %d trials (<6)', Ntr)); end
        continue;
    end

    % Optional per-trial, per-bin speed mask (from pos via ca_velocity)
    binMask = true(Nbins, Ntr);
    if useSpeedMask_day
        [tVel, vVel] = getSpeedFromPlace_via_ca_velocity(rat, day);
        if isempty(tVel)
            if opts.Debug, dbg(day,'UseSpeedMask requested but no speed/pos found; proceeding UNMASKED'); end
            useSpeedMask_day = false;
        else
            for tr = 1:Ntr
                edges = tEdges + usTimes(tr);
                mids  = (edges(1:end-1) + edges(2:end))/2;
                v_mid = interp1(tVel, vVel, mids, 'linear', 'extrap');
                binMask(:,tr) = v_mid(:) >= opts.SpeedThresh;
            end
            if opts.Debug
                dbg(day, sprintf('speed mask built: median frac bins passing per trial = %.2f', ...
                    median(mean(binMask,1))));
            end
        end
    end

    % Build counts: Ncells × Nbins × Ntr
    X = zeros(NcellsUsed, Nbins, Ntr, 'single');
    for tr = 1:Ntr
        relEdges = tEdges + usTimes(tr);
        for c = 1:NcellsUsed
            ev = evByCell{c};
            if ~isempty(ev)
                X(c,:,tr) = histcounts(ev, relEdges);
            end
        end
    end

    % Apply speed mask (trial-bin NaNs)
    if useSpeedMask_day
        for tr = 1:Ntr
            bad = ~binMask(:,tr);
            if any(bad), X(:,bad,tr) = NaN; end
        end
    end

    % ---------- per-cell preprocessing ----------
    switch lower(string(opts.Preprocess))
        case "none"
            % do nothing
        case "demean"
            mu = mean(X, [2 3], 'omitnan');          % cells×1×1
            X  = X - mu;                              % broadcast subtract
        case "zscore"
            mu = mean(X, [2 3], 'omitnan');
            sd = std(X, 0, [2 3], 'omitnan');
            sd(sd==0 | isnan(sd)) = eps;
            X  = (X - mu) ./ sd;
        case "meanrate"
          mu = mean(X, [2 3], 'omitnan');          % cells×1×1
          X  = X - mu;

        otherwise
            error('Preprocess must be ''none'', ''demean'', or ''zscore''.');
    end
    % -----------------------------------------------

    % Even/odd splits
    trIdx  = 1:Ntr;
    trEven = trIdx(mod(trIdx,2)==0);  trOdd = trIdx(mod(trIdx,2)==1);
    if isempty(trEven) || isempty(trOdd), continue; end

    % Mean PV per bin for each split
    PV_even = mean(X(:,:,trEven), 3, 'omitnan') / opts.BinSec;
    PV_odd  = mean(X(:,:,trOdd ), 3, 'omitnan') / opts.BinSec;

    % Corr matrix (bins × bins)
    C = nan(Nbins, Nbins, 'single');
    for i = 1:Nbins
        vi = PV_even(:,i); if all(isnan(vi)) || all(vi==0), continue; end
        for j = 1:Nbins
            vj = PV_odd(:,j);  if all(isnan(vj)) || all(vj==0), continue; end
            C(i,j) = corr(vi, vj, 'rows','complete','type','Pearson');
        end
    end
    diag_r = diag(C)';

    % Permutation null (shuffle even/odd labels)
    labels = zeros(1,Ntr); labels(trOdd)=1; labels(trEven)=2;
    diag_null = nan(opts.NShuff, Nbins, 'single');
    for s = 1:opts.NShuff
        sh = labels(randperm(Ntr));
        shOdd  = find(sh==1);
        shEven = find(sh==2);
        PV_e = mean(X(:,:,shEven), 3, 'omitnan') / opts.BinSec;
        PV_o = mean(X(:,:,shOdd ), 3, 'omitnan') / opts.BinSec;
        rs = nan(1,Nbins,'single');
        for k = 1:Nbins
            vi = PV_e(:,k); vj = PV_o(:,k);
            if all(isnan(vi)) || all(isnan(vj)) || all(vi==0) || all(vj==0), continue; end
            rs(k) = corr(vi, vj, 'rows','complete','type','Pearson');
        end
        diag_null(s,:) = rs;
    end

    % Per-bin two-sided p-values
    p_diag = nan(1,Nbins);
    for k = 1:Nbins
        r0 = diag_r(k);
        nullk = diag_null(:,k); nullk = nullk(~isnan(nullk));
        if isnan(r0) || isempty(nullk)
            p_diag(k) = NaN;
        else
            p_diag(k) = max(1/numel(nullk), mean(abs(nullk) >= abs(r0)));
        end
    end

    % FDR across bins (per day)
    pForFDR = replaceNaNwithOne(p_diag);
    [~, mask_fdr] = fdr_bh_local(pForFDR, opts.QFDR);

    % Segments
    [segments_s, duration_sig_s] = contiguousSegments(mask_fdr, tBins, opts.BinSec);

    % Store
    perDay(end+1).day           = day; %#ok<AGROW>
    perDay(end).C               = C;
    perDay(end).diag_r          = diag_r;
    perDay(end).p_diag          = p_diag;
    perDay(end).mask_fdr        = mask_fdr;
    perDay(end).segments_s      = segments_s;
    perDay(end).duration_sig_s  = duration_sig_s;
    perDay(end).Nbins           = Nbins;
    perDay(end).BinSec          = opts.BinSec;
    perDay(end).PostWin         = opts.PostWin;
    perDay(end).Ntrials         = Ntr;
    perDay(end).NcellsUsed      = NcellsUsed;

    % For pooled (within rat)
    zMats = cat(3, zMats, atanh(max(min(C,0.999999),-0.999999)));
    diagPs = [diagPs; p_diag]; %#ok<AGROW>
    wts    = [wts; sqrt(Ntr)]; %#ok<AGROW>
end

if isempty(perDay)
    if opts.Debug, fprintf('[%s] No valid days.\n', ratName); end
    R = []; return;
end

% Pooled within rat
W = reshape(wts,1,1,[]);
Z = sum(zMats .* W, 3, 'omitnan') ./ sum(W,3,'omitnan');
C_pooled = tanh(Z);
diag_pooled_r = diag(C_pooled)';

% Stouffer combine per-bin p across days (within rat)
pMat = diagPs;  pMat(isnan(pMat)) = 1;
Zs   = -norminv(pMat/2);
Wrep = repmat(wts,1,size(pMat,2));
Zcomb = sum(Wrep .* Zs, 1) ./ sqrt(sum(Wrep.^2,1));
p_pooled = 2*normcdf(-abs(Zcomb));

[~, mask_pool] = fdr_bh_local(replaceNaNwithOne(p_pooled), opts.QFDR);
[seg_pool_s, dur_pool_s] = contiguousSegments(mask_pool, tBins, opts.BinSec);

R.perDay = perDay;
R.pooled.C              = C_pooled;
R.pooled.diag_r         = diag_pooled_r;
R.pooled.p_diag         = p_pooled;
R.pooled.mask_fdr       = mask_pool;
R.pooled.segments_s     = seg_pool_s;
R.pooled.duration_sig_s = dur_pool_s;
R.pooled.tBins          = tBins;
R.params = struct('PostWin',opts.PostWin,'BinSec',opts.BinSec,'MinEvents',opts.MinEvents, ...
                  'QFDR',opts.QFDR,'NShuff',opts.NShuff,'Days',{daysToUse}, ...
                  'UseSpeedMask',opts.UseSpeedMask,'SpeedThresh',opts.SpeedThresh, ...
                  'SpeedMaskFrac',opts.SpeedMaskFrac,'Preprocess',opts.Preprocess);
end


% ===================== helpers =====================
function [us, why] = getUSTimes_guess(rat, day)
why = {}; us = [];
cands = candidateDayKeys(day);

% top-level rat.US_times.(dayKey)
if isfield(rat,'US_times')
    if isstruct(rat.US_times)
        for i = 1:numel(cands)
            dk = cands{i};
            if isfield(rat.US_times, dk)
                cand = rat.US_times.(dk);
                if isnumeric(cand) && ~isempty(cand)
                    us = cand(:).'; why{end+1}=sprintf('Found US in rat.US_times.%s',dk); return
                end
            end
        end
    elseif isnumeric(rat.US_times) && ~isempty(rat.US_times)
        us = rat.US_times(:).'; why{end+1}='Found numeric US in rat.US_times'; return
    end
end

parents = {'US','CSUS','Behavior','Events','Trials','Stim'};
keys    = {'US','us','US_onsets','USonset','us_onsets','US_times','USTimes'};
for p = 1:numel(parents)
    P = parents{p};
    if isfield(rat,P) && isstruct(rat.(P))
        for i = 1:numel(cands)
            dk = cands{i};
            if isfield(rat.(P), dk)
                blk = rat.(P).(dk);
                if isstruct(blk)
                    for k = 1:numel(keys)
                        K = keys{k};
                        if isfield(blk,K) && ~isempty(blk.(K)) && isnumeric(blk.(K))
                            us = blk.(K)(:).'; why{end+1}=sprintf('Found US in rat.%s.%s.%s',P,dk,K); return
                        end
                    end
                elseif isnumeric(blk) && ~isempty(blk)
                    us = blk(:).'; why{end+1}=sprintf('Found numeric US in rat.%s.%s',P,dk); return
                end
            end
        end
    end
end
why{end+1}='No explicit US field found.';

% fallback from CS (+0.75 s)
csKeys = {'CS','cs','CS_onsets','CSonset','cs_onsets','CS_times','CSTimes'};
if isfield(rat,'CS_times') && isstruct(rat.CS_times)
    for i=1:numel(cands)
        dk = cands{i};
        if isfield(rat.CS_times, dk)
            cs = rat.CS_times.(dk);
            if isnumeric(cs) && ~isempty(cs)
                us = cs(:).'+0.75; why{end+1}=sprintf('Derived US from rat.CS_times.%s (+0.75s)',dk); return
            end
        end
    end
end
for p = 1:numel(parents)
    P = parents{p};
    if isfield(rat,P) && isstruct(rat.(P))
        for i=1:numel(cands)
            dk = cands{i};
            if isfield(rat.(P), dk)
                blk = rat.(P).(dk);
                if isstruct(blk)
                    for k = 1:numel(csKeys)
                        K = csKeys{k};
                        if isfield(blk,K) && ~isempty(blk.(K)) && isnumeric(blk.(K))
                            cs = blk.(K)(:).'; us = cs+0.75;
                            why{end+1}=sprintf('Derived US from CS in rat.%s.%s.%s (+0.75s)',P,dk,K); return
                        end
                    end
                end
            end
        end
    end
end
why{end+1}='Could not derive US from CS either.';
end

function cands = candidateDayKeys(day)
S = string(day);
cand = unique(S);
cand = [cand, erase(S, ["CA_peaks_","Ca_peaks_","CA_Peaks_","peaks_","CS_","US_"])];
tmp = unique(cand);
cand = unique([cand, replace(tmp,"_","-"), replace(tmp,"-","_")]);
cand = unique([cand, "US_"+cand, "CS_"+cand]);
cands = cellstr(cand);
end

function [thr, mask] = fdr_bh_local(p, q)
p = p(:)'; p(isnan(p)) = 1; [ps, ~] = sort(p);
m = numel(ps); k = find(ps <= (1:m)/m*q, 1, 'last');
mask = false(size(p)); if ~isempty(k), crit = ps(k); mask = p <= crit; end
thr = [];
end

function [segments_s, totalDur_s] = contiguousSegments(mask, tBins, binSec)
segments_s = []; totalDur_s = 0;
if isempty(mask) || ~any(mask), return; end
idx = find(mask); gaps = find(diff(idx)>1); cuts = [0 gaps numel(idx)];
for i = 1:numel(cuts)-1
    chunk = idx(cuts(i)+1 : cuts(i+1));
    t0 = tBins(chunk(1)); t1 = tBins(chunk(end)) + binSec;
    segments_s = [segments_s; [t0 t1]]; %#ok<AGROW>
    totalDur_s = totalDur_s + (t1 - t0);
end
end

function whySkipped(day, reason)
fprintf('[skip] %s — %s\n', day, reason);
end

function [tVel, vVel] = getSpeedFromPlace_via_ca_velocity(rat, day)
% Your rule:
%   pos lives at rat.pos.pos_YYYY_MM_DD (numeric matrix)
%   ca_velocity(pos_data) -> 2xN, row1 = vel (cm/s), row2 = time (s)
tVel = []; vVel = [];
if ~isfield(rat,'pos') || ~isstruct(rat.pos), return; end
S = string(day);
base1 = erase(S, ["CA_peaks_","Ca_peaks_","CA_Peaks_","peaks_","CS_","US_"]);
base2 = replace(base1, "-", "_");
base3 = replace(base1, "_", "-");
m = regexp(char(base1), '(\d{4})[-_](\d{2})[-_](\d{2})', 'tokens', 'once');
if ~isempty(m)
    pure_us = sprintf('%s_%s_%s', m{1}, m{2}, m{3});
else
    pure_us = char(base2);
end
cands = unique({['pos_' char(base2)], ['pos_' pure_us], ['pos_' char(base1)], ['pos_' char(base3)], char(base2), char(base1)});
for i = 1:numel(cands)
    dk = cands{i};
    if isfield(rat.pos, dk)
        pos_data = rat.pos.(dk);
        try
            out = ca_velocity(pos_data);     % returns 2xN: [vel; time]
            if isnumeric(out) && size(out,1) >= 2 && ~isempty(out)
                vVel = out(1,:).';
                tVel = out(2,:).';
                return
            end
        catch
        end
    end
end
end

function dbg(day, msg)
fprintf('[info] %s — %s\n', day, msg);
end

function p = replaceNaNwithOne(p)
p(isnan(p)) = 1;
end
