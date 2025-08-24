function neuralTraceHeatmap(csusNum, varargin)

%%%% WORKING

% neuralTraceHeatmap(csusNum, 'RatNames', {'rat0222'}, 'Pre',0.5, 'Post',0.5, ...
%                    'K',5, 'ZscorePerCell',true, 'MakePlots',true, 'Sigma',0.3)
%
% Panels per RAT:
%  (1) Trial x Time heatmap of population mean trace (hot colormap)
%  (2) PV similarity (TT, TN, NN) across trials, mean ± SEM across last 3 CS days
%  (3) Early (first K) vs Late (last K) mean population trace, mean ± SEM across days
%
% Inputs (same as before).

% -------- inputs --------
p = inputParser;
p.addParameter('RatNames', {'rat0222'}, @(c)iscell(c)||isstring(c));
p.addParameter('Pre', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Post', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('K', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('ZscorePerCell', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('MakePlots', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('Sigma', 0.30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});

ratNames    = cellstr(p.Results.RatNames);
preWin      = p.Results.Pre;
postWin     = p.Results.Post;
K           = p.Results.K;
doZ         = p.Results.ZscorePerCell;
makePlots   = p.Results.MakePlots;
sigmaKernel = p.Results.Sigma;

SR = 7.5;   % sampling rate used for synthetic traces from peaks

for ii = 1:numel(ratNames)
    ratVar = ratNames{ii};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', ratVar))
        warning('Variable %s not found in base workspace. Skipping.', ratVar);
        continue;
    end
    rat = evalin('base', ratVar);
    fprintf('\n=== %s : within-session analysis (neural) ===\n', ratVar);

    % ---- choose last 3 CS days up to An ----
    [csKeys, csDates] = keepCSdays(getAllDayKeys(rat));
    anDate = datenum(strrep(rat.An,'_','-'));
    useMask = csDates <= anDate;
    csKeysUpToAn   = csKeys(useMask);
    csDatesUpToAn  = csDates(useMask);
    if numel(csKeysUpToAn) < 3
        warning('%s: Not enough CS days up to An (%s). Skipping.', ratVar, rat.An);
        continue;
    end
    [~, ord] = sort(csDatesUpToAn, 'ascend');
    daysToUse = csKeysUpToAn(ord(end-2:end));   % 3 days

    csusField = sprintf('csus%d', csusNum);
    if ~isfield(rat, csusField)
        warning('%s missing field %s. Skipping...', ratVar, csusField);
        continue;
    end

    % ----- collectors for plotting -----
    H_days = {}; tRel_days = {}; usLag_days = [];               % Panel 1
    TT_days = {}; TN_days = {}; NN_days = {};                   % Panel 2
    earlyTrace_days = {}; lateTrace_days = {};                  % Panel 3 (per day time-courses)

    for d = 1:numel(daysToUse)
        dayKey = daysToUse{d};              % 'CS_YYYY_MM_DD'
        usKey  = regexprep(dayKey,'^CS_','US_');
        tok    = regexp(dayKey,'(\d{4}_\d{2}_\d{2})','tokens','once');
        if isempty(tok), fprintf('  [%s] bad key — skipping.\n', dayKey); continue; end
        dateTok = tok{1};

        if ~isfield(rat.CS_times, dayKey) || ~isfield(rat.US_times, usKey)
            fprintf('  [%s] missing CS/US times — skipping.\n', dayKey);
            continue;
        end
        cs_on = valueVec(rat.CS_times.(dayKey));
        us_on = valueVec(rat.US_times.(usKey));
        if (contains(dayKey,'2022')==1)
          cs_on= cs_on(2:end);
          us_on = us_on(2:end);
        end

        nT    = min(numel(cs_on), numel(us_on));
        if nT < 6, fprintf('  [%s] only %d trials — skipping.\n', dayKey, nT); continue; end
        cs_on = cs_on(1:nT); us_on = us_on(1:nT);
        good  = (us_on - cs_on) > 0.1;
        cs_on = cs_on(good);  us_on = us_on(good);
        nT    = numel(cs_on);
        if nT < 6, fprintf('  [%s] <6 valid trials — skipping.\n', dayKey); continue; end

        % --- get Ca_traces (or synthesize from Ca_peaks) ---
        [X, tCa] = getDayTraceOrPeaks(rat, dateTok, SR, sigmaKernel);  % X: [cells x nTime], tCa: [nTime x 1]
        if isempty(X), fprintf('  [%s] no Ca trace/peaks for %s — skipping.\n', dayKey, dateTok); continue; end
        if doZ, X = zscorePerCellDay(X); end

        % ===== Panel 1: trial x time heatmap of population mean =====
        [H_day, tRel, usLag] = makeTrialHeatmap(X, tCa, cs_on, us_on, preWin, postWin, SR);
        H_days{end+1}   = H_day;      %#ok<AGROW>  trials x time
        tRel_days{end+1}= tRel;       %#ok<AGROW>
        usLag_days(end+1)= usLag;     %#ok<AGROW>

        % keep per-day early and late time-courses for Panel 3
        K_use = min(K, size(H_day,1));
        earlyTrace_days{end+1} = mean(H_day(1:K_use ,:), 1, 'omitnan'); %#ok<AGROW>
        lateTrace_days{end+1}  = mean(H_day(end-K_use+1:end,:), 1, 'omitnan'); %#ok<AGROW>

        % ===== Panel 2: PV similarity across trials (TT/TN/NN) =====
        csWin   = [0, 0.50];
        trWin   = [0.50, NaN]; %#ok<NASGU> (trace unused for PV vectors here)
        usWin   = [0.00, 0.15]; %#ok<NASGU>
        baseWin = [-1.0, 0.0];

        rate = tracesToRates(X, tCa, cs_on, us_on, csWin, trWin, usWin, baseWin);
        [TT_raw, TN_raw, NN_raw] = pvSimRaw(rate);

        % only keep if there is at least one finite value
        if any(isfinite(TT_raw)), TT_days{end+1} = TT_raw; end %#ok<AGROW>
        if any(isfinite(TN_raw)), TN_days{end+1} = TN_raw; end %#ok<AGROW>
        if any(isfinite(NN_raw)), NN_days{end+1} = NN_raw; end %#ok<AGROW>
    end

    if ~makePlots, continue; end
    if isempty(H_days)
        fprintf('  No valid days to plot for %s.\n', ratVar);
        continue;
    end

    % ===== FIGURE per RAT =====
    figure('Name', sprintf('%s | neural trace heatmap + PV-sim', ratVar), ...
           'Color','w', 'Position', [80 80 1200 420]);
    tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    % -------- Panel 1: Heatmap averaged across days --------
    nexttile; hold on;
    [H_mean, tRel_mean] = meanHeatmapAcrossDays(H_days, tRel_days);
    imagesc(tRel_mean, 1:size(H_mean,1), H_mean);
    set(gca,'YDir','normal');
    colormap(gca, hot); colorbar; caxis([nanmin(H_mean(:)), nanmax(H_mean(:))]);
    xline(0, '--k', 'CS', 'LabelVerticalAlignment','bottom');
    if ~isempty(usLag_days)
        xline(nanmean(usLag_days), ':k', 'US (mean)', 'LabelVerticalAlignment','bottom');
    end
    xlabel('Time from CS (s)'); ylabel('Trial #');
    title(sprintf('%s — mean pop. trace (z=%d)', ratVar, doZ));

    % -------- Panel 2: PV similarity vs trial (mean ± SEM across days) --------

%    nexttile; hold on; grid on;
    [TTm, TTsem, xTT] = meanCurveAcrossDays(TT_days);
    [TNm, TNsem, xTN] = meanCurveAcrossDays(TN_days);
    [NNm, NNsem, xNN] = meanCurveAcrossDays(NN_days);

%    co = lines(3);
%    ebfill(xTT, TTm, TTsem, co(1,:)); plot(xTT, TTm, 'o-', 'LineWidth',2, 'Color',co(1,:), 'DisplayName','TT');
%    ebfill(xTN, TNm, TNsem, co(2,:)); plot(xTN, TNm, 's-', 'LineWidth',2, 'Color',co(2,:), 'DisplayName','TN');
%    ebfill(xNN, NNm, NNsem, co(3,:)); plot(xNN, NNm, '^-', 'LineWidth',2, 'Color',co(3,:), 'DisplayName','NN');
%    xlabel('Trial (or pair index)'); ylabel('Cosine similarity');
%    title('PV similarity across trials (mean ± SEM across days)');
%    legend('Location','best');


    % -------- Panel 3: Early vs Late mean population trace (±SEM) --------
    nexttile; hold on; grid on;
    % --- Mean early and late traces as before ---
    [Emean, Esem] = meanTracesAcrossDays(earlyTrace_days, tRel_days, tRel_mean);
    [Lmean, Lsem] = meanTracesAcrossDays(lateTrace_days,  tRel_days, tRel_mean);


    % --- Also plot individual early trial traces (first K trials) ---
    K = min(5, min(cellfun(@(M) size(M,1), earlyTrace_days)));  % number of early trials
    if K > 0
        cEarly = [0 0.45 0.74];                 % same hue as the mean
      %  faint  = [0.7 0.7 0.9];                 % faint lines for individuals


        for et = 1:K
            % gather the et-th early trial from each day (skip days w/ < et trials)
            D = nan(numel(earlyTrace_days), numel(tRel_mean));   % day x time
            have = false(numel(earlyTrace_days),1);

            for d = 1:numel(earlyTrace_days)
                M = earlyTrace_days{d};     % trials x time for this day
                if size(M,1) >= et
                    % interpolate this trial to the common timebase
                    tr_d = M(et,:);                    % 1 x T_d
                    t_d  = tRel_days{d}(:)';           % 1 x T_d
                    if numel(t_d) >= 2 && numel(tr_d) >= 2
                        D(d,:) = interp1(t_d, tr_d, tRel_mean, 'linear', 'extrap');
                        have(d) = true;
                        % plot the individual day trace faintly
                        %plot(tRel_mean, D(d,:), '-','LineWidth', 1);
                    end
                end
            end

            % average across days that had this trial
            if any(have)
                mi = mean(D(have,:), 1, 'omitnan');
                plot(tRel_mean, mi, '-', 'Color', cEarly, 'LineWidth', 1.5, ...
                     'HandleVisibility','off');  % thin colored line up to the mean
            end
        end
    end


    % --- Population means (with SEM shading) ---
    c = [0 0.45 0.74];
    ebfill(tRel_mean, Emean, Esem, lighten(c,0.6));
    plot(tRel_mean, Emean, '--', 'LineWidth',2, 'Color',c, ...
         'DisplayName', sprintf('Early mean (first %d)',K));

    c = [0.45 0.74 0];
    ebfill(tRel_mean, Lmean, Lsem, lighten(c,0.4));
    plot(tRel_mean, Lmean, '-',  'LineWidth',2, 'Color',c, ...
         'DisplayName', sprintf('Late mean (last %d)',K));

    xline(0, '--k');
    if ~isempty(usLag_days), xline(nanmean(usLag_days), ':k'); end
    xlabel('Time from CS (s)'); ylabel('Mean dF/F');
    title('Population trace: early vs late (±SEM), with individual early trials');
    legend('Location','best');

end
end

% ===================== helpers =====================

function [X, tCa] = getDayTraceOrPeaks(rat, dateTok, SR, sigma)
    X = []; tCa = [];
    trFN = findFieldWithDate(rat.Ca_traces, dateTok);
    tsFN = findFieldWithDate(rat.Ca_ts,     dateTok);

    if ~isempty(trFN) && ~isempty(tsFN)
        % ---- real traces available ----

        Xraw = rat.Ca_traces.(trFN);   % [cells x nTime]
        tRaw = rat.Ca_ts.(tsFN);    % [nTime x 1]

        tRaw = tRaw(:,2)./1000;

        % enforce strictly increasing time and keep matching columns in X
        [tRaw, uniqIdx] = unique(tRaw, 'stable');


        % downsample 15 Hz -> 7.5 Hz by taking every other sample
        downsamp = 1:2:length(tRaw);
        tRaw = tRaw(downsamp);          % [nTime_ds x 1]

        tCa = tRaw;
        X   = Xraw;       % [cells x nTime_ds]

        if length(X)>length(tCa)
          X = X(1:length(tCa));
        elseif length(X)<length(tCa)
          tCa = tCa(1:length(X));
        end
        return;
    end
end

function Xz = zscorePerCellDay(X)
    mu = mean(X,2,'omitnan'); sd = std(X,0,2,'omitnan'); sd(sd==0) = 1;
    Xz = (X - mu) ./ sd;
end

function [H, tRel, usLag] = makeTrialHeatmap(X, tCa, cs_on, us_on, pre, post, SR)
    pop = mean(X,1,'omitnan'); pop = pop(:);
    usLag = nanmean(us_on - cs_on);
    T_end = max(us_on - cs_on) + post;
    tRel  = (-pre : 1/SR : T_end)';              % relative axis
    H     = nan(numel(cs_on), numel(tRel));
    % sanity: if lengths still differ, map pop to tCa via interp once
    if numel(tCa) ~= numel(pop)
        pop = interp1(linspace(tCa(1), tCa(end), numel(pop))', pop, tCa, 'linear','extrap');
    end
    for k = 1:numel(cs_on)
        tSeg = cs_on(k) + tRel;
        H(k,:) = interp1(tCa, pop, tSeg, 'linear', 'extrap');
    end
end

function rate = tracesToRates(X, tCa, cs_on, us_on, csWin, ~, ~, baseWin)
    nC = size(X,1); nT = numel(cs_on);
    rate.cs   = nan(nC,nT);
    rate.base = nan(nC,nT);
    tCa = tCa(:);
    if size(X,2) ~= numel(tCa)
        error('tracesToRates: length mismatch: size(X,2)=%d, numel(tCa)=%d', size(X,2), numel(tCa));
    end
    for t = 1:nT
        tCS = cs_on(t);
        wCS = tCS + csWin;                % [CS+a, CS+b]
        wBL = tCS + baseWin;              % [CS-1, CS)
        iCS = (tCa >= wCS(1)) & (tCa < wCS(2));
        iBL = (tCa >= wBL(1)) & (tCa < wBL(2));
        if ~any(iCS) || ~any(iBL), continue; end
        rate.cs(:,t)   = mean(X(:,iCS), 2, 'omitnan');
        rate.base(:,t) = mean(X(:,iBL), 2, 'omitnan');
    end
end

function [TT_raw, TN_raw, NN_raw] = pvSimRaw(rate)
    Vt = rate.cs';   Vt = Vt ./ max(vecnorm(Vt,2,2), eps);
    Vn = rate.base'; Vn = Vn ./ max(vecnorm(Vn,2,2), eps);
    nT = size(Vt,1);
    TT_raw = nan(nT-1,1);
    for i=1:nT-1, TT_raw(i) = dot(Vt(i,:), Vt(i+1,:)); end
    TN_raw = nan(nT,1);
    for i=1:nT,   TN_raw(i) = dot(Vt(i,:), Vn(i,:));   end
    NN_raw = nan(nT-1,1);
    for i=1:nT-1, NN_raw(i) = dot(Vn(i,:), Vn(i+1,:)); end
end

function [Hmean, tRel_mean] = meanHeatmapAcrossDays(H_days, tRel_days)
    [~,ix] = max(cellfun(@numel, tRel_days));
    tRel_mean = tRel_days{ix};
    maxT = max(cellfun(@(H) size(H,1), H_days));
    Hstack = nan(maxT, numel(tRel_mean), numel(H_days));
    for d=1:numel(H_days)
        H = H_days{d}; tRel = tRel_days{d};
        Hrg = nan(size(H,1), numel(tRel_mean));
        for tr=1:size(H,1)
            Hrg(tr,:) = interp1(tRel, H(tr,:), tRel_mean, 'linear','extrap');
        end
        Hstack(1:size(H,1),:,d) = Hrg;
    end
    Hmean = mean(Hstack, 3, 'omitnan');
end

function [m, s, x] = meanCurveAcrossDays(CURVES)
    if isempty(CURVES), m=[]; s=[]; x=[]; return; end
    Lmax = max(cellfun(@numel, CURVES));
    M = nan(numel(CURVES), Lmax);
    for d=1:numel(CURVES)
        v = CURVES{d}(:)'; M(d,1:numel(v)) = v;
    end
    m = mean(M,1,'omitnan')'; s = sem(M); x = (1:Lmax)';
end

function [m, s] = meanTracesAcrossDays(TRACES, tRel_days, tRef)
% Align each per-day time-course to tRef, then mean ± SEM across days
    if isempty(TRACES), m = []; s = []; return; end
    nD = numel(TRACES);
    A = nan(nD, numel(tRef));
    for d=1:nD
        y = TRACES{d};
        tr = tRel_days{d};
        yi = interp1(tr, y, tRef, 'linear', 'extrap');
        A(d,:) = yi;
    end
    m = mean(A,1,'omitnan'); s = sem(A);
end

function s = sem(X)
    if isvector(X)
        x = X(:); n = sum(isfinite(x)); s = std(x,'omitnan')/max(1,sqrt(n));
    else
        n = sum(isfinite(X),1);
        s = std(X,0,1,'omitnan') ./ max(1,sqrt(n)); s = s(:);
    end
end

function ebfill(x, m, se, col)
    if isempty(x) || isempty(m) || isempty(se), return; end
    xv = [x(:); flipud(x(:))];
    yv = [m(:)-se(:); flipud(m(:)+se(:))];
    patch('XData',xv,'YData',yv,'FaceColor',col,'FaceAlpha',0.15,'EdgeColor','none');
end

function c2 = lighten(c, frac)
% mix with white by 'frac'
    c2 = c + (1-c)*min(max(frac,0),1);
end

function out = getAllDayKeys(rat)
    buckets = {'CS_times','US_times','Ca_peaks','Ca_ts','Ca_traces','pos','csus15','csus30','csus45','csus60','csus90'};
    out = {};
    for b = 1:numel(buckets)
        if isfield(rat, buckets{b}) && isstruct(rat.(buckets{b}))
            out = [out; fieldnames(rat.(buckets{b}))]; %#ok<AGROW>
        end
    end
    out = unique(out);
end

function [keysCS, datesCS] = keepCSdays(keysIn)
    keysCS = {}; datesCS = [];
    for i = 1:numel(keysIn)
        k = keysIn{i};
        tok = regexp(k, '^CS_(\d{4})_(\d{2})_(\d{2})$', 'tokens', 'once');
        if ~isempty(tok)
            dnum = datenum(sprintf('%s-%s-%s', tok{1}, tok{2}, tok{3}));
            keysCS{end+1,1} = k; %#ok<AGROW>
            datesCS(end+1,1) = dnum; %#ok<AGROW>
        end
    end
end

function fn = findFieldWithDate(S, dateTok)
    fn = ''; if ~isstruct(S), return; end
    f = fieldnames(S); pat = strrep(dateTok,'_','[_-]?');
    for i = 1:numel(f)
        if ~isempty(regexp(f{i}, pat, 'once')), fn = f{i}; return; end
    end
end

function out = valueVec(x)
    if isnumeric(x), out = x(:)'; return; end
    if isstruct(x)
        f = fieldnames(x);
        for k = 1:numel(f)
            if isnumeric(x.(f{k})), out = x.(f{k})(: )'; return; end
        end
    end
    error('CS/US times are not numeric.');
end

function X = makeSyntheticTraceFromPeaks(S_pk, t, sigma)
    nC = size(S_pk,1); X  = zeros(nC, numel(t));
    for c=1:nC
        ev = S_pk(c,:); ev = ev(isfinite(ev) & ev>0);
        if isempty(ev), continue; end
        for e=1:numel(ev)
            X(c,:) = X(c,:) + exp(-0.5*((t - ev(e))/sigma).^2);
        end
    end
end
