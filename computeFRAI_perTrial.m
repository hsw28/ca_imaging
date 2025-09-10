function FRAI = computeFRAI_perTrial(ratName, varargin)
% computeFRAI_perTrial(ratName, 'Window', <'CSUS'|'CS'|'BASE'|'TRACE'|[a b]>, 'SkipFirst2022', true)
%
% Computes the firing-rate asymmetry index (FRAI) for each cell, per trial, on the
% last 3 CS days up to An for the rat struct living in the base workspace.
%
% FRAI definition (Skaggs-style):
%   If F1 and F2 are mean firing rates for the time during which the first and
%   second 50% of spikes were fired (within a trial window), then:
%       FRAI = (F1 - F2) / (F1 + F2)
% We handle small-spike cases robustly (see notes below).
%
% INPUTS
%   ratName  : char/string, e.g., 'rat0314'  (variable must exist in base)
%
% Name-Value options:
%   'Window'        : which within-trial window to use
%                     - 'CSUS'   (default): [CS_i, US_i)
%                     - 'CS'               : [CS_i, CS_i+0.5)
%                     - 'BASE'             : [CS_i-1.0, CS_i)
%                     - 'TRACE'            : [CS_i+0.5, US_i)
%                     - [a b] numeric      : [CS_i+a, CS_i+b)
%   'SkipFirst2022' : logical, default true. If true, drops first CS/US trial on 2022 dates.
%
% OUTPUT
%   FRAI : struct with fields:
%          - .byDay.(dayKey) : [nCells x nTrials] FRAI matrix for that day
%          - .meta.daysUsed  : cellstr of day keys
%          - .meta.window    : struct describing the window used
%          - .meta.ratName   : ratName
%          - .meta.cs_on/us_on.(dayKey) : vectors of CS/US used
%
% EXAMPLE
%   FRAI = computeFRAI_perTrial('rat0314','Window','TRACE');
%   imagesc(FRAI.byDay.CS_2023_05_04); colorbar; caxis([-1 1]);
%
% NOTES / Edge cases:
%   - If a trial window has < 2 spikes for a cell, FRAI is set to NaN (ill-defined halves).
%   - If the split time (when cumulative spikes reaches 50%) equals a boundary, we still
%     compute rates using exact durations on each half; a zero-duration half yields NaN.
%   - Ca_peaks is assumed to be [nCells x nEvents] with event times (seconds).
%   - CS/US times are pulled from rat.CS_times and rat.US_times as elsewhere in your code.
%
% Dependencies in this file:
%   local helpers: getAllDayKeys, keepCSdays, findFieldWithDate, valueVec,
%                  fraiFromSpikeTimes, buildTrialWindow

% ---------- parse inputs ----------
p = inputParser;
p.addRequired('ratName', @(s)ischar(s)||isstring(s));
p.addParameter('Window', 'CSUS', @(w)ischar(w)||isstring(w)|| ...
                                  (isnumeric(w)&&numel(w)==2&&isfinite(w(1))&&isfinite(w(2))));
p.addParameter('SkipFirst2022', true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(ratName, varargin{:});

ratName        = char(p.Results.ratName);
winOpt         = p.Results.Window;
skipFirst2022  = p.Results.SkipFirst2022;

% ---------- fetch rat struct ----------
if ~evalin('base', sprintf('exist(''%s'',''var'')', ratName))
    error('Variable %s not found in base workspace.', ratName);
end
rat = evalin('base', ratName);

% ---------- choose last 3 CS days up to An ----------
[csKeys, csDates] = keepCSdays(getAllDayKeys(rat));
anDate = datenum(strrep(rat.An,'_','-'));
useMask = csDates <= anDate;
csKeysUpToAn   = csKeys(useMask);
csDatesUpToAn  = csDates(useMask);

if numel(csKeysUpToAn) < 3
    error('%s: Not enough CS days up to An (%s).', ratName, rat.An);
end

[~, ord] = sort(csDatesUpToAn, 'ascend');
daysToUse = csKeysUpToAn(ord(end-2:end));   % last 3 CS days

% ---------- compute per-day FRAI matrices ----------
FRAI = struct();
FRAI.byDay = struct();
FRAI.meta.ratName = ratName;
FRAI.meta.daysUsed = daysToUse(:)';

for d = 1:numel(daysToUse)
    dayKey = daysToUse{d};               % 'CS_YYYY_MM_DD'
    usKey  = regexprep(dayKey,'^CS_','US_');

    % Pull CS/US times
    if ~isfield(rat.CS_times, dayKey) || ~isfield(rat.US_times, usKey)
        warning('[%s] missing CS/US times — skipping this day.', dayKey);
        continue;
    end
    cs_on = valueVec(rat.CS_times.(dayKey));
    us_on = valueVec(rat.US_times.(usKey));

    % Optionally drop the first trial on 2022 dates
    if skipFirst2022 && contains(dayKey, '2022') && numel(cs_on)>1 && numel(us_on)>1
        cs_on = cs_on(2:end);
        us_on = us_on(2:end);
    end

    % Basic validity and cleaning
    nT = min(numel(cs_on), numel(us_on));
    if nT < 1
        warning('[%s] no trials after sanity checks — skipping day.', dayKey);
        continue;
    end
    cs_on = cs_on(1:nT); us_on = us_on(1:nT);
    good  = (us_on - cs_on) > 0.01;
    cs_on = cs_on(good);  us_on = us_on(good);
    nT    = numel(cs_on);
    if nT < 1
        warning('[%s] no valid trials after filtering — skipping day.', dayKey);
        continue;
    end

    % Pull Ca_peaks for that date
    tok = regexp(dayKey,'(\d{4}_\d{2}_\d{2})','tokens','once');
    if isempty(tok)
        warning('[%s] could not parse date token — skipping day.', dayKey);
        continue;
    end
    dateTok = tok{1};
    pkFN = findFieldWithDate(rat.Ca_peaks, dateTok);
    if isempty(pkFN) || ~isfield(rat.Ca_peaks, pkFN)
        warning('[%s] no Ca_peaks found — skipping day.', dayKey);
        continue;
    end

    S_pk = rat.Ca_peaks.(pkFN);    % [nCells x nEvents] event times (s)
    if ~isnumeric(S_pk)
        warning('[%s] Ca_peaks is not numeric — skipping day.', dayKey);
        continue;
    end
    nCells = size(S_pk,1);

    % Compute FRAI per cell x trial
    FRAI_day = nan(nCells, nT);
    for c = 1:nCells
        spikes_c = S_pk(c,:);                         % 1 x nEvents
        spikes_c = spikes_c(isfinite(spikes_c) & spikes_c > 0);
        if isempty(spikes_c)
            continue;
        end
        for t = 1:nT
            [t0, t1] = buildTrialWindow(cs_on(t), us_on(t), winOpt);
            if ~(isfinite(t0)&&isfinite(t1)&&t1>t0)
                FRAI_day(c,t) = NaN; continue;
            end
            % spikes within [t0, t1)
            sp = spikes_c(spikes_c >= t0 & spikes_c < t1);
            FRAI_day(c,t) = fraiFromSpikeTimes(sp, t0, t1);
        end
    end

    FRAI.byDay.(dayKey) = FRAI_day;
    FRAI.meta.cs_on.(dayKey) = cs_on;
    FRAI.meta.us_on.(dayKey) = us_on;
end

% record window metadata
FRAI.meta.window = struct('option', winOpt);

end % main function

% ===================== helpers =====================

function [t0, t1] = buildTrialWindow(tCS, tUS, winOpt)
% Map a window option to absolute [t0,t1) for a given trial.
    if ischar(winOpt) || isstring(winOpt)
        switch upper(string(winOpt))
            case "CSUS"
                t0 = tCS;      t1 = tUS;
            case "CS"
                t0 = tCS;      t1 = tCS + 0.5;
            case "BASE"
                t0 = tCS - 1;  t1 = tCS;
            case "TRACE"
                t0 = tCS + 0.5; t1 = tUS;
            otherwise
                error('Unknown Window option: %s', string(winOpt));
        end
    elseif isnumeric(winOpt) && numel(winOpt)==2
        t0 = tCS + winOpt(1);
        t1 = tCS + winOpt(2);
    else
        error('Window must be ''CSUS'',''CS'',''BASE'',''TRACE'', or [a b].');
    end
end

function v = fraiFromSpikeTimes(spikeTimes, t0, t1)
% FRAI = (F1 - F2) / (F1 + F2), where F1/F2 are rates during the times that
% account for the first/second 50% of spikes. Robust to small N.
    v = NaN;
    spikeTimes = sort(spikeTimes(:));
    n = numel(spikeTimes);
    if n < 2
        % with <2 spikes, the half-split is ill-defined => leave NaN
        return;
    end
    % Half-count
    halfN = n/2;

    % Find the time where cumulative spike count reaches halfN (can be fractional between spikes)
    % We compute a "split time" tHalf by linear interpolation between spikes if n is odd.
    % Build a cumulative count step function:
    tAll = [t0; spikeTimes; t1];
    cAll = (0:numel(tAll)-1)';           % bogus, we will compute properly per spike
    % Better: directly compute tHalf
    if mod(n,2)==0
        % even: split exactly between spike #halfN and #halfN+1
        tHalf = (spikeTimes(halfN) + spikeTimes(halfN+1))/2;
    else
        % odd: the (halfN)-th is .5 fractional; place split at that spike's time
        tHalf = spikeTimes(ceil(halfN));
    end
    % clamp
    tHalf = max(min(tHalf, t1), t0);

    dur1 = max(tHalf - t0, 0);
    dur2 = max(t1 - tHalf, 0);
    if dur1<=0 || dur2<=0
        v = NaN; return;
    end

    % Count spikes occurring in each half (ties at tHalf go to second half by convention)
    n1 = sum(spikeTimes >= t0 & spikeTimes <  tHalf);
    n2 = sum(spikeTimes >= tHalf & spikeTimes < t1);

    F1 = n1 / dur1;
    F2 = n2 / dur2;
    denom = (F1 + F2);
    if denom<=0
        v = NaN;
    else
        v = (F1 - F2) / denom;
    end
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
    fn = '';
    if ~isstruct(S), return; end
    f = fieldnames(S);
    pat = strrep(dateTok,'_','[_-]?'); % tolerate '_' or '-'
    for i = 1:numel(f)
        if ~isempty(regexp(f{i}, pat, 'once'))
            fn = f{i}; return;
        end
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
