function run_MI_noCSUS_Shuffle(csusNum)
% run_MI_noCSUS_Shuffle  Wrapper to run MI shuffle excluding task periods (noCSUS)
% across the last 3 days for each rat.
%
%   run_MI_noCSUS_Shuffle(15) uses csus15 and MI_CSUS15 fields
%   run_MI_noCSUS_Shuffle(30) uses csus30 and MI_CSUS30 fields

if nargin < 1, csusNum = 15; end


ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};
Fs = 4;
dim = 2.5;
nShuff = 500;

for i = 1:numel(ratNames)
    ratVar = ratNames{i};
    rat = evalin('base', ratVar);

    % Get last 3 days
    dateList = autoDateList(rat);
    idx = find(strcmp(dateList, rat.An));
    if idx < 3
        warning('%s does not have enough days before An. Skipping...', ratVar);
        continue;
    end
    daysToUse = dateList(idx-2:idx);

    % Dynamically extract fields
    spike_struct = filterFieldsByDay(rat.Ca_peaks, daysToUse);
    pos_struct   = filterFieldsByDay(rat.pos, daysToUse);
    ts_struct    = filterFieldsByDay(rat.Ca_ts, daysToUse);
    csus_struct  = filterFieldsByDay(rat.(sprintf('csus%d', csusNum)), daysToUse);
    mi_struct    = filterFieldsByDay(rat.(sprintf('MI_noCSUS%d', csusNum)), daysToUse);

    % Run shuffle control
    MIshuff = mutualinfo_openfield_shuff_noCSUS( ...
        spike_struct, pos_struct, Fs, dim, ts_struct, csus_struct, mi_struct, nShuff);

    % Store into rat struct
    fieldname = sprintf('MI_noCSUS%d_shuff', csusNum);
    rat.(fieldname) = MIshuff;

    assignin('base', ratVar, rat);
end
end
