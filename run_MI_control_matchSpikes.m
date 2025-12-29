function run_MI_control_matchSpikes(csusNum)
% run_MI_control_matchSpikes  Wrapper to run MI control shuffle for all rats
%   run_MI_control_matchSpikes(15) uses CSUS15 fields
%   run_MI_control_matchSpikes(30) uses CSUS30 fields

if nargin < 1, csusNum = 15; end
%  ratNames = {'rat0816'};


ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};

Fs  = 7.5;
dim = 2.5;
nIter = 100;

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

    % Dynamically grab relevant fields
    spikes = filterFieldsByDay(rat.Ca_peaks, daysToUse);
    pos    = filterFieldsByDay(rat.pos, daysToUse);
    ts     = filterFieldsByDay(rat.Ca_ts, daysToUse);
    csus   = filterFieldsByDay(rat.(sprintf('csus%d', csusNum)), daysToUse);
    MI     = filterFieldsByDay(rat.(sprintf('MI_noCSUS%d', csusNum)), daysToUse);

    % Run control shuffle
    out = MI_control_matchSpikes(spikes, pos, 4, dim, ts, csus, MI, nIter);

    % Store in struct
    fieldname = sprintf('MI_noCSUS%d_control', csusNum);
    rat.(fieldname) = out;

    assignin('base', ratVar, rat);
end
end
