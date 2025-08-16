function run_MI_wCSUS_Shuffle
% Wrapper to run MI shuffle for only An and 2 days before for all 5 rats
%returns top 5%, mean, shuffled percent (want 0.95 or higher), and hertz during running periods


ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};
Fs = 4;
dim = 2.5;
nShuff = 500;

for i = 1:numel(ratNames)
    ratVar = ratNames{i};
    rat = evalin('base', ratVar);

    % Get target days
    dateList = autoDateList(rat);
    idx = find(strcmp(dateList, rat.An));
    if idx < 3
        warning('%s does not have enough days before An. Skipping...', ratVar);
        continue;
    end
    daysToUse = dateList(idx-2:idx);

    % Filter each structure to only include those days
    spike_struct = filterFieldsByDay(rat.Ca_peaks, daysToUse);
    pos_struct   = filterFieldsByDay(rat.pos, daysToUse);
    ts_struct    = filterFieldsByDay(rat.Ca_ts, daysToUse);

    csus_struct  = filterFieldsByDay(rat.csus90, daysToUse);
    mi_struct    = filterFieldsByDay(rat.MI_wCSUS, daysToUse);


    % Run the shuffle function


    MIshuff = mutualinfo_openfield_shuff_wCSUS(spike_struct, pos_struct, Fs, dim, ts_struct, csus_struct, mi_struct, nShuff);

    % Store back into rat structure
    if ~isfield(rat, 'MI_wCSUS_shuff')
        rat.MI_wCSUS_shuff = struct();
    end

    % Merge shuffled results in

    rat.MI_wCSUS_shuff = structmerge(rat.MI_wCSUS_shuff, MIshuff);


    assignin('base', ratVar, rat);
end
end


function merged = structmerge(a, b)
% Simple utility to merge two structs (b overrides a if duplicate)
    merged = a;
    fieldsB = fieldnames(b);
    for i = 1:numel(fieldsB)
        merged.(fieldsB{i}) = b.(fieldsB{i});
    end
end
