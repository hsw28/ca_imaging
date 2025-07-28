function run_bitsper_CSUS15_Shuffle
% Wrapper to run MI shuffle for only An and 2 days before for all 5 rats
%returns top 5%, mean, shuffled percent (want 0.95 or higher), and hertz during running periods

%IN PREP
%need to input ca_bitsper instead of MI

ratNames = {'rat0816'};

%ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};
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
    CSUS_structure   = filterFieldsByDay(rat.csus15, daysToUse);
    ts_struct    = filterFieldsByDay(rat.Ca_ts, daysToUse);
    bitsper_struct   = filterFieldsByDay(rat.bitsperCSUS, daysToUse)

    % Run the shuffle function


    [per_spike per_sec] = bitsper_CSUS15_shuff(spike_struct, CSUS_structure, 0, 500, bitsper_struct);

    % Store back into rat structure
    if ~isfield(rat, 'per_spike_shuffCSUS')
        rat.per_spike_shuffCSUS = struct();
    end
    if ~isfield(rat, 'per_sec_shuffCSUS')
        rat.per_sec_shuffCSUS = struct();
    end

    % Merge shuffled results in
    rat.per_spike_shuffCSUS = structmerge(rat.per_spike_shuffCSUS, per_spike);
    rat.per_sec_shuffCSUS = structmerge(rat.per_sec_shuffCSUS, per_sec);

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
