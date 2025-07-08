function outStruct = filterFieldsByDay(inStruct, keepDates)
% Keeps only the fields from inStruct whose suffix matches one of keepDates
    allFields = fieldnames(inStruct);
    outStruct = struct();
    for i = 1:numel(allFields)
        fld = allFields{i};
        if any(contains(fld, keepDates))
            outStruct.(fld) = inStruct.(fld);
        end
    end
end
