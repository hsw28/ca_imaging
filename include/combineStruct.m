function combinedStruct = combineStruct(struct1, struct2)
%combines structures and sorts them

    fields1 = fieldnames(struct1);
    fields2 = fieldnames(struct2);
    combinedStruct = struct1;

    for i = 1:length(fields2)
        field = fields2{i};
        combinedStruct.(field) = struct2.(field);
    end


combinedStruct = orderfields(combinedStruct);
