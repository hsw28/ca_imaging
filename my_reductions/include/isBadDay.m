function skip = isBadDay(animal, dateStr)
    csField = ['CS_' dateStr];
    skip = ~isfield(animal.CS_times, csField) || isempty(animal.CS_times.(csField)) || numel(animal.CS_times.(csField)) < 2 || all(isnan(animal.CS_times.(csField)));
end
