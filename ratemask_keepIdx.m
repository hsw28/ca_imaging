function keepIdx = ratemask_keepIdx(ratName, dayStr, nCells)
% Returns indices of cells to keep for this day (based on ratemask).
% If ratemask missing/mismatch -> [] (skip day).

keepIdx = [];

% rat struct from base workspace
if ~evalin('base', sprintf('exist(''%s'',''var'')', ratName))
    return
end
rat = evalin('base', ratName);

maskField = sprintf('ratemask_%s', dayStr);

if ~isfield(rat,'ratemask') || ~isfield(rat.ratemask, maskField)
    return
end

mv = rat.ratemask.(maskField);
mv = mv(:);  % force column

if numel(mv) ~= nCells
    return
end

keepIdx = find(mv == 1);
end
