function results = MI_findgood(MIshuff_structure)

%goes through every field of a shuffled mutual info structure and tells you the number of cells over 95% for each field


MI_shuff = MIshuff_structure;

% Get all field names of the MI_shuff structure
fieldNames = fieldnames(MI_shuff);

% Preallocate a vector to store the results
results = zeros(length(fieldNames), 1);

% Loop through each field name
for i = 1:length(fieldNames)
    % Current field name
    currentField = fieldNames{i};

    % Perform the operation on the current field
    currentData = MI_shuff.(currentField);

    % Check if the field has at least 3 columns
    if size(currentData, 2) >= 3
        % Calculate the ratio
        ratio = length(find(currentData(:,3) >= 0.95)) / length(currentData);
    else
        % Handle the case where there are less than 3 columns
        ratio = NaN;  % Assign NaN or some other value to indicate an error
        fprintf('Field %s does not have enough columns.\n', currentField);
    end

    % Store the result in the results vector
    results(i) = ratio;
end
