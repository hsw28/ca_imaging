function checkCSUS(A, B)
    % Get the field names of the structures
    fieldNamesA = fieldnames(A);
    fieldNamesB = fieldnames(B);

    % Ensure that both structures have the same fields
    if length(fieldNamesA) ~= length(fieldNamesB)
        error('The structures do not have the same number of fields.');
    end

    % Iterate over the fields
    for i = 1:length(fieldNamesA)
        fieldNameA = fieldNamesA{i};
        fieldNameB = fieldNamesB{i};
        vectorA = A.(fieldNameA);
        vectorB = B.(fieldNameB);

        if length(vectorB) ~= length(vectorA)
          fprintf('Fields have different numbers of values: %s\n', fieldNameA);
          continue
        end

        % Subtract the vectors
        diffVector = vectorB - vectorA;

        % Check for values outside the specified range
        if any(abs(diffVector) < .745 | abs(diffVector) > 0.755)
            fprintf('Field with out-of-range values: %s\n', fieldNameA);
        end
    end
end
