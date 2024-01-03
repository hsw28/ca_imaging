function concatenatedData = concatCSV

[fileNames, pathNames] = uigetfile('*.csv', 'Select CSV Files', 'MultiSelect', 'on');

% Check if any file is selected
if isequal(fileNames,0)
    disp('No files selected');
    return;
end

% Initialize an empty table for concatenation
concatenatedData = [];

% Ensure fileNames is a cell array, even for a single file selection
if ~iscell(fileNames)
    fileNames = {fileNames};
end

% Loop through selected files
for i = 1:length(fileNames)
    % Full path to the file
    fullPath = fullfile(pathNames, fileNames{i});

    % Read CSV file into a table
    tbl = readtable(fullPath);

    % Concatenate the table vertically
    concatenatedData = [concatenatedData; tbl];
end
