function concatenatedData = concatCSV

  % Ask the user to select CSV files
  [fileNames, pathNames] = uigetfile('*.csv', 'Select CSV Files', 'MultiSelect', 'on');

  % Check if any file is selected
  if isequal(fileNames,0)
      disp('No files selected');
      return;
  end

  % Ensure fileNames is a cell array, even for a single file selection
  if ~iscell(fileNames)
      fileNames = {fileNames};
  end

  % Initialize variables
  concatenatedData = [];
  firstFile = true;
  firstFileVariableNames = {};

  % Loop through selected files
  for i = 1:length(fileNames)
      % Full path to the file
      fullPath = fullfile(pathNames, fileNames{i});

      % Read CSV file into a table
      tbl = readtable(fullPath);

      % For the first file, keep the original variable names
      if firstFile
          firstFileVariableNames = tbl.Properties.VariableNames;
          firstFile = false;
      else
          % For subsequent files, standardize variable names to match the first file
          tbl.Properties.VariableNames = firstFileVariableNames;
      end

      % Concatenate the table vertically
      concatenatedData = [concatenatedData; tbl];
  end
