function processFolders_vid(currentDir, giantStructure)

  % Get a list of all directories in the current directory
     dirs = dir(currentDir);
     dirNames = {dirs([dirs.isdir]).name};
     dirNames = dirNames(~ismember(dirNames, {'.', '..'}));

     % Loop through each directory in the current directory
     for i = 1:numel(dirNames)
         subDir = fullfile(currentDir, dirNames{i});
         webCamDir = fullfile(subDir, 'My_WebCam');

         % Check if the current directory contains My_WebCam folder
         if isfolder(webCamDir)
             csvFiles = dir(fullfile(webCamDir, '*.csv'));

             % Check if My_WebCam folder contains two CSV files with timeStamps.csv
             if numel(csvFiles) == 2 && any(strcmpi({csvFiles.name}, 'timeStamps.csv'))
                 % Get the folder name three levels above My_WebCam
                 pos_folder = (fileparts(fileparts(webCamDir)))
                 % Split the text using '/' as the delimiter
                 pos_folder = strsplit(pos_folder, '/');
                 % Get the last part, which contains the date
                 pos_folder = pos_folder{end};

                 % Get the full file paths for pos and timestamps
                 posFile = fullfile(webCamDir, csvFiles(~strcmpi({csvFiles.name}, 'timeStamps.csv')).name);
                 timeStampFile = fullfile(webCamDir, 'timeStamps.csv');

                 % Read the data from the CSV files using readtable
                 posData = readtable(posFile);
                 timeStampData = readtable(timeStampFile);


                 % Call the function dlc_fixpos and store the output in the giant structure
                 output = dlc_fixpos(posData, timeStampData, 1);
                 test = (sprintf('pos_%s', pos_folder))
                 giantStructure.(sprintf('pos_%s', pos_folder)) = output;
             end
         end

         % If the subdirectory contains more directories, recursively process them
         if isfolder(subDir)
             processFolders_vid(subDir, giantStructure);
         end
     end
