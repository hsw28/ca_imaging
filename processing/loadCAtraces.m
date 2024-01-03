function [peak_stuct all_traces frame_timestamps] = loadCAtraces(currentDir)

%loads calcium transients and/or peaks for all sorted files in the selected folder

all_traces = struct();
peak_stuct = struct();
frame_timestamps = struct();
validFoldersList = getValidFolders(currentDir);

% Prompt the user to select which folders to process
if ~isempty(validFoldersList)
    validFoldersList
    selected_folders = input('Enter the number of the folders to process (e.g., [1, 3, 5]): ');

    for i = 1:numel(selected_folders)
        selected_folder_idx = selected_folders(i);
        if selected_folder_idx > 0 && selected_folder_idx <= numel(validFoldersList)
            selected_folder = validFoldersList{selected_folder_idx}



            [trace_file sorting timestamps] = processSelectedFolder(selected_folder);


            %this gets cells that were sorted as good
            good = sorting.validCNMFE;
            good = find(good==1);

            %this is what loads your traces
            traces = trace_file.cnmfeAnalysisOutput.extractedSignals;
            traces = traces(good,:);

            %gets peaks
            %[signalPeaks, signalPeaksArray] = computeSignalPeaks(traces,'numStdsForThresh',2,'minTimeBtEvents',0,'detectMethod','diff','medianFilterLength',225);

            peaks = trace_file.cnmfeAnalysisOutput.extractedPeaks;
            signalPeaks = peaks(good,:);
            %converts frames to timestamps


            [peaks times_sec] = converttotime(signalPeaks, timestamps);


            % Get the folder name three levels above My_WebCam
            folder = (fileparts(fileparts((selected_folder))));
            % Split the text using '/' as the delimiter
            folder_parts = strsplit(folder, '\');
            % Get the last part, which contains the date
            folder_date = folder_parts{end};

            folder = strsplit(selected_folder, '\');
            folder = folder{end};


            if strcmp(folder_date, 'extinction')==1
              folder = fileparts(fileparts(fileparts((selected_folder))));
              folder_parts = strsplit(folder, '\');
              folder_date = folder_parts{end};
              folder_date=sprintf('exinction_%s', folder_date);
            end




                            all_traces.(sprintf('CA_traces_%s', folder_date)) = traces;
                            peak_stuct.(sprintf('CA_peaks_%s', folder_date)) = peaks;
                            frame_timestamps.(sprintf('CA_frame_ts_%s', folder_date)) = times_sec;

                        else
                            fprintf('Invalid folder number: %d\n', selected_folder_idx);
                            %}
                        end
                    end
                else
                    fprintf('No valid folders found.\n');
                end

    peak_stuct = orderfields(peak_stuct);
    all_traces = orderfields(all_traces);

end





  function validFoldersList = getValidFolders(parentFolder)
        % Initialize an empty list to store valid folders
        validFoldersList = {};

        % Get a list of all subfolders within the parent folder
        subfolders = dir(parentFolder);
        subfolders = subfolders([subfolders.isdir]); % Keep only directories, remove files

        % Loop through each subfolder
        for i = 1:length(subfolders)
            folderName = subfolders(i).name;

            % Skip the current and parent folder entries
            if strcmp(folderName, '.') || strcmp(folderName, '..')
                continue;
            end

            % Check if the folder has a file whose name ends with '_cnmfeAnalysis.mat'
            eventFiles = dir(fullfile(parentFolder, folderName, '*cnmfeAnalysisSorted*'));

            % Initialize variables to keep track of the largest 'Events' file size and its index
            largestFileSize = 0;
            largestFileIndex = 0;

            % Loop through all 'Events' files and find the largest one
            for i = 1:length(eventFiles)
              fileSize = eventFiles(i).bytes;
              if fileSize > largestFileSize
                largestFileSize = fileSize;
                largestFileIndex = i;
              end
            end


            if ~isempty(eventFiles)
              %  if largestFileSize > 1e6  % 1 MB = 1,048,576 bytes
                % If there is at least large 'CNMFE' file, add the folder to the list
                validFoldersList{end+1} = fullfile(parentFolder, folderName);
              %end
            end


            % Recursively call the function for the current subfolder
            subValidFolders = getValidFolders(fullfile(parentFolder, folderName));

            % Merge the subfolder list with the main list
            validFoldersList = [validFoldersList, subValidFolders];
        end
end


    function [trace_file sorting timestamps] = processSelectedFolder(parentFolder)
      % Initialize variables to keep track of the largest file size and its path
     largestFileSize = 0;
     largestFilePath = '';

      % Check if the folder has a file whose name ends with '_cnmfeAnalysis.mat'
              eventFiles = dir(fullfile(parentFolder, '*cnmfeAnalysis.mat'));
              sorted = dir(fullfile(parentFolder, '*cnmfeAnalysisSorted*.mat'));
              timestamps = dir(fullfile(parentFolder, 'timeStamps.csv'));



            if ~isempty(eventFiles)
                % Loop through all 'Events' files in the folder
                for j = 1:length(eventFiles)
                    filePath = fullfile(parentFolder, eventFiles(j).name);

                    % Get the file size
                    fileInfo = dir(filePath);
                    fileSize = fileInfo.bytes;

                    % Convert 1 MB to bytes (1 MB = 1,048,576 bytes)


                    % Update largest file information if this file is larger
                    if fileSize > largestFileSize
                        largestFileSize = fileSize;
                        largestFilePath = filePath;

                    end
                end
            end


        % Check if any 'Events' file was found
        if isempty(largestFilePath)
             dir(fullfile(parentFolder, '*cnmfeAnalysis.mat'));

            disp('No file starting with ''Events'' found');
        else
            % Run the getRawTTLs function on the largest 'Events' file
            dir(fullfile(parentFolder, '*cnmfeAnalysis.mat'));
            trace_file = load(largestFilePath);

            filePath = fullfile(parentFolder, sorted.name);
            
            sorting = load(filePath);

            filePath = fullfile(parentFolder, timestamps.name);
            filePath
            timestamps = readtable(filePath);
        end
    end
