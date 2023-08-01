function [giantStructure_TS giantStructure_EMG] = processFolders_EMG(currentDir)
%makes a structure of neuralynx events for all your days in the specified dirctory
%currectDir should be the directory you wish to search from, ie './031423/'

    giantStructure_TS = struct();
    giantStructure_EMG = struct();
    validFoldersList = getValidFolders(currentDir);

    % Prompt the user to select which folders to process
    if ~isempty(validFoldersList)
        validFoldersList
        selected_folders = input('Enter the numbers of the folders to process (e.g., [1, 3, 5]): ');

        % Process the selected folders using dlc_fixpos
        for i = 1:numel(selected_folders)
            selected_folder_idx = selected_folders(i);
            if selected_folder_idx > 0 && selected_folder_idx <= numel(validFoldersList)
                selected_folder = validFoldersList{selected_folder_idx};
                [timestamps,dataSamples] = processSelectedFolder(selected_folder);


                % Get the folder name three levels above My_WebCam
                pos_folder = (selected_folder)
                % Find the position of the last slash '/'
                slashPos = find(pos_folder == '/', 1, 'last');
                % Extract the substring after the last slash, which contains the date
                dateSubstring = pos_folder(slashPos+1:end);
                % Find the position of the underscore '_'
                underscorePos = find(dateSubstring == '_', 1, 'first');
                % Extract the substring before the underscore, which is the date
                outputDate = dateSubstring(1:underscorePos-1);
                % Replace hyphens '-' with underscores '_'
                outputDate = strrep(outputDate, '-', '_');


                if ~isempty(strfind(selected_folder, 'extinction'))
                  outputDate=sprintf('exinction_%s', outputDate);
                end

                giantStructure_TS.(sprintf('EMGts_%s', outputDate)) = timestamps';
                giantStructure_EMG.(sprintf('EMG_%s', outputDate)) = dataSamples;
            else
                fprintf('Invalid folder number: %d\n', selected_folder_idx);
            end
        end
    else
        fprintf('No valid folders found.\n');
    end
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

        % Check if the folder has a file whose name starts with 'Events'
        eventFiles = dir(fullfile(parentFolder, folderName, 'DC2*'));
        % Initialize variables to keep track of the largest 'DC2' file size and its index
        largestFileSize = 0;
        largestFileIndex = 0;

        % Loop through all 'DC2' files and find the largest one
        for i = 1:length(eventFiles)
          fileSize = eventFiles(i).bytes;
          if fileSize > largestFileSize
            largestFileSize = fileSize;
            largestFileIndex = i;
          end
        end


        if ~isempty(eventFiles)
            if largestFileSize > 1e6  % 1 MB = 1,048,576 bytes
            % If there is at least one 'Events' file, add the folder to the list
            validFoldersList{end+1} = fullfile(parentFolder, folderName);
          end
        end


        % Recursively call the function for the current subfolder
        subValidFolders = getValidFolders(fullfile(parentFolder, folderName));

        % Merge the subfolder list with the main list
        validFoldersList = [validFoldersList, subValidFolders];
    end
end


function [timestamps,dataSamples] = processSelectedFolder(parentFolder)
  % Initialize variables to keep track of the largest file size and its path
 largestFileSize = 0;
 largestFilePath = '';

        % Check if the folder has a file whose name starts with 'Events'
        eventFiles = dir(fullfile(parentFolder, 'DC2*'));

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


    % Check if any 'DC2' file was found
    if isempty(largestFilePath)
        disp('No file starting with ''DC2'' found');
    else
        % Run the  function on the largest 'DC2' file
        [timestamps,dataSamples]  = getRawCSCData(largestFilePath, 1, 10000000000000000000000);
    end
end
