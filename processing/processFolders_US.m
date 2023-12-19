function [US_EMG, CS_EMG] = processFolders_US(currentDir)
%makes a structure of US times for all your days in the specified dirctory
%REMEMBER to STILL CORRECt tiMES USING fix_all_ts.m
%currectDir should be the directory you wish to search from, ie './031423/'


    US_EMG = struct();
    CS_EMG = struct();
    validFoldersList = getValidFolders(currentDir);

    % Prompt the user to select which folders to process
    if ~isempty(validFoldersList)
        validFoldersList
        selected_folders = input('Enter the numbers of the folders to process (e.g., [1, 3, 5]): ');

        for i = 1:numel(selected_folders)
            selected_folder_idx = selected_folders(i);
            if selected_folder_idx > 0 && selected_folder_idx <= numel(validFoldersList)
                selected_folder = validFoldersList{selected_folder_idx};
                [timestamps, dataSamplesUS, dataSamplesCS] = processSelectedFolder(selected_folder);


                % Get the folder name three levels above My_WebCam
                selected_folder
                pos_folder = (selected_folder);
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

                %upsample time
                Fs_original = 60;
                Fs_target = 30720;
                t_original = (0:numel(timestamps)-1) / Fs_original;
                t_target = (0:(numel(timestamps)*Fs_target/Fs_original)-1) / Fs_target;
                timestamps = interp1(t_original, timestamps, t_target, 'linear');


                %make a list of the indices for the peaks
                %using negative bc it actually dips before it rises
                [pks,locs_US] = findpeaks(dataSamplesUS*-1, 'MinPeakHeight', mean(dataSamplesUS)+(1*std(dataSamplesUS)), 'MinPeakDistance', 860160);
                [pks,locs_CS] = findpeaks(dataSamplesCS*-1, 'MinPeakHeight', mean(dataSamplesCS)+(1*std(dataSamplesCS)), 'MinPeakDistance', 860160);


                %relate that to the timestamp
                US_times = timestamps(locs_US);
                CS_times = timestamps(locs_CS);



                if ~isempty(strfind(selected_folder, 'extinction'))
                  outputDate=sprintf('exinction_%s', outputDate);
                end

                US_EMG.(sprintf('US_%s', outputDate)) = US_times;
                CS_EMG.(sprintf('CS_%s', outputDate)) = CS_times;
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

        % Check if the folder has a file whose name starts with 'DC3'
        eventFilesUS = dir(fullfile(parentFolder, folderName, 'DC3*'));
        eventFilesCS = dir(fullfile(parentFolder, folderName, 'DC1*'));

        % Initialize variables to keep track of the largest 'DC3' file size and its index
        largestFileSize = 0;
        largestFileIndex = 0;

        % Loop through all 'DC3' files and find the largest one
        for i = 1:length(eventFilesUS)
          fileSize = eventFilesUS(i).bytes;
          if fileSize > largestFileSize
            largestFileSize = fileSize;
            largestFileIndexUS = i;
          end
        end

        % Loop through all 'DC1' files and find the largest one
        for i = 1:length(eventFilesCS)
          fileSize = eventFilesCS(i).bytes;
          if fileSize > largestFileSize
            largestFileSize = fileSize;
            largestFileIndexCS = i;
          end
        end



        if ~isempty(eventFilesCS)
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


function [timestamps,dataSamplesUS, dataSamplesCS] = processSelectedFolder(parentFolder)
  % Initialize variables to keep track of the largest file size and its path
 largestFileSize_US = 0;
 largestFilePath_US = '';
 largestFileSize_CS = 0;
 largestFilePath_CS = '';

        % Check if the folder has a file whose name starts with 'Events'
        eventFilesUS = dir(fullfile(parentFolder, 'DC3*'));
        eventFilesCS = dir(fullfile(parentFolder, 'DC1*'));

        if ~isempty(eventFilesUS)
            % Loop through all 'Events' files in the folder
            for j = 1:length(eventFilesUS)
                filePath = fullfile(parentFolder, eventFilesUS(j).name);

                % Get the file size
                fileInfo = dir(filePath);
                fileSize = fileInfo.bytes;

                % Convert 1 MB to bytes (1 MB = 1,048,576 bytes)


                % Update largest file information if this file is larger
                if fileSize > largestFileSize_US
                    largestFileSize_US = fileSize;
                    largestFilePath_US = filePath;

                end
            end
        end

        if ~isempty(eventFilesCS)
            % Loop through all 'Events' files in the folder
            for j = 1:length(eventFilesCS)
                filePath = fullfile(parentFolder, eventFilesCS(j).name);

                % Get the file size
                fileInfo = dir(filePath);
                fileSize = fileInfo.bytes;

                % Convert 1 MB to bytes (1 MB = 1,048,576 bytes)


                % Update largest file information if this file is larger
                if fileSize > largestFileSize_CS
                    largestFileSize_CS = fileSize;
                    largestFilePath_CS = filePath;

                end
            end
        end


    % Check if any 'DC3' file was found
    if isempty(largestFilePath_CS)
        disp('No file starting with ''DC3'' found');
    else
        % Run the  function on the largest 'DC3' file
        [timestamps,dataSamplesUS] = getRawCSCData(largestFilePath_US, 1, 10000000000000000000000);
        [timestamps,dataSamplesCS] = getRawCSCData(largestFilePath_CS, 1, 10000000000000000000000);


    end
end
