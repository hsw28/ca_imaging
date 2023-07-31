function giantStructure = processFolders_Pos(currentDir, type1_for_cage_2_for_oval)
%makes a position structure with all your run days for one animal using dlc_fixpos
%%FIXES POS
%input 1 for the cage and 2 for the oval
%currectDir should be the directory you wish to search from, ie './031423/'

    giantStructure = struct()
    validFolders = getValidFolders(currentDir);

    % Prompt the user to select which folders to process
    if ~isempty(validFolders)
        selected_folders = input('Enter the numbers of the folders to process (e.g., [1, 3, 5]): ');

        % Process the selected folders using dlc_fixpos
        for i = 1:numel(selected_folders)
            selected_folder_idx = selected_folders(i);
            if selected_folder_idx > 0 && selected_folder_idx <= numel(validFolders)
                selected_folder = validFolders{selected_folder_idx};
                output = processSelectedFolder(selected_folder, type1_for_cage_2_for_oval);


                % Get the folder name three levels above My_WebCam
                pos_folder = (fileparts(fileparts((selected_folder))));
                % Split the text using '/' as the delimiter
                pos_folder_parts = strsplit(pos_folder, '/');
                % Get the last part, which contains the date
                pos_folder_date = pos_folder_parts{end};

                pos_folder = strsplit(selected_folder, '/');
                pos_folder = pos_folder{end};


                if strcmp(pos_folder_date, 'extinction')==1
                  pos_folder = fileparts(fileparts(fileparts((selected_folder))));
                  pos_folder_parts = strsplit(pos_folder, '/');
                  pos_folder_date = pos_folder_parts{end};
                  pos_folder_date=sprintf('exinction_%s', pos_folder_date);
                end

                giantStructure.(sprintf('pos_%s', pos_folder_date)) = output;
            else
                fprintf('Invalid folder number: %d\n', selected_folder_idx);
            end
        end
    else
        fprintf('No valid folders found.\n');
    end
end

function validFolders = getValidFolders(currentDir)
    validFolders = {};

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
                pos_folder = (fileparts(fileparts(fileparts(webCamDir))));

                % Split the text using '/' as the delimiter
                pos_folder_parts = strsplit(pos_folder, '/');

                % Get the last part, which contains the date
                pos_folder_date = pos_folder_parts{end};

                % Print the folder and date information
                fprintf('Folder: %s, Date: %s\n', webCamDir, pos_folder_date);

                % Add the folder to the list of valid folders
                validFolders{end + 1} = webCamDir;
            end
        end

        % If the subdirectory contains more directories, recursively process them
        if isfolder(subDir)
            % Get the valid folders from the recursive call
            validFoldersSubDir = getValidFolders(subDir);
            % Append the valid folders from the subdirectory to the list
            validFolders = [validFolders, validFoldersSubDir];
        end
    end
end

function output = processSelectedFolder(folder, type1_for_cage_2_for_oval)

    % Get the full file path for timestamps
    %webCamDir = fullfile(folder, 'My_WebCam');
    timeStampFile = fullfile(folder, 'timeStamps.csv');

    % Read the timestamps data from the CSV file using readtable
    timeStampData = readtable(timeStampFile);

    % Find the pos file in the My_WebCam folder
    posFiles = dir(fullfile(folder, '*.csv'));

    % Check if any file in the My_WebCam folder is not timeStamps.csv
    posFileNames = {posFiles.name};
    posFileName = posFileNames{~strcmpi(posFileNames, 'timeStamps.csv')};

    % Get the full file path for the pos file
    posFile = fullfile(folder, posFileName)

    % Read the data from the pos CSV file using readtable
    posData = readtable(posFile);

    % Call the function dlc_fixpos and store the output in the giant structure
    output = dlc_fixpos(posData, timeStampData, type1_for_cage_2_for_oval);

end
