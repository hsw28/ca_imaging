function giantStructure = processFolders_CaTime(currentDir, type1_for_cage_2_for_oval)
%makes a position structure with all your run days for one animal using dlc_fixposor
%should auto detect if oval or rectangle
%%FIXES POS USING dlc_fixpos
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

                % Get the folder name three levels above My_WebCam
                ca_folder = (fileparts(fileparts((selected_folder))));
                % Split the text using '/' as the delimiter
                ca_folder_parts = strsplit(ca_folder, '/');
                % Get the last part, which contains the date
                ca_folder_date = ca_folder_parts{end};

                ca_folder = strsplit(selected_folder, '/');
                ca_folder = ca_folder{end};


                if strcmp(ca_folder_date, 'extinction')==1
                  ca_folder = fileparts(fileparts(fileparts((selected_folder))));
                  ca_folder_parts = strsplit(ca_folder, '/');
                  ca_folder_date = ca_folder_parts{end};
                  ca_folder_date=sprintf('exinction_%s', ca_folder_date);
                en
                  giantStructure.(sprintf('pos_%s', ca_folder_date)) = output;

            else
                fprintf('Invalid folder number: %d\n', selected_folder_idx);
            end
        end
    else
        fprintf('No valid folders found.\n');
    end

    giantStructure = orderfields(giantStructure);
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
        webCamDir = fullfile(subDir, 'My_V4_Miniscope');

        % Check if the current directory contains My_WebCam folder
        if isfolder(webCamDir)
            csvFiles = dir(fullfile(webCamDir, '*.csv'));

            % Check if My_WebCam folder contains two CSV files with timeStamps.csv
            if numel(csvFiles) == 2 && any(strcmpi({csvFiles.name}, 'timeStamps.csv'))
                % Get the folder name three levels above My_WebCam
                ca_folder = (fileparts(fileparts(fileparts(webCamDir))));

                % Split the text using '/' as the delimiter
                ca_folder_parts = strsplit(ca_folder, '/');

                % Get the last part, which contains the date
                ca_folder_date = ca_folder_parts{end};

                % Print the folder and date information
                fprintf('Folder: %s, Date: %s\n', webCamDir, ca_folder_date);

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

function [output type1_for_cage_2_for_oval] = processSelectedFolder(folder)

    % Get the full file path for timestamps
    webCamDir = fullfile(folder, 'My_V4_Miniscope');
    timeStampFile = fullfile(folder, 'timeStamps.csv');

    % Read the timestamps data from the CSV file using readtable
    timeStampData = readtable(timeStampFile);
    timeStampData = table2array(timeStampData);

    output = dlc_fixpos(timeStampData);

end
