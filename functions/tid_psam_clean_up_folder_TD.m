function tid_psam_clean_up_folder_TD(folderPath)
% tid_psam_clean_up_folder_TD.m
% 
% Moves files to a timestamped archive subfolder within an "archive" folder
%
% Usage:
%   tid_psam_clean_up_folder_TD('Path')
%
% Input:
%   folderPath - The path of the folder to check and clean up.
%
% Description:
%   This function checks if a specified folder is empty. If the folder
%   contains any files, it ensures an "archive" folder exists within the
%   specified folder, creates a new subfolder named archive_DD_MM_YYYY_HH-MM-SS
%   inside the archive folder, and moves all files from the original folder
%   to this new archive subfolder.
%
% Tim Dressler, 12.11.24

% Check if the folder exists
if ~isfolder(folderPath)
    error('Specified folder does not exist');
end

% Get list of files in the folder (excluding '.' and '..')
files = dir(folderPath);
files = files(~ismember({files.name}, {'.', '..', 'archive'}));

% If folder is not empty, proceed with archiving
if ~isempty(files)
    % Define the archive folder path
    archiveFolderPath = fullfile(folderPath, 'archive');

    % Create the archive folder if it does not exist
    if ~exist(archiveFolderPath, 'dir')
        mkdir(archiveFolderPath);
    end

    % Get current date and time in DD_MM_YYYY_TIME format
    currentDateTime = datestr(now, 'dd_mm_yyyy_HH-MM-SS');

    % Greate the timestamped archive subfolder within the archive folder
    timestampedArchiveFolder = fullfile(archiveFolderPath, ['archive_' currentDateTime]);

    % Create the timestamped archive folder
    mkdir(timestampedArchiveFolder);

    % Move files to the timestamped archive folder
    for k = 1:length(files)
        sourceFile = fullfile(folderPath, files(k).name);
        destinationFile = fullfile(timestampedArchiveFolder, files(k).name);

        % Move the file
        movefile(sourceFile, destinationFile);
    end

    fprintf('Files moved to: %s\n', timestampedArchiveFolder);
else
    fprintf('The folder is empty. No action taken.\n');
end
end
