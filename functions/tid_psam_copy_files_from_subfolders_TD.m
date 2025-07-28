function tid_psam_copy_files_from_subfolders_TD(baseFolder, folderPattern, filePattern, destinationFolder)
% tid_psam_copy_files_from_subfolders_TD - Copies specific files from subfolders to a destination folder
%
% Usage:
%   tid_psam_copy_files_from_subfolders_TD(baseFolder, folderPattern, filePattern, destinationFolder)
%
% Inputs:
%   baseFolder        - Path to the folder containing the subfolders
%   folderPattern     - Pattern to match subfolders (e.g., 'sub-*')
%   filePattern       - Pattern to match files inside each subfolder (e.g., 'sub-*_accuracy_heatmaps.*')
%   destinationFolder - Name of the destination folder to create within baseFolder
%
% Description:
%   This function searches for subfolders within a specified parent folder that match
%   a given pattern. For each matching subfolder, it copies all files that match a 
%   specified pattern into a new destination folder within the parent directory.
%
%   This is useful for collecting specific files (e.g., figures, logs, metrics) from 
%   multiple subject or session subfolders into one central location for review or analysis.
%
% Example:
%   tid_psam_copy_files_from_subfolders_TD(...
%       '/Users/you/Projects/data', ...
%       'sub-*', ...
%       'sub-*_accuracy_heatmaps.*', ...
%       'individual_accuracy_heatmaps');
%
% Tim Dressler, 28.07.25

% Check if the base folder exists
if ~isfolder(baseFolder)
    error('Specified base folder does not exist');
end

% Full path to the destination folder
targetFolder = fullfile(baseFolder, destinationFolder);

% Create destination folder if it does not exist
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

% Get list of subfolders matching the pattern
subFolders = dir(fullfile(baseFolder, folderPattern));
subFolders = subFolders([subFolders.isdir]);

% Loop through each matching subfolder
for i = 1:length(subFolders)
    subFolderPath = fullfile(baseFolder, subFolders(i).name);

    % Get files matching the file pattern
    files = dir(fullfile(subFolderPath, filePattern));

    for j = 1:length(files)
        sourceFile = fullfile(subFolderPath, files(j).name);
        destinationFile = fullfile(targetFolder, files(j).name);

        % Copy file
        copyfile(sourceFile, destinationFile);
    end
end

fprintf('Matching files copied to: %s\n', targetFolder);
end
