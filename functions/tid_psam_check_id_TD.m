function tid_psam_check_id_TD(varargin)
% tid_psam_check_folder_TD checks if the specified subject ID path already exist.
% If a subject ID path is already present, the function warns the user but performs no action(!).
%
% Usage:
%   tid_psam_check_id_TD('Path1', 'Path2', ...)
%
% Inputs:
%   varargin - any number of folder paths as strings
%
% Tim Dressler, 18.03.2024

existingFolders = {};

for i = 1:numel(varargin)
    folderPath = varargin{i};
    if ~exist(folderPath)
        continue
    else
        existingFolders{end+1} = folderPath;
    end
end


if ~isempty(existingFolders)
    folderList = strjoin(existingFolders, '\n');
    message = sprintf('The following folders exist:\n%s', folderList);

    choice = questdlg(message, ...
        'Warning', ...
        'IGNORE WARNING', 'ARCHIVE PREVIOUS FOLDERS', 'IGNORE WARNING');

    switch choice
        case 'IGNORE WARNING'
            fprintf('User chose to ignore the warning.\n');
        case 'ARCHIVE PREVIOUS FOLDERS'
            fprintf('User chose to archive previous folders.\n');
            %%archiveFolders(existingFolders);
    end
else
    fprintf('No existing folders detected - SubjectID OK \n');
end

