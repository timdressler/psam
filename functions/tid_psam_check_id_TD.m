function tid_psam_check_id_TD(varargin)
% tid_psam_check_folder_TD checks if the specified subject ID path already exist.
% If a subject ID path is already present, the function warns the user and gives 
% the option to either ignore the warning (which potentially leads to files being
% overwritten) or to move the files to an archive (see tid_psam_archive_subj).
% No files are deleted, only moved.
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
            fprintf('Warning ignored. Potentially overwriting files...\n');
        case 'ARCHIVE PREVIOUS FOLDERS'
            fprintf('Archiving files... \n');
            subj = regexp(existingFolders{1}, 'sub-(\d+)', 'tokens');
            subj = subj{1}{:};
            subj = str2num(subj);
            tid_psam_archive_subj_TD(subj);
    end
else
    fprintf('No existing folders detected - SubjectID OK \n');
end

