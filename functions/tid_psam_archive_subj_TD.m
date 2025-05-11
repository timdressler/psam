function tid_psam_archive_subj_TD(subj_number)
% tid_psam_archive_subj_TD - Archive all files for a given subject ID
%
% Usage:
%   tid_psam_archive_subj_TD(subjectID)
%
% Input:
%   subjectID - The number of a subject (only the number, without sub-)
%
% Description:
%   This function moves all files (eeg, beh and stimuli) to a archive
%   folder outside of the BIDS structure. It relies on the current folder
%   structure implemented for the tid_psam project. 
%   No files are deleted, only moved.
%
% Tim Dressler, 19.03.2025

% Format subject ID
if ~isempty(subj_number)
    subj = sprintf('%02d', subj_number); % Ensure two-digit format
    subj = ['sub-' subj]; % Construct final subject ID
end

% Set up general paths
SCRIPTPATH = fileparts(mfilename('fullpath'));
if regexp(SCRIPTPATH, regexptranslate('wildcard','*psam\functions')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, '\functions');

% Set up subject paths
STIMULIPATH = fullfile(MAINPATH, ['data/BIDS/stimuli/' subj]);
SUBJPATH = fullfile(MAINPATH, ['data/BIDS/' subj]);

% Checks whetehr folders exist
if ~isfolder(STIMULIPATH) ||  ~isfolder(SUBJPATH)
    error('Subject-Folder(s) not found.')
end

currentDateTime = datestr(now, 'yyyymmdd_HHMMSS'); % Get current date and time
ARCHIVEPATH = fullfile(MAINPATH, 'data/archive/', [subj '_archive_' currentDateTime]);
archive_STIMULIPATH = fullfile(ARCHIVEPATH, 'stimuli/');
archive_SUBJPATH = ARCHIVEPATH;

tid_psam_check_folder_TD(archive_STIMULIPATH, archive_SUBJPATH)

% Function to move files and folders 
moveToArchive = @(srcPath, destPath) movefile(srcPath, destPath);

% Move STIMULIPATH contents to archive
if isfolder(STIMULIPATH)
    moveToArchive(STIMULIPATH, archive_STIMULIPATH);
end

% Move SUBJPATH contents to archive
if isfolder(SUBJPATH)
    moveToArchive(SUBJPATH, archive_SUBJPATH);
end

clc
fprintf('Files have been archived. \n')
