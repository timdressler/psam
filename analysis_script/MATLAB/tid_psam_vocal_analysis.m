% tid_psam_vocal_analysis.m
%
% Performs analysis of vocal data.
% Has to be executed AFTER tid_psam_vocal_analysis.praat.
%
% Tim Dressler, 04.04.2025

clear
close all
clc

set(0,'DefaultTextInterpreter','none')

% Set up paths
SCRIPTPATH = cd;
if regexp(SCRIPTPATH, regexptranslate('wildcard','*psam\analysis_script\MATLAB')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, '\analysis_script\MATLAB');
INPATH = fullfile(MAINPATH, 'data\BIDS\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\vocal_analysis');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit


% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*'));

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
%%wb = waitbar(0,'starting tid_psam_vocal_analysis.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    %%waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_vocal_analysis.m'])

    % Get log file & Sanity Check: One log file per subject
    subj_log_filename = dir(fullfile(INPATH, [subj '\beh\*.csv']));
    if numel(subj_log_filename) == 1
        subj_log = readtable(fullfile(subj_log_filename.folder, subj_log_filename.name),'VariableNamingRule', 'preserve');
    else
        %%error('No matching CSV files found.');
    end


    


end
