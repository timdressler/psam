% tid_psam_exclude_trials.m
%
%
% Tim Dressler, 04.04.2025

% Add RT calculation rt_tab - (go_onset-mic.started)

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
INPATH_ERP = fullfile(MAINPATH, 'data\processed_data\erp_preprocessed\');
INPATH_BEH = fullfile(MAINPATH, 'data\processed_data\beh_preprocessed_2\');
OUTPATH_ERP = fullfile(MAINPATH, 'data\processed_data\erp_preprocessed_clean\');
OUTPATH_BEH = fullfile(MAINPATH, 'data\processed_data\beh_preprocessed_clean\');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH_ERP, OUTPATH_BEH)
tid_psam_clean_up_folder_TD(OUTPATH_ERP, OUTPATH_BEH)

% Variables to edit
n_trials = 960;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*'));

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
%%wb = waitbar(0,'starting tid_psam_vocal_analysis.m');

clear subj_idx n_exluded_trials
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    %%waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_vocal_analysis.m'])

    tic;
    % Load vocal data

    % Load ERP data

    % Sanity Check: Same number of trials (=960)

    % Merge to be excluded trials

    % Get number of excldued trials
    % Stop processing for subjects with too little remaining trials

    % Remove to be excluded trials from vocal data

    % Remove to be excluded trials from ERP data


    % Update Protocol
    subj_time = toc;
    protocol{subj_idx,1} = subj;
    protocol{subj_idx,2} = subj_time;
    if any(strcmp(marked_subj, subj), 'all')
        protocol{subj_idx,3} = 'MARKED';
    else
        protocol{subj_idx,3} = 'OK';
    end

end

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'tid_psam_vocal_preprocessing_2_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_vocal_preprocessing_2_marked_subj.xlsx'))
end

%%close(wb)

check_done = 'tid_psam_vvocal_preprocessing_2_DONE'
