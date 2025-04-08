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
OUTPATH = fullfile(MAINPATH, 'data\processed_data\exclude_trials\');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_ERP, INPATH_BEH, OUTPATH_ERP, OUTPATH_BEH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH_ERP)
tid_psam_clean_up_folder_TD(OUTPATH_BEH)
tid_psam_clean_up_folder_TD(OUTPATH)


% Variables to edit
n_trials = 960;

% Get directory content
dircont_subj_erp = dir(fullfile(INPATH_ERP, 'sub-*.set'));
dircont_subj_beh = dir(fullfile(INPATH_BEH, 'sub-*.xlsx'));

% Sanity Check: Same length of directory contents
if length(dircont_subj_erp) == length(dircont_subj_beh)
else
    error('Different number of files')
end

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
%%wb = waitbar(0,'starting tid_psam_exclude_trials.m');

clear subj_idx 
for subj_idx= 1:length(dircont_subj_erp)

    % Get current ID & Sanity Check: Same ID for ERP data and behavioural data
    subj_erp = dircont_subj_erp(subj_idx).name;
    subj_erp = regexp(subj_erp, 'sub-\d+', 'match', 'once');
    subj_beh = dircont_subj_erp(subj_idx).name;
    subj_beh = regexp(subj_beh, 'sub-\d+', 'match', 'once');

    if strcmp(subj_erp,subj_beh)
        subj = subj_erp;
    else 
        error('Subject IDs do not match')
    end

    % Update progress bar
    %%waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_exclude_trials.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load behavioural data
    beh = readtable(fullfile(INPATH_BEH, [subj '_beh_preprocessed.xlsx']));
    load(fullfile(INPATH_BEH, [subj '_excluded_trials_beh_preprocessing.mat'])) % Loads variable excluded_trials_beh 

    % Load ERP data
    EEG = pop_loadset('filename',[subj '_erp_preprocessed.set'],'filepath',INPATH_ERP);

    % Sanity Check: Same number of trials (=960)
    if length(excluded_trials_beh) == n_trials && size(EEG.data,3)  == n_trials && height(beh)  == n_trials
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'number_of_trials';
    end

    % Merge to be excluded trials
    excluded_trials_all = excluded_trials_beh' | EEG.reject.rejglobal;

    % Get number of excluded trials
    n_excluded = sum(excluded_trials_all);

    % Remove to be excluded trials from vocal data

    % Remove to be excluded trials from ERP data

    % Sanity Check: Number of trials in each condition


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
writetable(protocol,fullfile(OUTPATH, 'tid_psam_exclude_trials_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_exclude_trials_marked_subj.xlsx'))
end

%%close(wb)

check_done = 'tid_psam_exclude_trials_DONE'
