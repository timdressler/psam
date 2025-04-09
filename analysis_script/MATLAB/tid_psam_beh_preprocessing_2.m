% tid_psam_beh_preprocessing_2.m
%
% Performs preprocessing of behavioural vocal data.
%
% Preprocessing includes the following steps
%
% Loads log file and removes not needed columns
% Loads vocal data (preprocessed in Praat, see tid_psam_beh_preproecessing_1.praat)
% Merges data
%
% Marks trials for exclusion based on the following criteria
    % If a trial's F0 differs more than 3 SDs from the mean F0
    % If a 'Passive' trial includes a reponses
    % If a 'Active' trial includes no response
%
% Stores dataset
%
% Note. Trials are only marked for exclusion but not excluded yet!
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
INPATH = fullfile(MAINPATH, 'data\BIDS\');
INPATH_PRAAT = fullfile(MAINPATH, 'data\processed_data\beh_preprocessed_1\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\beh_preprocessed_2\');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
n_trials = 960;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*'));

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
%%wb = waitbar(0,'starting tid_psam_beh_preprocessing_2.m');

clear subj_idx 
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    %%waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_beh_preprocessing_2.m'])

    tic;

    % Get log file & Sanity Check: One log file per subject
    subj_log_filename = dir(fullfile(INPATH, [subj '\beh\*.csv']));
    if numel(subj_log_filename) == 1
        subj_log = readtable(fullfile(subj_log_filename.folder, subj_log_filename.name));
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'number_of_log_files';
    end

    % Remove all rows from log file which don not resembel experimental trials (e.g. Instruction trials)
    subj_log_cleaned = subj_log(~isnan(subj_log.("mic_started")), :);
    % Only included needed columns
    subj_log_cleaned = subj_log_cleaned(:, {'ISI_started', 'ISI_stopped', 'go_stim_started', 'conditions_file', ...
        'date', 'expName', 'go_port_started', 'mic_clip', ...
        'mic_started', 'mic_stopped', 'probe', 'probe_duration', 'probe_intensity', 'probe_marker', ...
        'probe_marker_port_started', 'probe_onset', 'probe_onset_cat', 'probe_type', ...
        'psychopyVersion', 'rec_duration', 'stim_file', 'subj', 'task', 'task_marker', ...
         'trial_started'});

    % Sanity Check: Equal db for unaltered and altered probes
    if height(subj_log_cleaned) == n_trials
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'log_dimensions';
    end
    % Extract filename from full path for each vocal recording in order to match it with the table produced in Praat
    subj_log_cleaned.filename_tab = cellfun(@(x) strrep(regexp(x, 'recording_mic_.*(?=\.wav)', 'match', 'once'), '.', '_'), subj_log_cleaned.("mic_clip"), 'UniformOutput', false);

    % Get probe properties file & Sanity Check: Correct dimensions
    subj_probe_properties = readtable(fullfile(INPATH, ['stimuli\' subj '\' subj '_probe_properties.xlsx']));
    if height(subj_probe_properties) == 1 && width(subj_probe_properties) == 7
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'probe_properties_dimensions';
    end
    % Sanity Check: Equal db for unaltered and altered probes
    if round(subj_probe_properties.db_tab_normal) == round(subj_probe_properties.db_tab_pitched)
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'probe_properties_db';
    end

    % Get F0 & RT Table & Sanity Check: Correct dimensions
    subj_f0_rt = readtable(fullfile(INPATH_PRAAT, [subj '_f0_rt_table.csv']), Delimiter=',');
    subj_f0_rt = standardizeMissing(subj_f0_rt,9999);

    if height(subj_f0_rt) == n_trials && width(subj_f0_rt) == 14
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'f0_rt_dimensions';
    end

    % Concatinate log file, probe properties file and F0 & RT Table
    subj_full = join(subj_f0_rt, subj_log_cleaned, 'Keys', 'filename_tab');
    % Exclude 'filename_tab' column and repeat the remaining columns for each row of subj_full
    subj_full = [subj_full, repmat(subj_probe_properties(:, setdiff(1:width(subj_probe_properties), find(strcmp(subj_probe_properties.Properties.VariableNames, 'filename_tab')))), height(subj_full), 1)];

    % Exlcude trials with outlier F0
    % Get threshold
    mean_f0 = mean(subj_full.f0_tab, 'omitnan');
    std_f0 = std(subj_full.f0_tab, 'omitnan');
    threshold_upper = mean_f0 + 3 * std_f0;
    threshold_lower = mean_f0 - 3 * std_f0;
    % Add a column indicating F0 outliers
    subj_full.f0_outlier = (subj_full.f0_tab > threshold_upper) | (subj_full.f0_tab < threshold_lower);
    % Excluded trials based on F0 outliers
    excluded_trials_outliers = (subj_full.f0_tab > threshold_upper) | (subj_full.f0_tab < threshold_lower);

    % Exlude trial with incorrect responses (e.g. vocalizing in 'Passive' trials)
    % Add a column indicating (in)correct responses based on classification in Praat (see vocal_analysis.praat)
    pas_correct = strcmp(subj_full.condition_tab, 'pas') & subj_full.vocal_response_tab == 0;
    act_correct = strcmp(subj_full.condition_tab, 'act') & subj_full.vocal_response_tab == 1;
    subj_full.correct_vocal_response(pas_correct | act_correct) = 1;
    % Excluded trials based on (in)correct vocal responses
    excluded_trials_responses = ~subj_full.correct_vocal_response;

    % Combine excluded trials
    excluded_trials_beh = excluded_trials_outliers | excluded_trials_responses;

    % Only included needed columns
    subj_full_cleaned = subj_full(:, {'mic_clip', ...
        'mic_started', 'mic_stopped', 'probe', ...
        'probe_onset', 'probe_onset_cat', 'probe_type', ...
        'subj', 'task', 'task_marker', 'f0_tab', 'rt_tab', 'condition_tab', ...
        'min_intensity_tab', 'max_intensity_tab', 'vocal_response_tab', ...
        'duration_vocal_tab','onset_vocal_tab', 'filename_tab', 'f0_tab_normal', ...
        'f0_tab_pitched', 'loudness', 'change_attempts', 'correct_vocal_response', ...
        'f0_outlier'});

    writetable(subj_full_cleaned,fullfile(OUTPATH, [subj '_beh_preprocessed.xlsx']))

    save(fullfile(OUTPATH, [subj '_excluded_trials_beh_preprocessing.mat']),"excluded_trials_beh")

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
writetable(protocol,fullfile(OUTPATH, 'tid_psam_beh_preprocessing_2_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_beh_preprocessing_2_marked_subj.xlsx'))
end

%%close(wb)

check_done = 'tid_psam_beh_preprocessing_2_DONE'
