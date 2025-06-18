% tid_psam_beh_preprocessing_2.m
%
% Performs preprocessing of behavioural vocal data.
%
% Preprocessing includes the following steps
%
    % Loads log file and removes not needed columns
    % Loads vocal data (preprocessed in Praat, see tid_psam_beh_preproecessing_1.praat)
    % Merges data
    % Calculates reaction time relative to the 'go-signal' based on the start
    %   time of the recording
%
    %  Marks trials for exclusion based on the following criteria
    %   If a trial's onset time differs more than 3 SDs from the mean onset time
    %   If a trial's F0 differs more than 3 SDs from the mean F0
    %   If a 'Passive' trial includes a reponses
    %   If a 'Active' trial includes no response
%
    % Saves data
%
% Note. Trials are only marked for exclusion but not excluded yet!
%
% Tim Dressler, 04.04.2025

clear
close all
clc

rng(123)
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
wb = waitbar(0,'starting tid_psam_beh_preprocessing_2.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_beh_preprocessing_2.m'])

    tic;

    % Sanity Check: One log file per subject
    subj_log_filename = dir(fullfile(INPATH, [subj '\beh\s*.csv']));
    if numel(subj_log_filename) == 1
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'number_of_log_files';
    end

    % Get log file
    subj_log = readtable(fullfile(subj_log_filename.folder, subj_log_filename.name));

    % Re-code probe_type and task variable
    subj_log.probe_type(strcmp(subj_log.probe_type, 'Normal')) = {'Unaltered'};
    subj_log.probe_type(strcmp(subj_log.probe_type, 'Pitch')) = {'Altered'};

    subj_log.task(strcmp(subj_log.task, '/ga/')) = {'Active'};
    subj_log.task(strcmp(subj_log.task, '/xx/')) = {'Passive'};

    % Remove all rows from log file which don not resembel experimental trials (e.g. Instruction trials)
    subj_log_clean = subj_log(~isnan(subj_log.("mic_started")), :);
    % Only included needed columns
    subj_log_clean = subj_log_clean(:, {'ISI_started', 'ISI_stopped', 'go_stim_started', 'conditions_file', ...
        'date', 'expName', 'go_port_started', 'mic_clip', ...
        'mic_started', 'mic_stopped', 'probe', 'probe_duration', 'probe_intensity', 'probe_marker', ...
        'probe_marker_port_started', 'probe_onset', 'probe_onset_cat', 'probe_type', ...
        'psychopyVersion', 'rec_duration', 'stim_file', 'subj', 'task', 'task_marker', ...
        'trial_started', 'conditions_file'});

    % Rename columns
    renamed_column = {
        'ISI_started'                 'ISI_started';
        'ISI_stopped'                 'ISI_stopped';
        'go_stim_started'             'go_stim_started';
        'conditions_file'             'conditions_file';
        'date'                        'date';
        'expName'                     'experiment_name';
        'go_port_started'             'go_port_started';
        'mic_clip'                    'recording_path';
        'mic_started'                 'mic_started';
        'mic_stopped'                 'mic_stopped';
        'probe'                       'probe';
        'probe_duration'              'probe_duration';
        'probe_intensity'             'probe_intensity';
        'probe_marker'                'probe_marker';
        'probe_marker_port_started'   'probe_marker_port_started';
        'probe_onset'                 'probe_onset';
        'probe_onset_cat'             'probe_onset_cat';
        'probe_type'                  'probe_type';
        'psychopyVersion'             'psychopyVersion';
        'rec_duration'                'recording_duration';
        'stim_file'                   'probe_file';
        'subj'                        'subj';
        'task'                        'task_instruction';
        'task_marker'                 'task_marker';
        'trial_started'               'trial_started';
        'conditions_file'             'block';  
        };

    old_columns = renamed_column(:, 1);
    new_columns = renamed_column(:, 2);

    subj_log_clean = renamevars(subj_log_clean, old_columns, new_columns);

    % Extract block number
    subj_log_clean.block = cellfun(@(x) str2double(regexp(x, '_([0-9]+)\.xlsx$', 'tokens', 'once')), subj_log_clean.block);

     % Sanity Check: Plausible blocks extracted
    if all(ismember( subj_log_clean.block, [1, 2, 3, 4, 5, 6, 7, 8]))
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'wrong_block_number';
    end

    % Sanity Check: Equal db for unaltered and altered probes
    if height(subj_log_clean) == n_trials
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'log_dimensions';
    end

    % Get start latency of go stimulus relative to trial start
    subj_log_clean.go_stim_started_trial_start = subj_log_clean.go_stim_started - subj_log_clean.trial_started;

    % Extract filename from full path for each vocal recording in order to match it with the table produced in Praat
    subj_log_clean.recording_file = cellfun(@(x) strrep(regexp(x, 'recording_mic_.*(?=\.wav)', 'match', 'once'), '.', '_'), subj_log_clean.("recording_path"), 'UniformOutput', false);

    % Get probe properties file
    subj_probe_properties = readtable(fullfile(INPATH, ['stimuli\' subj '\' subj '_probe_properties.xlsx']));

    % Rename columns
    renamed_column = {
        'f0_tab_normal'         'probe_unaltered_f0';
        'f0_tab_pitched'        'probe_altered_f0';
        'db_tab_normal'         'probe_unaltered_db';
        'db_tab_pitched'        'probe_altered_db';
        'filename_tab'          'probe_file';
        'change_attempts'       'probe_change_attempts';
        'loudness'              'probe_loudness_attenuation';
        };

    old_columns = renamed_column(:, 1);
    new_columns = renamed_column(:, 2);

    subj_probe_properties = renamevars(subj_probe_properties, old_columns, new_columns);

    % Sanity Check: Correct dimensions
    if height(subj_probe_properties) == 1 && width(subj_probe_properties) == 7
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'probe_properties_dimensions';
    end
    % Sanity Check: Equal db for unaltered and altered probes
    if round(subj_probe_properties.probe_unaltered_db) == round(subj_probe_properties.probe_altered_db)
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'probe_properties_db';
    end

    % Get F0 & vocal onset Table
    subj_f0_rt = readtable(fullfile(INPATH_PRAAT, [subj '_f0_rt_table.csv']), Delimiter=',');
    subj_f0_rt = standardizeMissing(subj_f0_rt,9999);

    % Rename columns
    renamed_column = {
        'subj_tab'                      'subj';
        'f0_tab'                        'recording_f0';
        'rt_tab'                        'recording_vot';
        'duration_file_tab'             'recoding_duration';
        'condition_tab'                 'task';
        'min_intensity_tab'             'recording_min_intensity';
        'max_intensity_tab'             'recording_max_intensity';
        'ratio_intensity_tab'           'recording_ratio_intensity';
        'vocal_response_tab'            'vocal_response';
        'n_intervals_tab'               'recording_n_intervals';
        'duration_vocal_tab'            'recording_duration_vocal';
        'onset_vocal_tab'               'recording_onset_vocal';
        'offset_vocal_tab'              'recording_offset_vocal';
        'filename_tab'                  'recording_file';
        };

    old_columns = renamed_column(:, 1);
    new_columns = renamed_column(:, 2);

    subj_f0_rt = renamevars(subj_f0_rt, old_columns, new_columns);

    % Sanity Check: Correct dimensions
    if height(subj_f0_rt) == n_trials && width(subj_f0_rt) == 14
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'f0_rt_dimensions';
    end

    % Concatinate log file, probe properties file and F0 & vocal onset table
    subj_full = join(subj_f0_rt, subj_log_clean, 'Keys', {'subj', 'recording_file'});
    % Exclude 'filename_tab' column and repeat the remaining columns for each row of subj_full
    subj_full = [subj_full, repmat(subj_probe_properties(:, setdiff(1:width(subj_probe_properties), find(strcmp(subj_probe_properties.Properties.VariableNames, 'probe_file')))), height(subj_full), 1)];

    % Get reaction time relative to go stimulus onset
    subj_full.recording_vot = subj_full.recording_onset_vocal - (subj_full.go_stim_started_trial_start- subj_full.mic_started);

    % Exlcude trials with outlier vocal onset time
    % Get threshold
    mean_rt = mean(subj_full.recording_vot, 'omitnan');
    std_rt = std(subj_full.recording_vot, 'omitnan');
    threshold_upper = mean_rt + 3 * std_rt;
    threshold_lower = mean_rt - 3 * std_rt;
    % Add a column indicating F0 outliers
    subj_full.recording_vot_outlier = (subj_full.recording_vot > threshold_upper) | (subj_full.recording_vot < threshold_lower);
    % Excluded trials based on F0 outliers
    excluded_trials_rt_outliers = (subj_full.recording_vot > threshold_upper) | (subj_full.recording_vot < threshold_lower);

    % Exlcude trials with outlier F0
    % Get threshold
    mean_f0 = mean(subj_full.recording_f0, 'omitnan');
    std_f0 = std(subj_full.recording_f0, 'omitnan');
    threshold_upper = mean_f0 + 3 * std_f0;
    threshold_lower = mean_f0 - 3 * std_f0;
    % Add a column indicating F0 outliers
    subj_full.recording_f0_outlier = (subj_full.recording_f0 > threshold_upper) | (subj_full.recording_f0 < threshold_lower);
    % Excluded trials based on F0 outliers
    excluded_trials_f0_outliers = (subj_full.recording_f0 > threshold_upper) | (subj_full.recording_f0 < threshold_lower);

    % Exlude trial with incorrect responses (e.g. vocalizing in 'Passive' trials)
    % Add a column indicating (in)correct responses based on classification in Praat (see vocal_analysis.praat)
    pas_correct = strcmp(subj_full.task, 'pas') & subj_full.vocal_response == 0;
    act_correct = strcmp(subj_full.task, 'act') & subj_full.vocal_response == 1;
    subj_full.correct_vocal_response(pas_correct | act_correct) = 1;
    % Excluded trials based on (in)correct vocal responses
    excluded_trials_responses = ~subj_full.correct_vocal_response;

    % Combine excluded trials
    excluded_trials_beh = excluded_trials_rt_outliers | excluded_trials_responses;
    excluded_trials_beh = excluded_trials_beh | excluded_trials_f0_outliers;

    % Only included needed columns
    subj_full_clean = subj_full(:, {'recording_path', ...
        'mic_started', 'mic_stopped', 'probe', ...
        'probe_onset', 'probe_onset_cat', 'probe_type', ...
        'subj', 'task_instruction', 'task_marker', 'recording_f0', 'recording_vot', ...
        'recording_min_intensity', 'recording_max_intensity', 'vocal_response', ...
        'recording_duration_vocal', 'recording_onset_vocal', 'recording_file', ...
        'probe_unaltered_f0', 'probe_altered_f0', 'probe_loudness_attenuation', ...
        'probe_change_attempts', 'correct_vocal_response', 'recording_f0_outlier', 'task', ...
        'trial_started', 'go_stim_started', 'go_stim_started_trial_start', 'block'});

    writetable(subj_full_clean,fullfile(OUTPATH, [subj '_beh_preprocessed.xlsx']))

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

check_done = 'tid_psam_beh_preprocessing_2_DONE'

delete(wb)
