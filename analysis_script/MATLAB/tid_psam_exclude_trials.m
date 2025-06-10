% tid_psam_exclude_trials.m
%
% Excludes invaled trials from behavioural, SVM and ERP data and excludes
%   subjects with too many excluded trials.
%
% The proccessing steps include
%
    % Loading behavioural and EEG-ERP data
    % Matches invalid trials identified in EEG-ERP data and behavioural data
    % Removes tagged trials from both EEG-ERP data and behavioural data
    % Performs multiple transformations using on the included trials
    % Checks number of valid trials per condition
    % Loading behavioural and EEG-SVM data
    % Matches invalid trials identified in EEG-SVM data and behavioural data
    % Removes tagged trials only from the EEG-SVM data
    % Excludes subjects if the number of trials in below threshold for
    %   either condition
%
% Saves data
%
% Concatinates all individual behavioural data sets to one large one and stores it.
%
% Note. Since multiple data types (EEG data and vocal data) and different EEG preprocessing pipelines
% (ERP-specific preprocessing and single-trial-specific preprocessing) will be used (see below), the flagged
% trials will need to be merged in a systematic way. Because the ERP analysis and the behavioural analysis are linked,
% in the sense that the vocal responses made are thought to influence the ERP analysis, the flagged trials will be merged.
% In other words, if a trial is flagged based on the behavioural data, it will also be removed from the ERP analysis.
% Likewise, if a trial is flagged based on the ERP-specific preprocessing, it will also be removed from the behavioural analysis.
% Thus, the ERP analysis and the behavioural analysis will include the exact same trials. Since the classification analysis
% does not rely on any quantities associated with the vocal responses themselves, a different approach will be used.
% A trial will be excluded from the classification analysis if it is flagged based on the single-trial-specific preprocessing
% and/or the behavioural preprocessing. However, the single-trial-specific preprocessing will not affect the trials included in
% the behavioural analysis. Therefore, if a trial is flagged based on the behavioural preprocessing, it will also be excluded
% from the classification analysis. Conversely, if a trial is (only) flagged based on the single-trial-specific preprocessing,
% it will not be excluded from the behavioural analysis. This approach is used to maximize the number of usable trials in each analysis.
% If the flagged trials were matched across all preprocessing pipelines, this could lead to trials being unnecessarily excluded.
% For example, as the epochs for the EEG analyses will be extracted based on different events, an artefact could be present
% at one time-point but not another. Thus, matching flagged trials could result in the loss of valid data.
% Finally, participants will only be included if more than 30 trials per condition remain for the ERP analysis,
% and more than 100 trials per condition remain for the classification analysis. Participants who do not meet these
% criteria will be excluded from all analyses. Furthermore, if data from any dataset is missing, the
% participant will also be excluded from all analyses. Thus, while individual trials will be matched between the
% ERP and behavioural analyses—but not with the classification analysis—the same set of participants will be included across all analyses.
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
INPATH_ERP = fullfile(MAINPATH, 'data\processed_data\erp_preprocessed\');
INPATH_SVM = fullfile(MAINPATH, 'data\processed_data\svm_preprocessed\');
INPATH_BEH = fullfile(MAINPATH, 'data\processed_data\beh_preprocessed_2\');
INPATH_QUEST = fullfile(MAINPATH, 'data\questionnaire_data\');
OUTPATH_ERP = fullfile(MAINPATH, 'data\processed_data\erp_preprocessed_clean\');
OUTPATH_SVM = fullfile(MAINPATH, 'data\processed_data\svm_preprocessed_clean\');
OUTPATH_BEH = fullfile(MAINPATH, 'data\processed_data\beh_preprocessed_clean\');
OUTPATH_QUEST = fullfile(MAINPATH, 'data\questionnaire_data_clean\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\exclude_trials\');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_ERP, INPATH_SVM, INPATH_BEH, OUTPATH_ERP, OUTPATH_SVM,OUTPATH_BEH,OUTPATH_QUEST ,OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH_ERP)
tid_psam_clean_up_folder_TD(OUTPATH_SVM)
tid_psam_clean_up_folder_TD(OUTPATH_BEH)
tid_psam_clean_up_folder_TD(OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH_QUEST)

% Variables to edit
n_trials = 960;
n_trials_no_probe = 480;
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};
N_TRIALS_THRESH = 25; % More than 25 trials trials per condition required for ERP analysis (excluded)
N_TRIALS_NO_PROBE_THRESH = 100; % More than 100 trials trials required per condition for SVM analysis (excluded)
N_TRIALS_BLOCK_THRESH = 60; % More than 60 trials trials required per block (only marked, not excluded)
SUBJ_TO_EXCLUDE = {};

% Get directory content
dircont_subj_erp = dir(fullfile(INPATH_ERP, 'sub-*.set'));
dircont_subj_svm = dir(fullfile(INPATH_SVM, 'sub-*.set'));
dircont_subj_beh = dir(fullfile(INPATH_BEH, 'sub-*.xlsx'));

% Sanity Check: Same length of directory contents
if length(dircont_subj_erp) == length(dircont_subj_beh) && length(dircont_subj_beh) == length(dircont_subj_svm)
else
    error('Different number of files')
end

% Initialize sanity check variables
excluded_subj = {};
marked_subj = {};
protocol = {};

% Initialize empty table used to get data from all subjects
all_subj_beh_clean = table;

% Setup progress bar
wb = waitbar(0,'starting tid_psam_exclude_trials.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj_erp)

    % Get current ID & Sanity Check: Same ID for ERP data and behavioural data
    subj_erp = dircont_subj_erp(subj_idx).name;
    subj_erp = regexp(subj_erp, 'sub-\d+', 'match', 'once');
    subj_svm = dircont_subj_svm(subj_idx).name;
    subj_svm = regexp(subj_svm, 'sub-\d+', 'match', 'once');
    subj_beh = dircont_subj_beh(subj_idx).name;
    subj_beh = regexp(subj_beh, 'sub-\d+', 'match', 'once');

    if all(strcmp(subj_erp, {subj_beh, subj_svm}))
        subj = subj_erp;
    else
        error('Subject IDs do not match')
    end

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj_erp),wb, [subj ' tid_psam_exclude_trials.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load behavioural data
    beh = readtable(fullfile(INPATH_BEH, [subj '_beh_preprocessed.xlsx']));
    load(fullfile(INPATH_BEH, [subj '_excluded_trials_beh_preprocessing.mat'])) % Loads variable excluded_trials_beh

    % Load EEG data (ERP preprocessing)
    EEG = pop_loadset('filename',[subj '_erp_preprocessed.set'],'filepath',INPATH_ERP);

    % Sanity Check: Same number of trials (=960)
    if length(excluded_trials_beh) == n_trials && size(EEG.data,3)  == n_trials && height(beh)  == n_trials
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'number_of_trials';
        error('Invalid number of trials')
    end

    % Merge to be excluded trials
    excluded_trials_all = excluded_trials_beh' | EEG.reject.rejglobal;

    % Remove to be excluded trials from vocal data
    beh_clean = beh(~excluded_trials_all',:);

    % Remove to be excluded trials from ERP data
    EEG.reject.rejglobal = excluded_trials_all;
    EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);
    EEG.setname = [subj '_erp_clean'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Apply transformations
    % Z-transform F0 values
    f0_mean = mean(beh_clean.recording_f0, 'omitmissing');
    f0_sd = std(beh_clean.recording_f0, 'omitmissing');
    beh_clean.recording_f0_z = (beh_clean.recording_f0 - f0_mean) ./ f0_sd;
    % Z-transform probe F0 values (relative to vocal responses)
    beh_clean.probe_f0_unaltered_z = (beh_clean.probe_unaltered_f0 - f0_mean) ./ f0_sd; % Unaltered Probe
    beh_clean.probe_f0_altered_z = (beh_clean.probe_altered_f0 - f0_mean) ./ f0_sd; % Altered Probe
    % Z-transform vocal onset times
    vot_mean = mean(beh_clean.recording_vot, 'omitmissing');
    vot_sd = std(beh_clean.recording_vot, 'omitmissing');
    beh_clean.recording_vot_z = (beh_clean.recording_vot - vot_mean) ./ vot_sd;

    % Sanity Check: Number of trials in each condition (ERP)
    for cond = 1:length(EVENTS)
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off'); % Create dataset for each condition
        EEG.setname = [subj '_' EVENTS{cond}];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        n_trials_condition_erp(cond).condition = EVENTS{cond};
        n_trials_condition_erp(cond).n_trials = size(EEG.data,3);
    end

    % Sanity Check: Number of trials in each condition (BEH)
    validRows = ~strcmp(beh_clean.probe_onset_cat, 'None') & ...
        ~strcmp(beh_clean.probe_type, 'None');
    filteredRows = beh_clean(validRows, :);
    % Group by the three factors
    [G, condTab, onsetCat, probeType] = findgroups( ...
        filteredRows.task, ...
        filteredRows.probe_onset_cat, ...
        filteredRows.probe_type);
    % Count number of trials per group
    counts = splitapply(@numel, filteredRows.probe_type, G);
    % Map probe_type to 'alt' / 'unalt'
    probeMap = containers.Map({'Altered', 'Unaltered'}, {'alt', 'unalt'});
    probeStr = values(probeMap, probeType);
    % Build condition names using cellfun to match ERP data
    conditionNames = strcat(condTab, '_', lower(onsetCat), '_', probeStr(:));
    n_trials_condition_beh = struct('condition', conditionNames, 'n_trials', num2cell(counts))';

    % Sanity Check: Matching number of trials per condition in ERP data and behavioural data
    % Extract condition names and trial counts
    behConditions = {n_trials_condition_beh.condition};
    behCounts = [n_trials_condition_beh.n_trials];
    erpConditions = {n_trials_condition_erp.condition};
    erpCounts = [n_trials_condition_erp.n_trials];
    % Find matching conditions between the two structs
    [commonConds, behIdx, erpIdx] = intersect(behConditions, erpConditions, 'stable');
    if all(behCounts(behIdx) == erpCounts(erpIdx),'all')
    else
        error('Number of trials per conndition do not match')
    end

    % Sanity Check: Enough trials per block
    for block_num = 1:8
        if sum([beh_clean.block] == block_num) <= N_TRIALS_BLOCK_THRESH
            marked_subj{end+1,1} = subj;
            marked_subj{end,2} = ['number_of_trials_block_' num2str(block_num)];
        end
    end

    % Exclude subjects if not enough trials per condition are left
    for cond_num = 1:12
        if n_trials_condition_erp(cond_num).n_trials <= N_TRIALS_THRESH
            excluded_subj{end+1,1} = subj;
            excluded_subj{end,2} = ['number_of_trials_' n_trials_condition_erp(cond_num).condition '_' num2str(n_trials_condition_erp(cond_num).n_trials)];
        end
    end

    % Exclude subject if manual exclusion was defined
    if any(strcmp(subj, SUBJ_TO_EXCLUDE))
        excluded_subj{end+1,1} = subj;
        excluded_subj{end,2} = 'manually_excluded';
    end

    % Load SVM data and merge to be excluded trials with behavioural data
    % Note. ERP preprocessed excluded and SVM preprocessed excluded trials are not merged
    % Update excluded subjects
    if ~any(strcmp(excluded_subj, subj), 'all') % Only load data for subjects, which are not already excluded
        % Load EEG data (SVM preprocessing)
        EEG = pop_loadset('filename',[subj '_svm_preprocessed.set'],'filepath',INPATH_SVM);

        % Extract no-probe trials from behavioural data
        excluded_trials_beh_no_probe = excluded_trials_beh(strcmp(beh.probe, 'No'));

        % Sanity Check: Correct number of no-probe trials (=480)
        if length(excluded_trials_beh_no_probe) == n_trials_no_probe && size(EEG.data,3) == n_trials_no_probe
        else
            marked_subj{end+1,1} = subj;
            marked_subj{end,2} = 'number_of_no_probe_trials';
        end

        % Merge to be excluded trials
        excluded_trials_all_no_probe = excluded_trials_beh_no_probe' | EEG.reject.rejglobal;

        % Remove to be excluded trials from SVM data
        EEG.reject.rejglobal = excluded_trials_all_no_probe;
        EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);
        EEG.setname = [subj '_svm_clean'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        % Exclude subjects if not enough trials per condition are left
        if sum(strcmp({EEG.event.type}, 'go_act')) > N_TRIALS_NO_PROBE_THRESH && sum(strcmp({EEG.event.type}, 'go_pas')) > N_TRIALS_NO_PROBE_THRESH
        else
            excluded_subj{end+1,1} = subj;
            excluded_subj{end,2} = 'number_of_no_probe_trials';
        end
    end

    % Store and save datasets for non-excluded subjects
    if ~any(strcmp(excluded_subj, subj), 'all')
        % EEG data (ERP preprocessing)
        EEG = ALLEEG(1);
        CURRENTSET = 1;
        EEG = pop_saveset(EEG, 'filename',[subj '_erp_preprocessed_clean.set'],'filepath', OUTPATH_ERP);

        % EEG data (SVM preprocessing)
        EEG = ALLEEG(14);
        CURRENTSET = 14;
        EEG = pop_saveset(EEG, 'filename',[subj '_svm_preprocessed_clean.set'],'filepath', OUTPATH_SVM);

        % Behavioural data
        writetable(beh_clean,fullfile(OUTPATH_BEH, [subj '_beh_preprocessed_clean.xlsx']))

        % Get behavioral data for all subjects
        all_subj_beh_clean = vertcat(all_subj_beh_clean, beh_clean);

        % Save to-be-exlcuded trials for use in tid_psam_hilbert_preparation
        save(fullfile(OUTPATH, [subj '_excluded_no_probe_trials.mat']),'excluded_trials_all_no_probe')
    end

    % Update Protocol
    subj_time = toc;
    protocol{subj_idx,1} = subj;
    protocol{subj_idx,2} = subj_time;
    if any(strcmp(excluded_subj, subj), 'all')
        protocol{subj_idx,3} = 'EXCLUDED';
    elseif any(strcmp(marked_subj, subj), 'all') && ~any(strcmp(excluded_subj, subj), 'all')
        protocol{subj_idx,3} = 'MARKED';
    else
        protocol{subj_idx,3} = 'OK';
    end

end

% Remove excluded subjects from questionnaire data
fal_data = readtable(fullfile(INPATH_QUEST, 'fal_data.xlsx')); % Load data
nasatlx_data = readtable(fullfile(INPATH_QUEST, 'nasatlx_data.xlsx')); % Load data
sam_data = readtable(fullfile(INPATH_QUEST, 'sam_data.xlsx')); % Load data

if ~isempty(excluded_subj)
    excluded_subj_unique = unique({excluded_subj{:,1}}); % Get IDs of excluded subjects

    fal_data_clean = fal_data(~strcmp(fal_data.subj, excluded_subj_unique' ),:); % Remove excluded subjects
    nasatlx_data_clean = nasatlx_data(~strcmp(nasatlx_data.subj, excluded_subj_unique' ),:); % Remove excluded subjects
    sam_data_clean = sam_data(~strcmp(sam_data.subj, excluded_subj_unique' ),:); % Remove excluded subjects
else
    fal_data_clean = fal_data;
    nasatlx_data_clean = nasatlx_data;
    sam_data_clean = sam_data;
end

writetable(fal_data_clean,fullfile(OUTPATH_QUEST, 'fal_data_clean.xlsx')); % Store data
writetable(nasatlx_data_clean,fullfile(OUTPATH_QUEST, 'nasatlx_data_clean.xlsx')); % Store data
writetable(sam_data_clean,fullfile(OUTPATH_QUEST, 'sam_data_clean.xlsx')); % Store data

% Sanity Check: Correct number of saved datasets and rows in questionnaire data
n_saved_erp = length(dir(fullfile(OUTPATH_ERP, 'sub-*.set')));
n_saved_svm = length(dir(fullfile(OUTPATH_SVM, 'sub-*.set')));
n_saved_beh = length(dir(fullfile(OUTPATH_BEH, 'sub-*.xlsx')));
n_rows_fal = size(fal_data_clean,1);
n_rows_nasatlx = size(nasatlx_data_clean,1);
n_rows_sam = size(sam_data_clean,1);
if isempty(excluded_subj)
    n_excluded = 0;
else
    n_excluded = length(unique({excluded_subj{:,1}}));
end
n_subjects_all = length(dircont_subj_erp);

if n_saved_erp == n_saved_beh && n_saved_erp == n_saved_svm && n_saved_erp == (n_subjects_all - n_excluded) && n_saved_erp == n_rows_fal && n_rows_fal == n_rows_nasatlx && n_rows_fal == n_rows_sam
else
    error('Number of saved files incorrect')
end

% Sanity Check: Same subjects in behavioural, EEG (ERP), EEG (SVM) and questionnaire data
subj_list_erp = regexp({dir(fullfile(OUTPATH_ERP, 'sub-*.set')).name}, 'sub-\d+', 'match', 'once')';
subj_list_svm = regexp({dir(fullfile(OUTPATH_SVM, 'sub-*.set')).name}, 'sub-\d+', 'match', 'once')';
subj_list_beh = regexp({dir(fullfile(OUTPATH_BEH, 'sub-*.xlsx')).name}, 'sub-\d+', 'match', 'once')'; % Individual data
subj_list_beh2 = unique(all_subj_beh_clean.subj); % Merged data containing all subjects
subj_list_fal = fal_data_clean.subj;
subj_list_nasatlx = nasatlx_data_clean.subj;
subj_list_sam = sam_data_clean.subj;

subj_list_all = {subj_list_erp, subj_list_svm, subj_list_beh, subj_list_beh2, subj_list_fal, subj_list_nasatlx, subj_list_sam};
if ~all(cellfun(@(x) isequal(x, subj_list_all{1}), subj_list_all))
    error('Different subjects included in data')
end

% Store behavioral data from all subjects
writetable(all_subj_beh_clean,fullfile(OUTPATH_BEH, 'all_subj_beh_preprocessed_clean.xlsx'))

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'tid_psam_exclude_trials_protocol.xlsx'))

if ~isempty(excluded_subj)
    excluded_subj = cell2table(excluded_subj, 'VariableNames',{'subj','issue'})
    writetable(excluded_subj,fullfile(OUTPATH, 'tid_psam_exclude_trials_excluded_subj.xlsx'))
end

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_exclude_trials_marked_subj.xlsx'))
end

check_done = 'tid_psam_exclude_trials_DONE'

delete(wb)
