% tid_psam_exclude_trials.m
%
% Excludes invalied trials from behavioural, SVM and ERP data and excludes
%   subjects with too many excluded trials.
%
% The proccessing steps include:
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
    % Stores datasets
%
% Concatinates all individual behavioural data sets to one large one and stores it.
%
% Note. The ERP analysis and the behavioural analysis include the exact
%   same trials and participants after the described procedures. The SVM
%   analysis might include trials that are excluded for the ERP analysis and
%   vice versa. Further, the behavioural analysis might include trials which
%   are excluded from the SVM analysis. To sum up, the trials for the ERP
%   analysis and the behavoural analysis are perfectly matched. The SVM
%   analysis is somewhat independent, potentially including different trials.
%   The rationale behind this approach is to maximize the usable trials
%   for each analysis. Since the ERP analysis relates to quantities like
%   F0 (differences), it makes sense to match this analysis with the
%   behavioural one. On the other hand, the SVM analysis is not related to
%   such quantities, making a more liberal approach possible. 
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
OUTPATH_ERP = fullfile(MAINPATH, 'data\processed_data\erp_preprocessed_clean\');
OUTPATH_SVM = fullfile(MAINPATH, 'data\processed_data\svm_preprocessed_clean\');
OUTPATH_BEH = fullfile(MAINPATH, 'data\processed_data\beh_preprocessed_clean\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\exclude_trials\');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_ERP, INPATH_SVM, INPATH_BEH, OUTPATH_ERP, OUTPATH_SVM,OUTPATH_BEH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH_ERP)
tid_psam_clean_up_folder_TD(OUTPATH_SVM)
tid_psam_clean_up_folder_TD(OUTPATH_BEH)
tid_psam_clean_up_folder_TD(OUTPATH)


% Variables to edit
n_trials = 960;
n_trials_no_probe = 480;
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};
N_TRIALS_THRESH = 30; % More than 30 trials trials per condition required for ERP analysis (excluded)
N_TRIALS_NO_PROBE_THRESH = 100; % More than 100 trials trials required per condition for SVM analysis (excluded)
N_TRIALS_BLOCK_THRESH = 60; % More than 60 trials trials required per block (only marked, not excluded)

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
    if any([n_trials_condition_erp.n_trials] <= N_TRIALS_THRESH)
        excluded_subj{end+1,1} = subj;
        excluded_subj{end,2} = 'number_of_trials';
    end

    % Load SVM data and merge to be excluded trials with behavioural data
    % Note. ERP preprocessed excluded and SVM preprocessed excluded trials are not merged
    % Update excluded subjects
    if ~any(strcmp(excluded_subj, subj), 'all') % Only load data for subjects, which are not already excluded
        % Load EEG data (SVM preprocessing)
        EEG = pop_loadset('filename',[subj '_svm_preprocessed.set'],'filepath',INPATH_SVM);

        % Extract no-probe trials from behavioural data
        excluded_trials_beh_no_probe = excluded_trials_beh(strcmp(beh.probe, 'No'));

        % % a = find(strcmp(beh.probe, 'No'));
        % % b = zeros(1,960);
        % % b(a) = EEG.reject.rejglobal;
        % % d = b(strcmp(beh.probe, 'No'));

        % Sanity Check: Correct number of no-probe trials (=480)
        if length(excluded_trials_beh_no_probe) == n_trials_no_probe && size(EEG.data,3) == n_trials_no_probe
        else
            marked_subj{end+1,1} = subj;
            marked_subj{end,2} = 'number_of_no_probe_trials';
        end

        % Merge to be excluded trials
        excluded_trials_all_no_probe = excluded_trials_beh_no_probe' | EEG.reject.rejglobal;

        % Remove to be excluded trials from SVM data
        EEG.reject.rejglobal = excluded_trials_beh_no_probe;
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

    % Sanity Check: Correct number of saved datasets
    n_saved_erp = length(dir(fullfile(OUTPATH_ERP, 'sub-*.set')));
    n_saved_beh = length(dir(fullfile(OUTPATH_BEH, 'sub-*.xlsx')));
    n_excluded = size(excluded_subj,1);
    n_subjects_all = length(dircont_subj_erp);

    if n_saved_erp == n_saved_beh && n_saved_erp == (n_subjects_all - n_excluded)
    else
        error('Number of saved files incorrect')
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
