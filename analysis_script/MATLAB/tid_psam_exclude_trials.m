% tid_psam_exclude_trials.m
%
% Excludes invalied trials from behavioural and ERP data and excludes
%   subjects with too many excluded trials.
%
% The proccessing steps include:
    % Loading behavioural and ERP data
    % Matches invalid trials identified in either dataset
    % Performs multiple transformations using on the included trials
    % Checks number of valid trials per condition
    % Excludes subjects if the number of trials in below threshold for
    %   either condition
    % Stores datasets 
%
% Concatinates all individual data sets to one large one and stores it.
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
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};
N_TRIALS_THRESH = 30; % Number of minimal required trials per condition (excluded)
N_TRIALS_BLOCK_THRESH = 60; % Number of minimal required trials per block (only marked, not excluded)

% Get directory content
dircont_subj_erp = dir(fullfile(INPATH_ERP, 'sub-*.set'));
dircont_subj_beh = dir(fullfile(INPATH_BEH, 'sub-*.xlsx'));

% Sanity Check: Same length of directory contents
if length(dircont_subj_erp) == length(dircont_subj_beh)
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
    subj_beh = dircont_subj_erp(subj_idx).name;
    subj_beh = regexp(subj_beh, 'sub-\d+', 'match', 'once');

    if strcmp(subj_erp,subj_beh)
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

    % Load EEG data
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
    EEG.setname = [subj '_clean'];
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
            'deleteevents','off','deleteepochs','on','invertepochs','off'); %create dataset for each condition
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

    % Sanity Check: Matching number of trials per condition in EPR data and behavioural data
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
    if sum([beh_clean.block] == block_num) < N_TRIALS_BLOCK_THRESH 
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = ['number_of_trials_block_' num2str(block_num)];
    end
    end

    % Exclude subjects if not enough trials per condition are left
    if any([n_trials_condition_beh.n_trials] < N_TRIALS_THRESH)
        excluded_subj{end+1,1} = subj;
        excluded_subj{end,2} = 'number_of_trials';
    end

    % Store and save datasets for non-excluded subjects
    if ~any(strcmp(excluded_subj, subj), 'all')
        % EEG data
        EEG = ALLEEG(1);
        CURRENTSET = 1;
        EEG = pop_saveset(EEG, 'filename',[subj '_erp_preprocessed_clean.set'],'filepath', OUTPATH_ERP);

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
