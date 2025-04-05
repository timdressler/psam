% tid_psam_create_conditions_file.m
%
% Creates BIDS-inspired folder structure and conditions file.
% Has to be executed BEFORE tid_psam_create_conditions_file.m,
%   tid_psam_stimuli_recording_adapted.py and tid_psam_prepare_stimuli.praat.
%
% The pseudorandomization ensures that task conditions are shuffled while
% maintaining constraints on consecutive repetitions.
% Note. The number of trials (n_trials) has to be dividable by 16 due to
% the used miniblock system.
%
% Randomization Process
%   Miniblock Setup: Each iteration of the loop creates a miniblock of 16 trials, with 8 "Yes" (probe) trials and 8 "No" (no probe) trials
%   Shuffling with Constraints:
%       A temporary copy of the task list is made
%       The script selects a random valid trial from the remaining list, ensuring that no task type ("Yes" or "No") repeats more than max_repeats (4 times) in a row
%       If no valid choice is available, the process resets and retries to generate a valid sequence
%       Task Assignment: Once a valid randomized sequence is generated:
%           Half of the "Yes" and "No" trials are assigned to "Active" and the other half to "Passive" conditions
%           "Yes" (probe) trials are further divided into "Early" and "Late" onset probes, with equal distribution of "Normal" and "Pitch" probe types
%       Finalizing the Block:
%           Probe onset times, durations, and intensities are assigned
%           The subject ID is added, and the miniblock is stored
%           The process repeats for all required iterations to build the full experiment trial list.
%
% Tim Dressler, 07.03.2025

clear
close all
clc

% Get subject ID
prompt = {'Enter subject number:'};
dlgtitle = 'Subject Number Input';
dims = [1 35];

answer = inputdlg(prompt, dlgtitle, dims);

% Format subject ID
if ~isempty(answer)
    num = str2double(answer{1}); % Convert input to a number
    formattedNum = sprintf('%02d', num); % Ensure two-digit format
    subj = ['sub-' formattedNum]; % Construct final subject ID
end

% Set up paths
SCRIPTPATH = cd;
if regexp(SCRIPTPATH, regexptranslate('wildcard','*psam\experiment_script')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, '\experiment_script');
STIMULIPATH = fullfile(MAINPATH, ['data/BIDS/stimuli/' subj]);
STIMULIPATH_Normal = fullfile(STIMULIPATH,'/all_normal');
STIMULIPATH_Pitch = fullfile(STIMULIPATH,'/all_pitch');
STIMULIPATH_Raw = fullfile(STIMULIPATH,'/all_raw');
SUBJPATH = fullfile(MAINPATH,['data/BIDS/' subj]);
EEGPATH = fullfile(SUBJPATH, '/eeg');
BEHAVIORALPATH = fullfile(SUBJPATH, '/beh');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_id_TD(STIMULIPATH, SUBJPATH)
tid_psam_check_folder_TD(MAINPATH,STIMULIPATH_Normal,STIMULIPATH_Pitch,STIMULIPATH_Raw, EEGPATH, BEHAVIORALPATH)

% Set up experiment paramaters
early_onset = 2.8; % relative to trial onset
late_onset = 2.9; % relative to trial onset
n_trials = 960; % has to be dividable by 16
max_repeats = 4; % how often probe (or no probe) trials can be repeated
active_recording_duration = 1; % how long the audio recording is
passive_recording_duration = 0.1; % how long the audio recording is

% Checks whether the number of trials is diviable by 16
if ~mod(n_trials, 16) == 0
    error('trial number has to be divideable by 16')
else
    disp('Trial number OK')
end

% Get needed number of miniblocks
num_iterations = n_trials/16;

% Set up filenames for probes
path_normal_probe = fullfile(STIMULIPATH, [subj '_normal_probe.wav']);
path_pitch_probe = fullfile(STIMULIPATH, [subj '_pitch_probe.wav']);

% Start randomization process
all_trials = {};
for iter = 1:num_iterations
    % Set up miniblock
    miniblock_task = {'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', ...
        'No', 'No', 'No', 'No', 'No', 'No', 'No', 'No'};
    success = false;
    while ~success
        try
            % Shuffle while maintaining max repetition constraint
            miniblock_task_rando = {};
            temp_task = miniblock_task; % copy to avoid modifying original

            while ~isempty(temp_task)
                valid_indices = [];

                % Find indices that won't break max_repeats constraint
                for i = 1:length(temp_task)
                    if length(miniblock_task_rando) < max_repeats || ~all(strcmp(miniblock_task_rando(end-(max_repeats-1):end), temp_task{i}))
                        valid_indices(end+1) = i;
                    end
                end

                % If no valid indices are found, trigger an error
                if isempty(valid_indices)
                    error('Could not generate sequence with the given constraints. Retrying...');
                end

                % Randomly pick from valid indices
                random_index = valid_indices(randi(length(valid_indices)));
                miniblock_task_rando{end+1} = temp_task{random_index};
                temp_task(random_index) = [];
            end

            success = true; % if successful, exit loop
        catch
            % If an error occurs, retry
            %%warning('Sequence generation failed. Retrying...');
        end
    end

    % Get indices of probe and no proeb trials
    miniblock_task_rando = miniblock_task_rando';
    yes_idx = find(strcmp(miniblock_task_rando, 'Yes'));
    no_idx = find(strcmp(miniblock_task_rando, 'No'));

    % Assign 50% of each type (probe + no probe) to Passive and Active
    yes_active = yes_idx(randperm(length(yes_idx), 4));
    yes_passive = setdiff(yes_idx, yes_active);

    no_active = no_idx(randperm(length(no_idx), 4));
    no_passive = setdiff(no_idx, no_active);

    % Create task column
    miniblock_task_rando(yes_active,2) = {'Active'};
    miniblock_task_rando(yes_active,11) = {21};
    miniblock_task_rando(no_active,2) = {'Active'};
    miniblock_task_rando(no_active,11) = {21};
    miniblock_task_rando(yes_passive,2) = {'Passive'};
    miniblock_task_rando(yes_passive,11) = {22};
    miniblock_task_rando(no_passive,2) = {'Passive'};
    miniblock_task_rando(no_passive,11) = {22};

    % Assign recording duration
    miniblock_task_rando(yes_active,12) = {active_recording_duration};
    miniblock_task_rando(no_active,12) = {active_recording_duration};
    miniblock_task_rando(yes_passive,12) = {passive_recording_duration};
    miniblock_task_rando(no_passive,12) = {passive_recording_duration};

    % Assign probe_type, probe_onset, probe_file and probe_onset_cat
    probe_type_active_idx = yes_active(randperm(length(yes_active)));
    miniblock_task_rando(probe_type_active_idx(1),3) = {'Early'};
    miniblock_task_rando(probe_type_active_idx(1),5) = {early_onset};
    miniblock_task_rando(probe_type_active_idx(1),4) = {'Normal'};
    miniblock_task_rando(probe_type_active_idx(1),7) = {path_normal_probe};
    miniblock_task_rando(probe_type_active_idx(1),10) = {31};
    miniblock_task_rando(probe_type_active_idx(2),3) = {'Early'};
    miniblock_task_rando(probe_type_active_idx(2),5) = {early_onset};
    miniblock_task_rando(probe_type_active_idx(2),4) = {'Pitch'};
    miniblock_task_rando(probe_type_active_idx(2),7) = {path_pitch_probe};
    miniblock_task_rando(probe_type_active_idx(2),10) = {32};
    miniblock_task_rando(probe_type_active_idx(3),3) = {'Late'};
    miniblock_task_rando(probe_type_active_idx(3),5) = {late_onset};
    miniblock_task_rando(probe_type_active_idx(3),4) = {'Normal'};
    miniblock_task_rando(probe_type_active_idx(3),7) = {path_normal_probe};
    miniblock_task_rando(probe_type_active_idx(3),10) = {33};
    miniblock_task_rando(probe_type_active_idx(4),3) = {'Late'};
    miniblock_task_rando(probe_type_active_idx(4),5) = {late_onset};
    miniblock_task_rando(probe_type_active_idx(4),4) = {'Pitch'};
    miniblock_task_rando(probe_type_active_idx(4),7) = {path_pitch_probe};
    miniblock_task_rando(probe_type_active_idx(4),10) = {34};

    probe_type_passive_idx = yes_passive(randperm(length(yes_passive)));
    miniblock_task_rando(probe_type_passive_idx(1),3) = {'Early'};
    miniblock_task_rando(probe_type_passive_idx(1),5) = {early_onset};
    miniblock_task_rando(probe_type_passive_idx(1),4) = {'Normal'};
    miniblock_task_rando(probe_type_passive_idx(1),7) = {path_normal_probe};
    miniblock_task_rando(probe_type_passive_idx(1),10) = {41};
    miniblock_task_rando(probe_type_passive_idx(2),3) = {'Early'};
    miniblock_task_rando(probe_type_passive_idx(2),5) = {early_onset};
    miniblock_task_rando(probe_type_passive_idx(2),4) = {'Pitch'};
    miniblock_task_rando(probe_type_passive_idx(2),7) = {path_pitch_probe};
    miniblock_task_rando(probe_type_passive_idx(2),10) = {42};
    miniblock_task_rando(probe_type_passive_idx(3),3) = {'Late'};
    miniblock_task_rando(probe_type_passive_idx(3),5) = {late_onset};
    miniblock_task_rando(probe_type_passive_idx(3),4) = {'Normal'};
    miniblock_task_rando(probe_type_passive_idx(3),7) = {path_normal_probe};
    miniblock_task_rando(probe_type_passive_idx(3),10) = {43};
    miniblock_task_rando(probe_type_passive_idx(4),3) = {'Late'};
    miniblock_task_rando(probe_type_passive_idx(4),5) = {late_onset};
    miniblock_task_rando(probe_type_passive_idx(4),4) = {'Pitch'};
    miniblock_task_rando(probe_type_passive_idx(4),7) = {path_pitch_probe};
    miniblock_task_rando(probe_type_passive_idx(4),10) = {44};

    % Assign probe_onset
    miniblock_task_rando(no_idx,5) = {0};

    % Assign probe_marker for no-probe trials
    miniblock_task_rando(no_idx,10) = {99};

    % Assign probe_duration
    miniblock_task_rando(no_idx,8) = {0};
    miniblock_task_rando(yes_idx,8) = {0.08};

    % Assign probe_intensity
    miniblock_task_rando(no_idx,6) = {0};
    miniblock_task_rando(yes_idx,6) = {1};

    % Add subject ID
    miniblock_task_rando(:,9) = {subj};

    % Store current block in overall array
    all_trials = vertcat(all_trials, miniblock_task_rando);

end

% Create and export file
conditions_table = cell2table(all_trials, "VariableNames",{'probe', 'task', 'probe_onset_cat', 'probe_type', 'probe_onset', 'probe_intensity', 'stim_file', 'probe_duration' ,'subj', 'probe_marker', 'task_marker', 'rec_duration'});

% Split table into 4 tables
if mod(height(conditions_table),4) == 0
    block_size = height(conditions_table)/4;
else
    warning('Table-Splitting not possible! Cannot divide by 4!')
end

if mod(block_size,16) == 0
    conditions_table_1 = conditions_table(1:block_size, :);
    conditions_table_2 = conditions_table(block_size+1:2*block_size, :);
    conditions_table_3 = conditions_table(2*block_size+1:3*block_size, :);
    conditions_table_4 = conditions_table(3*block_size+1:end, :);
else
    warning('Table-Splitting not possible! Minimal Blocksize!')
end

% Add filenames & meta-table
conditions_table_1_filename = fullfile(STIMULIPATH, [subj '_conditions_1.xlsx']);
conditions_table_2_filename = fullfile(STIMULIPATH, [subj '_conditions_2.xlsx']);
conditions_table_3_filename = fullfile(STIMULIPATH, [subj '_conditions_3.xlsx']);
conditions_table_4_filename = fullfile(STIMULIPATH, [subj '_conditions_4.xlsx']);

conditions_table_filenames = {conditions_table_1_filename; conditions_table_2_filename; ...
    conditions_table_3_filename; conditions_table_4_filename};

condition_table_meta = cell2table(conditions_table_filenames, "VariableNames",{'conditions_file'});


% Export tables
writetable(conditions_table, fullfile(STIMULIPATH, [subj '_conditions.xlsx']));
writetable(conditions_table_1, conditions_table_1_filename);
writetable(conditions_table_2, conditions_table_2_filename);
writetable(conditions_table_3, conditions_table_3_filename);
writetable(conditions_table_4, conditions_table_4_filename);
writetable(condition_table_meta, fullfile(STIMULIPATH, [subj '_conditions_meta.xlsx']));


% End of processing
disp('Conditions file saved!');
disp('READY FOR STIMULI RECORDINGS')