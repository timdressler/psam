clear
close all
clc

subj = 'sub-89';

% Set up paths
SCRIPTPATH = cd;
if regexp(SCRIPTPATH, regexptranslate('wildcard','*psam\experiment_script')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, '\experiment_script');
STIMULIPATH = fullfile(MAINPATH, ['BIDS/' subj '/stimuli']);
STIMULIPATH_Normal = fullfile(MAINPATH, ['BIDS/' subj '/stimuli/all_normal']);
STIMULIPATH_Pitch = fullfile(MAINPATH, ['BIDS/' subj '/stimuli/all_pitch']);
STIMULIPATH_Raw = fullfile(MAINPATH, ['BIDS/' subj '/stimuli/all_raw']);
EEGPATH = fullfile(MAINPATH, ['BIDS/' subj '/eeg']);
BEHAVIORALPATH = fullfile(MAINPATH, ['BIDS/' subj '/beh']);

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);


mt_check_folder_TD(MAINPATH, STIMULIPATH,STIMULIPATH_Normal,STIMULIPATH_Pitch,STIMULIPATH_Raw, EEGPATH, BEHAVIORALPATH)

early_onset = 2.8;
late_onset = 2.9;
n_trials = 960;
if ~mod(n_trials, 16) == 0
    error('trial number has to be divideable by 16')
else
    disp('Trial number OK')
end

max_repeats = 4; 
num_iterations = n_trials/16; 

path_normal_probe = fullfile(STIMULIPATH, [subj '_normal_probe.wav']);
path_pitch_probe = fullfile(STIMULIPATH, [subj '_pitch_probe.wav']);

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
            temp_task = miniblock_task; % Copy to avoid modifying original
            
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
            
            success = true; % If successful, exit loop
        catch
            % If an error occurs, retry
            %%warning('Sequence generation failed. Retrying...');
        end
    end

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
    miniblock_task_rando(no_active,2) = {'Active'};
    miniblock_task_rando(yes_passive,2) = {'Passive'};
    miniblock_task_rando(no_passive,2) = {'Passive'};

    % Assign probe_type, probe_onset, probe_file and probe_onset_cat
    probe_type_active_idx = yes_active(randperm(length(yes_active)));
    miniblock_task_rando(probe_type_active_idx(1),3) = {'Early'};
    miniblock_task_rando(probe_type_active_idx(1),5) = {early_onset};
    miniblock_task_rando(probe_type_active_idx(1),4) = {'Normal'};
    miniblock_task_rando(probe_type_active_idx(1),7) = {path_normal_probe};
    miniblock_task_rando(probe_type_active_idx(2),3) = {'Early'};
    miniblock_task_rando(probe_type_active_idx(2),5) = {early_onset};
    miniblock_task_rando(probe_type_active_idx(2),4) = {'Pitch'};
    miniblock_task_rando(probe_type_active_idx(2),7) = {path_pitch_probe};
    miniblock_task_rando(probe_type_active_idx(3),3) = {'Late'};
    miniblock_task_rando(probe_type_active_idx(3),5) = {late_onset};
    miniblock_task_rando(probe_type_active_idx(3),4) = {'Normal'};
    miniblock_task_rando(probe_type_active_idx(3),7) = {path_normal_probe};
    miniblock_task_rando(probe_type_active_idx(4),3) = {'Late'};
    miniblock_task_rando(probe_type_active_idx(4),5) = {late_onset};
    miniblock_task_rando(probe_type_active_idx(4),4) = {'Pitch'};
    miniblock_task_rando(probe_type_active_idx(4),7) = {path_pitch_probe};


    probe_type_passive_idx = yes_passive(randperm(length(yes_passive)));
    miniblock_task_rando(probe_type_passive_idx(1),3) = {'Early'};
    miniblock_task_rando(probe_type_passive_idx(1),5) = {early_onset};
    miniblock_task_rando(probe_type_passive_idx(1),4) = {'Normal'};
    miniblock_task_rando(probe_type_passive_idx(1),7) = {path_normal_probe};
    miniblock_task_rando(probe_type_passive_idx(2),3) = {'Early'};
    miniblock_task_rando(probe_type_passive_idx(2),5) = {early_onset};
    miniblock_task_rando(probe_type_passive_idx(2),4) = {'Pitch'};
    miniblock_task_rando(probe_type_passive_idx(2),7) = {path_pitch_probe};
    miniblock_task_rando(probe_type_passive_idx(3),3) = {'Late'};
    miniblock_task_rando(probe_type_passive_idx(3),5) = {late_onset};
    miniblock_task_rando(probe_type_passive_idx(3),4) = {'Normal'};
    miniblock_task_rando(probe_type_passive_idx(3),7) = {path_normal_probe};
    miniblock_task_rando(probe_type_passive_idx(4),3) = {'Late'};
    miniblock_task_rando(probe_type_passive_idx(4),5) = {late_onset};
    miniblock_task_rando(probe_type_passive_idx(4),4) = {'Pitch'};
    miniblock_task_rando(probe_type_passive_idx(4),7) = {path_pitch_probe};

    % Assign probe_onset
    miniblock_task_rando(no_idx,5) = {0};

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

conditions_table = cell2table(all_trials, "VariableNames",{'probe', 'task', 'probe_onset_cat', 'probe_type', 'probe_onset', 'probe_intensity', 'stim_file', 'probe_duration' ,'subj'});

writetable(conditions_table, fullfile(STIMULIPATH, [subj '_conditions.xlsx']));

