% tid_psam_select_stimuli.m
%
% Selects stimuli for main experiment.
% Has to be executed AFTER tid_psam_create_conditions_file.m,
%   tid_psam_stimuli_recording_adapted.py and tid_psam_prepare_stimuli.praat.
%
% Selection Process:
%   Selects probe with median F0
%   Presents selected probe 5 times   
%   Experimenter can choose whether to listen again, change to another
%       probe or keep the probe
%   In case a change is requested, the probe 'next to the median probe' is 
%       selected (first the one below the median, then the one above)
%   If the probe is still not usable the process is repeated with the
%       probes that 'lie 2 steps from the median one'. 
%   After that, no probe change is possible anymore and the recoring has to
%   b   e repeated
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

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

% Load table and determine median index
f0_table = readtable(fullfile(STIMULIPATH_Raw, [subj '_f0_table.csv']), 'Delimiter', ',');
sorted_f0_table = sortrows(f0_table, 1);
height_sorted_f0_table = height(sorted_f0_table);
median_index = (height_sorted_f0_table + 1) / 2;

% Initialize variables
change_attempts = 0;  % track number of changes requested by the experimenter
shift_values = [-1, 1, -2, 2];  % first select alternative file one below the median, then one above, ...
selection_done = false;  % set up flag to stop selection process

% Function to select file details
select_file = @(idx) table2cell(sorted_f0_table(idx, ["f0_tab_normal", "f0_tab_pitched", "db_tab_normal", "db_tab_pitched", "filename_tab"]));

% Start with median index
current_index = median_index;
final_probe_properties = select_file(current_index);

while ~selection_done
    % Define source file paths
    normal_probe_file = fullfile(STIMULIPATH_Normal, [final_probe_properties{1,5} '_normal.wav']);
    pitch_probe_file = fullfile(STIMULIPATH_Pitch, [final_probe_properties{1,5} '_pitch.wav']);

    % Load and plot
    [normal_probe_data, FS_normal] = audioread(normal_probe_file);
    [pitch_probe_data, FS_pitch] = audioread(pitch_probe_file);

    % Experimenter decision dialog
    change_file_yes_no = 'listen_again';
    while strcmp(change_file_yes_no, 'listen_again')
        change_file_yes_no = 'no';  % Default value

        figure;
        subplot(121)
        plot(normal_probe_data);
        title('Normal Probe');
        subplot(122)
        plot(pitch_probe_data);
        title('Pitch Probe');
        sgtitle(['Change: ' num2str(change_attempts) ' (max. 4 allowed)'])

        % Play sounds
        for i = 1:5
            pause(1)
            sound(normal_probe_data, FS_normal); 
            pause(1);
            sound(pitch_probe_data, FS_pitch); 
            pause(1);
        end

        d = dialog('Position', [500 300 350 150], 'Name', 'File Selection');

        uicontrol('Parent', d, 'Style', 'text', 'Position', [50 90 250 40], ...
            'String', 'Do you want to keep or change the file?', 'FontSize', 10);

        uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [20 30 80 30], ...
            'String', 'KEEP FILE', 'Callback', 'change_file_yes_no = ''no''; delete(gcbf);');

        uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [120 30 100 30], ...
            'String', 'CHANGE FILE', 'Callback', 'change_file_yes_no = ''yes''; delete(gcbf);');

        uicontrol('Parent', d, 'Style', 'pushbutton', 'Position', [240 30 100 30], ...
            'String', 'LISTEN AGAIN', 'Callback', 'change_file_yes_no = ''listen_again''; delete(gcbf);');

        uiwait(d);
    end

    % Handle experimenter selection
    if strcmp(change_file_yes_no, 'no')
        selection_done = true;  % Stop selection loop

    elseif strcmp(change_file_yes_no, 'yes')
        if change_attempts < 4
            change_attempts = change_attempts + 1;
            % Apply shift pattern for requested changes
            current_index = median_index + shift_values(change_attempts);
            final_probe_properties = select_file(current_index);
        else
            close all
            error('Max number of file changes reached. No further changes allowed. NO PROBE SAVED! REPEAT STIMULI RECORDING!');
            selection_done = true;  % stop selection loop
        end
    end
end

% Check Loudness
if round(final_probe_properties{1,3},4) == round(final_probe_properties{1,4},4)
    disp('Loudness OK')
else
    warning('Loudness not OK')
end

% Define final destination paths
destinationFile_normal = fullfile(STIMULIPATH, [subj '_normal_probe.wav']);
destinationFile_pitch = fullfile(STIMULIPATH, [subj '_pitch_probe.wav']);

% Copy selected files to final location
copyfile(fullfile(STIMULIPATH_Normal, [final_probe_properties{1,5} '_normal.wav']), destinationFile_normal);
copyfile(fullfile(STIMULIPATH_Pitch, [final_probe_properties{1,5} '_pitch.wav']), destinationFile_pitch);

% Export final probe properties 
exportFile = fullfile(STIMULIPATH, [subj '_probe_properties.xlsx']);
final_probe_table = cell2table(final_probe_properties, ...
    'VariableNames', {'f0_tab_normal', 'f0_tab_pitched', 'db_tab_normal', 'db_tab_pitched', 'filename_tab'});
final_probe_table.change_attempts = change_attempts;
writetable(final_probe_table, exportFile, 'FileType', 'spreadsheet');

% Display successful performance
disp(['Final probe properties saved as: ' exportFile]);
disp('READY FOR MAIN EXPERIMENT')

close all

