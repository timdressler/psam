clear
close all
clc

subj = 'sub-99'

% Set up paths
SCRIPTPATH = cd;
if regexp(SCRIPTPATH, regexptranslate('wildcard','*psam\experiment_script')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, '\experiment_script');
STIMULIPATH = fullfile(MAINPATH, ['BIDS/' subj '/stimuli']);
STIMULIPATH_Normal = fullfile(STIMULIPATH, 'all_normal');
STIMULIPATH_Pitch = fullfile(STIMULIPATH, 'all_pitch');
STIMULIPATH_Raw = fullfile(STIMULIPATH, 'all_raw');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

% Load table and determine median index
f0_table = readtable(fullfile(STIMULIPATH_Raw, [subj '_f0_table.csv']), 'Delimiter', ',');
sorted_f0_table = sortrows(f0_table, 1);
height_sorted_f0_table = height(sorted_f0_table);
median_index = (height_sorted_f0_table + 1) / 2;

% Initialize variables
change_attempts = 0;  % Track number of changes
shift_values = [-1, 1, -2, 2];  % Change sequence pattern
selection_done = false;  % Flag to stop selection process

% Function to select file details
select_file = @(idx) table2cell(sorted_f0_table(idx, ["f0_tab_normal", "f0_tab_pitched", "filename_tab"]));

% Start with median index
current_index = median_index;
final_probe_properties = select_file(current_index);

while ~selection_done
    % Define source file paths
    normal_probe_file = fullfile(STIMULIPATH_Normal, [final_probe_properties{1,3} '_normal.wav']);
    pitch_probe_file = fullfile(STIMULIPATH_Pitch, [final_probe_properties{1,3} '_pitch.wav']);

    % Load and plot
    [normal_probe_data, FS_normal] = audioread(normal_probe_file);
    [pitch_probe_data, FS_pitch] = audioread(pitch_probe_file);





    % User decision dialog
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
            sound(normal_probe_data, FS_normal); pause();
            sound(pitch_probe_data, FS_pitch); pause();
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

    % Handle user selection
    if strcmp(change_file_yes_no, 'no')
        selection_done = true;  % Stop selection loop

    elseif strcmp(change_file_yes_no, 'yes')
        if change_attempts < 4
            change_attempts = change_attempts + 1;
            % Apply shift pattern for changes
            current_index = median_index + shift_values(change_attempts);
            final_probe_properties = select_file(current_index);

        else
            close all
            error('Max number of file changes reached. No further changes allowed. NO PROBE SAVED! REPEAT STIMULI RECORDING!');
            selection_done = true;  % Stop selection loop
        end
    end
end

% Define final destination paths
destinationFile_normal = fullfile(STIMULIPATH, [subj '_normal_probe.wav']);
destinationFile_pitch = fullfile(STIMULIPATH, [subj '_pitch_probe.wav']);

% Copy selected files to final location
copyfile(fullfile(STIMULIPATH_Normal, [final_probe_properties{1,3} '_normal.wav']), destinationFile_normal);
copyfile(fullfile(STIMULIPATH_Pitch, [final_probe_properties{1,3} '_pitch.wav']), destinationFile_pitch);

% Export final probe properties as Excel file
exportFile = fullfile(STIMULIPATH, [subj '_probe_properties.xlsx']);
final_probe_table = cell2table(final_probe_properties, ...
    'VariableNames', {'f0_tab_normal', 'f0_tab_pitched', 'filename_tab'});
final_probe_table.change_attemps = change_attempts;

writetable(final_probe_table, exportFile, 'FileType', 'spreadsheet');
disp(['Final probe properties saved as: ' exportFile]);
disp('READY FOR MAIN EXPERIMENT')

close all
