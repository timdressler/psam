% tid_psam_determine_loudness.m
%
% Plays selected stimuli to determine a pleasant loudness.
% Has to be executed AFTER tid_psam_create_conditions_file.m,
%   tid_psam_stimuli_recording_adapted.py, tid_psam_prepare_stimuli.praat
%   and tid_psam_select_stimuli.m.
%
% Tim Dressler, 30.03.2025

clear
close all
clc

set(0,'DefaultTextInterpreter','none')

% Set up paths
SCRIPTPATH = cd;
if regexp(SCRIPTPATH, regexptranslate('wildcard','*psam\experiment_script')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, '\experiment_script');
STIMULIPATH = fullfile(MAINPATH,'data/BIDS/stimuli');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH,STIMULIPATH)

% Get subject ID
prompt = {'Enter subject number:'};
dlgtitle = 'Subject Number Input';
dims = [1 35];
answer = inputdlg(prompt, dlgtitle, dims);

% Format subject ID
if ~isempty(answer)
    num = str2double(answer{1});
    formattedNum = sprintf('%02d', num);
    subj = ['sub-' formattedNum];
end

% Get filename
normal_probe_file = fullfile(STIMULIPATH, subj, [subj  '_normal_probe.wav']);

% Load audio
[normal_probe_data, FS] = audioread(normal_probe_file);

% Extract the right channel (assuming 2-channel stereo)
if size(normal_probe_data, 2) == 2
    normal_probe_data = normal_probe_data(:,2);
end

% Initial playback
for i = 1:10
    sound(normal_probe_data, FS);
    pause(1.5)
end

% Allow repeated playback
while true
    choice = questdlg('Do you want to listen again?', 'Listen again?', 'No, exit', 'Yes', 'No, exit');
    if strcmp(choice, 'Yes')
        for i = 1:10
            sound(normal_probe_data, FS);
            pause(1.5)
        end
    else
        break;
    end
end

% Ask for loudness selection
loudness_input = inputdlg('Enter selected loudness:', 'Loudness Selection', [1 35]);
if isempty(loudness_input)
    error('No loudness selected. Exiting.');
end
selected_loudness = str2double(loudness_input{1});

% Add selected loudness to probe_properties.xlsx
prope_properties_file = fullfile(STIMULIPATH,subj,[subj '_probe_properties.xlsx']);
prop_properties = readtable(prope_properties_file);
prop_properties.loudness = selected_loudness;
writetable(prop_properties, prope_properties_file);

% End of processing
disp('probe_properties.xlsx updated')
disp('READY FOR MAIN EXPERIMENT')
