% tid_psam_markertest.m
%
% Gets marker latencies.
%
% Tim Dressler, 07.03.2025

clear 
close all
clc

% Start eeglab
eeglab;

% Load data
EEG = pop_loadbv('C:\Users\timdr\OneDrive\Uni_Oldenburg\4_Semester\Master_Thesis\Analysis_Experiment\psam\testing\timingtest\timingtest_240325_4_1024\sub-99\eeg\', 'timingtest_24032025_4_1024.vhdr', [], []);
EEG = eeg_checkset( EEG );

% Extract event latencies and types
event_types = {EEG.event.type}; 
event_latencies = [EEG.event.latency]; 

% Define expected event markers
expected_events = {'S 21', 'S 22', 'S 31', 'S 32', 'S 33', 'S 34', 'S 41', 'S 42', 'S 43', 'S 44'};

% Check for missing events
missing_events = setdiff(expected_events, unique(event_types));

if ~isempty(missing_events)
    warning('The following expected events are missing from the data: %s', strjoin(missing_events, ', '));
end

% Set up variables
time_diffs_31_32 = [];
time_diffs_33_34 = [];
time_diffs_41_42 = [];
time_diffs_43_44 = [];

% Loop through events and find pairs
for i = 1:length(event_types)-1
    if strcmp(event_types{i}, 'S 21')
        if any(strcmp(event_types{i+1}, {'S 31', 'S 32'}))
            delay = (event_latencies(i+1) - event_latencies(i)) / EEG.srate;
            time_diffs_31_32 = [time_diffs_31_32, delay];
        elseif any(strcmp(event_types{i+1}, {'S 33', 'S 34'}))
            delay = (event_latencies(i+1) - event_latencies(i)) / EEG.srate;
            time_diffs_33_34 = [time_diffs_33_34, delay];
        end
    elseif strcmp(event_types{i}, 'S 22')
        if any(strcmp(event_types{i+1}, {'S 41', 'S 42'}))
            delay = (event_latencies(i+1) - event_latencies(i)) / EEG.srate;
            time_diffs_41_42 = [time_diffs_41_42, delay];
        elseif any(strcmp(event_types{i+1}, {'S 43', 'S 44'}))
            delay = (event_latencies(i+1) - event_latencies(i)) / EEG.srate;
            time_diffs_43_44 = [time_diffs_43_44, delay];
        end
    end
end

% Compute mean and standard deviation for each group
mean_delay_31_32 = mean(time_diffs_31_32);
sd_delay_31_32 = std(time_diffs_31_32);

mean_delay_33_34 = mean(time_diffs_33_34);
sd_delay_33_34 = std(time_diffs_33_34);

mean_delay_41_42 = mean(time_diffs_41_42);
sd_delay_41_42 = std(time_diffs_41_42);

mean_delay_43_44 = mean(time_diffs_43_44);
sd_delay_43_44 = std(time_diffs_43_44);

% Display results
disp(['Mean delay for S 31/S 32: ', num2str(mean_delay_31_32), ' seconds']);
disp(['Standard deviation for S 31/S 32: ', num2str(sd_delay_31_32), ' seconds']);

disp(['Mean delay for S 33/S 34: ', num2str(mean_delay_33_34), ' seconds']);
disp(['Standard deviation for S 33/S 34: ', num2str(sd_delay_33_34), ' seconds']);

disp(['Mean delay for S 41/S 42: ', num2str(mean_delay_41_42), ' seconds']);
disp(['Standard deviation for S 41/S 42: ', num2str(sd_delay_41_42), ' seconds']);

disp(['Mean delay for S 43/S 44: ', num2str(mean_delay_43_44), ' seconds']);
disp(['Standard deviation for S 43/S 44: ', num2str(sd_delay_43_44), ' seconds']);
