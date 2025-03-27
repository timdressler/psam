% tid_psam_set_markers.m
%
% Set and rename event markers based on detected peaks
%
% Tim Dressler, 27.03.2025

clear
close all
clc

% Start EEGLAB
eeglab;
tic;

% Load EEG data
EEG = pop_loadbv('C:\Users\timdr\OneDrive\Uni_Oldenburg\4_Semester\Master_Thesis\Analysis_Experiment\psam\testing\timingtest\timingtest_240325\sub-99\eeg\', ...
                 'timingtest_24032025.vhdr', [], []);
EEG = eeg_checkset(EEG);

% Detect peaks in Channel 1 with thresholding and refractory period
threshold = 1000;  % Peak detection threshold
min_interval = 100; % Minimum time (ms) between detected peaks to avoid duplicates

peak_latencies = [];
last_peak_time = -inf; % Track last detected peak time

for i = 1:length(EEG.times)
    if abs(EEG.data(1, i)) > threshold % Check if value exceeds threshold
        peak_time = EEG.times(i);
        if peak_time - last_peak_time > min_interval % Avoid duplicate detections
            peak_latencies = [peak_latencies; peak_time]; % Store valid peaks
            last_peak_time = peak_time;
        end
    end
end

% Display number of detected peaks
num_peaks = length(peak_latencies);
disp(['Number of detected sound onsets: ' num2str(num_peaks)]);

% Extract existing event latencies (converted to ms)
event_latencies = [EEG.event.latency] / EEG.srate * 1000;
event_types = {EEG.event.type}; % Get event types

% Find closest markers and add new markers
marker_window = 200; % ±200ms search window
new_marker_count = 0; % Counter for new markers

for peak_time = peak_latencies'
    % Find closest event within ±200 ms
    [min_diff, closest_idx] = min(abs(event_latencies - peak_time));
    
    if min_diff <= marker_window && closest_idx > 0
        old_marker = event_types{closest_idx}; % Get closest event marker

        if startsWith(old_marker, 'S ')
            old_marker_num = str2double(strtrim(extractAfter(old_marker, 2))); % Extract marker number
            if ~isnan(old_marker_num)
                new_marker_num = 900 + old_marker_num; % Generate new marker name
                new_marker_name = ['S ' num2str(new_marker_num)];

                % Add new event directly to EEG.event (faster than pop_editeventvals)
                EEG.event(end+1).latency = peak_time / 1000 * EEG.srate; % Convert ms to samples
                EEG.event(end).type = new_marker_name;
                EEG.event(end).duration = 0.001;
                EEG.event(end).channel = 0;
                EEG.event(end).code = 'Stimulus';

                new_marker_count = new_marker_count + 1; % Increment counter
            end
        end
    end
end

% Update EEG structure
EEG = eeg_checkset(EEG);

% Save updated dataset
pop_saveset(EEG, 'filename', 'updated_EEG.set', 'filepath', ...
    'C:\Users\timdr\OneDrive\Uni_Oldenburg\4_Semester\Master_Thesis\Analysis_Experiment\psam\testing\timingtest\timingtest_240325\sub-99\eeg\');

% Display number of new markers added
disp(['Number of new markers added: ' num2str(new_marker_count)]);
disp('Processing complete.');

toc
