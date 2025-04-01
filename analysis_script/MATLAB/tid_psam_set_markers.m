% tid_psam_set_markers.m
%
% Detects sound onsets and modifies EEG.event structure accordingly.
%
% Tim Dressler, 27.03.2025

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
INPATH = fullfile(MAINPATH,'data/BIDS/');
OUTPATH = fullfile(MAINPATH,'data\processed_data\cleaned\');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH,INPATH,OUTPATH)

% Variables to edit
THRESHOLD = 1000;  % Peak detection THRESHOLD
MIN_INTERVAL = 100; % Minimum time (ms) between detected peaks to avoid duplicates
MARKER_CHAN = 31;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub*'));

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_set_markers.m');

clear subj_idx
for subj_idx = 1:length(dircont_subj)

    % Get subject ID
    subj = dircont_subj(subj_idx).name;

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_set_markers.m'])

    % Sanity Check: Number of EEG files == 1
    if length(dir(fullfile(INPATH, [subj '\eeg\*.vhdr']))) == 1

    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'num_of_eeg_files';
    end

    % Get filename & path
    subj_file = [subj '.vhdr'];
    subj_path = fullfile(INPATH, [subj '\eeg\']);

    % Start EEGLAB
    eeglab;
    tic;

    % Load EEG data
    EEG = pop_loadbv(subj_path, subj_file, [], []);
    EEG = eeg_checkset(EEG);

    peak_latencies = [];
    last_peak_time = -inf; % Track last detected peak time

    for i = 1:length(EEG.times)
        if abs(EEG.data(1, i)) > THRESHOLD % Check if value exceeds THRESHOLD
            peak_time = EEG.times(i);
            if peak_time - last_peak_time > MIN_INTERVAL % Avoid duplicate detections
                peak_latencies = [peak_latencies; peak_time]; % Store valid peaks
                last_peak_time = peak_time;
            end
        end
    end

    % Sanity Check: Number of detected peaks == 480
    if length(peak_latencies) == 480

    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'num_of_detected_peaks';
    end

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

                    % Add new event directly to EEG.event
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

    % Sanity Check: Number of added marker == 480
    if new_marker_count == 480

    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'num_of_new_markers';
    end

    % Update EEG structure
    EEG = eeg_checkset(EEG);

    % Save updated dataset
    pop_saveset(EEG, 'filename', [subj '_cleaned.set'], 'filepath', OUTPATH);

    % Update Protocol
    subj_time = toc;
    protocol{subj_idx,1} = subj;
    protocol{subj_idx,2} = subj_time;
    if any(strcmp(marked_subj, subj), 'all')
        protocol{subj_idx,3} = 'MARKED';
    else
        protocol{subj_idx,3} = 'OK';
    end

end

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'tid_psam_cleaned_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_marked_subj.xlsx'))
end

close(wb)

check_done = 'tid_psam_set_markers_DONE'