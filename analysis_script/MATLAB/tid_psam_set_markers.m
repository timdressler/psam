% tid_psam_set_markers.m
%
% Detects sound onsets and modifies EEG.event structure accordingly & includes control markers.
% The control markers are placed in no-probe trials at the timepoints were
% probes whould have been presented in case of active trials. To
% increase the SNR, only early and late conditions are reflected. In
% other words, instead of classifying no-probe trials as Control -
% Active - Early - Unaltered; Control - Active - Early - Altered;
% Control - Active - Late - Unaltered; ..., no-probe trial are
% classified as Control - Active - Early; Control - Active - Late;
%
% Conditions
% S 931 - Active - Early Probe - Unaltered
% S 932 - Active - Early Probe - Altered
% S 933 - Active - Late Probe - Unaltered
% S 934 - Active - Late Probe - Altered
%
% S 941 - Passive - Early Probe - Unaltered
% S 942 - Passive - Early Probe - Altered
% S 943 - Passive - Late Probe - Unaltered
% S 944 - Passive - Late Probe - Altered
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
OUTPATH = fullfile(MAINPATH,'data\processed_data\markers_included\');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH,INPATH,OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
THRESHOLD = 1000;  % Peak detection THRESHOLD
MIN_INTERVAL = 100; % Minimum time (ms) between detected peaks to avoid duplicates
MARKER_CHAN = 31;
EARLY_ONSET = 2.6;
LATE_ONSET = 2.8;

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
    subj_file = ['tid_psam_' subj '.vhdr'];
    subj_path = fullfile(INPATH, [subj '\eeg\']);

    % Start EEGLAB
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    tic;

    % Load EEG data
    EEG = pop_loadbv(subj_path, subj_file, [], []);
    EEG = eeg_checkset(EEG);

    % Delete not needed markers
    needed_markers = {  'S  1'  'S  5'  'S 21'  'S 22'  'S 31'  'S 32'  'S 33'  'S 34'  'S 41'  'S 42'  'S 43'  'S 44'  'boundary'  };
    all_markers = {EEG.event.type};
    not_needed_markers_indices = find(~ismember(all_markers, needed_markers));
    EEG = pop_editeventvals(EEG,'delete',not_needed_markers_indices);

    % Set control markers
    con_act_conditions_type = ["con_act_early", "con_act_late"];
    con_act_conditions = repelem(con_act_conditions_type, 120); % Create an array with 120 repetitions of each category
    con_act_conditions = con_act_conditions(randperm(length(con_act_conditions)));

    con_pas_conditions_type = ["con_pas_early", "con_pas_late"];
    con_pas_conditions = repelem(con_pas_conditions_type, 120); % Create an array with 120 repetitions of each category
    con_pas_conditions = con_pas_conditions(randperm(length(con_pas_conditions)));


    new_con_act_marker_count = 0;
    new_con_pas_marker_count = 0;
    for event = 1:length(EEG.event)-1
        if strcmp(EEG.event(event).type, 'S 21') && strcmp(EEG.event(event+1).type, 'S  5')
            new_con_act_marker_count = new_con_act_marker_count + 1;

            if strcmp(con_act_conditions(new_con_act_marker_count), 'con_act_early')
                new_con_act_latency = EEG.event(event).latency + (EARLY_ONSET .* EEG.srate);
            elseif strcmp(con_act_conditions(new_con_act_marker_count), 'con_act_late')
                new_con_act_latency = EEG.event(event).latency + (LATE_ONSET .* EEG.srate);
            end

            % Add new event
            EEG.event(end+1).latency = new_con_act_latency;
            EEG.event(end).type = con_act_conditions(new_con_act_marker_count);
            EEG.event(end).duration = 1;
            EEG.event(end).channel = 0;
            EEG.event(end).code = 'Stimulus';

        elseif strcmp(EEG.event(event).type, 'S 22') && strcmp(EEG.event(event+1).type, 'S  5')
            new_con_pas_marker_count = new_con_pas_marker_count + 1;

            if strcmp(con_pas_conditions(new_con_pas_marker_count), 'con_pas_early')
                new_con_pas_latency = EEG.event(event).latency + (EARLY_ONSET .* EEG.srate);
            elseif strcmp(con_pas_conditions(new_con_pas_marker_count), 'con_pas_late')
                new_con_pas_latency = EEG.event(event).latency + (LATE_ONSET .* EEG.srate);
            end

            % Add new event
            EEG.event(end+1).latency = new_con_pas_latency;
            EEG.event(end).type = con_pas_conditions(new_con_pas_marker_count);
            EEG.event(end).duration = 1;
            EEG.event(end).channel = 0;
            EEG.event(end).code = 'Stimulus';
        end
    end

    % Sanity Check: Number of added control markers
    if new_con_act_marker_count == 240 && new_con_pas_marker_count == 240

    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'num_of_added_control_markers';
    end

    % Detect peaks in Marker channel
    peak_latencies = [];
    last_peak_time = -inf; % Track last detected peak time

    for i = 1:length(EEG.times)
        if abs(EEG.data(MARKER_CHAN, i)) > THRESHOLD % Check if value exceeds THRESHOLD
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
    wrong_marker_count = 0; % Counter for wrong markers
    new_marker_latency = []; % Stores latencies of added markers

    for peak_time = peak_latencies'
        % Find closest event within ±200 ms
        [min_diff, closest_idx] = min(abs(event_latencies - peak_time));

        if min_diff <= marker_window && closest_idx > 0
            old_marker = event_types{closest_idx}; % Get closest event marker

            if startsWith(old_marker, 'S ')
                old_marker_num = str2double(strtrim(extractAfter(old_marker, 2))); % Extract marker number
                if ~isnan(old_marker_num) && any(old_marker_num == [31:34 41:44])
                    new_marker_num = 900 + old_marker_num; % Generate new marker name
                    new_marker_name = ['S ' num2str(new_marker_num)];

                    % Add new event
                    EEG.event(end+1).latency = peak_time / 1000 * EEG.srate; % Convert ms to samples
                    EEG.event(end).type = new_marker_name;
                    EEG.event(end).duration = 1;
                    EEG.event(end).channel = 0;
                    EEG.event(end).code = 'Stimulus';

                    new_marker_count = new_marker_count + 1;

                    new_marker_latency = [new_marker_latency; peak_time];
                else
                    warning('Invalid Marker')
                    wrong_marker_count = wrong_marker_count + 1;
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

    % Sanity Check: No wrong markers
    if wrong_marker_count ~= 0
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'wrong_marker';
    end

    % Update EEG structure
    EEG = eeg_checkset(EEG);

    % Sort events
    [~, sortIdx] = sort([EEG.event.latency]);
    EEG.event = EEG.event(sortIdx);
    EEG = eeg_checkset(EEG);

    % Save updated dataset
    pop_saveset(EEG, 'filename', [subj '_markers_inlcuded.set'], 'filepath', OUTPATH);

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
writetable(protocol,fullfile(OUTPATH, 'tid_psam_markers_included_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_markers_included_marked_subj.xlsx'))
end

check_done = 'tid_psam_set_markers_DONE'

delete(wb)