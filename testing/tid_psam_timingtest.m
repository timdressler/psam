% tid_psam_timingtest.m
%
% Performs timing test.
%
% Tim Dressler, 07.03.2025

clear 
close all
clc

% Variables to edit
MARKER_CHAN = 31;

% Start eeglab
eeglab;

% Load data
EEG = pop_loadbv('C:\Users\timdr\OneDrive\Uni_Oldenburg\4_Semester\Master_Thesis\Analysis_Experiment\psam\testing\timingtest\timingtest_240325\sub-99\eeg\', 'timingtest_24032025.vhdr', [], []);
EEG = eeg_checkset( EEG );

% Epoch and BL-Correction
%%EEG2 = pop_epoch( EEG, {'S 31', 'S 32','S 33','S 34','S 41','S 42','S 43','S 44'},[-0.2  0.2], 'epochinfo', 'yes');
EEG2 = pop_epoch( EEG, {'S 931', 'S 932','S 933','S 934','S 941','S 942','S 943','S 944'},[-0.2  0.2], 'epochinfo', 'yes');
%%EEG2 = pop_epoch( EEG, {  'S  1'  }, [-0.2         0.2], 'epochinfo', 'yes');
EEG2 = pop_rmbase( EEG2, [-200 0] ,[]);

% Plot: Sound 
fig=figure;
set(fig,'defaultTextFontSize',14);
set(fig,'defaultAxesFontSize',14);
subplot(1,2,1)
for k=1:size(EEG2.data,3)
plot(EEG2.times,(EEG2.data(MARKER_CHAN,:,k)))
hold on
end
xlabel('time (ms)')
title('Audio marker')

% Find sound onsets
for k=1:size(EEG2.data,3)
   a=find(abs(squeeze(EEG2.data(MARKER_CHAN,1:end,k))')>1000);  
   sound_onsets(k)=EEG2.times(a(1));
   
end

% Plot: Latency
subplot(1,2,2)
plot(sound_onsets,'bo')
xlabel('trial')
ylabel('sound onset relative to trigger (ms)')

% Display results
disp(['Mean Latency: ' num2str(mean(sound_onsets))])
disp(['SD Latency: ' num2str(std(sound_onsets))])



%archive

% % clear
% % close all
% % clc
% % 
% % % Start eeglab
% % eeglab;
% % tic;
% % % Load data
% % EEG = pop_loadbv('C:\Users\timdr\OneDrive\Uni_Oldenburg\4_Semester\Master_Thesis\Analysis_Experiment\psam\testing\timingtest\timingtest_240325\sub-99\eeg\', 'timingtest_24032025.vhdr', [], []);
% % EEG = eeg_checkset(EEG);
% % 
% % % Detect peaks in Channel 1 with thresholding and refractory period
% % threshold = 1000;  % Peak detection threshold
% % min_interval = 100; % Minimum time (ms) between detected peaks to avoid duplicates
% % peak_latencies = []; % Store detected peak latencies
% % 
% % last_peak_time = -inf; % Track last detected peak time
% % 
% % for i = 1:length(EEG.times)
% %     if abs(EEG.data(1, i)) > threshold % Check if the value exceeds threshold
% %         peak_time = EEG.times(i);
% %         if peak_time - last_peak_time > min_interval % Avoid double detections
% %             peak_latencies = [peak_latencies peak_time]; % Store valid peaks
% %             last_peak_time = peak_time;
% %         end
% %     end
% % end
% % 
% % % Display the number of detected peaks
% % num_peaks = length(peak_latencies);
% % disp(['Number of detected sound onsets: ' num2str(num_peaks)]);
% % 
% % % Extract existing event latencies (converted to ms)
% % event_latencies = [EEG.event.latency] / EEG.srate * 1000;
% % event_types = {EEG.event.type}; % Get event types
% % 
% % % Find closest markers and add new markers
% % marker_window = 200; % ±200ms search window
% % new_marker_count = 0; % Counter for new markers
% % 
% % for peak_idx = 1:num_peaks
% %     peak_time = peak_latencies(peak_idx);
% % 
% %     % Find indices of events within ±200 ms of the peak
% %     nearby_idx = find(abs(event_latencies - peak_time) <= marker_window);
% % 
% %     if ~isempty(nearby_idx) % If a nearby event exists
% %         % Find the closest marker
% %         [~, closest_idx] = min(abs(event_latencies(nearby_idx) - peak_time));
% %         closest_marker_idx = nearby_idx(closest_idx);
% %         old_marker = event_types{closest_marker_idx}; % Get closest event marker
% % 
% %         if startsWith(old_marker, 'S ')
% %             old_marker_num = str2double(strtrim(extractAfter(old_marker, 2))); % Extract marker number
% %             if ~isnan(old_marker_num)
% %                 new_marker_num = 900 + old_marker_num; % Generate new marker name
% %                 new_marker_name = ['S ' num2str(new_marker_num)];
% % 
% %                 % Convert latency from ms to samples
% %                 new_latency = peak_time ./ 1000;
% % 
% %                 % Insert new event using pop_editeventvals
% %                 evalc('EEG = pop_editeventvals(EEG, ''insert'', {1,[],[],[],[],[],[],[],[],[]}, ''changefield'', {1, ''latency'', new_latency}, ''changefield'', {1, ''duration'', 0.001}, ''changefield'', {1, ''channel'', 0}, ''changefield'', {1, ''type'', new_marker_name}, ''changefield'', {1, ''code'', ''Stimulus''})');
% % 
% % 
% %                 new_marker_count = new_marker_count + 1; % Increment counter
% %             end
% %         end
% %     end
% % end
% % 
% % % Update EEG structure
% % EEG = eeg_checkset(EEG);
% % 
% % % Save updated dataset
% % pop_saveset(EEG, 'filename', 'updated_EEG.set', 'filepath', 'C:\Users\timdr\OneDrive\Uni_Oldenburg\4_Semester\Master_Thesis\Analysis_Experiment\psam\testing\timingtest\timingtest_240325\sub-99\eeg\');
% % 
% % % Plot detected peaks
% % figure;
% % plot(EEG.times, EEG.data(1, :)); hold on;
% % scatter(peak_latencies, ones(size(peak_latencies)) * threshold, 'ro'); % Mark detected peaks
% % xlabel('Time (ms)');
% % ylabel('EEG Signal (µV)');
% % title('Detected Peaks in Channel 1');
% % hold off;
% % 
% % % Display the number of new markers added
% % disp(['Number of new markers added: ' num2str(new_marker_count)]);
% % disp('Processing complete.');
% % 
% % toc