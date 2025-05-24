% tid_psam_SVM_preparation.m
%
% Performs preparations for SVM analysis and saves data.
%
% Preparation includes the following steps
% Removes EOG electrodes
% Extracts early and late time windows relative to the 'go-signal'
% Extract the following features for each epoch, time window and channel
%   Mean amplitude, RMS, standard deviation of the amplitude, maximum and minimum amplitude, kurtosis, skewness, and zero-crossing rate
%   Bandpower for alpha (8-13 Hz), beta (13 - 30 Hz) and (low) gamma (30 -37 Hz)
%   Hjorth Parameters (Activity, Mobility and Complexity)
%
% Saves data
%
% Also see: tid_psam_hjorth_activity_TD, tid_psam_hjorth_mobility_TD, tid_psam_hjorth_complexity_TD.
%
% Tim Dressler, 05.05.2025

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
INPATH = fullfile(MAINPATH, 'data\processed_data\svm_preprocessed_clean\');
INPATH_HILBERT = fullfile(MAINPATH, 'data\processed_data\hilbert_prepared_clean');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\svm_prepared_clean');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EOG_CHAN = {'E29','E30'}; % Labels of EOG electrodes
EVENTS = {'go_act', 'go_pas'};
CHANI = 1; % Channel to plot sanity check ERP
WIN_EARLY_FROM = -700;
WIN_EARLY_TILL = -401;
WIN_LATE_FROM = -400;
WIN_LATE_TILL = -101;

% Sanity Check: Same window length
win_early_length = abs(WIN_EARLY_FROM - WIN_EARLY_TILL)+1;
win_late_length = abs(WIN_LATE_FROM - WIN_LATE_TILL)+1;
if win_early_length ~= win_late_length
    error('Different window lengths!')
end

% Set colors
main_blue = '#004F9F';
main_blue = sscanf(main_blue(2:end),'%2x%2x%2x',[1 3])/255;
main_red = '#D53D0E';
main_red = sscanf(main_red(2:end),'%2x%2x%2x',[1 3])/255;
main_green = '#00786B';
main_green = sscanf(main_green(2:end),'%2x%2x%2x',[1 3])/255;
light_blue = '#5BC5F2';
light_blue = sscanf(light_blue(2:end),'%2x%2x%2x',[1 3])/255;
main_yellow = '#FDC300';
main_yellow = sscanf(main_yellow(2:end),'%2x%2x%2x',[1 3])/255;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*.set'));

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_svm_preparation.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_svm_preparation.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load epoched data
    EEG = pop_loadset('filename',[subj '_svm_preprocessed_clean.set'],'filepath',INPATH);

    % Load epoched Hilbert transformed data
    load(fullfile(INPATH_HILBERT, [subj '_hilbert_preparation_clean.mat'])) % loads variable 'data_hilbert' and 'freq_bands'

    % Remove EOG channels as they are not used for classification
    EEG = pop_select( EEG, 'rmchannel',EOG_CHAN);

    % Get channel labels
    chan_labels = {EEG.chanlocs.labels};

    % Get epoch labels
    labels = cellfun(@(x) x{2}, {EEG.epoch.eventtype}, 'UniformOutput', false);

    % Define time windows
    windows = struct( ...
        'early', struct('from', WIN_EARLY_FROM, 'till', WIN_EARLY_TILL), ...
        'late', struct('from', WIN_LATE_FROM, 'till', WIN_LATE_TILL) ...
        );

    % Initialize feature struct
    features = struct();

    % Loop through windows
    for win_label = ["early", "late"]
        label = win_label;

        % Get time indices
        [~,start_idx] = min(abs(EEG.times - windows.(label).from));
        [~,end_idx] = min(abs(EEG.times - windows.(label).till));

        % Get window data
        data_win = EEG.data(:, start_idx:end_idx, :);
        features_data.(label).data_win = data_win;

        % Sanity Check: Correct dimensions
        if size(data_win,1) ~= EEG.nbchan || size(data_win,2) ~= win_early_length || size(data_win,3) == 1
            marked_subj{end+1,1} = subj;
            marked_subj{end,2} = 'wrong_dimensions_data';
        end

        % Time-domain features
        features.(label).mean = squeeze(mean(data_win, 2));
        features.(label).rms = squeeze(rms(data_win, 2));
        features.(label).sd = squeeze(std(data_win, [], 2));
        features.(label).min = squeeze(min(data_win, [], 2));
        features.(label).max = squeeze(max(data_win, [], 2));
        features.(label).kurtosis = squeeze(kurtosis(data_win, [], 2));
        features.(label).skewness = squeeze(skewness(data_win, [], 2));

        for i = 1:size(data_win, 3)
            features.(label).zerocrossing(:, i) = zerocrossrate(data_win(:,:,i)')';
        end

        % Frequency-domain features
        for freq_band_num = 1:length(data_hilbert)
            % Get window Hilbert transformed data
            data_hilbert_win = data_hilbert(freq_band_num).hilbert_transformed_epochs_clean(:, start_idx:end_idx, :);
            data_hilbert_label = sprintf('data_hilbert_%d_%d', freq_bands{freq_band_num});
            features_data.(label).(data_hilbert_label) = data_hilbert_win;

            % Sanity Check: Correct dimensions
            if size(data_hilbert_win,1) ~= EEG.nbchan || size(data_hilbert_win,2) ~= win_early_length || size(data_hilbert_win,3) == 1
                marked_subj{end+1,1} = subj;
                marked_subj{end,2} = 'wrong_dimensions_hilbert_data';
            end

            % Extract bandpower
            bandpower_label = sprintf('band_power_%d_%d', freq_bands{freq_band_num});
            features.(label).(bandpower_label) = squeeze(mean(data_hilbert_win, 2));
        end

        % Hjorth parameters
        % Get Activity
        features.(label).activity = tid_psam_hjorth_activity_TD(data_win);
        % Get Mobility
        features.(label).mobility = tid_psam_hjorth_mobility_TD(data_win);
        % Get Complexity
        features.(label).complexity = tid_psam_hjorth_complexity_TD(data_win);

        % Bring features into the format epochs x features
        feature_fields = fieldnames(features.(label));
        field_tables = cell(1, length(feature_fields));  % Preallocate

        for f = 1:length(feature_fields)
            field_name = feature_fields{f};
            field_data = features.(label).(field_name);  % channels x epochs

            % Transpose to epochs x channels
            field_data = field_data';

            % Create column names
            n_chans = size(field_data, 2);
            % Sanity Check: Matching number of channels
            if ~(n_chans == length(chan_labels))
                error('number of channels does not match')
            end
            col_names = arrayfun(@(i) sprintf('%s_%s', field_name, chan_labels{i}), 1:n_chans, 'UniformOutput', false);

            % Convert to table
            field_tables{f} = array2table(field_data, 'VariableNames', col_names);
        end

        % Concatenate and store result
        final_table = cat(2, field_tables{:});
        features.(label).feature_table = final_table;

        % Sanity Check: Correct dimensions of the feature table
        if ~(size(features.(label).feature_table,2) == length(chan_labels)*length(feature_fields) && size(features.(label).feature_table,1) == length(labels))
            error('dimensions of final feature table')
        end

        % Attach labels
        features.(label).feature_table.labels = labels';

        % Save features
        writetable(features.(label).feature_table,fullfile(OUTPATH, [subj '_features_' char(label) '.csv']))
    end

    % Sanity Check: Plot ERPs
    EEG.setname = [subj '_all_conds'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    y_lim_lower = -8;
    y_lim_upper = 8;
    figure;
    for cond = 1:length(EVENTS)
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off');
        erp = mean(EEG.data, 3);
        subplot(1,2,cond)
        plot(EEG.times, erp(CHANI,:), 'LineWidth', 1.5, 'Color',main_blue)
        xlim([EEG.times(1) EEG.times(end)])
        ylim([y_lim_lower y_lim_upper])
        xlabel('Time [ms]')
        ylabel('Amplitude [µV]')
        title(['ERPs for ' EVENTS{cond} ' condition'])
        hold on
        fill([WIN_LATE_FROM WIN_LATE_TILL WIN_LATE_TILL WIN_LATE_FROM], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], main_yellow, 'FaceAlpha',0.1, 'EdgeColor','none');
        fill([WIN_EARLY_FROM WIN_EARLY_TILL WIN_EARLY_TILL WIN_EARLY_FROM], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], main_red, 'FaceAlpha',0.1, 'EdgeColor','none');
        hold off
    end
    sgtitle(['Sanity Check ERPs for ' subj])
    saveas(gcf,fullfile(OUTPATH, [subj '_sanity_erp.png']));

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
writetable(protocol,fullfile(OUTPATH, 'tid_psam_svm_preparation_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_svm_preparation_marked_subj.xlsx'))
end

check_done = 'tid_psam_svm_preparation_DONE'

delete(wb)
