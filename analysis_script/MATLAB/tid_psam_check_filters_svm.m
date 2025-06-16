% tid_psam_check_filters_svm.m
%
% Creates plots to investigate the effects of filterimg based on the SVM-specific preprocessing (see tid_psam_svm_preprocessing.m, tid_psam_svm_preparation.m).
%
% Conditions
% S 931 = Active - Early Probe - Unaltered
% S 932 = Active - Early Probe - Altered
% S 933 = Active - Late Probe - Unaltered
% S 934 = Active - Late Probe - Altered
%
% S 941 = Passive - Early Probe - Unaltered
% S 942 = Passive - Early Probe - Altered
% S 943 = Passive - Late Probe - Unaltered
% S 944 = Passive - Late Probe - Altered
%
% Processing includes the following steps
%
%   Perform preprocesing once with and once without filtering (see tid_psam_svm_preprocessing.m, tid_psam_svm_preparation.m)
%   Excludes markerd trials based on tid_psam_exclude_trials.m
%   Create plots to illustarte effects of filtering.
%
% Saves plots
%
% Tim Dressler, 16.06.2025

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
INPATH_RAW = fullfile(MAINPATH, 'data\processed_data\markers_included\');
INPATH_ICA = fullfile(MAINPATH, 'data\processed_data\ica_preprocessed\');
INPATH_HILBERT = fullfile(MAINPATH, 'data\processed_data\hilbert_prepared_clean');
INPATH_EXCLUDED = fullfile(MAINPATH, 'data\processed_data\exclude_trials\');

OUTPATH = fullfile(MAINPATH, 'data\analysis_data\check_filters_svm');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_RAW, INPATH_ICA, INPATH_EXCLUDED, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
%   Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LCF = 0.1;
HCF = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EOG_CHAN = {'E29','E30'}; % Labels of EOG electrodes
EPO_FROM = -1;
EPO_TILL = 0.1;
BL_FROM = -1000;
BL_TILL = -800;
SD_PROB = 3;
SD_PROB_ICA = 3;
EVENTS = {'go_act', 'go_pas'};
%   Analysis
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
dircont_subj = dir(fullfile(INPATH_RAW, 'sub-*.set'));

% Sanity Check: Same number of files for raw data and ICA data
if length(dir(fullfile(INPATH_RAW, 'sub-*.set'))) == length(dir(fullfile(INPATH_ICA, 'sub-*.set'))) && length(dir(fullfile(INPATH_RAW, 'sub-*.set'))) == length(dir(fullfile(INPATH_EXCLUDED, 'sub-*.mat')))
else
    error('Number of raw data files, exclusion files and ICA data files does not match')
end

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_check_filters_svm.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    %% Preprocessing
    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_check_filters_svm.m'])
    tic;

    % Get to-be excluded trials based on tid_psam_exclude_trials.m
    exclusion_filename = fullfile(INPATH_EXCLUDED,[subj '_excluded_trials.mat']);
    load(exclusion_filename) % Loads variables 'excluded_trials_erp_beh' (not used here) and 'excluded_trials_svm' (used here)

    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load ICA data
    EEG = pop_loadset('filename',[subj '_ica_weights.set'],'filepath',INPATH_ICA);

    % Get flagged components as identified in tid_psam_ica_preprocessing.m
    flagged_comps = EEG.reject.gcompreject;

    EEG.setname = [subj '_ICA_weights'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Get bad channels as identified in tid_psam_ica_preprocessing.m
    chans_to_interp = EEG.badchans;

    % Load raw data
    EEG = pop_loadset('filename',[subj '_markers_inlcuded.set'],'filepath',INPATH_RAW);

    % Remove Marker-Channel
    EEG = pop_select( EEG, 'rmchannel',{'M'});

    % Add channel locations
    EEG.chanlocs = readlocs( fullfile(MAINPATH,'\config\elec_96ch_adapted.elp'));
    EEG = eeg_checkset( EEG );

    % Add type = EOG for EOG electrodes
    eog_chani = find(ismember({EEG.chanlocs.labels}, EOG_CHAN));
    [EEG.chanlocs(eog_chani).type] = deal('EOG');
    EEG = eeg_checkset( EEG );

    % Store channel locations in another field
    EEG.urchanlocs = EEG.chanlocs;

    % Rename events
    clear event
    for event = 1:length(EEG.event)
        if strcmp(EEG.event(event).type, 'S  5') && strncmp(EEG.event(event-1).type, 'con', 3) % Get all go-signals (S  5) during no-probe trials
            if strcmp(EEG.event(event-2).type, 'S 21')
                EEG.event(event).type = 'go_act';
            elseif strcmp(EEG.event(event-2).type, 'S 22')
                EEG.event(event).type = 'go_pas';
            else
                error('Unkown Marker!')
            end
        end
    end

    EEG.setname = [subj '_ready_for_preprocessing'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply Filter
    % Highpass-Filter
    %%LCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth(LCF)); % Get filter order (also see pop_firwsord.m, tid_psam_get_transition_bandwidth.m)
    %%EEG = pop_firws(EEG, 'fcutoff', LCF, 'ftype', 'highpass', 'wtype', 'hamming', 'forder',LCF_ord, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
    % Lowpass-Filter
    HCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth(HCF)); % Get filter order (also see pop_firwsord.m, tid_psam_get_transition_bandwidth.m)
    EEG = pop_firws(EEG, 'fcutoff', HCF, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder',HCF_ord, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

    %%EEG = pop_eegfiltnew(EEG, 'locutoff',LCF,'hicutoff',0);
    %%EEG = pop_eegfiltnew(EEG, 'hicutoff',HCF,'plotfreqz',0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Remove bad channels as identified in tid_psam_ica_preprocessing.m (see above)
    EEG.badchans = chans_to_interp;
    EEG = pop_select(EEG,'nochannel', EEG.badchans);

    EEG.setname = [subj '_ready_for_ICA_weights'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Attach ICA weight to main data
    EEG = pop_editset(EEG,'run', [], 'icaweights','ALLEEG(1).icaweights', 'icasphere','ALLEEG(1).icasphere');
    % Label ICA components with IC Label Plugin (Pion-Tonachini et al., 2019)
    %%EEG = pop_iclabel(EEG, 'default');
    %%EEG = pop_icflag(EEG, [0 0;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1]);
    % Sanity Check: Plot flagged ICs
    %%tid_psam_plot_flagged_ICs_TD(EEG,['ICs for ' subj], 'SavePath' ,fullfile(OUTPATH, [subj '_ic_topo.png']), 'PlotOn', false)
    % Remove previously flagged ICs (see tid_psam_ica_preprocessing.m)
    EEG.reject.gcompreject = flagged_comps;
    EEG = pop_subcomp( EEG, [], 0);

    % Interpolate bad channels
    if ~isempty(EEG.badchans)
        EEG = pop_interp(EEG, EEG.urchanlocs , 'spherical'); % and interpolate them using urchanlocs
    end

    % Sanity Check: Plot RMS in 10s bins for each electode (not used here)
    %%tid_psam_plot_rms_bins_TD(EEG, [subj ' RMS bins'], 'SavePath', fullfile(OUTPATH, [subj '_channel_rms.png']), 'PlotOn', false)

    % Epoching
    EEG = pop_epoch( EEG, EVENTS, [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');

    % Baseline-Removal
    EEG = pop_rmbase( EEG, [BL_FROM BL_TILL] ,[]);

    % Remove to be excluded trials from ERP data
    EEG.reject.rejglobal = excluded_trials_svm;
    EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);
    EEG.setname = [subj '_erp_clean'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Save dataset
    EEG.setname = [subj '_preprocessed'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    %% Analysis/Plots

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

    % Sanity Check: Plot ERPs & Topographies
    EEG.setname = [subj '_all_conds'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    [~,win_early_start] = min(abs(EEG.times-WIN_EARLY_FROM));
    [~,win_early_end] = min(abs(EEG.times-WIN_EARLY_TILL));

    [~,win_late_start] = min(abs(EEG.times-WIN_LATE_FROM));
    [~,win_late_end] = min(abs(EEG.times-WIN_LATE_TILL));

    [~,t_zero] = min(abs(EEG.times));

    y_lim_lower = -8;
    y_lim_upper = 8;
    cb_lim_lower = -5;
    cb_lim_upper = 2;

    figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.7, 0.7]);
    for cond = 1:2
        EEG = pop_selectevent( ALLEEG(6), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off');
        erp = mean(EEG.data, 3);
        erp_se = std(EEG.data, [], 3) / sqrt(size(EEG.data,3));

        % Store ERP
        if strcmp('go_act', EVENTS{cond})
            all_erp_go_act(:,:,subj_idx) = erp;
        elseif strcmp('go_pas', EVENTS{cond})
            all_erp_go_pas(:,:,subj_idx) = erp;
        end

        subplot(2, 4, 2*(cond-1) + [1 2])
        shadedErrorBar(EEG.times, erp(CHANI,:), erp_se(CHANI,:), ...
            'lineprops', '-b', 'transparent', true);
        xlim([EEG.times(1) EEG.times(end)])
        ylim([y_lim_lower y_lim_upper])
        xlabel('Time [ms]')
        ylabel('Amplitude [µV]')
        title(['ERP for ' EVENTS{cond}])
        hold on
        fill([WIN_LATE_FROM WIN_LATE_TILL WIN_LATE_TILL WIN_LATE_FROM], ...
            [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], main_yellow, 'FaceAlpha',0.1, 'EdgeColor','none');
        fill([WIN_EARLY_FROM WIN_EARLY_TILL WIN_EARLY_TILL WIN_EARLY_FROM], ...
            [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], main_red, 'FaceAlpha',0.1, 'EdgeColor','none');
        hold off

        subplot(2, 4, 4 + (cond-1)*2 + 1)
        topoplot(mean(erp(:,win_early_start:win_early_end),2), EEG.chanlocs, ...
            'emarker2', {CHANI,'o','r',5,1});
        colormap("parula")
        cb = colorbar;
        title(cb, 'Amplitude [µV]')
        clim([cb_lim_lower cb_lim_upper])
        title('Early Window')

        subplot(2, 4, 4 + (cond-1)*2 + 2)
        topoplot(mean(erp(:,win_late_start:win_late_end),2), EEG.chanlocs, ...
            'emarker2', {CHANI,'o','r',5,1});
        colormap("parula")
        cb = colorbar;
        title(cb, 'Amplitude [µV]')
        clim([cb_lim_lower cb_lim_upper])
        title('Late Window')
    end

    sgtitle(['Sanity Check ERPs for ' subj ' (including SE across trials)'])
    saveas(gcf, fullfile(OUTPATH, [subj '_sanity_erp_topo.png']));

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

% Get grandaverage ERPs
grandaverage_erp_act = mean(all_erp_go_act,3);
grandaverage_erp_pas = mean(all_erp_go_pas,3);

% Sanity Check: Plot grandaverage ERPs
figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.7, 0.7]);
subplot(2,4,1:2)
plot(EEG.times, grandaverage_erp_act(CHANI,:,:), 'LineWidth', 1.5, 'Color',main_blue)
xlim([EEG.times(1) EEG.times(end)])
ylim([y_lim_lower y_lim_upper])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
title('Grandaverage ERP for go_act condition')
hold on
fill([WIN_LATE_FROM WIN_LATE_TILL WIN_LATE_TILL WIN_LATE_FROM], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], main_yellow, 'FaceAlpha',0.1, 'EdgeColor','none');
fill([WIN_EARLY_FROM WIN_EARLY_TILL WIN_EARLY_TILL WIN_EARLY_FROM], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], main_red, 'FaceAlpha',0.1, 'EdgeColor','none');
hold off

subplot(2,4,3:4)
plot(EEG.times, grandaverage_erp_pas(CHANI,:,:), 'LineWidth', 1.5, 'Color',main_blue)
xlim([EEG.times(1) EEG.times(end)])
ylim([y_lim_lower y_lim_upper])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
title('Grandaverage ERP for go_pas condition')
hold on
fill([WIN_LATE_FROM WIN_LATE_TILL WIN_LATE_TILL WIN_LATE_FROM], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], main_yellow, 'FaceAlpha',0.1, 'EdgeColor','none');
fill([WIN_EARLY_FROM WIN_EARLY_TILL WIN_EARLY_TILL WIN_EARLY_FROM], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], main_red, 'FaceAlpha',0.1, 'EdgeColor','none');
hold off

subplot(2, 4, 5)
topoplot(mean(grandaverage_erp_act(:,win_early_start:win_early_end),2), EEG.chanlocs, ...
    'emarker2', {CHANI,'o','r',5,1});
colormap("parula")
cb = colorbar;
title(cb, 'Amplitude [µV]')
clim([cb_lim_lower cb_lim_upper])
title('Early Window')

subplot(2, 4, 6)
topoplot(mean(grandaverage_erp_act(:,win_late_start:win_late_end),2), EEG.chanlocs, ...
    'emarker2', {CHANI,'o','r',5,1});
colormap("parula")
cb = colorbar;
title(cb, 'Amplitude [µV]')
clim([cb_lim_lower cb_lim_upper])
title('Late Window')

subplot(2, 4, 7)
topoplot(mean(grandaverage_erp_pas(:,win_early_start:win_early_end),2), EEG.chanlocs, ...
    'emarker2', {CHANI,'o','r',5,1});
colormap("parula")
cb = colorbar;
title(cb, 'Amplitude [µV]')
clim([cb_lim_lower cb_lim_upper])
title('Early Window')

subplot(2, 4, 8)
topoplot(mean(all_erp_go_pas(:,win_late_start:win_late_end),2), EEG.chanlocs, ...
    'emarker2', {CHANI,'o','r',5,1});
colormap("parula")
cb = colorbar;
title(cb, 'Amplitude [µV]')
clim([cb_lim_lower cb_lim_upper])
title('Late Window')

sgtitle('Sanity Check ERPs (grandaverage)')
saveas(gcf,fullfile(OUTPATH, 'grandaverage_sanity_erp_topo.png'));

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'tid_psam_check_filters_svm_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_check_filters_svm_marked_subj.xlsx'))
end

check_done = 'tid_psam_check_filters_svm_DONE'

close all; delete(wb)