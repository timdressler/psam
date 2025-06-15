% tid_psam_check_filters_erp.m
%
% Creates plots to investigate the effects of filterimg based on the ERP-specific prerpocessing (see tid_psam_erp_preprocessing.m, tid_psam_erp_analysis.m).
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
%   Perform preprocesing once with and once without filtering (see tid_psam_erp_preprocessing.m)
%   Excludes markerd trials based on tid_psam_exclude_trials.m
%   Create plots to illustarte effects of filtering
%
% Saves data
%
% Tim Dressler, 15.06.2025

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
INPATH_EXCLUDED = fullfile(MAINPATH, 'data\processed_data\exclude_trials\');
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\check_filters_erp');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_RAW, INPATH_ICA, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
%   Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LCF = 0.3;
HCF = 30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EOG_CHAN = {'E29','E30'}; % Labels of EOG electrodes
EPO_FROM = -0.2;
EPO_TILL = 0.400;
BL_FROM = -200;
SD_PROB = 3;
SD_PROB_ICA = 3;
ALL_EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};

%   Anaylsis/Plots
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ... % Real events
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt'};
CON_EVENTS = {'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};
CHANI = 1; % Channel to plot ERP from
ERP_FROM = 70;
ERP_TILL = 130;

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
% Concatinate colors
colors = {
    main_yellow;
    main_red;
    main_green;
    main_blue;
    };

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

% Initialize overall data storage for filtered and unfiltered
all_erp_filtered = cell(0,3); % Initialized as empty cell matrix with 3 columns
all_con_erp_filtered = cell(0,3); % Initialized as empty cell matrix with 3 columns
all_cor_erp_filtered = cell(0,3); % Initialized as empty cell matrix with 3 columns
all_cor_erp_data_filtered = cell(0,9); % Initialized as empty cell matrix with 9 columns

all_erp_unfiltered = cell(0,3); % Initialized as empty cell matrix with 3 columns
all_con_erp_unfiltered = cell(0,3); % Initialized as empty cell matrix with 3 columns
all_cor_erp_unfiltered = cell(0,3); % Initialized as empty cell matrix with 3 columns
all_cor_erp_data_unfiltered = cell(0,9); % Initialized as empty cell matrix with 9 columns

% Setup progress bar
wb = waitbar(0,'starting tid_psam_check_filters_erp.m');

% Store EEG info for plotting (times and chanlocs from a processed EEG)
EEG_times_for_plotting = [];
EEG_chanlocs_for_plotting = [];


for subj_idx= 1:length(dircont_subj)

    %% Preprocessing
    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_check_filters_erp.m'])

    tic;

    % Store subject-specific ERP data temporarily for individual plots
    % This will hold the data for the current subject for all conditions, for both filtered and unfiltered states
    current_subj_data = struct();

    % Loop for filtered vs. non-filtered processing
    filter_settings = {'filtered', 'unfiltered'}; % First iteration: filtered, Second: unfiltered

    for f_idx = 1:length(filter_settings)
        current_filter_setting = filter_settings{f_idx};
        filter_on = strcmp(current_filter_setting, 'filtered'); % true for 'filtered', false for 'unfiltered'

        % Get file of to-be excluded trials based on tid_psam_exclude_trials.m
        exclusion_filename = fullfile(INPATH_EXCLUDED,[subj '_excluded_trials.mat']);
        load(exclusion_filename) % loads variables 'excluded_trials_erp_beh' (used here) and 'excluded_trials_svm' (not used here)

        % Start eeglab (clears ALLEEG, EEG, CURRENTSET, ALLCOM)
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
        for event = 1:length(EEG.event)
            switch EEG.event(event).type
                case 'S 931'
                    EEG.event(event).type = 'act_early_unalt';
                case 'S 932'
                    EEG.event(event).type = 'act_early_alt';
                case 'S 933'
                    EEG.event(event).type = 'act_late_unalt';
                case 'S 934'
                    EEG.event(event).type = 'act_late_alt';
                case 'S 941'
                    EEG.event(event).type = 'pas_early_unalt';
                case 'S 942'
                    EEG.event(event).type = 'pas_early_alt';
                case 'S 943'
                    EEG.event(event).type = 'pas_late_unalt';
                case 'S 944'
                    EEG.event(event).type = 'pas_late_alt';
            end
        end

        EEG.setname = [subj '_ready_for_preprocessing'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Apply Filter (Conditional based on filter_on)
        if filter_on
            % Highpass-Filter
            LCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth(LCF)); % Get filter order (also see pop_firwsord.m, tid_psam_get_transition_bandwidth.m)
            EEG = pop_firws(EEG, 'fcutoff', LCF, 'ftype', 'highpass', 'wtype', 'hamming', 'forder',LCF_ord, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
            % Lowpass-Filter
            HCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth(HCF)); % Get filter order (also see pop_firwsord.m, tid_psam_get_transition_bandwidth.m)
            EEG = pop_firws(EEG, 'fcutoff', HCF, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder',HCF_ord, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

            %%EEG = pop_firws(EEG, 'fcutoff', LCF, 'ftype', 'highpass', 'wtype', 'hamming', 'forder',(2 * ceil((3*(EEG.srate/LCF)) / 2)), 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
            %%EEG = pop_firws(EEG, 'fcutoff', HCF, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder',(2 * ceil((3*(EEG.srate/HCF)) / 2)), 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

            %%EEG = pop_eegfiltnew(EEG, 'locutoff',LCF,'hicutoff',0);
            %%EEG = pop_eegfiltnew(EEG, 'hicutoff',HCF,'plotfreqz',0);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        % Sanity Check: Plot RMS in 10s bins for each electode
        tid_psam_plot_rms_bins_TD(EEG, [subj ' RMS bins'], 'SavePath', fullfile(OUTPATH, [subj '_channel_rms.png']), 'PlotOn', false)

        % Epoching
        EEG = pop_epoch( EEG, ALL_EVENTS, [EPO_FROM         EPO_TILL], 'epochinfo', 'yes');

        % Baseline-Removal
        EEG = pop_rmbase( EEG, [BL_FROM 0] ,[]);

        % Remove to be excluded trials from ERP data
        EEG.reject.rejglobal = excluded_trials_erp_beh;
        EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);
        EEG.setname = [subj '_erp_clean'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        % Remove EOG channels
        EEG = pop_select( EEG, 'rmchannel',{'E29','E30'});

        % Store data (final EEG for current filter setting)
        EEG.setname = [subj '_' current_filter_setting '_all_conds'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        % Store EEG times and channel locations for plotting (only need to do this once)
        if isempty(EEG_times_for_plotting)
            EEG_times_for_plotting = EEG.times;
            EEG_chanlocs_for_plotting = EEG.chanlocs;
        end

        %% Analysis/Plots for current filter setting
        % Corrected initialization for temp_erp_store and temp_cor_erp_store
        temp_erp_store = cell(length(EVENTS), 3); % Holds erp, subj, event
        temp_con_erp_store = cell(length(CON_EVENTS), 3); % Holds con_erp, subj, event
        temp_cor_erp_store = cell(length(EVENTS), 3); % Holds cor_erp, subj, event
        temp_cor_erp_data_store = cell(length(EVENTS), 9);

        temp_con_counter_inner = 1; % For unique control ERPs within this filter setting

        % Get ERPs for each condition
        for cond = 1:length(EVENTS)
            EEG_cond = pop_selectevent( ALLEEG(end), 'latency','-2<=2','type',EVENTS(cond), ...
                'deleteevents','off','deleteepochs','on','invertepochs','off');

            % Get ERP
            erp = mean(EEG_cond.data, 3);

            % Correct ERP by subtracting the ERP from the equivalent no-probe condition (see Daliri & Max, 2016)
            % Get matching control condition
            con_cond_idx = find(strcmp(erase(EVENTS{cond}, {'_alt', '_unalt'}), erase(string(CON_EVENTS), 'con_')));
            % Get control ERP
            EEG_con = pop_selectevent( ALLEEG(end), 'latency','-2<=2','type',CON_EVENTS(con_cond_idx), ...
                'deleteevents','off','deleteepochs','on','invertepochs','off');
            con_erp = mean(EEG_con.data, 3);
            % Get corrected ERP
            cor_erp = erp - con_erp;

            % Setup ERP analysis
            [~,win_start] = min(abs(EEG.times-ERP_FROM));
            [~,win_end] = min(abs(EEG.times-ERP_TILL));
            %[~,t_zero] = min(abs(EEG.times)); % not used

            % ERP analysis (control)
            % Store ERP (control)
            % Only store control ERP once per unique control event per filter setting
            % Check if this control condition has already been stored for this filter setting and subject
            is_stored = false;
            for k = 1:temp_con_counter_inner-1
                if strcmp(temp_con_erp_store{k,3}, CON_EVENTS{con_cond_idx})
                    is_stored = true;
                    break;
                end
            end

            if ~is_stored
                temp_con_erp_store{temp_con_counter_inner,1} = con_erp;
                temp_con_erp_store{temp_con_counter_inner,2} = subj;
                temp_con_erp_store{temp_con_counter_inner,3} = CON_EVENTS{con_cond_idx};
                temp_con_counter_inner = temp_con_counter_inner+1;
            end


            % ERP analysis (uncorrected)
            % Get N100 amplitude (uncorrected)
            erp_amp = min(erp(CHANI,win_start:win_end));
            % Get N100 latency (uncorrected)
            erp_sam = find(erp(CHANI,:) == erp_amp);
            erp_lat = EEG.times(erp_sam(1)); % Take first if multiple samples have the same min amp
            % Store ERP
            temp_erp_store{cond,1} = erp;
            temp_erp_store{cond,2} = subj;
            temp_erp_store{cond,3} = EVENTS{cond};

            % ERP analysis (corrected)
            % Get N100 amplitude (corrected)
            cor_erp_amp = min(cor_erp(CHANI,win_start:win_end));
            % Get N100 latency (corrected)
            cor_erp_sam = find(cor_erp(CHANI,:) == cor_erp_amp);
            cor_erp_lat = EEG.times(cor_erp_sam(1)); % Take first if multiple samples have the same min amp
            % Store ERP (corrected)
            temp_cor_erp_store{cond,1} = cor_erp;
            temp_cor_erp_store{cond,2} = subj;
            temp_cor_erp_store{cond,3} = EVENTS{cond};
            % Store ERP data (corrected)
            temp_cor_erp_data_store{cond,1} = subj;
            temp_cor_erp_data_store{cond,2} = ERP_FROM;
            temp_cor_erp_data_store{cond,3} = ERP_TILL;
            temp_cor_erp_data_store{cond,4} = cor_erp_amp;
            temp_cor_erp_data_store{cond,5} = cor_erp_lat;
            temp_cor_erp_data_store{cond,6} = EVENTS{cond};

            % Split condition labels into multiple columns (corrected)
            cond_parts = strsplit(EVENTS{cond}, '_');
            % Map and rename task condition
            switch cond_parts{1}
                case 'act'
                    task_label = 'Active';
                case 'pas'
                    task_label = 'Passive';
                otherwise
                    error('Invalid label')
            end
            % Map and rename probe-onset condition (corrected)
            switch cond_parts{2}
                case 'early'
                    probe_onset_label = 'Early';
                case 'late'
                    probe_onset_label = 'Late';
                otherwise
                    error('Invalid label')
            end
            % Map and rename probe-type condition (corrected)
            switch cond_parts{3}
                case 'alt'
                    probe_type_label = 'Altered';
                case 'unalt'
                    probe_type_label = 'Unaltered';
                otherwise
                    error('Invalid label')
            end

            % Store mapped and renamed values (corrected)
            temp_cor_erp_data_store{cond, 7} = task_label;
            temp_cor_erp_data_store{cond, 8} = probe_onset_label;
            temp_cor_erp_data_store{cond, 9} = probe_type_label;

        end % End of cond loop for data extraction

        % Store the current subject's data for this filter setting
        if filter_on
            current_subj_data.filtered.erp = temp_erp_store;
            current_subj_data.filtered.con_erp = temp_con_erp_store(1:temp_con_counter_inner-1,:); % Only store actual entries
            current_subj_data.filtered.cor_erp = temp_cor_erp_store;
            current_subj_data.filtered.cor_erp_data = temp_cor_erp_data_store;

            % Append to global filtered storage
            all_erp_filtered = [all_erp_filtered; temp_erp_store];
            all_con_erp_filtered = [all_con_erp_filtered; temp_con_erp_store(1:temp_con_counter_inner-1,:)];
            all_cor_erp_filtered = [all_cor_erp_filtered; temp_cor_erp_store];
            all_cor_erp_data_filtered = [all_cor_erp_data_filtered; temp_cor_erp_data_store];

        else
            current_subj_data.unfiltered.erp = temp_erp_store;
            current_subj_data.unfiltered.con_erp = temp_con_erp_store(1:temp_con_counter_inner-1,:);
            current_subj_data.unfiltered.cor_erp = temp_cor_erp_store;
            current_subj_data.unfiltered.cor_erp_data = temp_cor_erp_data_store;

            % Append to global unfiltered storage
            all_erp_unfiltered = [all_erp_unfiltered; temp_erp_store];
            all_con_erp_unfiltered = [all_con_erp_unfiltered; temp_con_erp_store(1:temp_con_counter_inner-1,:)];
            all_cor_erp_unfiltered = [all_cor_erp_unfiltered; temp_cor_erp_store];
            all_cor_erp_data_unfiltered = [all_cor_erp_data_unfiltered; temp_cor_erp_data_store];
        end

    end % End of filter_settings loop

    % Individual plots (moved out of initial conds loop, now using stored subj data)
    % Define limits
    ylim_upper = 25;
    ylim_lower = -15;
    cb_lim_upper = 5;
    cb_lim_lower = -5;

    individual_plot = figure('Units', 'normalized', 'Position', [0 0 1 1]);

    % Plot rows
    r_erp_unc = 1; r_topo_unc = 2;
    r_erp_con = 3; r_topo_con = 4;
    r_erp_cor = 5; r_topo_cor = 6;

    % New conds loop for individual plots
    for cond_plot_idx = 1:length(EVENTS)

        % Retrieve data for current subject and condition for both filtered and unfiltered
        % Filtered data
        erp_filtered = current_subj_data.filtered.erp{cond_plot_idx, 1};
        % Need to find the correct control ERP based on the condition index
        con_cond_idx_for_plot = find(strcmp(erase(EVENTS{cond_plot_idx}, {'_alt', '_unalt'}), erase(string(CON_EVENTS), 'con_')));
        con_erp_filtered_row_idx = cellfun(@(x) strcmp(x, CON_EVENTS{con_cond_idx_for_plot}), current_subj_data.filtered.con_erp(:,3));
        con_erp_filtered = current_subj_data.filtered.con_erp{con_erp_filtered_row_idx,1};
        cor_erp_filtered = current_subj_data.filtered.cor_erp{cond_plot_idx, 1};

        erp_amp_filtered = current_subj_data.filtered.cor_erp_data{cond_plot_idx, 4};
        erp_lat_filtered = current_subj_data.filtered.cor_erp_data{cond_plot_idx, 5};
        cor_erp_amp_filtered = current_subj_data.filtered.cor_erp_data{cond_plot_idx, 4};
        cor_erp_lat_filtered = current_subj_data.filtered.cor_erp_data{cond_plot_idx, 5};

        % Unfiltered data
        erp_unfiltered = current_subj_data.unfiltered.erp{cond_plot_idx, 1};
        con_erp_unfiltered_row_idx = cellfun(@(x) strcmp(x, CON_EVENTS{con_cond_idx_for_plot}), current_subj_data.unfiltered.con_erp(:,3));
        con_erp_unfiltered = current_subj_data.unfiltered.con_erp{con_erp_unfiltered_row_idx,1};
        cor_erp_unfiltered = current_subj_data.unfiltered.cor_erp{cond_plot_idx, 1};

        erp_amp_unfiltered = current_subj_data.unfiltered.cor_erp_data{cond_plot_idx, 4};
        erp_lat_unfiltered = current_subj_data.unfiltered.cor_erp_data{cond_plot_idx, 5};
        cor_erp_amp_unfiltered = current_subj_data.unfiltered.cor_erp_data{cond_plot_idx, 4};
        cor_erp_lat_unfiltered = current_subj_data.unfiltered.cor_erp_data{cond_plot_idx, 5};

        % Recalculate win_start and win_end as EEG might be cleared, but EEG_times_for_plotting is set
        [~,win_start] = min(abs(EEG_times_for_plotting-ERP_FROM));
        [~,win_end] = min(abs(EEG_times_for_plotting-ERP_TILL));


        % Plot: Uncorrected ERP and Topoplot
        subplot(6, 8, (r_erp_unc - 1) * 8 + cond_plot_idx)
        plot(EEG_times_for_plotting, erp_filtered(CHANI,:), 'LineWidth', 1.5, 'Color', colors{4}); % Filtered
        hold on
        plot(EEG_times_for_plotting, erp_unfiltered(CHANI,:), '--', 'LineWidth', 1.5, 'Color', colors{4}); % Unfiltered (dashed)
        % Removed scatter points for N100
        % scatter(erp_lat_filtered, erp_amp_filtered, 'filled', 'MarkerEdgeColor', colors{4}, 'MarkerFaceColor', colors{4});
        % scatter(erp_lat_unfiltered, erp_amp_unfiltered, 'o', 'MarkerEdgeColor', colors{4}, 'MarkerFaceColor', 'none');
        hold off
        title(EVENTS{cond_plot_idx}, 'Interpreter', 'none')
        ylim([ylim_lower ylim_upper]); xlim([-200 400])
        if cond_plot_idx == 1
            ylabel('Uncorrected');
            % Add a legend for the first subplot
            leg = legend({'Filtered', 'Unfiltered'}, 'Location', 'northwest');
            set(leg, 'Interpreter', 'none');
        end


        subplot(6, 8, (r_topo_unc - 1) * 8 + cond_plot_idx)
        topoplot(mean(erp_filtered(:,win_start:win_end),2), EEG_chanlocs_for_plotting, ...
            'emarker2', {CHANI,'o','r',2,2},'emarker', {'.','k',0.1,1});
        colormap("parula")
        colorbar;
        clim([cb_lim_lower cb_lim_upper])

        % Plot: Control ERP and Topoplot
        subplot(6, 8, (r_erp_con - 1) * 8 + cond_plot_idx)
        plot(EEG_times_for_plotting, con_erp_filtered(CHANI,:), 'LineWidth', 1.5, 'Color', colors{4});
        hold on
        plot(EEG_times_for_plotting, con_erp_unfiltered(CHANI,:), '--', 'LineWidth', 1.5, 'Color', colors{4});
        hold off
        ylim([ylim_lower ylim_upper]); xlim([-200 400])
        if cond_plot_idx == 1
            ylabel('Control');
            % Add a legend for the first subplot
            leg = legend({'Filtered', 'Unfiltered'}, 'Location', 'northwest');
            set(leg, 'Interpreter', 'none');
        end


        subplot(6, 8, (r_topo_con - 1) * 8 + cond_plot_idx)
        topoplot(mean(con_erp_filtered(:,win_start:win_end),2), EEG_chanlocs_for_plotting, ...
            'emarker2', {CHANI,'o','r',2,2},'emarker', {'.','k',0.1,1});
        colormap("parula")
        colorbar;
        clim([cb_lim_lower cb_lim_upper])

        % Plot: Corrected ERP and Topoplot
        subplot(6, 8, (r_erp_cor - 1) * 8 + cond_plot_idx)
        plot(EEG_times_for_plotting, cor_erp_filtered(CHANI,:), 'LineWidth', 1.5, 'Color', colors{4});
        hold on
        plot(EEG_times_for_plotting, cor_erp_unfiltered(CHANI,:), '--', 'LineWidth', 1.5, 'Color', colors{4});
        % Removed scatter points for N100
        % scatter(cor_erp_lat_filtered, cor_erp_amp_filtered, 'filled', 'MarkerEdgeColor', colors{4}, 'MarkerFaceColor', colors{4});
        % scatter(cor_erp_lat_unfiltered, cor_erp_amp_unfiltered, 'o', 'MarkerEdgeColor', colors{4}, 'MarkerFaceColor', 'none');
        hold off
        ylim([ylim_lower ylim_upper]); xlim([-200 400])
        if cond_plot_idx == 1
            ylabel('Corrected');
            % Add a legend for the first subplot
            leg = legend({'Filtered', 'Unfiltered'}, 'Location', 'northwest');
            set(leg, 'Interpreter', 'none');
        end


        subplot(6, 8, (r_topo_cor - 1) * 8 + cond_plot_idx)
        topoplot(mean(cor_erp_filtered(:,win_start:win_end),2), EEG_chanlocs_for_plotting, ...
            'emarker2', {CHANI,'o','r',2,2},'emarker', {'.','k',0.1,1});
        colormap("parula")
        colorbar;
        clim([cb_lim_lower cb_lim_upper])

    end % End of individual plot loop (new conds loop)

    sgtitle(['Sanity Check ERPs and Topoplots for ' subj], 'Interpreter', 'none');
    saveas(individual_plot, fullfile(OUTPATH, [subj '_sanity_erp_topo.png']));
    close(individual_plot)

    % Update Protocol
    subj_time = toc;
    protocol{subj_idx,1} = subj;
    protocol{subj_idx,2} = subj_time;
    if any(strcmp(marked_subj, subj), 'all')
        protocol{subj_idx,3} = 'MARKED';
    else
        protocol{subj_idx,3} = 'OK';
    end

end % End of subj_idx loop

% Get grandaverage ERPs for filtered data
erp_data_col = 1;
cond_col = 3;

% Ensure conditions are character vectors before passing to unique
conditions_to_check_filtered = all_erp_filtered(:, cond_col);
is_char_cell_filtered = cellfun(@ischar, conditions_to_check_filtered);
clean_conditions_filtered = conditions_to_check_filtered(is_char_cell_filtered);
[unique_conditions_filtered, ~, condition_idx_filtered] = unique(clean_conditions_filtered);

n_conds_filtered = numel(unique_conditions_filtered);
grandaverage_erp_filtered = cell(n_conds_filtered, 2);
for i = 1:n_conds_filtered
    erp_group = all_erp_filtered(condition_idx_filtered == i, erp_data_col);
    erp_stack = cat(3, erp_group{:});
    grandaverage_erp_filtered{i, 1} = mean(erp_stack, 3);
    grandaverage_erp_filtered{i, 2} = unique_conditions_filtered{i};
end

% Get grandaverage control ERP for each condition (filtered)
conditions_to_check_con_filtered = all_con_erp_filtered(:, cond_col);
is_char_cell_con_filtered = cellfun(@ischar, conditions_to_check_con_filtered);
clean_conditions_con_filtered = conditions_to_check_con_filtered(is_char_cell_con_filtered);
[unique_conditions_con_filtered, ~, condition_idx_con_filtered] = unique(clean_conditions_con_filtered);

n_conds_con_filtered = numel(unique_conditions_con_filtered);
grandaverage_erp_con_filtered = cell(n_conds_con_filtered, 2);
for i = 1:n_conds_con_filtered
    erp_group = all_con_erp_filtered(condition_idx_con_filtered == i, erp_data_col);
    erp_stack = cat(3, erp_group{:});
    grandaverage_erp_con_filtered{i, 1} = mean(erp_stack, 3);
    grandaverage_erp_con_filtered{i, 2} = unique_conditions_con_filtered{i};
end

% Get grandaverage corrected ERP for each condition (filtered)
conditions_to_check_cor_filtered = all_cor_erp_filtered(:, cond_col);
is_char_cell_cor_filtered = cellfun(@ischar, conditions_to_check_cor_filtered);
clean_conditions_cor_filtered = conditions_to_check_cor_filtered(is_char_cell_cor_filtered);
[unique_conditions_cor_filtered, ~, condition_idx_cor_filtered] = unique(clean_conditions_cor_filtered);

n_conds_cor_filtered = numel(unique_conditions_cor_filtered);
grandaverage_erp_cor_filtered = cell(n_conds_cor_filtered, 2);
for i = 1:n_conds_cor_filtered
    erp_group = all_cor_erp_filtered(condition_idx_cor_filtered == i, erp_data_col);
    erp_stack = cat(3, erp_group{:});
    grandaverage_erp_cor_filtered{i, 1} = mean(erp_stack, 3);
    grandaverage_erp_cor_filtered{i, 2} = unique_conditions_cor_filtered{i};
end

% Get grandaverage ERPs for unfiltered data
% Get grandaverage uncorrected ERP for each condition (unfiltered)
conditions_to_check_unfiltered = all_erp_unfiltered(:, cond_col);
is_char_cell_unfiltered = cellfun(@ischar, conditions_to_check_unfiltered);
clean_conditions_unfiltered = conditions_to_check_unfiltered(is_char_cell_unfiltered);
[unique_conditions_unfiltered, ~, condition_idx_unfiltered] = unique(clean_conditions_unfiltered);

n_conds_unfiltered = numel(unique_conditions_unfiltered);
grandaverage_erp_unfiltered = cell(n_conds_unfiltered, 2);
for i = 1:n_conds_unfiltered
    erp_group = all_erp_unfiltered(condition_idx_unfiltered == i, erp_data_col);
    erp_stack = cat(3, erp_group{:});
    grandaverage_erp_unfiltered{i, 1} = mean(erp_stack, 3);
    grandaverage_erp_unfiltered{i, 2} = unique_conditions_unfiltered{i};
end

% Get grandaverage control ERP for each condition (unfiltered)
conditions_to_check_con_unfiltered = all_con_erp_unfiltered(:, cond_col);
is_char_cell_con_unfiltered = cellfun(@ischar, conditions_to_check_con_unfiltered);
clean_conditions_con_unfiltered = conditions_to_check_con_unfiltered(is_char_cell_con_unfiltered);
[unique_conditions_con_unfiltered, ~, condition_idx_con_unfiltered] = unique(clean_conditions_con_unfiltered);

n_conds_con_unfiltered = numel(unique_conditions_con_unfiltered);
grandaverage_erp_con_unfiltered = cell(n_conds_con_unfiltered, 2);
for i = 1:n_conds_con_unfiltered
    erp_group = all_con_erp_unfiltered(condition_idx_con_unfiltered == i, erp_data_col);
    erp_stack = cat(3, erp_group{:});
    grandaverage_erp_con_unfiltered{i, 1} = mean(erp_stack, 3);
    grandaverage_erp_con_unfiltered{i, 2} = unique_conditions_con_unfiltered{i};
end

% Get grandaverage corrected ERP for each condition (unfiltered)
conditions_to_check_cor_unfiltered = all_cor_erp_unfiltered(:, cond_col);
is_char_cell_cor_unfiltered = cellfun(@ischar, conditions_to_check_cor_unfiltered);
clean_conditions_cor_unfiltered = conditions_to_check_cor_unfiltered(is_char_cell_cor_unfiltered);
[unique_conditions_cor_unfiltered, ~, condition_idx_cor_unfiltered] = unique(clean_conditions_cor_unfiltered);

n_conds_cor_unfiltered = numel(unique_conditions_cor_unfiltered);
grandaverage_erp_cor_unfiltered = cell(n_conds_cor_unfiltered, 2);
for i = 1:n_conds_cor_unfiltered
    erp_group = all_cor_erp_unfiltered(condition_idx_cor_unfiltered == i, erp_data_col);
    erp_stack = cat(3, erp_group{:});
    grandaverage_erp_cor_unfiltered{i, 1} = mean(erp_stack, 3);
    grandaverage_erp_cor_unfiltered{i, 2} = unique_conditions_cor_unfiltered{i};
end


% Save ERPs (control, uncorrected, corrected)
% Save filtered and unfiltered separately, or together in a struct
save(fullfile(OUTPATH, 'all_subj_erp_filtered.mat'),'all_cor_erp_filtered')
save(fullfile(OUTPATH, 'all_subj_erp_unfiltered.mat'),'all_cor_erp_unfiltered')


% Save ERP data
all_cor_erp_data_filtered = cell2table(all_cor_erp_data_filtered, 'VariableNames',{'subj','erp_from', 'erp_till', 'erp_amp', 'erp_lat', 'condition_full', 'task_instruction', 'probe_onset_cat', 'probe_type'});
writetable(all_cor_erp_data_filtered, fullfile(OUTPATH, 'all_subj_cor_erp_data_filtered.xlsx'))
all_cor_erp_data_unfiltered = cell2table(all_cor_erp_data_unfiltered, 'VariableNames',{'subj','erp_from', 'erp_till', 'erp_amp', 'erp_lat', 'condition_full', 'task_instruction', 'probe_onset_cat', 'probe_type'});
writetable(all_cor_erp_data_unfiltered, fullfile(OUTPATH, 'all_subj_cor_erp_data_unfiltered.xlsx'))


% Save grandaverage ERP
save(fullfile(OUTPATH, 'grandaverage_cor_erp_filtered.mat'),'grandaverage_erp_cor_filtered')
save(fullfile(OUTPATH, 'grandaverage_cor_erp_unfiltered.mat'),'grandaverage_erp_cor_unfiltered')


% Preparation for plots
% Get channel ID
chani = CHANI;
% Get values for dynamic plot limits
y_lim_lower = min(min(cellfun(@(x) min(x(:)), grandaverage_erp_cor_filtered(:,1))), min(cellfun(@(x) min(x(:)), grandaverage_erp_cor_unfiltered(:,1))))-1;
y_lim_upper = max(max(cellfun(@(x) max(x(:)), grandaverage_erp_cor_filtered(:,1))), max(cellfun(@(x) max(x(:)), grandaverage_erp_cor_unfiltered(:,1))))+1;

% Rename condition labels
rename_conditions_map = {
    'act_early_unalt', 'Active - Early - Unaltered';
    'act_early_alt',   'Active - Early - Altered';
    'act_late_unalt',  'Active - Late - Unaltered';
    'act_late_alt',    'Active - Late - Altered';
    'pas_early_unalt', 'Passive - Early - Unaltered';
    'pas_early_alt',   'Passive - Early - Altered';
    'pas_late_unalt',  'Passive - Late - Unaltered';
    'pas_late_alt',    'Passive - Late - Altered';
    };

for i = 1:size(grandaverage_erp_filtered,1) % for uncorrected ERP filtered
    old_label = grandaverage_erp_filtered{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp_filtered{i,2} = new_label{1};
end
for i = 1:size(grandaverage_erp_cor_filtered,1) % for corrected ERP filtered
    old_label = grandaverage_erp_cor_filtered{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp_cor_filtered{i,2} = new_label{1};
end

for i = 1:size(grandaverage_erp_unfiltered,1) % for uncorrected ERP unfiltered
    old_label = grandaverage_erp_unfiltered{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp_unfiltered{i,2} = new_label{1};
end
for i = 1:size(grandaverage_erp_cor_unfiltered,1) % for corrected ERP unfiltered
    old_label = grandaverage_erp_cor_unfiltered{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp_cor_unfiltered{i,2} = new_label{1};
end


% Split data into early and late conditions for filtered and unfiltered
grandaverage_erp_cor_early_filtered = grandaverage_erp_cor_filtered(contains(grandaverage_erp_cor_filtered(:,2), 'Early'), :);
grandaverage_erp_cor_late_filtered = grandaverage_erp_cor_filtered(contains(grandaverage_erp_cor_filtered(:,2), 'Late'), :);

grandaverage_erp_cor_early_unfiltered = grandaverage_erp_cor_unfiltered(contains(grandaverage_erp_cor_unfiltered(:,2), 'Early'), :);
grandaverage_erp_cor_late_unfiltered = grandaverage_erp_cor_unfiltered(contains(grandaverage_erp_cor_unfiltered(:,2), 'Late'), :);


% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'});
writetable(protocol,fullfile(OUTPATH, 'tid_psam_check_filters_erp_protocol.xlsx'));

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'});
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_check_filters_erp_marked_subj.xlsx'));
end

check_done = 'tid_psam_check_filters_erp_DONE';

close all; delete(wb);

%% Plots
close all

% Plot: ERPs for all conditions with Topoplots (Combined filtered and unfiltered)
% Create ERP plot
figure('Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
sgtitle('Corrected ERPs and Topoplots for all Conditions (Filtered vs. Unfiltered)')
subplot(3,4,5:6)
for cond_early_num = 1:length(grandaverage_erp_cor_early_filtered)
    erp_2_plot_filtered = grandaverage_erp_cor_early_filtered{cond_early_num,1};
    erp_2_plot_unfiltered = grandaverage_erp_cor_early_unfiltered{cond_early_num,1};

    plot(EEG_times_for_plotting, erp_2_plot_filtered(chani,:), 'LineWidth', 1.5, ...
        'Color', colors{cond_early_num}, 'DisplayName', [grandaverage_erp_cor_early_filtered{cond_early_num,2} ' (Filtered)']);
    hold on
    plot(EEG_times_for_plotting, erp_2_plot_unfiltered(chani,:), '--', 'LineWidth', 1.5, ...
        'Color', colors{cond_early_num}, 'DisplayName', [grandaverage_erp_cor_early_unfiltered{cond_early_num,2} ' (Unfiltered)']);
end
xlim([-200 400])
ylim([y_lim_lower y_lim_upper])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
title('ERPs for early probes')
fill([70 130 130 70], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], light_blue, 'FaceAlpha',0.1, 'EdgeColor','none');
l = legend('Location', 'northwest'); % legend should dynamically pick up DisplayNames
fontsize(l,10,'points')
hold off

subplot(3,4,7:8)
for cond_late_num = 1:length(grandaverage_erp_cor_late_filtered)
    erp_cor_2_plot_filtered = grandaverage_erp_cor_late_filtered{cond_late_num,1};
    erp_cor_2_plot_unfiltered = grandaverage_erp_cor_late_unfiltered{cond_late_num,1};

    plot(EEG_times_for_plotting, erp_cor_2_plot_filtered(chani,:), 'LineWidth', 1.5, ...
        'Color', colors{cond_late_num}, 'DisplayName', [grandaverage_erp_cor_late_filtered{cond_late_num,2} ' (Filtered)']);
    hold on
    plot(EEG_times_for_plotting, erp_cor_2_plot_unfiltered(chani,:), '--', 'LineWidth', 1.5, ...
        'Color', colors{cond_late_num}, 'DisplayName', [grandaverage_erp_cor_late_unfiltered{cond_late_num,2} ' (Unfiltered)']);
end
xlim([-200 400])
ylim([y_lim_lower y_lim_upper])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
title('ERPs for late probes')
fill([70 130 130 70], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], light_blue, 'FaceAlpha',0.1, 'EdgeColor','none');
l = legend('Location', 'northwest');
fontsize(l,10,'points')
hold off

% Create Topoplots (using filtered data for topo, as unfiltered would be difficult to represent comparatively on one topo)
% Get limit values
% Combine filtered and unfiltered topo data to get overall min/max for colorbar
all_topo_2_plot_combined = [];
% Recalculate win_start and win_end here as well to be safe
[~,win_start] = min(abs(EEG_times_for_plotting-ERP_FROM));
[~,win_end] = min(abs(EEG_times_for_plotting-ERP_TILL));

for cond_num = 1:length(grandaverage_erp_cor_filtered)
    topo_2_plot_filtered = grandaverage_erp_cor_filtered{cond_num,1};
    all_topo_2_plot_combined(:,end+1) = mean(topo_2_plot_filtered(:,win_start:win_end),2);

    topo_2_plot_unfiltered = grandaverage_erp_cor_unfiltered{cond_num,1};
    all_topo_2_plot_combined(:,end+1) = mean(topo_2_plot_unfiltered(:,win_start:win_end),2);
end
cb_lim_lower = min(all_topo_2_plot_combined, [], 'all');
cb_lim_upper = max(all_topo_2_plot_combined, [], 'all');


% Set up indices for Topoplots in Subplot
topo_idx = [1 2 3 4 9 10 11 12];

% Create Topoplot (Displaying filtered data, as per original script, no direct way to show both on one topo)
for cond_num = 1:length(grandaverage_erp_cor_filtered)
    subplot(3,4,topo_idx(cond_num))
    topo_2_plot = grandaverage_erp_cor_filtered{cond_num,1}; % Using filtered data for topoplot
    topoplot(mean(topo_2_plot(:,win_start:win_end),2), EEG_chanlocs_for_plotting, 'emarker2', {chani,'o','r',5,1});
    title(grandaverage_erp_cor_filtered{cond_num,2})
    colormap("parula")
    cb = colorbar;
    clim([cb_lim_lower cb_lim_upper]);
    cb.Label.String = 'Amplitude [µV]';
    cb.Label.Position = [2.2, -1, 0];
    cb.Label.HorizontalAlignment = 'center';
end

% Save plot
saveas(gcf,fullfile(OUTPATH, 'tid_psam_plot_grandaverage_erp_topo_combined.png'))

% Plot: ERP, control ERP and correct ERP for each condition (Combined filtered and unfiltered)
% Concatinate colors for the three ERP types (Uncorrected, Corrected, Control)
colors_erp_types = {
    main_blue; % Uncorrected (filtered)
    main_red;  % Corrected (filtered)
    main_green; % Control (filtered)
    };

figure('Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
sgtitle('Uncorrected, Corrected and Control ERPs for all Conditions (Filtered vs. Unfiltered)')

con_cond_num_filtered = 1;
con_cond_num_unfiltered = 1;

for cond_num = 1:length(grandaverage_erp_cor_filtered) % Loop through conditions
    % Filtered ERPs
    erp_2_plot_filtered = grandaverage_erp_filtered{cond_num,1};
    erp_cor_2_plot_filtered = grandaverage_erp_cor_filtered{cond_num,1};
    erp_con_2_plot_filtered = grandaverage_erp_con_filtered{con_cond_num_filtered,1}; % Using specific counter for controls

    % Unfiltered ERPs
    erp_2_plot_unfiltered = grandaverage_erp_unfiltered{cond_num,1};
    erp_cor_2_plot_unfiltered = grandaverage_erp_cor_unfiltered{cond_num,1};
    erp_con_2_plot_unfiltered = grandaverage_erp_con_unfiltered{con_cond_num_unfiltered,1}; % Using specific counter for controls

    subplot(2,4,cond_num)
    % Plot filtered data
    plot(EEG_times_for_plotting, erp_2_plot_filtered(chani,:), 'LineWidth', 1.5, 'Color', colors_erp_types{1}, 'DisplayName', 'Uncorrected (Filtered)')
    hold on
    plot(EEG_times_for_plotting, erp_cor_2_plot_filtered(chani,:), 'LineWidth', 1.5, 'Color', colors_erp_types{2}, 'DisplayName', 'Corrected (Filtered)')
    plot(EEG_times_for_plotting, erp_con_2_plot_filtered(chani,:), 'LineWidth', 1.5, 'Color', colors_erp_types{3}, 'DisplayName', 'Control (Filtered)')

    % Plot unfiltered data (dashed)
    plot(EEG_times_for_plotting, erp_2_plot_unfiltered(chani,:), '--', 'LineWidth', 1.5, 'Color', colors_erp_types{1}, 'DisplayName', 'Uncorrected (Unfiltered)')
    plot(EEG_times_for_plotting, erp_cor_2_plot_unfiltered(chani,:), '--', 'LineWidth', 1.5, 'Color', colors_erp_types{2}, 'DisplayName', 'Corrected (Unfiltered)')
    plot(EEG_times_for_plotting, erp_con_2_plot_unfiltered(chani,:), '--', 'LineWidth', 1.5, 'Color', colors_erp_types{3}, 'DisplayName', 'Control (Unfiltered)')

    xlim([-200 400])
    ylim([y_lim_lower-1 y_lim_upper+1])
    xlabel('Time [ms]')
    ylabel('Amplitude [µV]')
    title(grandaverage_erp_cor_filtered{cond_num,2}) % Title based on filtered condition name
    l = legend('Location', 'northwest'); % legend should dynamically pick up DisplayNames
    fontsize(l,8,'points')
    hold off

    % Increment control condition counters appropriately
    % There are 4 unique control conditions and 8 main conditions (2 main conditions for each control)
    % So, increment every two main conditions.
    if mod(cond_num,2) == 0
        con_cond_num_filtered = con_cond_num_filtered +1;
        con_cond_num_unfiltered = con_cond_num_unfiltered +1;
    end
end

% Save plot
saveas(gcf,fullfile(OUTPATH, 'tid_psam_plot_grandaverage_erp_con_cor_combined.png'))
