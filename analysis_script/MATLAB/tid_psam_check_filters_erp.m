% tid_psam_check_filters_erp.m
%
% Creates plots to investigate the effects of filterimg based on the ERP-specific preprocessing (see tid_psam_erp_preprocessing.m, tid_psam_erp_analysis.m).
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
%   Create plots to illustarte effects of filtering (see tid_psam_erp_analsis.m).
%
% Saves plots
%
% Tim Dressler, 15.06.2025

clear
close all
clc

rng(123)
set(0,'DefaultTextInterpreter','none')

% Set up paths
SCRIPTPATH = cd;
normalizedPath = strrep(SCRIPTPATH, filesep, '/');
expectedSubpath = 'psam/analysis_script/MATLAB';

if contains(normalizedPath, expectedSubpath)
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = strrep(SCRIPTPATH, fullfile('analysis_script', 'MATLAB'), '');
INPATH_RAW = fullfile(MAINPATH, 'data', 'processed_data', 'markers_included');
INPATH_CLEAN = fullfile(MAINPATH, 'data', 'processed_data', 'erp_preprocessed_clean');
INPATH_ICA = fullfile(MAINPATH, 'data', 'processed_data', 'ica_preprocessed');
INPATH_EXCLUDED = fullfile(MAINPATH, 'data', 'processed_data', 'exclude_trials');
OUTPATH = fullfile(MAINPATH, 'data', 'analysis_data', 'check_filters_erp');
FUNPATH = fullfile(MAINPATH, 'functions');

addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_RAW, INPATH_ICA, INPATH_EXCLUDED, OUTPATH, INPATH_CLEAN)
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
ERP_FROM = 90;
ERP_TILL = 140;

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
dircont_subj = dir(fullfile(INPATH_CLEAN, 'sub-*.set')); % Only runs through non-excluded subjects

% Sanity Check: Same number of files for raw data and ICA data
if length(dir(fullfile(INPATH_RAW, 'sub-*.set'))) == length(dir(fullfile(INPATH_ICA, 'sub-*.set'))) 
else
    error('Number of raw data files, exclusion files and ICA data files does not match')
end

% Initialize variables
all_erp_filtered = cell(0,3); 
all_con_erp_filtered = cell(0,3); 
all_cor_erp_filtered = cell(0,3); 
all_erp_unfiltered = cell(0,3); 
all_con_erp_unfiltered = cell(0,3); 
all_cor_erp_unfiltered = cell(0,3); 

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_check_filters_erp.m');

for subj_idx= 1:length(dircont_subj)

    %% Preprocessing
    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_check_filters_erp.m'])
    tic;

    % Store subject-specific ERP data temporarily for individual plots
    current_subj_data = struct(); % Holds the data for the current subject for all conditions, for both filtered and unfiltered states

    % Loop-Variable for filtered vs. non-filtered processing
    filter_settings = {'filtered', 'unfiltered'}; 

    for f_idx = 1:length(filter_settings)
        current_filter_setting = filter_settings{f_idx};
        use_filters = strcmp(current_filter_setting, 'filtered'); % true for 'filtered', false for 'unfiltered'

        % Get to-be excluded trials based on tid_psam_exclude_trials.m
        exclusion_filename = fullfile(INPATH_EXCLUDED,[subj '_excluded_trials.mat']);
        load(exclusion_filename) % Loads variables 'excluded_trials_erp_beh' (used here) and 'excluded_trials_svm' (not used here)

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
        % Apply Filter (Conditional based on use_filters)
        if use_filters
            % Highpass-Filter
            LCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth_TD(LCF)); % Get filter order (also see pop_firwsord.m, tid_psam_get_transition_bandwidth_TD.m)
            %%EEG = pop_firws(EEG, 'fcutoff', LCF, 'ftype', 'highpass', 'wtype', 'hamming', 'forder',LCF_ord, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
            % Lowpass-Filter
            HCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth_TD(HCF)); % Get filter order (also see pop_firwsord.m, tid_psam_get_transition_bandwidth_TD.m)
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

        % Sanity Check: Plot RMS in 10s bins for each electode (not used here)
        %%tid_psam_plot_rms_bins_TD(EEG, [subj ' RMS bins'], 'SavePath', fullfile(OUTPATH, [subj '_channel_rms.png']), 'PlotOn', false)

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

        % Store data 
        EEG.setname = [subj '_' current_filter_setting '_all_conds'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %% Analysis/Plots 
        % Initialite temporary variables
        temp_erp_store = cell(length(EVENTS), 3); % Holds erp, subj, event
        temp_con_erp_store = cell(0, 3); % Changed to dynamic growth
        temp_cor_erp_store = cell(length(EVENTS), 3); % Holds cor_erp, subj, event
        
        temp_con_counter_inner = 1; % For unique control ERPs within this filter setting

        % Get ERPs for each condition
        for cond = 1:length(EVENTS)
            EEG_cond = pop_selectevent( ALLEEG(5), 'latency','-2<=2','type',EVENTS(cond), ...
                'deleteevents','off','deleteepochs','on','invertepochs','off');

            % Get ERP
            erp = mean(EEG_cond.data, 3);
            % Get number of trials (uncorrected)
            n_trials_uncor = size(EEG_cond.data,3);

            % Correct ERP by subtracting the ERP from the equivalent no-probe condition (see Daliri & Max, 2016)
            % Get matching control condition
            con_cond_idx = find(strcmp(erase(EVENTS{cond}, {'_alt', '_unalt'}), erase(string(CON_EVENTS), 'con_')));
            % Get control ERP
            EEG_con = pop_selectevent( ALLEEG(end), 'latency','-2<=2','type',CON_EVENTS(con_cond_idx), ...
                'deleteevents','off','deleteepochs','on','invertepochs','off');
            con_erp = mean(EEG_con.data, 3);
            % Get number of trials (control)
            n_trials_con = size(EEG_con.data,3);

            % Get corrected ERP
            cor_erp = erp - con_erp;

            % Setup ERP analysis
            [~,win_start] = min(abs(EEG.times-ERP_FROM));
            [~,win_end] = min(abs(EEG.times-ERP_TILL));
            [~,t_zero] = min(abs(EEG.times)); 
            
            % ERP analysis (control)
            % Store ERP (control)
            temp_con_erp_store{temp_con_counter_inner,1} = con_erp;
            temp_con_erp_store{temp_con_counter_inner,2} = subj;
            temp_con_erp_store{temp_con_counter_inner,3} = CON_EVENTS{con_cond_idx};
            temp_con_counter_inner = temp_con_counter_inner+1;


            % ERP analysis (uncorrected)
            % Store ERP (uncorrected)
            temp_erp_store{cond,1} = erp;
            temp_erp_store{cond,2} = subj;
            temp_erp_store{cond,3} = EVENTS{cond};

            % ERP analysis (corrected)
            % Store ERP (corrected)
            temp_cor_erp_store{cond,1} = cor_erp;
            temp_cor_erp_store{cond,2} = subj;
            temp_cor_erp_store{cond,3} = EVENTS{cond};
            
        end 

        % Store the current subject's data for this filter setting
        if use_filters
            current_subj_data.filtered.erp = temp_erp_store;
            current_subj_data.filtered.con_erp = temp_con_erp_store(1:temp_con_counter_inner-1,:); 
            current_subj_data.filtered.cor_erp = temp_cor_erp_store;
            
            % Append to global filtered storage
            all_erp_filtered = [all_erp_filtered; temp_erp_store];
            all_con_erp_filtered = [all_con_erp_filtered; temp_con_erp_store(1:temp_con_counter_inner-1,:)];
            all_cor_erp_filtered = [all_cor_erp_filtered; temp_cor_erp_store];
            
        else
            current_subj_data.unfiltered.erp = temp_erp_store;
            current_subj_data.unfiltered.con_erp = temp_con_erp_store(1:temp_con_counter_inner-1,:);
            current_subj_data.unfiltered.cor_erp = temp_cor_erp_store;
            
            % Append to global unfiltered storage
            all_erp_unfiltered = [all_erp_unfiltered; temp_erp_store];
            all_con_erp_unfiltered = [all_con_erp_unfiltered; temp_con_erp_store(1:temp_con_counter_inner-1,:)];
            all_cor_erp_unfiltered = [all_cor_erp_unfiltered; temp_cor_erp_store];
        end

    end 

    % Individual plots 
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

    for cond_plot_idx = 1:length(EVENTS)

        % Get data for current subject and condition for both filtered and unfiltered
        % Filtered data
        erp_filtered = current_subj_data.filtered.erp{cond_plot_idx, 1};
        % Get the correct control ERP based on the condition index
        con_cond_idx_for_plot = find(strcmp(erase(EVENTS{cond_plot_idx}, {'_alt', '_unalt'}), erase(string(CON_EVENTS), 'con_')));
        con_erp_filtered_row_idx = cellfun(@(x) strcmp(x, CON_EVENTS{con_cond_idx_for_plot}), current_subj_data.filtered.con_erp(:,3));
        con_erp_filtered = current_subj_data.filtered.con_erp{con_erp_filtered_row_idx,1};
        cor_erp_filtered = current_subj_data.filtered.cor_erp{cond_plot_idx, 1};

        % Unfiltered data
        erp_unfiltered = current_subj_data.unfiltered.erp{cond_plot_idx, 1};
        con_erp_unfiltered_row_idx = cellfun(@(x) strcmp(x, CON_EVENTS{con_cond_idx_for_plot}), current_subj_data.unfiltered.con_erp(:,3));
        con_erp_unfiltered = current_subj_data.unfiltered.con_erp{con_erp_unfiltered_row_idx,1};
        cor_erp_unfiltered = current_subj_data.unfiltered.cor_erp{cond_plot_idx, 1};

        % Plot: Uncorrected ERP 
        subplot(6, 8, (r_erp_unc - 1) * 8 + cond_plot_idx)
        % Plot filtered data (dashed)
        plot(EEG.times, erp_filtered(CHANI,:), ':', 'LineWidth', 1.5, 'Color', colors{4}); % Filtered data
        hold on
        % Plot unfiltered data 
        plot(EEG.times, erp_unfiltered(CHANI,:), 'LineWidth', 1.5, 'Color', colors{4}); % Unfiltered data
        hold off
        title(EVENTS{cond_plot_idx}, 'Interpreter', 'none')
        ylim([ylim_lower ylim_upper]); xlim([-200 400])
        if cond_plot_idx == 1
            ylabel('Uncorrected');
        end

        subplot(6, 8, (r_topo_unc - 1) * 8 + cond_plot_idx)
        topoplot(mean(erp_filtered(:,win_start:win_end),2), EEG.chanlocs, ... 
            'emarker2', {CHANI,'o','r',2,2},'emarker', {'.','k',0.1,1});
        colormap("parula")
        colorbar;
        clim([cb_lim_lower cb_lim_upper])

        % Plot: Control ERP 
        subplot(6, 8, (r_erp_con - 1) * 8 + cond_plot_idx)
        % Plot filtered data (dashed)
        plot(EEG.times, con_erp_filtered(CHANI,:), ':', 'LineWidth', 1.5, 'Color', colors{4}); 
        hold on
        % Plot unfiltered data 
        plot(EEG.times, con_erp_unfiltered(CHANI,:), 'LineWidth', 1.5, 'Color', colors{4}); 
        hold off
        ylim([ylim_lower ylim_upper]); xlim([-200 400])
        if cond_plot_idx == 1
            ylabel('Control');
        end

        subplot(6, 8, (r_topo_con - 1) * 8 + cond_plot_idx)
        topoplot(mean(con_erp_filtered(:,win_start:win_end),2), EEG.chanlocs, ... 
            'emarker2', {CHANI,'o','r',2,2},'emarker', {'.','k',0.1,1});
        colormap("parula")
        colorbar;
        clim([cb_lim_lower cb_lim_upper])

        % Plot: Corrected ERP 
        subplot(6, 8, (r_erp_cor - 1) * 8 + cond_plot_idx)
        % Plot filtered data (dashed)
        plot(EEG.times, cor_erp_filtered(CHANI,:), ':', 'LineWidth', 1.5, 'Color', colors{4}); 
        hold on
        % Plot unfiltered data 
        plot(EEG.times, cor_erp_unfiltered(CHANI,:), 'LineWidth', 1.5, 'Color', colors{4}); 
        hold off
        ylim([ylim_lower ylim_upper]); xlim([-200 400])
        if cond_plot_idx == 1
            ylabel('Corrected');
        end

        subplot(6, 8, (r_topo_cor - 1) * 8 + cond_plot_idx)
        topoplot(mean(cor_erp_filtered(:,win_start:win_end),2), EEG.chanlocs, ... 
            'emarker2', {CHANI,'o','r',2,2},'emarker', {'.','k',0.1,1});
        colormap("parula")
        colorbar;
        clim([cb_lim_lower cb_lim_upper])

    end 

    sgtitle(['Filter Check ERPs and Topoplots for ' subj ' (Filtered = Dashed Line, Unfiltered = Solid Line) (Topoplots filtered)'], 'Interpreter', 'none');
    saveas(individual_plot, fullfile(OUTPATH, [subj '_filter_check_erp_topo.png']));
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

end 

% Get grandaverage ERPs (filtered)
erp_data_col = 1;
cond_col = 3;
% Get condition labels
conditions_to_check_filtered = all_erp_filtered(:, cond_col);
%%is_char_cell_filtered = cellfun(@ischar, conditions_to_check_filtered);
%%clean_conditions_filtered = conditions_to_check_filtered(is_char_cell_filtered);
[unique_conditions_filtered, ~, condition_idx_filtered] = unique(conditions_to_check_filtered);

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
%%is_char_cell_con_filtered = cellfun(@ischar, conditions_to_check_con_filtered);
%%clean_conditions_con_filtered = conditions_to_check_con_filtered(is_char_cell_con_filtered);
[unique_conditions_con_filtered, ~, condition_idx_con_filtered] = unique(conditions_to_check_con_filtered);

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
%%is_char_cell_cor_filtered = cellfun(@ischar, conditions_to_check_cor_filtered);
%%clean_conditions_cor_filtered = conditions_to_check_cor_filtered(is_char_cell_cor_filtered);
[unique_conditions_cor_filtered, ~, condition_idx_cor_filtered] = unique(conditions_to_check_cor_filtered);

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
%%is_char_cell_unfiltered = cellfun(@ischar, conditions_to_check_unfiltered);
%%clean_conditions_unfiltered = conditions_to_check_unfiltered(is_char_cell_unfiltered);
[unique_conditions_unfiltered, ~, condition_idx_unfiltered] = unique(conditions_to_check_unfiltered);

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
%%is_char_cell_con_unfiltered = cellfun(@ischar, conditions_to_check_con_unfiltered);
%%clean_conditions_con_unfiltered = conditions_to_check_con_unfiltered(is_char_cell_con_unfiltered);
[unique_conditions_con_unfiltered, ~, condition_idx_con_unfiltered] = unique(conditions_to_check_con_unfiltered);

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
%%is_char_cell_cor_unfiltered = cellfun(@ischar, conditions_to_check_cor_unfiltered);
%%clean_conditions_cor_unfiltered = conditions_to_check_cor_unfiltered(is_char_cell_cor_unfiltered);
[unique_conditions_cor_unfiltered, ~, condition_idx_cor_unfiltered] = unique(conditions_to_check_cor_unfiltered);

n_conds_cor_unfiltered = numel(unique_conditions_cor_unfiltered);
grandaverage_erp_cor_unfiltered = cell(n_conds_cor_unfiltered, 2);
for i = 1:n_conds_cor_unfiltered
    erp_group = all_cor_erp_unfiltered(condition_idx_cor_unfiltered == i, erp_data_col);
    erp_stack = cat(3, erp_group{:});
    grandaverage_erp_cor_unfiltered{i, 1} = mean(erp_stack, 3);
    grandaverage_erp_cor_unfiltered{i, 2} = unique_conditions_cor_unfiltered{i};
end

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

for i = 1:size(grandaverage_erp_filtered,1) % for uncorrected ERP (filtered)
    old_label = grandaverage_erp_filtered{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp_filtered{i,2} = new_label{1};
end
for i = 1:size(grandaverage_erp_cor_filtered,1) % for corrected ERP (filtered)
    old_label = grandaverage_erp_cor_filtered{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp_cor_filtered{i,2} = new_label{1};
end

for i = 1:size(grandaverage_erp_unfiltered,1) % for uncorrected ERP (unfiltered)
    old_label = grandaverage_erp_unfiltered{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp_unfiltered{i,2} = new_label{1};
end
for i = 1:size(grandaverage_erp_cor_unfiltered,1) % for corrected ERP (unfiltered)
    old_label = grandaverage_erp_cor_unfiltered{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp_cor_unfiltered{i,2} = new_label{1};
end


% Split data into early and late conditions for filtered and unfiltered data
grandaverage_erp_cor_early_filtered = grandaverage_erp_cor_filtered(contains(grandaverage_erp_cor_filtered(:,2), 'Early'), :);
grandaverage_erp_cor_late_filtered = grandaverage_erp_cor_filtered(contains(grandaverage_erp_cor_filtered(:,2), 'Late'), :);

grandaverage_erp_cor_early_unfiltered = grandaverage_erp_cor_unfiltered(contains(grandaverage_erp_cor_unfiltered(:,2), 'Early'), :);
grandaverage_erp_cor_late_unfiltered = grandaverage_erp_cor_unfiltered(contains(grandaverage_erp_cor_unfiltered(:,2), 'Late'), :);

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'tid_psam_check_filters_erp_protocol.xlsx'))

close all; delete(wb);

%% Grandaverage plots
close all

% Plot: ERPs for all conditions (filtered and unfiltered) 
figure('Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
sgtitle('Corrected ERPs for all Conditions (Filtered = Dashed Line, Unfiltered = Solid Line)')

% Plot early probes 
subplot(2,1,1)
hold on; 
plot_handles_early = gobjects(1, length(grandaverage_erp_cor_early_filtered) * 2); % Preallocate for plot handles
idx_handle = 1;
for cond_early_num = 1:length(grandaverage_erp_cor_early_filtered)
    erp_2_plot_filtered = grandaverage_erp_cor_early_filtered{cond_early_num,1};
    erp_2_plot_unfiltered = grandaverage_erp_cor_early_unfiltered{cond_early_num,1};
    % Plot filtered data (dashed)
    p1 = plot(EEG.times, erp_2_plot_filtered(chani,:), ':', 'LineWidth', 1.5, ...
        'Color', colors{cond_early_num}, 'DisplayName', [grandaverage_erp_cor_early_filtered{cond_early_num,2}], 'HandleVisibility', 'off');
        % Plot unfiltered data 
    p2 = plot(EEG.times, erp_2_plot_unfiltered(chani,:), 'LineWidth', 1.5, ...
        'Color', colors{cond_early_num}, 'DisplayName', [grandaverage_erp_cor_early_unfiltered{cond_early_num,2}]);
    
    plot_handles_early(idx_handle) = p1;
    plot_handles_early(idx_handle+1) = p2;
    idx_handle = idx_handle + 2;
end
xlim([-200 400])
ylim([y_lim_lower y_lim_upper])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
title('ERPs for early probes')
legend('Location','northwest')
fill([70 130 130 70], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], light_blue, 'FaceAlpha',0.1, 'EdgeColor','none', 'HandleVisibility', 'off');
hold off

% Plot late probes 
subplot(2,1,2)
hold on; 
plot_handles_late = gobjects(1, length(grandaverage_erp_cor_late_filtered) * 2); % Preallocate for plot handles
idx_handle = 1;
for cond_late_num = 1:length(grandaverage_erp_cor_late_filtered)
    erp_cor_2_plot_filtered = grandaverage_erp_cor_late_filtered{cond_late_num,1};
    erp_cor_2_plot_unfiltered = grandaverage_erp_cor_late_unfiltered{cond_late_num,1};

    p1 = plot(EEG.times, erp_cor_2_plot_filtered(chani,:), ':', 'LineWidth', 1.5, ...
        'Color', colors{cond_late_num}, 'DisplayName', [grandaverage_erp_cor_late_filtered{cond_late_num,2}], 'HandleVisibility', 'off');
    p2 = plot(EEG.times, erp_cor_2_plot_unfiltered(chani,:), 'LineWidth', 1.5, ...
        'Color', colors{cond_late_num}, 'DisplayName', [grandaverage_erp_cor_late_unfiltered{cond_late_num,2}]);
    
    plot_handles_late(idx_handle) = p1;
    plot_handles_late(idx_handle+1) = p2;
    idx_handle = idx_handle + 2;
end
xlim([-200 400])
ylim([y_lim_lower y_lim_upper])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
title('ERPs for late probes')
legend('Location','northwest')
fill([70 130 130 70], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], light_blue, 'FaceAlpha',0.1, 'EdgeColor','none', 'HandleVisibility', 'off');
hold off

% Save plot
print(gcf, fullfile(OUTPATH, 'grandaverage_filter_check_erp.png'), '-dpng', '-r900')


% Plot: ERP, control ERP and correct ERP for each condition (filtered and unfiltered)
% Concatinate colors for the three ERP types (Uncorrected, Corrected, Control)
colors_erp_types = {
    main_blue; % Uncorrected (filtered)
    main_red;  % Corrected (filtered)
    main_green; % Control (filtered)
    };

figure('Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
sgtitle('Uncorrected, Corrected and Control ERPs for all Conditions (Filtered = Dashed Line, Unfiltered = Solid Line)')

con_cond_num_filtered = 1;
con_cond_num_unfiltered = 1;

for cond_num = 1:length(grandaverage_erp_cor_filtered) 
    subplot(2,4,cond_num)
    hold on; 
    % Filtered ERPs
    erp_2_plot_filtered = grandaverage_erp_filtered{cond_num,1};
    erp_cor_2_plot_filtered = grandaverage_erp_cor_filtered{cond_num,1};
    erp_con_2_plot_filtered = grandaverage_erp_con_filtered{con_cond_num_filtered,1}; 

    % Unfiltered ERPs
    erp_2_plot_unfiltered = grandaverage_erp_unfiltered{cond_num,1};
    erp_cor_2_plot_unfiltered = grandaverage_erp_cor_unfiltered{cond_num,1};
    erp_con_2_plot_unfiltered = grandaverage_erp_con_unfiltered{con_cond_num_unfiltered,1}; 

    % Plot filtered data (dashed)
    p1 = plot(EEG.times, erp_2_plot_filtered(chani,:), ':', 'LineWidth', 1.5, 'Color', colors_erp_types{1}, 'HandleVisibility', 'off');
    p2 = plot(EEG.times, erp_cor_2_plot_filtered(chani,:), ':', 'LineWidth', 1.5, 'Color', colors_erp_types{2}, 'HandleVisibility', 'off');
    p3 = plot(EEG.times, erp_con_2_plot_filtered(chani,:), ':', 'LineWidth', 1.5, 'Color', colors_erp_types{3}, 'HandleVisibility', 'off');

    % Plot unfiltered data 
    p4 = plot(EEG.times, erp_2_plot_unfiltered(chani,:), 'LineWidth', 1.5, 'Color', colors_erp_types{1}, 'DisplayName', 'Uncorrected');
    p5 = plot(EEG.times, erp_cor_2_plot_unfiltered(chani,:), 'LineWidth', 1.5, 'Color', colors_erp_types{2}, 'DisplayName', 'Corrected');
    p6 = plot(EEG.times, erp_con_2_plot_unfiltered(chani,:),  'LineWidth', 1.5, 'Color', colors_erp_types{3}, 'DisplayName', 'Control');
    
    xlim([-200 400])
    ylim([y_lim_lower-1 y_lim_upper+1])
    xlabel('Time [ms]')
    ylabel('Amplitude [µV]')
    title(grandaverage_erp_cor_filtered{cond_num,2})
    legend('Location','northwest')
    hold off

    % Increment control condition counters 
    if mod(cond_num,2) == 0
        con_cond_num_filtered = con_cond_num_filtered +1;
        con_cond_num_unfiltered = con_cond_num_unfiltered +1;
    end
end

% Save plot
print(gcf, fullfile(OUTPATH, 'grandaverage_filter_check_erp_con_cor.png'), '-dpng', '-r900')


close all
