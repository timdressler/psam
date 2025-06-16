% tid_psam_erp_analysis.m
%
% Performs ERP analysis, creates plots and exports data for further analysis.
%
% Processing steps include
%
% Extract ERP and control ERP for each condition
%   Control ERPs are extracted from no-probe trials and are thought to
%   only reflect non-auditory processes (e.g. cognitive ones). The
%   uncorrected ERPs are thought to reflect both auditory and
%   non-auditory processes. The corrected ERPs are thought to relfect
%   mainly auditory processes (see Daliri & Max, 2016)
% Get corrected ERP
% Extract ERP component infortmation
% Plot uncorrected, corrected and control ERPs for each condition
% Plot corrected ERPs and Topoplots for each condition
%
% Literature
% Daliri, A., & Max, L. (2016). Modulation of auditory responses to speech vs.
%   nonspeech stimuli during speech movement planning.
%   Frontiers in human neuroscience, 10. https://doi.org/10.3389/fnhum.2016.00234
%
% Tim Dressler, 17.04.2025

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
INPATH = fullfile(MAINPATH, 'data\processed_data\erp_preprocessed_clean\');
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\erp_analysis');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ... % Real events
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt'};
CON_EVENTS = {'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};
CHANI = 1; % Channel to plot ERP from
ERP_FROM = 70;
ERP_TILL = 130;
INDIVIDUAL_PLOTS = true; % Whether or not to create individual ERP plots and Topoplots

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
dircont_subj = dir(fullfile(INPATH, 'sub-*_clean.set'));

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_erp_analysis.m');

clear subj_idx
counter = 1;
cor_counter = 1;
con_counter = 1;
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_erp_analysis.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load data
    EEG = pop_loadset('filename',[subj '_erp_preprocessed_clean.set'],'filepath',INPATH);

    % Remove EOG channels
    EEG = pop_select( EEG, 'rmchannel',{'E29','E30'});

    % Store data
    EEG.setname = [subj '_all_conds'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Get ERPs for each condition
    for cond = 1:length(EVENTS)
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off');

        % Get ERP
        erp = mean(EEG.data, 3);

        % Get number of trials (uncorrected)
        n_trials_uncor = size(EEG.data,3);

        % Correct ERP by subtracting the ERP from the equivalent no-probe condition (see Daliri & Max, 2016)
        % Get matching control condition
        con_cond = find(strcmp(erase(EVENTS{cond}, {'_alt', '_unalt'}), erase(string(CON_EVENTS), 'con_')));
        % Get control ERP
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',CON_EVENTS(con_cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off');
        con_erp = mean(EEG.data, 3);
        % Get number of trials (control)
        n_trials_con = size(EEG.data,3);
        % Get corrected ERP
        cor_erp = erp - con_erp;

        % Setup ERP analysis
        [~,win_start] = min(abs(EEG.times-ERP_FROM));
        [~,win_end] = min(abs(EEG.times-ERP_TILL));
        [~,t_zero] = min(abs(EEG.times));

        % ERP analysis (control)
        % Store ERP (control)
        all_con_erp{con_counter,1} = con_erp;
        all_con_erp{con_counter,2} = subj;
        all_con_erp{con_counter,3} = CON_EVENTS{con_cond};
        con_counter = con_counter+1;

        % ERP analysis (uncorrected)
        % Get N100 amplitude (uncorrected)
        erp_amp = min(erp(CHANI,win_start:win_end));
        % Get N100 latency (uncorrected)
        erp_sam = find(erp(CHANI,:) == erp_amp);
        erp_lat = EEG.times(erp_sam);
        % Store ERP (uncorrected)
        all_erp{counter,1} = erp;
        all_erp{counter,2} = subj;
        all_erp{counter,3} = EVENTS{cond};
        counter = counter+1;

        % ERP analysis (corrected)
        % Get N100 amplitude (corrected)
        cor_erp_amp = min(cor_erp(CHANI,win_start:win_end));
        % Get N100 latency (corrected)
        cor_erp_sam = find(cor_erp(CHANI,:) == cor_erp_amp);
        cor_erp_lat = EEG.times(cor_erp_sam);
        % Store ERP (corrected)
        all_cor_erp{cor_counter,1} = cor_erp;
        all_cor_erp{cor_counter,2} = subj;
        all_cor_erp{cor_counter,3} = EVENTS{cond};
        % Store ERP data (corrected)
        all_cor_erp_data{cor_counter,1} = subj;
        all_cor_erp_data{cor_counter,2} = ERP_FROM;
        all_cor_erp_data{cor_counter,3} = ERP_TILL;
        all_cor_erp_data{cor_counter,4} = cor_erp_amp;
        all_cor_erp_data{cor_counter,5} = cor_erp_lat;
        all_cor_erp_data{cor_counter,6} = EVENTS{cond};
        all_cor_erp_data{cor_counter,10} = n_trials_uncor;
        all_cor_erp_data{cor_counter,11} = n_trials_con;

        % Split condition labels into multiple columns (corrected)
        cond_parts = strsplit(EVENTS{cond}, '_');
        % Map and rename task condition (corrected)
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
        all_cor_erp_data{cor_counter, 7} = task_label;
        all_cor_erp_data{cor_counter, 8} = probe_onset_label;
        all_cor_erp_data{cor_counter, 9} = probe_type_label;
        cor_counter = cor_counter+1;

        % Individual plots
        if INDIVIDUAL_PLOTS
            % Define limits
            ylim_upper = 25;
            ylim_lower = -15;
            cb_lim_upper = 5;
            cb_lim_lower = -5;

            % Create figure once
            if cond == 1
                individual_plot = figure('Units', 'normalized', 'Position', [0 0 1 1]);
            end

            % Plot rows
            r_erp_unc = 1; r_topo_unc = 2;
            r_erp_con = 3; r_topo_con = 4;
            r_erp_cor = 5; r_topo_cor = 6;

            % Plot: Uncorrected ERP and Topoplot
            subplot(6, 8, (r_erp_unc - 1) * 8 + cond)
            plot(EEG.times, erp(CHANI,:), 'LineWidth', 1.5, 'Color', colors{4});
            hold on
            scatter(erp_lat, erp_amp, 'filled')
            hold off
            title(EVENTS{cond}, 'Interpreter', 'none')
            ylim([ylim_lower ylim_upper]); xlim([-200 400])
            if cond == 1
                ylabel('Uncorrected'); 
            end

            subplot(6, 8, (r_topo_unc - 1) * 8 + cond)
            topoplot(mean(erp(:,win_start:win_end),2), EEG.chanlocs, ...
                'emarker2', {CHANI,'o','r',2,2},'emarker', {'.','k',0.1,1});            
            colormap("parula")
            colorbar;
            clim([cb_lim_lower cb_lim_upper])

            % Plot: Control ERP and Topoplot
            subplot(6, 8, (r_erp_con - 1) * 8 + cond)
            plot(EEG.times, con_erp(CHANI,:), 'LineWidth', 1.5, 'Color', colors{4});
            ylim([ylim_lower ylim_upper]); xlim([-200 400])
            if cond == 1
                ylabel('Control'); 
            end

            subplot(6, 8, (r_topo_con - 1) * 8 + cond)
            topoplot(mean(con_erp(:,win_start:win_end),2), EEG.chanlocs, ...
                'emarker2', {CHANI,'o','r',2,2},'emarker', {'.','k',0.1,1});
            colormap("parula")
            colorbar;
            clim([cb_lim_lower cb_lim_upper])

            % Plot: Corrected ERP and Topoplot
            subplot(6, 8, (r_erp_cor - 1) * 8 + cond)
            plot(EEG.times, cor_erp(CHANI,:), 'LineWidth', 1.5, 'Color', colors{4});
            hold on
            scatter(cor_erp_lat, cor_erp_amp, 'filled')
            hold off
            ylim([ylim_lower ylim_upper]); xlim([-200 400])
            if cond == 1
                ylabel('Corrected'); 
            end

            subplot(6, 8, (r_topo_cor - 1) * 8 + cond)
            topoplot(mean(cor_erp(:,win_start:win_end),2), EEG.chanlocs, ...
                'emarker2', {CHANI,'o','r',2,2},'emarker', {'.','k',0.1,1});            
            colormap("parula")
            colorbar;
            clim([cb_lim_lower cb_lim_upper])


            % Save after last condition
            if cond == length(EVENTS)
                sgtitle(['Sanity Check ERPs and Topoplots for ' subj], 'Interpreter', 'none');
                saveas(individual_plot, fullfile(OUTPATH, [subj '_sanity_erp_topo.png']));
                close(individual_plot)
            end
        end
    end

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
erp_data_col = 1;
cond_col = 3;
% Get grandaverage uncorrected ERP for each condition
% Extract unique conditions
[unique_conditions, ~, condition_idx] = unique(all_erp(:, cond_col));
n_conds = numel(unique_conditions);
grandaverage_erp = cell(n_conds, 2); % Preallocate
% Group and compute grand average
for i = 1:n_conds
    % Get all ERP matrices for this condition
    erp_group = all_erp(condition_idx == i, erp_data_col);

    % Convert to 3D array (chan x time x subjects)
    erp_stack = cat(3, erp_group{:});

    % Get grand average
    grandaverage_erp{i, 1} = mean(erp_stack, 3);
    grandaverage_erp{i, 2} = unique_conditions{i};
end

% Get grandaverage control ERP for each condition
% Extract unique conditions
[unique_conditions_con, ~, condition_idx] = unique(all_con_erp(:, cond_col));
n_conds = numel(unique_conditions_con);
grandaverage_erp_con = cell(n_conds, 2); % Preallocate
% Group and compute grand average
for i = 1:n_conds
    % Get all ERP matrices for this condition
    erp_group = all_con_erp(condition_idx == i, erp_data_col);

    % Convert to 3D array (chan x time x subjects)
    erp_stack = cat(3, erp_group{:});

    % Get grand average
    grandaverage_erp_con{i, 1} = mean(erp_stack, 3);
    grandaverage_erp_con{i, 2} = unique_conditions_con{i};
end

% Get grandaverage corrected ERP for each condition
% Extract unique conditions
[unique_conditions_cor, ~, condition_idx] = unique(all_cor_erp(:, cond_col));
n_conds = numel(unique_conditions_cor);
grandaverage_erp_cor = cell(n_conds, 2); % Preallocate
% Group and compute grand average
for i = 1:n_conds
    % Get all ERP matrices for this condition
    erp_group = all_cor_erp(condition_idx == i, erp_data_col);

    % Convert to 3D array (chan x time x subjects)
    erp_stack = cat(3, erp_group{:});

    % Get grand average
    grandaverage_erp_cor{i, 1} = mean(erp_stack, 3);
    grandaverage_erp_cor{i, 2} = unique_conditions_cor{i};
end

% Save ERPs (control, uncorrected, corrected)
save(fullfile(OUTPATH, 'all_subj_erp.mat'),'all_cor_erp')

% Save ERP data
all_cor_erp_data = cell2table(all_cor_erp_data, 'VariableNames',{'subj','erp_from', 'erp_till', 'erp_amp', 'erp_lat', 'condition_full', 'task_instruction', 'probe_onset_cat', 'probe_type', 'n_trials_uncorrected', 'n_trials_control'});
writetable(all_cor_erp_data, fullfile(OUTPATH, 'all_subj_cor_erp_data.xlsx'))

% Save grandaverage ERP
save(fullfile(OUTPATH, 'grandaverage_cor_erp.mat'),'grandaverage_erp_cor')

% Preparation for plots
% Get channel ID
chani = CHANI;
% Get values for dynamic plot limits
y_lim_lower = min(cellfun(@(x) min(x(:)), grandaverage_erp_cor(:,1)))-1;
y_lim_upper = max(cellfun(@(x) max(x(:)), grandaverage_erp_cor(:,1)))+1;
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
for i = 1:size(grandaverage_erp,1) % for uncorrected ERP
    old_label = grandaverage_erp{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp{i,2} = new_label{1};
end
for i = 1:size(grandaverage_erp_cor,1) % for corrected ERP
    old_label = grandaverage_erp_cor{i,2};
    new_label = rename_conditions_map(strcmp(rename_conditions_map(:,1), old_label), 2);
    grandaverage_erp_cor{i,2} = new_label{1};
end
% Split data into early and late conditions
grandaverage_erp_cor_early = grandaverage_erp_cor(contains(grandaverage_erp_cor(:,2), 'Early'), :);
grandaverage_erp_cor_late  = grandaverage_erp_cor(contains(grandaverage_erp_cor(:,2), 'Late'), :);

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'erp_analysis_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_erp_analysis_marked_subj.xlsx'))
end

check_done = 'tid_psam_erp_analysis_DONE'

delete(wb); close all;

%% Plots
close all

% Plot: ERPs for all conditions with Topoplots
% Create ERP plot
figure('Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
sgtitle('Corrected ERPs and Topoplots for all Conditions')
subplot(3,4,5:6)
for cond_early_num = 1:length(grandaverage_erp_cor_early)
    erp_2_plot = grandaverage_erp_cor_early{cond_early_num,1};

    plot(EEG.times, erp_2_plot(chani,:), 'LineWidth', 1.5, ...
        'Color', colors{cond_early_num});
    hold on
end
xlim([-200 400])
ylim([y_lim_lower y_lim_upper])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
title('ERPs for early probes')
fill([70 130 130 70], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], light_blue, 'FaceAlpha',0.1, 'EdgeColor','none');
l = legend([grandaverage_erp_cor_early(:,2)]', 'Interpreter','none');
l.Location = 'northwest';
fontsize(l,10,'points')
hold off

subplot(3,4,7:8)
for cond_late_num = 1:length(grandaverage_erp_cor_late)
    erp_cor_2_plot = grandaverage_erp_cor_late{cond_late_num,1};

    plot(EEG.times, erp_cor_2_plot(chani,:), 'LineWidth', 1.5, ...
        'Color', colors{cond_late_num});
    hold on
end
xlim([-200 400])
ylim([y_lim_lower y_lim_upper])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
title('ERPs for late probes')
fill([70 130 130 70], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], light_blue, 'FaceAlpha',0.1, 'EdgeColor','none');
l = legend([grandaverage_erp_cor_late(:,2)]', 'Interpreter','none');
l.Location = 'northwest';
fontsize(l,10,'points')
hold off

% Create Topoplots
% Get limit values
for cond_num = 1:length(grandaverage_erp_cor)
    topo_2_plot = grandaverage_erp_cor{cond_num,1};
    all_topo_2_plot(:,cond_num) = mean(topo_2_plot(:,win_start:win_end),2);
    cb_lim_lower = min(all_topo_2_plot, [], 'all');
    cb_lim_upper = max(all_topo_2_plot, [], 'all');
end

% Set up indices for Topoplots in Subplot
topo_idx = [1 2 3 4 9 10 11 12];

% Create Topoplot
for cond_num = 1:length(grandaverage_erp_cor)
    subplot(3,4,topo_idx(cond_num))
    topo_2_plot = grandaverage_erp_cor{cond_num,1};
    topoplot(mean(topo_2_plot(:,win_start:win_end),2), EEG.chanlocs, 'emarker2', {chani,'o','r',5,1});
    title(grandaverage_erp_cor{cond_num,2})
    colormap("parula")
    cb = colorbar;
    clim([cb_lim_lower cb_lim_upper]);
    cb.Label.String = 'Amplitude [µV]';
    cb.Label.Position = [2.2, -1, 0];
    cb.Label.HorizontalAlignment = 'center';
end

% Save plot
saveas(gcf,fullfile(OUTPATH, 'tid_psam_plot_grandaverage_erp_topo.png'))

% Plot: ERP, control ERP and correct ERP for each condition
% Concatinate colors
colors = {
    main_blue;
    main_red;
    main_green;
    };
% Create plot
figure('Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
sgtitle('Uncorrected, Corrected and Control ERPs for all Conditions')

con_cond_num = 1;
for cond_num = 1:length(grandaverage_erp_cor)
    erp_2_plot = grandaverage_erp{cond_num,1};
    erp_cor_2_plot = grandaverage_erp_cor{cond_num,1};
    erp_con_2_plot = grandaverage_erp_con{con_cond_num,1};
    subplot(2,4,cond_num)
    plot(EEG.times, erp_2_plot(chani,:), 'LineWidth', 1.5, 'Color', colors{1})
    hold on
    plot(EEG.times, erp_cor_2_plot(chani,:), 'LineWidth', 1.5, 'Color', colors{2})
    plot(EEG.times, erp_con_2_plot(chani,:), 'LineWidth', 1.5, 'Color', colors{3})
    xlim([-200 400])
    ylim([y_lim_lower-1 y_lim_upper+1])
    xlabel('Time [ms]')
    ylabel('Amplitude [µV]')
    title(grandaverage_erp_cor{cond_num,2})
    l = legend({'Uncorrected', 'Corrected', 'Control'}, 'Interpreter','none');
    l.Location = 'northwest';
    fontsize(l,10,'points')
    hold off
    if mod(cond_num,2) == 0 % There is no difference being made between unaltered control trials and altered control trials,
        %   thus every control condition is used twice (e.g. both Active - Late - Unaltered and Active - Late - Altered are corrected
        %   using the same control ERP (Control - Active - Late)
        con_cond_num = con_cond_num +1;
    end
end

% Save plot
saveas(gcf,fullfile(OUTPATH, 'tid_psam_plot_grandaverage_erp_con_cor.png'))




















