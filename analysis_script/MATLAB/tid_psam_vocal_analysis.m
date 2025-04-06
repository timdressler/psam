% tid_psam_vocal_analysis.m
%
% Performs analysis of vocal data.
% Has to be executed AFTER tid_psam_vocal_analysis.praat.
%
% Tim Dressler, 04.04.2025

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
INPATH = fullfile(MAINPATH, 'data\BIDS\');
INPATH_PRAAT = fullfile(MAINPATH, 'data\analysis_data\vocal_analysis_preparation\');
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\vocal_analysis');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
n_trials = 960;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*'));

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
%%wb = waitbar(0,'starting tid_psam_vocal_analysis.m');

clear subj_idx n_exluded_trials
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    %%waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_vocal_analysis.m'])

    tic;

    % Get log file & Sanity Check: One log file per subject
    subj_log_filename = dir(fullfile(INPATH, [subj '\beh\*.csv']));
    if numel(subj_log_filename) == 1
        subj_log = readtable(fullfile(subj_log_filename.folder, subj_log_filename.name),'VariableNamingRule', 'preserve');
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'number_of_log_files';
    end

    % Remove all rows from log file which don not resembel experimental trials (e.g. Instruction trials)
    subj_log_cleaned = subj_log(~isnan(subj_log.("mic.started")), :);
    % Only included needed columns
    subj_log_cleaned = subj_log_cleaned(:, {'ISI.started', 'ISI.stopped', 'aa_green_onset', 'conditions_file', ...
        'cue_stim.started', 'cue_stim.stopped', 'date', 'expName', 'go_port.started', 'mic.clip', ...
        'mic.started', 'mic.stopped', 'probe', 'probe_duration', 'probe_intensity', 'probe_marker', ...
        'probe_marker_port.started', 'probe_onset', 'probe_onset_cat', 'probe_stim.started', 'probe_type', ...
        'psychopyVersion', 'rec_duration', 'stim_file', 'subj', 'task', 'task_marker', ...
        'target_stim.started', 'trial.started'});

    % Sanity Check: Equal db for unaltered and altered probes
    if height(subj_log_cleaned) == n_trials
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'log_dimensions';
    end
    % Extract filename from full path for each vocal recording in order to match it with the table produced in Praat
    subj_log_cleaned.filename_tab = cellfun(@(x) strrep(regexp(x, 'recording_mic_.*(?=\.wav)', 'match', 'once'), '.', '_'), subj_log_cleaned.("mic.clip"), 'UniformOutput', false);

    % Get probe properties file & Sanity Check: Correct dimensions
    subj_probe_properties = readtable(fullfile(INPATH, ['stimuli\' subj '\' subj '_probe_properties.xlsx']),'VariableNamingRule', 'preserve');
    if height(subj_probe_properties) == 1 && width(subj_probe_properties) == 7
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'probe_properties_dimensions';
    end
    % Sanity Check: Equal db for unaltered and altered probes
    if round(subj_probe_properties.db_tab_normal) == round(subj_probe_properties.db_tab_pitched)
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'probe_properties_db';
    end

    % Get F0 & RT Table & Sanity Check: Correct dimensions
    subj_f0_rt = readtable(fullfile(INPATH_PRAAT, [subj '_f0_rt_table.csv']),'VariableNamingRule', 'preserve', Delimiter=',');
    subj_f0_rt = standardizeMissing(subj_f0_rt,9999);

    if height(subj_f0_rt) == n_trials && width(subj_f0_rt) == 14
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'f0_rt_dimensions';
    end

    % Concatinate log file, probe properties file and F0 & RT Table
    subj_full = join(subj_f0_rt, subj_log_cleaned, 'Keys', 'filename_tab');
    % Exclude 'filename_tab' column and repeat the remaining columns for each row of subj_full
    subj_full = [subj_full, repmat(subj_probe_properties(:, setdiff(1:width(subj_probe_properties), find(strcmp(subj_probe_properties.Properties.VariableNames, 'filename_tab')))), height(subj_full), 1)];

    % Exlcude trials with outlier F0
    % Get threshold
    mean_f0 = mean(subj_full.f0_tab, 'omitnan');
    std_f0 = std(subj_full.f0_tab, 'omitnan');
    threshold_upper = mean_f0 + 3 * std_f0;
    threshold_lower = mean_f0 - 3 * std_f0;
    % Update excluded trials based on F0 outliers
    excluded_idx = (subj_full.f0_tab > threshold_upper) | (subj_full.f0_tab < threshold_lower);
    % Exclude trials based on F0 outliers
    subj_full = subj_full(~excluded_idx, :);

    % Store number of exluded trials
    n_exluded_trials = sum(excluded_idx);

    % Only included needed columns
    subj_full_cleaned = subj_full(:, {'mic.clip', ...
        'mic.started', 'mic.stopped', 'probe', ...
        'probe_onset', 'probe_onset_cat', 'probe_type', ...
        'subj', 'task', 'task_marker', 'f0_tab', 'rt_tab', 'condition_tab', ...
        'min_intensity_tab', 'max_intensity_tab', 'vocal_response_tab', ...
        'duration_vocal_tab','onset_vocal_tab', 'filename_tab', 'f0_tab_normal', ...
        'f0_tab_pitched', 'loudness', 'change_attempts'});

    % Compare F0 distribution with F0 of unaltered and altered probes
    [h_normal, p_normal, ci_normal, stats_normal] = ttest(subj_full_cleaned.f0_tab, subj_full_cleaned.f0_tab_normal(1));
    cohens_d_normal = stats_normal.tstat / sqrt(stats_normal.df + 1);
    [h_pitched, p_pitched, ci_pitched, stats_pitched] = ttest(subj_full_cleaned.f0_tab, subj_full_cleaned.f0_tab_pitched(1));
    cohens_d_pitched = stats_pitched.tstat / sqrt(stats_pitched.df + 1);

    % Plot: F0 Distribution including probe F0's
    figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
    [f, xi] = ksdensity(subj_full_cleaned.f0_tab); % Density of F0 vocal responses
    plot(xi, f, 'LineWidth', 2, 'Color','k');
    hold on;
    xline(subj_full_cleaned.f0_tab_pitched(1), 'r--', 'LineWidth', 2); % F0 altered probe
    xline(subj_full_cleaned.f0_tab_normal(1), 'b--', 'LineWidth', 2); % F0 unaltered probe

    ylims = ylim; % Get current y-axis limits
    text(subj_full_cleaned.f0_tab_pitched(1)+0.25, ylims(2)*0.9, ... % Test-Statistic for altered probe
        sprintf('t=%.2f, p=%.3f, d=%.2f', stats_pitched.tstat, p_pitched, cohens_d_pitched), ...
        'Color', 'red', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    text(subj_full_cleaned.f0_tab_normal(1)+0.25, ylims(2)*0.8, ... % Test-Statistic for unaltered probe
        sprintf('t=%.2f, p=%.3f, d=%.2f', stats_normal.tstat, p_normal, cohens_d_normal), ...
        'Color', 'blue', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    xlabel('F0 [Hz]');
    ylabel('Density');
    title(['F0 Density of Vocal Responses and F0 Values for Probes for ' subj]);

    legend('F0 Distribution Vocal Responses', 'F0 Altered Probe', 'F0 Unaltered Probe', 'Location', 'southwest', 'Interpreter', 'none');
    hold off;

    exportgraphics(gcf,fullfile(OUTPATH, ['tid_psam_f0_distribution_' subj '.png']),'Resolution',1000)

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
writetable(protocol,fullfile(OUTPATH, 'tid_psam_vocal_analysis_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_vocal_analysis_marked_subj.xlsx'))
end

%%close(wb)

check_done = 'tid_psam_vocal_analysis_DONE'
