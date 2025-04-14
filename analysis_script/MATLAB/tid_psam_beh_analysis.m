% tid_psam_beh_analysis.m
%
% Performs analysis of vocal data.
%
% Tim Dressler, 07.04.2025

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
INPATH = fullfile(MAINPATH, 'data\processed_data\beh_preprocessed_clean\');
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\beh_analysis\');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*clean.xlsx'));

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_vocal_analysis.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_vocal_analysis.m'])

    tic;
    % Load data
    beh_clean = readtable(fullfile(INPATH,[subj '_beh_preprocessed_clean.xlsx']));

    %% F0 Analysis

    % Z-Transformation F0 data
    % Z-transform F0's of responses and probes based on the mean and SD of the responses)
    f0_mean = mean(beh_clean.recording_f0, 'omitmissing');
    f0_sd = std(beh_clean.recording_f0, 'omitmissing');

    beh_clean.recording_f0_z = (beh_clean.recording_f0 - f0_mean) ./ f0_sd; % Responses
    beh_clean.recording_f0_unaltered_z = (beh_clean.probe_unaltered_f0 - f0_mean) ./ f0_sd; % Unaltered Probe
    beh_clean.recording_f0_altered_z = (beh_clean.probe_altered_f0 - f0_mean) ./ f0_sd; % Altered Probe

    recording_f0_z = beh_clean.recording_f0_z(~isnan(beh_clean.recording_f0_z)); % Responses without NaNs

    recording_f0_z_no_probe = beh_clean.recording_f0_z(strcmp(beh_clean.probe_type, 'None')); % Responses during no-probe trials
    recording_f0_z_unaltered_probe = beh_clean.recording_f0_z(strcmp(beh_clean.probe_type, 'Normal')); % Responses during no-probe trials
    recording_f0_z_altered_probe = beh_clean.recording_f0_z(strcmp(beh_clean.probe_type, 'Pitch')); % Responses during no-probe trials

    %% Plots

    % Plot: Z-transformed F0 Distribution including Probe F0's
    figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
    [f, xi] = ksdensity(recording_f0_z); % Density of F0 vocal responses
    plot(xi, f, 'LineWidth', 2, 'Color','b');
    hold on;
    xline(beh_clean.recording_f0_unaltered_z(1), 'k--', 'LineWidth', 2); % F0 unaltered probe
    xline(beh_clean.recording_f0_altered_z(1), 'r--', 'LineWidth', 2); % F0 altered probe

    ylims = ylim; % Get current y-axis limits
    xlabel('F0 [Z-Transformed]');
    ylabel('Density');
    title(['Z-transformed F0 Density of Vocal Responses and F0 Values for Probes for ' subj]);

    legend('F0 Distribution Vocal Responses', 'F0 Unaltered Probe', 'F0 Altered Probe', 'Location', 'northwest', 'Interpreter', 'none');
    hold off;

    exportgraphics(gcf,fullfile(OUTPATH, ['tid_psam_z_f0_distribution_' subj '.png']),'Resolution',1000)

    % Plot: Z-transformed F0 Distribution including Probe F0's (binnedin one plot)
    % Prepare data
    n_bins = 10;
    min_trials_per_bin = 25;
    n_data = numel(recording_f0_z);
    bin_edges = round(linspace(1, n_data + 1, n_bins + 1));
    bin_counts = diff(bin_edges);

    % Prepare figure
    figure('Units', 'normalized', 'Position', [0.5, 0.5, 1.2, 1.2]);
    colors = parula(n_bins);
    hold on;

    % Plot density for each bin
    for bin_idx = 1:n_bins
        bin_data = recording_f0_z(bin_edges(bin_idx):bin_edges(bin_idx+1)-1);
        [f, xi] = ksdensity(bin_data);
        plot(xi, f, 'LineWidth', 2, 'Color', colors(bin_idx, :));
    end

    xline(beh_clean.recording_f0_unaltered_z(1), 'k--', 'LineWidth', 2); % F0 unaltered probe
    xline(beh_clean.recording_f0_altered_z(1), 'r--', 'LineWidth', 2); % F0 altered probe

    xlabel('F0 [Z-Transformed]');
    ylabel('Density');
    title(['Z-transformed F0 Distribution across Bins for ' subj]);
    legend_entries = arrayfun(@(i) sprintf('Bin %d', i), 1:n_bins, 'UniformOutput', false);
    legend([legend_entries, {'F0 Unaltered Probe', 'F0 Altered Probe'}], ...
        'Location', 'northeast', 'Interpreter', 'none');
    grid on;
    box on;
    hold off;

    exportgraphics(gcf, fullfile(OUTPATH, ['tid_psam_z_f0_distribution_binned_' subj '.png']), 'Resolution', 1000);

    % Plot: Z-transformed F0 Distribution including Probe F0's (binned in subplot)
    % Prepare data
    n_bins = 10;
    min_trials_per_bin = 25;
    n_data = numel(recording_f0_z);
    bin_edges = round(linspace(1, n_data + 1, n_bins + 1));
    bin_counts = diff(bin_edges);

    all_xi = [];
    for bin_idx = 1:n_bins
        bin_data = recording_f0_z(bin_edges(bin_idx):bin_edges(bin_idx+1)-1);
        [f, xi] = ksdensity(bin_data);
        all_xi = [all_xi, xi];  % Collect all xi values to determine the range
    end
    x_limits = [min(all_xi), max(all_xi)];  % Set consistent x-axis limits

    % Prepare figure
    figure('Units', 'normalized', 'Position', [0.5, 0.5, 1.2, 1.2]);
    colors = parula(n_bins);
    hold on;

    % Plot density for each bin
    for bin_idx = 1:n_bins
        bin_data = recording_f0_z(bin_edges(bin_idx):bin_edges(bin_idx+1)-1);
        [f, xi] = ksdensity(bin_data);

        subplot(n_bins, 1, bin_idx); 
        plot(xi, f, 'LineWidth', 2, 'Color', colors(bin_idx, :));

        xline(beh_clean.recording_f0_unaltered_z(1), 'k--', 'LineWidth', 2); % F0 unaltered probe
        xline(beh_clean.recording_f0_altered_z(1), 'r--', 'LineWidth', 2); % F0 altered probe

        xlabel('F0 [Z-Transformed]');
        ylabel('Density');
        title(['Z-transformed F0 Distribution for Bin ' num2str(bin_idx) ' n_trials = ' num2str(length(bin_data))]);
        grid on;
        box on;
        xlim(x_limits);
    end

    % Overall figure title
    sgtitle(['Z-transformed F0 Distribution across Bins for ' subj]);
    hold off;

    exportgraphics(gcf, fullfile(OUTPATH, ['tid_psam_z_f0_distribution_binned_subplots_' subj '.png']), 'Resolution', 1000);

    % Plot: Z-transformed F0 Distribution including Probe F0's (for no-probe, unaltered probe and altered probe trials)
    figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
    [f, xi] = ksdensity(recording_f0_z_no_probe); % Density of F0 vocal responses after no-probe trials
    plot(xi, f, 'LineWidth', 2, 'Color',[.7 .7 .7]);
    hold on;
    [f, xi] = ksdensity(recording_f0_z_unaltered_probe); % Density of F0 vocal responses after unaltered probe trials
    plot(xi, f, 'LineWidth', 2, 'Color','k');
    [f, xi] = ksdensity(recording_f0_z_altered_probe); % Density of F0 vocal responses after altered probe trials
    plot(xi, f, 'LineWidth', 2, 'Color','r');
    xline(beh_clean.recording_f0_unaltered_z(1), 'k--', 'LineWidth', 2); % F0 unaltered probe
    xline(beh_clean.recording_f0_altered_z(1), 'r--', 'LineWidth', 2); % F0 altered probe

    ylims = ylim; % Get current y-axis limits
    xlabel('F0 [Z-Transformed]');
    ylabel('Density');
    title(['Z-transformed F0 Density of Vocal Responses after different Probe Types and F0 Values for Probes for ' subj]);

    legend('No Probe','Unaltered Probe','Altered Probe', 'F0 Unaltered Probe', 'F0 Altered Probe', 'Location', 'northwest', 'Interpreter', 'none');
    hold off;

    exportgraphics(gcf,fullfile(OUTPATH, ['tid_psam_z_f0_distribution_probe_types_' subj '.png']),'Resolution',1000)

    %% Reaction-Time Analysis


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

check_done = 'tid_psam_vocal_analysis_DONE'

delete(wb)




% % [h_unaltered, p_unaltered, ci_unaltered, stats_unaltered] = ttest(subj_full_cleaned.recording_f0, subj_full_cleaned.recording_f0_unaltered(1));
% % cohens_d_unaltered = stats_unaltered.tstat / sqrt(stats_unaltered.df + 1);
% % [h_pitched, p_pitched, ci_pitched, stats_pitched] = ttest(subj_full_cleaned.recording_f0, subj_full_cleaned.recording_f0_pitched(1));
% % cohens_d_pitched = stats_pitched.tstat / sqrt(stats_pitched.df + 1);


