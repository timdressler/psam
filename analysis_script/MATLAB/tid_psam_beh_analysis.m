% tid_psam_beh_analysis.m
%
% Performs analysis of vocal data.
%
% Processing includes the following steps 
%
    % Plots each subjects' z-transformed F0 values of vocal responses relative 
    %   to the F0 values of each subjects' probes
    % Plots each subjects' z-transformed F0 values of vocal responses relative 
    %   to the F0 values of each subjects' probes for each block (in one plot 
    %   and in a subplot)
    % Plots each subjects' z-transformed F0 values of vocal responses for
    %   trials including no probe, an unaltered probe and an altered probe
    % Plots each subject vocal onset time for no-probe trials and probe
    %   trials, for trials including early and late probes and for all
    %   trials
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
N_BLOCKS = 8; % Shoud be left at 8

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

    % F0 Analysis

    recording_f0_z = beh_clean.recording_f0_z(~isnan(beh_clean.recording_f0_z)); % Responses without NaNs
    recording_f0_z_no_probe = beh_clean.recording_f0_z(strcmp(beh_clean.probe_type, 'None')); % Responses during no-probe trials
    recording_f0_z_unaltered_probe = beh_clean.recording_f0_z(strcmp(beh_clean.probe_type, 'Unaltered')); % Responses during unaltered probe trials
    recording_f0_z_altered_probe = beh_clean.recording_f0_z(strcmp(beh_clean.probe_type, 'Altered')); % Responses during altered probe trials

    % Plots

    % Plot: Z-transformed F0 Distribution including Probe F0's
    figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
    [f, xi] = ksdensity(recording_f0_z); % Density of F0 vocal responses
    plot(xi, f, 'LineWidth', 2, 'Color','b');
    hold on;
    xline(beh_clean.probe_f0_unaltered_z(1), 'k--', 'LineWidth', 2); % F0 unaltered probe
    xline(beh_clean.probe_f0_altered_z(1), 'r--', 'LineWidth', 2); % F0 altered probe

    xlabel('F0 [Z-Transformed]');
    ylabel('Density');
    title(['Z-transformed F0 Density of Vocal Responses and F0 Values for Probes for ' subj]);

    legend('F0 Distribution Vocal Responses', 'F0 Unaltered Probe', 'F0 Altered Probe', 'Location', 'northwest', 'Interpreter', 'none');
    hold off;

    saveas(gcf,fullfile(OUTPATH, ['tid_psam_z_f0_distribution_' subj '.png']))

    % Plot: Z-transformed F0 Distribution including Probe F0's (binnedin one plot)
    % Prepare figure
    figure('Units', 'normalized', 'Position', [0.5, 0.5, 1.2, 1.2]);
    colors = parula(N_BLOCKS);
    hold on;

    % Plot density for each bin
    for block_idx = 1:N_BLOCKS
        block_data = beh_clean.recording_f0_z(beh_clean.block == block_idx);
        [f, xi] = ksdensity(block_data);
        plot(xi, f, 'LineWidth', 2, 'Color', colors(block_idx, :));
    end

    xline(beh_clean.probe_f0_unaltered_z(1), 'k--', 'LineWidth', 2); % F0 unaltered probe
    xline(beh_clean.probe_f0_altered_z(1), 'r--', 'LineWidth', 2); % F0 altered probe

    xlabel('F0 [Z-Transformed]');
    ylabel('Density');
    title(['Z-transformed F0 Distribution across Bins for ' subj]);
    legend_entries = arrayfun(@(i) sprintf('Block %d', i), 1:N_BLOCKS, 'UniformOutput', false);
    legend([legend_entries, {'F0 Unaltered Probe', 'F0 Altered Probe'}], ...
        'Location', 'northeast', 'Interpreter', 'none');
    grid on;
    box on;
    hold off;

    saveas(gcf, fullfile(OUTPATH, ['tid_psam_z_f0_distribution_blocks_' subj '.png']));

    % Plot: Z-transformed F0 Distribution including Probe F0's (binned in subplots)

    % Prepare data
    all_xi = [];
    for block_idx = 1:N_BLOCKS
        block_data = beh_clean.recording_f0_z(beh_clean.block == block_idx);
        [f, xi] = ksdensity(block_data);
        all_xi = [all_xi, xi];  % Collect all xi values to determine the range
    end
    all_xi = [all_xi, beh_clean.probe_f0_unaltered_z(1), beh_clean.probe_f0_altered_z(1)]; % Add x value of the probes
    x_limits = [min(all_xi), max(all_xi)];  % Set consistent x-axis limits

    % Prepare figure
    figure('Units', 'normalized', 'Position', [0.5, 0.5, 1.2, 1.2]);
    colors = parula(N_BLOCKS);
    hold on;

    % Plot density for each bin
    for block_idx = 1:N_BLOCKS
        block_data = beh_clean.recording_f0_z(beh_clean.block == block_idx);
        [f, xi] = ksdensity(block_data);

        subplot(N_BLOCKS, 1, block_idx);
        plot(xi, f, 'LineWidth', 2, 'Color', colors(block_idx, :));

        xline(beh_clean.probe_f0_unaltered_z(1), 'k--', 'LineWidth', 2); % F0 unaltered probe
        xline(beh_clean.probe_f0_altered_z(1), 'r--', 'LineWidth', 2); % F0 altered probe

        xlabel('F0 [Z-Transformed]');
        ylabel('Density');
        title(['Z-transformed F0 Distribution for Block ' num2str(block_idx) ' n_trials = ' num2str(length(block_data))]);
        grid on;
        box on;
        xlim(x_limits);
    end

    % Overall figure title
    sgtitle(['Z-transformed F0 Distribution across Bins for ' subj]);
    hold off;

    saveas(gcf, fullfile(OUTPATH, ['tid_psam_z_f0_distribution_blocks_subplots_' subj '.png']));

    % Plot: Z-transformed F0 Distribution including Probe F0's (for no-probe, unaltered probe and altered probe trials)
    figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
    [f, xi] = ksdensity(recording_f0_z_no_probe); % Density of F0 vocal responses after no-probe trials
    plot(xi, f, 'LineWidth', 2, 'Color',[.7 .7 .7]);
    hold on;
    [f, xi] = ksdensity(recording_f0_z_unaltered_probe); % Density of F0 vocal responses after unaltered probe trials
    plot(xi, f, 'LineWidth', 2, 'Color','k');
    [f, xi] = ksdensity(recording_f0_z_altered_probe); % Density of F0 vocal responses after altered probe trials
    plot(xi, f, 'LineWidth', 2, 'Color','r');
    xline(beh_clean.probe_f0_unaltered_z(1), 'k--', 'LineWidth', 2); % F0 unaltered probe
    xline(beh_clean.probe_f0_altered_z(1), 'r--', 'LineWidth', 2); % F0 altered probe

    xlabel('F0 [Z-Transformed]');
    ylabel('Density');
    title(['Z-transformed F0 Distribution of Vocal Responses after different Probe Types and F0 Values for Probes for ' subj]);

    legend('No Probe','Unaltered Probe','Altered Probe', 'F0 Unaltered Probe', 'F0 Altered Probe', 'Location', 'northwest', 'Interpreter', 'none');
    hold off;

    saveas(gcf,fullfile(OUTPATH, ['tid_psam_z_f0_distribution_probe_types_' subj '.png']))

    % Vocal Onset Time Analysis

    recording_vot_early_probe = beh_clean.recording_vot(strcmp(beh_clean.probe_onset_cat, 'Early')); % Responses during early probe trials
    recording_vot_late_probe = beh_clean.recording_vot(strcmp(beh_clean.probe_onset_cat, 'Late')); % Responses during late probe trials
    recording_vot_no_probe = beh_clean.recording_vot(strcmp(beh_clean.probe, 'No')); % Responses during no-probe trials
    recording_vot_probe = beh_clean.recording_vot(strcmp(beh_clean.probe, 'Yes')); % Responses during probe trials

    % Plots

    % Plot: Vocal onset time distribution for all data, all probe-trials and no-probe trials and for late- and early probe-trials
    figure('Units', 'normalized', 'Position', [0.5, 0.5, 1.2, 1.2]);

    % Subplot: Vocal onset time distribution for late- and early probe-trials
    subplot(2,2,1);
    hold on
    [f_early, xi_early] = ksdensity(recording_vot_early_probe);
    [f_late, xi_late] = ksdensity(recording_vot_late_probe);
    plot(xi_early, f_early, 'Color', 'g', 'LineWidth', 2);
    plot(xi_late, f_late, 'Color', 'm', 'LineWidth', 2);
    xline(mean(recording_vot_early_probe, 'omitnan'), 'g--', 'LineWidth', 2); % RT early probe
    xline(mean(recording_vot_late_probe, 'omitnan'), 'm--', 'LineWidth', 2); % RT late probe

    legend({'Early', 'Late', 'Mean Early', 'Mean Late'});
    xlabel('RT [ms]');  
    ylabel('Density');
    title('Distribution of Recording VOT for early and late probe trials');
    box off
    hold off

    % Subplot: Vocal onset time distribution for no-probe and probe trials
    subplot(2,2,2);
    hold on
    [f_probe, xi_probe] = ksdensity(recording_vot_probe);
    [f_no_probe, xi_no_probe] = ksdensity(recording_vot_no_probe);
    plot(xi_probe, f_probe, 'Color', 'k', 'LineWidth', 2);
    plot(xi_no_probe, f_no_probe, 'Color', [.7 .7 .7], 'LineWidth', 2);
    xline(mean(recording_vot_probe, 'omitnan'), 'k--', 'LineWidth', 2); % RT probe trials
    xline(mean(recording_vot_no_probe, 'omitnan'), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2); % RT no-probe trials

    legend({'Probe', 'No-Probe', 'Mean Probe', 'Mean No-Probe'});
    xlabel('RT [ms]');  
    ylabel('Density');
    title('Distribution of Recording VOT for probe and no-probe trials');
    box off
    hold off

    % Subplot: Vocal onset time distribution for all trials
    subplot(2,2,[3 4]);
    hold on;
    [f_early, xi_early] = ksdensity(beh_clean.recording_vot);
    plot(xi_early, f_early, 'Color', 'k', 'LineWidth', 2);
    xline(mean(beh_clean.recording_vot, 'omitnan'), 'k--', 'LineWidth', 2); % RT 

    legend({'All Trials', 'Mean'});
    xlabel('RT [ms]'); 
    ylabel('Density');
    title('Distribution of Recording RT for all trials');
    box off;
    hold off;

    sgtitle(['VOT Distributions for ' subj]);

    saveas(gcf,fullfile(OUTPATH, ['tid_psam_rt_distributions_' subj '.png']))

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

delete(wb); close all;





