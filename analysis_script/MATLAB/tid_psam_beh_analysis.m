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
    % Plots F0 distribution and probe F0 values across all participants
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

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*clean.xlsx'));

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_beh_analysis.m');

clear subj_idx
all_recording_f0_z = {};
all_probe_f0_z = {};
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Get subject outpath
    subj_outpath = fullfile(OUTPATH,subj);
    mkdir(subj_outpath);

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_beh_analysis.m'])

    tic;

    % Load data
    beh_clean = readtable(fullfile(INPATH,[subj '_beh_preprocessed_clean.xlsx']));

    % F0 Analysis
    recording_f0_z = beh_clean.recording_f0_z(~isnan(beh_clean.recording_f0_z)); % Responses without NaNs
    recording_f0_z_no_probe = beh_clean.recording_f0_z(strcmp(beh_clean.probe_type, 'None')); % Responses during no-probe trials
    recording_f0_z_unaltered_probe = beh_clean.recording_f0_z(strcmp(beh_clean.probe_type, 'Unaltered')); % Responses during unaltered probe trials
    recording_f0_z_altered_probe = beh_clean.recording_f0_z(strcmp(beh_clean.probe_type, 'Altered')); % Responses during altered probe trials

    % Store values
    all_recording_f0_z{subj_idx,1} = subj;
    all_recording_f0_z{subj_idx,2} = recording_f0_z;
    all_probe_f0_z{subj_idx,1} = subj;
    all_probe_f0_z{subj_idx,2} = beh_clean.probe_f0_unaltered_z(1); % unaltered probe
    all_probe_f0_z{subj_idx,3} = beh_clean.probe_f0_altered_z(1); % altered probe

    % Plots
    if INDIVIDUAL_PLOTS
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

        saveas(gcf,fullfile(subj_outpath, [subj '_z_f0_distribution.png']))

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

        saveas(gcf, fullfile(subj_outpath, [subj '_z_f0_distribution_blocks.png']));

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

        saveas(gcf, fullfile(subj_outpath, [subj '_z_f0_distribution_blocks_subplots.png']));


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

        saveas(gcf,fullfile(subj_outpath, [subj '_z_f0_distribution_probe_types.png']))

    end

    % Vocal Onset Time Analysis

    recording_vot_early_probe = beh_clean.recording_vot(strcmp(beh_clean.probe_onset_cat, 'Early')); % Responses during early probe trials
    recording_vot_late_probe = beh_clean.recording_vot(strcmp(beh_clean.probe_onset_cat, 'Late')); % Responses during late probe trials
    recording_vot_no_probe = beh_clean.recording_vot(strcmp(beh_clean.probe, 'No')); % Responses during no-probe trials
    recording_vot_probe = beh_clean.recording_vot(strcmp(beh_clean.probe, 'Yes')); % Responses during probe trials

    % Plots
    if INDIVIDUAL_PLOTS
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

        saveas(gcf,fullfile(subj_outpath, [subj '_rt_distributions.png']))
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

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'tid_psam_beh_analysis_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_beh_analysis_marked_subj.xlsx'))
end

check_done = 'tid_psam_beh_analysis_DONE'

delete(wb); close all;

%% Plots
close all

% Plot: Z-transformed F0 Distribution including Probe F0's for all subjects
% Concatinate colors
colors = {
    main_yellow;
    main_red;
    main_green;
    main_blue;
    };

% Create F0 plot
figure('Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
x = linspace(-6, 6, 1000);  % Common F0 range for all subjects
hold on;
view(3);  % 3D view

for subj_idx= 1:length(dircont_subj)
    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Sanity Check: Matching F0 data and probe data
    subj_idx_f0 = find(strcmp(all_recording_f0_z(:,1), subj));
    subj_idx_probe = find(strcmp(all_probe_f0_z(:,1), subj));
    if ~(subj_idx_f0 == subj_idx_probe)
        error('F0 data and probe data do not match')
    end

    % Z-transformed F0 data
    f0_data = all_recording_f0_z{subj_idx_f0,2}';
    [density, xi] = ksdensity(f0_data, x);

    % Z-transformed Probe values 
    unaltered_probe = all_probe_f0_z{subj_idx_probe,2};
    altered_probe = all_probe_f0_z{subj_idx_probe,3};

    % Plot F0 distribution 
    y = subj_idx * ones(size(xi)); % Create Y-axis vector (same for whole curve)
    fill3([xi, fliplr(xi)], ...
          [y, fliplr(y)], ...
          [density, zeros(size(density))], ...
          main_blue, 'FaceAlpha', 0.5, 'EdgeColor', main_blue); 
    plot3(xi, y, density,'Color' ,main_blue, 'LineWidth', 1.5); 

    % Plot probe F0 values
    plot3([unaltered_probe unaltered_probe], [subj_idx subj_idx], [0 max(density)], '--', 'Color','k', 'LineWidth', 1.2);
    plot3([altered_probe altered_probe], [subj_idx subj_idx], [0 max(density)],'--','Color', main_red, ...
        'LineWidth', 1.2);
end

xlabel('F0 [Z-Transformed]');
%%ylabel('Subject');
zlabel('Density');
title('F0 Distributions and Probe F0s Across Subjects');
grid on;
view(-12.6193,77.8152) % for fixed 3D view

% Adjust y-axis ticks and labels
yticks(1:length(dircont_subj));  % only integer ticks for subjects
subj_labels = cell(length(dircont_subj),1);
for idx = 1:length(dircont_subj)
    s = dircont_subj(idx).name;
    s = regexp(s, 'sub-\d+', 'match', 'once');
    subj_labels{idx} = s;
end
yticklabels(subj_labels);

% Save plot
saveas(gcf,fullfile(OUTPATH, 'tid_psam_all_z_f0_distribution.png'))

