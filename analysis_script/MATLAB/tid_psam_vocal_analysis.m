% tid_psam_vocal_analysis.m
%
% Performs analysis of vocal data.
%
% Tim Dressler, 07.04.2025

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
    % % text(subj_full_cleaned.f0_tab_pitched(1)+0.25, ylims(2)*0.9, ... % Test-Statistic for altered probe
    % %     sprintf('t=%.2f, p=%.3f, d=%.2f', stats_pitched.tstat, p_pitched, cohens_d_pitched), ...
    % %     'Color', 'red', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    % % 
    % % text(subj_full_cleaned.f0_tab_normal(1)+0.25, ylims(2)*0.8, ... % Test-Statistic for unaltered probe
    % %     sprintf('t=%.2f, p=%.3f, d=%.2f', stats_normal.tstat, p_normal, cohens_d_normal), ...
    % %     'Color', 'blue', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    xlabel('F0 [Hz]');
    ylabel('Density');
    title(['F0 Density of Vocal Responses and F0 Values for Probes for ' subj]);

    legend('F0 Distribution Vocal Responses', 'F0 Altered Probe', 'F0 Unaltered Probe', 'Location', 'southwest', 'Interpreter', 'none');
    hold off;

    exportgraphics(gcf,fullfile(OUTPATH, ['tid_psam_f0_distribution_' subj '.png']),'Resolution',1000)