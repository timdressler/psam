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
%%wb = waitbar(0,'starting tid_psam_vocal_analysis.m');

clear subj_idx 
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    %%waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_vocal_analysis.m'])

    tic;
    % Load data
    beh_clean = readtable(fullfile(INPATH,[subj '_beh_preprocessed_clean.xlsx']));

    %% F0 Analysis

    % Z-Transformation F0 data 
    % Z-transform F0's of responses and probes based on the mean and SD of the responses)
    f0_mean = mean(beh_clean.recording_f0, 'omitmissing');
    f0_sd = std(beh_clean.recording_f0, 'omitmissing');
    beh_clean.recording_f0_z = (beh_clean.recording_f0 - f0_mean) ./ f0_sd; % Responses
    beh_clean.recording_f0_normal_z = (beh_clean.probe_unaltered_f0 - f0_mean) ./ f0_sd; % Unaltered Probe
    beh_clean.recording_f0_altered_z = (beh_clean.probe_altered_f0 - f0_mean) ./ f0_sd; % Altered Probe

    % Plot: Z-transformed F0 Distribution including Probe F0's
    figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
    [f, xi] = ksdensity(beh_clean.recording_f0_z); % Density of F0 vocal responses
    plot(xi, f, 'LineWidth', 2, 'Color','k');
    hold on;
    xline(beh_clean.recording_f0_normal_z(1), 'r--', 'LineWidth', 2); % F0 altered probe
    xline(beh_clean.recording_f0_altered_z(1), 'b--', 'LineWidth', 2); % F0 unaltered probe

    ylims = ylim; % Get current y-axis limits
    xlabel('F0 [Z-Transformed]');
    ylabel('Density');
    title(['Z-transformed F0 Density of Vocal Responses and F0 Values for Probes for ' subj]);

    legend('F0 Distribution Vocal Responses', 'F0 Unaltered Probe', 'F0 Altered Probe', 'Location', 'northwest', 'Interpreter', 'none');
    hold off;

    exportgraphics(gcf,fullfile(OUTPATH, ['tid_psam_z_f0_distribution_' subj '.png']),'Resolution',1000)

    % Save z-transformed F0 values for unaltered and altered probes






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




% % [h_normal, p_normal, ci_normal, stats_normal] = ttest(subj_full_cleaned.recording_f0, subj_full_cleaned.recording_f0_normal(1));
% % cohens_d_normal = stats_normal.tstat / sqrt(stats_normal.df + 1);
% % [h_pitched, p_pitched, ci_pitched, stats_pitched] = ttest(subj_full_cleaned.recording_f0, subj_full_cleaned.recording_f0_pitched(1));
% % cohens_d_pitched = stats_pitched.tstat / sqrt(stats_pitched.df + 1);


