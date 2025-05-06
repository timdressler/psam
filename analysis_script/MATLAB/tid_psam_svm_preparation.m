% tid_psam_SVM_preparation.m
%
% Performs preparations for SVM analysis and saves data.
%
% Preparation includes the following steps
%
% Load data and data set containing ICA weights (see
%
% Stores data
%
% Note. Trials are only marked for exclusion but not excluded yet!
%
% Literature
% Pion-Tonachini, L., Kreutz-Delgado, K., & Makeig, S. (2019).
%   ICLabel: An automated electroencephalographic independent component classifier, dataset, and website.
%   NeuroImage, 198, 181-197.
%
% Tim Dressler, 05.05.2025

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
INPATH = fullfile(MAINPATH, 'data\processed_data\svm_preprocessed_clean\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\svm_prepared_clean');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
WIN_EARLY_FROM = -600;
WIN_EARLY_TILL = -401;
WIN_LATE_FROM = -400;
WIN_LATE_TILL = -201;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*.set'));

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_svm_preparation.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_svm_preparation.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load data
    EEG = pop_loadset('filename',[subj '_svm_preprocessed_clean.set'],'filepath',INPATH);

    % Extract data from early and late time window
    [~,win_early_start] = min(abs(EEG.times-WIN_EARLY_FROM)); % Start of the early window
    [~,win_early_end] = min(abs(EEG.times-WIN_EARLY_TILL)); % End of the early window
    [~,win_late_start] = min(abs(EEG.times-WIN_LATE_FROM)); % Start of the late window
    [~,win_late_end] = min(abs(EEG.times-WIN_LATE_TILL)); % End of the late window

    win_early_data = EEG.data(:,win_early_start:win_early_end,:);
    win_late_data = EEG.data(:,win_late_start:win_late_end,:);

    % Extract the features
    % Time-domain features
    % Mean amplitude
    win_early_mean = squeeze(mean(win_early_data, 2));

    % RMS
    win_early_rms = squeeze(rms(win_early_data, 2));

    % SD of the amplitude
    win_early_sd = squeeze(std(win_early_data, [], 2));

    % Maximum amplitude
    win_early_min = squeeze(min(win_early_data, [], 2));

    % Minimum amplitude
    win_early_max = squeeze(max(win_early_data, [], 2));

    % Kurtosis
    win_early_kurtosis = squeeze(kurtosis(win_early_data, [], 2));

    % Skewness
    win_early_skewness = squeeze(skewness(win_early_data, [], 2));

    % Zero-crossing rate
    for i = 1:size(win_early_data,3)
        win_early_zerocrossing(:,i) = zerocrossrate(win_early_data(:,:,i)')';
    end

    % Frequency-domain features
    % Preparations

    % Bandpower for alpha, beta and (low) gamma bands

    % Hjorth Parameters
    % Activity
    win_early_activity = tid_psam_hjorth_activity_TD(win_early_data);

    % Mobility
    win_early_mobility = tid_psam_hjorth_mobility_TD(win_early_data);

    % Complexity
    win_early_complexity = tid_psam_hjorth_complexity_TD(win_early_data);








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
writetable(protocol,fullfile(OUTPATH, 'tid_psam_svm_preparation.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_svm_preparation_marked_subj.xlsx'))
end

check_done = 'tid_psam_svm_preparation_DONE'

delete(wb)


