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
EOG_CHAN = {'E29','E30'}; % Labels of EOG electrodes
WIN_EARLY_FROM = -600;
WIN_EARLY_TILL = -401;
WIN_LATE_FROM = -400;
WIN_LATE_TILL = -201;
ZEROPADDING = 1000; % Zeropadding to 1 s (1000 Samples due to 1000 Hz SR)
APLHA_FROM = 8;
APLHA_TILL = 13;
BETA_FROM = 14;
BETA_TILL = 30;
GAMMA_FROM = 31;
GAMMA_TILL = 37;

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

    % Remove EOG channels as they are not used for classification
    EEG = pop_select( EEG, 'rmchannel',EOG_CHAN);

    % Get epoch labels
    labels = cellfun(@(x) x{2}, {EEG.epoch.eventtype}, 'UniformOutput', false);

    % Define time windows
    windows = struct( ...
        'early', struct('from', WIN_EARLY_FROM, 'till', WIN_EARLY_TILL), ...
        'late', struct('from', WIN_LATE_FROM, 'till', WIN_LATE_TILL) ...
        );

    % Initialize feature struct
    features = struct();

    % Loop through windows
    for win_label = ["early", "late"]
        label = char(win_label);

        % Get time indices
        [~,start_idx] = min(abs(EEG.times - windows.(label).from));
        [~,end_idx] = min(abs(EEG.times - windows.(label).till));

        % Get window data
        data = EEG.data(:, start_idx:end_idx, :);
        features.(label).data = data;

        % Sanity Check: Correct dimensions
        if size(data,1) ~= EEG.nbchan || size(data,2) ~= 200 || size(data,3) == 1
            marked_subj{end+1,1} = subj;
            marked_subj{end,2} = 'wrong_dimensions';
        end

        % Time-domain features
        features.(label).mean = squeeze(mean(data, 2));
        features.(label).rms = squeeze(rms(data, 2));
        features.(label).sd = squeeze(std(data, [], 2));
        features.(label).min = squeeze(min(data, [], 2));
        features.(label).max = squeeze(max(data, [], 2));
        features.(label).kurtosis = squeeze(kurtosis(data, [], 2));
        features.(label).skewness = squeeze(skewness(data, [], 2));

        for i = 1:size(data, 3)
            features.(label).zerocrossing(:, i) = zerocrossrate(data(:,:,i)')';
        end

        % Frequency-domain features
        % Preparation
        % Apply Hamming window to each epoch for each electrode
        ham = hamming(size(data,2));

        for e = 1:size(data,3)
            for ch = 1:size(data,1)
                data(ch, :, e) = squeeze(data(ch, :, e)) .* ham';
            end
        end

        % Apply zeropadding to 1 s (i.e., 1000 Samples due to SR = 1000 Hz)
        padding_needed = ZEROPADDING - size(data,2);
        data_padded = padarray(data, [0, padding_needed, 0], 0, 'post');

        % Sanity check: Correct dimensions after padding
        if size(data_padded,1) ~= EEG.nbchan || size(data_padded,2) ~= ZEROPADDING || size(data_padded,3) == 1
            marked_subj{end+1,1} = subj;
            marked_subj{end,2} = 'wrong_dimensions_after_padding';
        end
        if ~all(data_padded(201:end,:) == 0, 'all')
            marked_subj{end+1,1} = subj;
            marked_subj{end,2} = 'wrong_values_after_padding';
        end

        % Apply FFT
        fft_data = abs(fft(data_padded, [], 2)) / size(data, 2);
        fft_data = fft_data(:, 1:end/2, :);
        fft_data(:, 2:end, :) = fft_data(:, 2:end, :) * 2;

        % Get bandpower for alpha, beta and gamma bands
        freq_vec = 0:1/(ZEROPADDING/1000):EEG.srate/2 - (1/(ZEROPADDING/1000));
        [~,a_start] = min(abs(freq_vec - APLHA_FROM));
        [~,a_end] = min(abs(freq_vec - APLHA_TILL));
        [~,b_start] = min(abs(freq_vec - BETA_FROM));
        [~,b_end] = min(abs(freq_vec - BETA_TILL));
        [~,g_start] = min(abs(freq_vec - GAMMA_FROM));
        [~,g_end] = min(abs(freq_vec - GAMMA_TILL));

        features.(label).alpha = squeeze(mean(fft_data(:, a_start:a_end, :), 2));
        features.(label).beta = squeeze(mean(fft_data(:, b_start:b_end, :), 2));
        features.(label).gamma = squeeze(mean(fft_data(:, g_start:g_end, :), 2));

        % Hjorth parameters
        % Get Activity
        features.(label).activity = tid_psam_hjorth_activity_TD(data);
        % Get Mobility
        features.(label).mobility = tid_psam_hjorth_mobility_TD(data);
        % Get Complexity
        features.(label).complexity = tid_psam_hjorth_complexity_TD(data);
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
writetable(protocol,fullfile(OUTPATH, 'tid_psam_svm_preparation.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_svm_preparation_marked_subj.xlsx'))
end

check_done = 'tid_psam_svm_preparation_DONE'

delete(wb)

%% 
plot(freq_vec,fft_data(1,:,1))