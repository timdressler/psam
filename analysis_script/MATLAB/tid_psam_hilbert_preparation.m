% tid_psam_hilbert_preparation.m
%
% Performs a hilbert transformation later used for the SVM analysis and saves data.
%
% Preparation includes the following steps
    % Removes EOG electrodes
%
% Stores data
%
%
% Tim Dressler, 14.05.2025

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
IDPATH = fullfile(MAINPATH, 'data\processed_data\svm_preprocessed_clean\'); % Only used to get subject IDs of non-excluded participants, data is loaded from INPATH
INPATH = fullfile(MAINPATH, 'data\processed_data\svm_preprocessed_clean\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\hilbert_prepared_clean');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EOG_CHAN = {'E29','E30'}; % Labels of EOG electrodes
EVENTS = {'go_act', 'go_pas'};
APLHA_FROM = 8;
APLHA_TILL = 13;
BETA_FROM = 14;
BETA_TILL = 30;
GAMMA_FROM = 31;
GAMMA_TILL = 37;

% Get directory content
dircont_subj = dir(fullfile(IDPATH, 'sub-*.set'));

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_hilbert_preparation.m');

% Aggregate frequency bands
freq_bands = {[APLHA_FROM APLHA_TILL], [BETA_FROM BETA_TILL], [GAMMA_FROM GAMMA_TILL]};

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_hilbert_preparation.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load epoched data
    EEG = pop_loadset('filename',[subj '_svm_preprocessed_clean.set'],'filepath',INPATH);

    % Remove EOG channels as they are not used for classification
    EEG = pop_select( EEG, 'rmchannel',EOG_CHAN);

    % Save data set

    % Loop over frequency bands
    for freq_band_num = 1:length(freq_bands)

        % Apply BP-Filter

        % Hilbert transform

        % Epoch

        % Reject bad epochs
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

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'tid_psam_hilbert_preparation.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_hilbert_preparation_marked_subj.xlsx'))
end

check_done = 'tid_psam_svm_preparation_DONE'

delete(wb)