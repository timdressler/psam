% tid_psam_ica_preprocessing.m
%
% Performs ICA preprocessing and saves datasets.
%
% Preprocessing includes the following steps
%
%   Apply a 1 Hz HP-Filter
%   Idenitify and remove bad channels using the
%       bemobil_detect_bad_channels() funtion
%   Create regular 1s epochs
%   Remove bad epochs based on probability
%   Calculate ICA weights
%   Label bad components using the ICLabel Plugin (Pion-Tonachini et al., 2019)
%
% Saves data
%
% Literature
% Pion-Tonachini, L., Kreutz-Delgado, K., & Makeig, S. (2019).
%   ICLabel: An automated electroencephalographic independent component classifier, dataset, and website.
%   NeuroImage, 198, 181-197.
%
% Tim Dressler, 04.04.2025

clear
close all
clc

rng(123)
set(0,'DefaultTextInterpreter','none')

% Set up paths
SCRIPTPATH = cd;
if regexp(SCRIPTPATH, regexptranslate('wildcard','*psam\analysis_script\MATLAB')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, '\analysis_script\MATLAB');
INPATH = fullfile(MAINPATH, 'data\processed_data\markers_included\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\ica_preprocessed');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
%%tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EOG_CHAN = {'E29','E30'}; % Labels of EOG electrodes
EPO_FROM = -0.2;
EPO_TILL = 0.400;
LCF_ICA = 1;
BL_FROM = -200;
THRESH = 75;
SD_PROB_ICA = 3;
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*.set'));

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_ica_preprocessing.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_ica_preprocessing.m'])

    % Check if subject was already processed
    subj_file = fullfile(OUTPATH, [subj '_ica_weights.set']);
    if exist(subj_file, 'file')
        disp(['Skipping ' subj ' (already run)'])
        protocol{subj_idx,1} = subj;
        protocol{subj_idx,2} = NaN;
        protocol{subj_idx,3} = 'SKIPPED';
        continue
    end

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % ICA-specific preprocessing
    % Load raw data
    EEG = pop_loadset('filename',[subj '_markers_inlcuded.set'],'filepath',INPATH);

    % Remove Marker-Channel
    EEG = pop_select( EEG, 'rmchannel',{'M'});

    % Add channel locations
    EEG.chanlocs = readlocs( fullfile(MAINPATH,'\config\elec_96ch_adapted.elp'));
    EEG = eeg_checkset( EEG );

    % Add type = EOG for EOG electrodes for bemobil_detect_bad_channels to ignore them
    eog_chani = find(ismember({EEG.chanlocs.labels}, EOG_CHAN));
    [EEG.chanlocs(eog_chani).type] = deal('EOG');
    EEG = eeg_checkset( EEG );

    EEG.setname = [subj '_ready_for_ICA_preprocessing'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Highpass-Filter
    LCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth_TD(LCF_ICA)); % Get filter order (also see pop_firwsord.m, tid_psam_get_transition_bandwidth_TD.m)
    EEG = pop_firws(EEG, 'fcutoff', LCF_ICA, 'ftype', 'highpass', 'wtype', 'hamming', 'forder',LCF_ord, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

    % Remove bad channels
    EEG.badchans = [];
    [chans_to_interp, chan_detected_fraction_threshold, detected_bad_channels] = bemobil_detect_bad_channels(EEG, ALLEEG, 1, [], [], 5);
    EEG.badchans = chans_to_interp;
    EEG = pop_select(EEG,'nochannel', EEG.badchans);

    % Create 1s epochs & remove bad ones
    EEG = eeg_regepochs(EEG);

    EEG = pop_jointprob(EEG,1,[1:EEG.nbchan] ,SD_PROB_ICA,0,0,0,[],0);
    EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,SD_PROB_ICA,0,0,0,[],0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);

    EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);

    % Run ICA
    rng(123) % Set seed again
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');

    % Label ICA components with IC Label Plugin (Pion-Tonachini et al., 2019)
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, [0 0;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1]);
    % Sanity Check: Plot flagged ICs
    tid_psam_plot_flagged_ICs_TD(EEG,['ICs for ' subj], 'SavePath' ,fullfile(OUTPATH, [subj '_ic_topo.png']), 'PlotOn', false)

    % Save dataset
    EEG.setname = [subj '__ICA_weights'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_saveset(EEG, 'filename',[subj '_ica_weights.set'],'filepath', OUTPATH);

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
writetable(protocol,fullfile(OUTPATH, 'tid_psam_ica_preprocessing_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_ica_preprocessing_marked_subj.xlsx'))
end

check_done = 'tid_psam_ica_preprocessing_DONE'

delete(wb)