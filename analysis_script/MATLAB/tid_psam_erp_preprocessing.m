% tid_psam_erp_preprocessing.m
%
% Performs ERP preprocessing and saves datasets.
%
% Conditions
% S 931 - Active - Early Probe - Unaltered
% S 932 - Active - Early Probe - Altered
% S 933 - Active - Late Probe - Unaltered
% S 934 - Active - Late Probe - Altered
%
% S 941 - Passive - Early Probe - Unaltered
% S 942 - Passive - Early Probe - Altered
% S 943 - Passive - Late Probe - Unaltered
% S 944 - Passive - Late Probe - Altered
%
% Preprocessing includes the following steps
%
% Loads raw data
% Renames the events
% Applies band-pass filter
% Performs PREP pipeline (Bigdely-Shamlo et al., 2015)
%
% Performs ICA-specific processing
% Loads raw datasets for talk & listen condition
% Merges datasets
% Applies band-pass filter
% Resamples dataset
% Calculates ICA weights
%
% Applies ICA to original data
% Marks and rejects bad components using the IC Label Plugin (Pion-Tonachini et al., 2019)
% Applies band-pass filter
% Epochs data
% Performs baseline correction
% Reject bad epochs using threshold and probability
%
% Stores dataset
%
%
% Literature
% Bigdely-Shamlo N, Mullen T, Kothe C, Su K-M and Robbins KA (2015).
% The PREP pipeline: standardized preprocessing for large-scale EEG analysis
% Front. Neuroinform. 9:16. doi: 10.3389/fninf.2015.00016
% Pion-Tonachini, L., Kreutz-Delgado, K., & Makeig, S. (2019).
% ICLabel: An automated electroencephalographic independent component classifier, dataset, and website.
% NeuroImage, 198, 181-197.
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
INPATH = fullfile(MAINPATH, 'data\processed_data\markers_included\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\erp_preprocessed');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EPO_FROM = -0.2;
EPO_TILL = 0.400;
LCF = 1;
HCF = 60;
LCF_2 = 1;
HCF_2 = 25;
LCF_ICA = 1;
HCF_ICA = 30;
BL_FROM = -200;
THRESH = 75;
SD_PROB = 3;
RESAM_ICA = 250;
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*.set'));

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_erp_preprocessing.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_erp_preprocessing.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load data
    EEG = pop_loadset('filename',[subj '_markers_inlcuded.set'],'filepath',INPATH);

    % Remove Marker-Channel
    EEG = pop_select( EEG, 'rmchannel',{'M'});

    % Rename events
    clear event
    for event = 1:length(EEG.event)
        switch EEG.event(event).type
            case 'S 931'
                EEG.event(event).type = 'act_early_unalt';
            case 'S 932'
                EEG.event(event).type = 'act_early_alt';
            case 'S 933'
                EEG.event(event).type = 'act_late_unalt';
            case 'S 934'
                EEG.event(event).type = 'act_late_alt';
            case 'S 941'
                EEG.event(event).type = 'pas_early_unalt';
            case 'S 942'
                EEG.event(event).type = 'pas_early_alt';
            case 'S 943'
                EEG.event(event).type = 'pas_late_unalt';
            case 'S 944'
                EEG.event(event).type = 'pas_late_alt';
        end
    end

    % Bandpass-Filter
    EEG = pop_eegfiltnew(EEG, 'locutoff',LCF,'hicutoff',HCF,'plotfreqz',0);

    % Run PREP pipeline (Bigdely-Shamlo et al., 2015)
    % Settings
    params = struct();
    params.lineFrequencies = 50:50:EEG.srate/2-50; % Set line noise to 50Hz
    params.ignoreBoundaryEvents = true;  % Ignore boundary events
    params.reportMode = 'skipReport';  % Suppress report
    params.keepFiltered = true; % Remove trend
    % Run
    EEG = pop_prepPipeline(EEG, params);
    % Store dataset
    EEG.setname = [subj '_after_PREP'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);


    % ICA
    % Reload and merge raw data
    EEG = pop_loadset('filename',[subj '_markers_inlcuded.set'],'filepath',INPATH);
    % Remove Marker-Channel
    EEG = pop_select( EEG, 'rmchannel',{'M'});
    % Bandpass-Filter
    EEG = pop_eegfiltnew(EEG, 'locutoff',LCF_ICA,'hicutoff',HCF_ICA,'plotfreqz',0);
    % Resample
    EEG = pop_resample( EEG, RESAM_ICA);
    % Run ICA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
    EEG.setname = [subj '_ICA_weights'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Attach ICA weight to main data
    EEG = ALLEEG(1);
    CURRENTSET = 1;
    EEG = pop_editset(EEG,'run', [], 'icaweights','ALLEEG(2).icaweights', 'icasphere','ALLEEG(2).icasphere');
    EEG.setname = [subj '_after_PREP_ICA'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Label ICA components with IC Label Plugin (Pion-Tonachini et al., 2019)
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, [0 0.2;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);
    EEG = pop_subcomp( EEG, [], 0);
    % Bandpass-Filter
    EEG = pop_eegfiltnew(EEG, 'locutoff',LCF_2,'hicutoff',HCF_2,'plotfreqz',0);

    % Epoching
    EEG = pop_epoch( EEG, EVENTS, [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');

    % Baseline-Removal
    EEG = pop_rmbase( EEG, [BL_FROM 0] ,[]);

    % Threshold removal
    EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan] ,-THRESH,THRESH,EPO_FROM,EPO_TILL,0,0);
    % Probability-based removal
    EEG = pop_jointprob(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
    EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);

    % Save dataset
    EEG.setname = [subj '_preprocessed'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_saveset(EEG, 'filename',[SUBJ '_erp_preprocessed.set'],'filepath', OUTPATH);

end



% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'erp_pilot_protocol.xlsx'))

close(wb)

check_done = 'tid_psam_erp_preprocessing_DONE'












