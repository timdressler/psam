% tid_psam_erp_preprocessing.m
%
% Performs ERP preprocessing and saves datasets.
%
% Conditions
% S 931 = Active - Early Probe - Unaltered
% S 932 = Active - Early Probe - Altered
% S 933 = Active - Late Probe - Unaltered
% S 934 = Active - Late Probe - Altered
%
% S 941 = Passive - Early Probe - Unaltered
% S 942 = Passive - Early Probe - Altered
% S 943 = Passive - Late Probe - Unaltered
% S 944 = Passive - Late Probe - Altered
%
% Preprocessing includes the following steps
%
    % Load data and data set containing ICA weights (see
    %   tid_psam_ica_preprocessing.m)
    % Rename events
    % Apply a 0.3 Hz HP-Filter
    % Apply a 30 Hz LP-Filter
    % Remove bad channels as identified in tid_psam_ica_preprocessing.m
    % Attach ICA weights and remove bad components using the ICLabel Plugin
    %   (Pion-Tonachini et al., 2019)
    % Interpolate bad (and removed) channels
    % Epoch data around probe events (see above)
    % Apply baseline correction
    % Mark bad epochs based on probability
%
% Saves data
%
% Note. Trials are only marked for exclusion but not excluded yet!
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

set(0,'DefaultTextInterpreter','none')

% Set up paths
SCRIPTPATH = cd;
if regexp(SCRIPTPATH, regexptranslate('wildcard','*psam\analysis_script\MATLAB')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, '\analysis_script\MATLAB');
INPATH_RAW = fullfile(MAINPATH, 'data\processed_data\markers_included\');
INPATH_ICA = fullfile(MAINPATH, 'data\processed_data\ica_preprocessed\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\erp_preprocessed');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_RAW, INPATH_ICA, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EOG_CHAN = {'E29','E30'}; % Labels of EOG electrodes
EPO_FROM = -0.2;
EPO_TILL = 0.400;
LCF = 0.3;
HCF = 30;
BL_FROM = -200;
SD_PROB = 3;
SD_PROB_ICA = 3;
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};

% Get directory content
dircont_subj = dir(fullfile(INPATH_RAW, 'sub-*.set'));

% Sanity Check: Same number of files for raw data and ICA data
if length(dir(fullfile(INPATH_RAW, 'sub-*.set'))) == length(dir(fullfile(INPATH_ICA, 'sub-*.set')))
else
    error('Number of raw data files and ICA data files does not match')
end

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

    % Load ICA data
    EEG = pop_loadset('filename',[subj '_ica_weights.set'],'filepath',INPATH_ICA);
    EEG.setname = [subj '_ICA_weights'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Get bad channels as identified in tid_psam_ica_preprocessing.m
    chans_to_interp = EEG.badchans;

    % Load raw data
    EEG = pop_loadset('filename',[subj '_markers_inlcuded.set'],'filepath',INPATH_RAW);

    % Remove Marker-Channel
    EEG = pop_select( EEG, 'rmchannel',{'M'});

    % Add channel locations
    EEG.chanlocs = readlocs( fullfile(MAINPATH,'\config\elec_96ch_adapted.elp'));
    EEG = eeg_checkset( EEG );

    % Add type = EOG for EOG electrodes 
    eog_chani = find(ismember({EEG.chanlocs.labels}, EOG_CHAN));
    [EEG.chanlocs(eog_chani).type] = deal('EOG');
    EEG = eeg_checkset( EEG );

    % Store channel locations in another field
    EEG.urchanlocs = EEG.chanlocs;

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

    EEG.setname = [subj '_ready_for_preprocessing'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Apply Filter
    %%EEG = pop_eegfiltnew(EEG, 'locutoff',LCF,'hicutoff',0);
    %%EEG = pop_eegfiltnew(EEG, 'hicutoff',HCF,'plotfreqz',0);
    EEG = pop_firws(EEG, 'fcutoff', LCF, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', (2 * ceil((3*(EEG.srate/LCF)) / 2)), 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
    EEG = pop_firws(EEG, 'fcutoff', HCF, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder',(2 * ceil((3*(EEG.srate/HCF)) / 2)), 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

    % Remove bad channels as identified in tid_psam_ica_preprocessing.m (see above)
    EEG.badchans = chans_to_interp;
    EEG = pop_select(EEG,'nochannel', EEG.badchans);

    EEG.setname = [subj '_ready_for_ICA_weights'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Attach ICA weight to main data
    EEG = pop_editset(EEG,'run', [], 'icaweights','ALLEEG(1).icaweights', 'icasphere','ALLEEG(1).icasphere');
    % Label ICA components with IC Label Plugin (Pion-Tonachini et al., 2019)
    %%EEG = pop_iclabel(EEG, 'default');
    %%EEG = pop_icflag(EEG, [0 0;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1]);
    % Sanity Check: Plot flagged ICs
    %%tid_psam_plot_flagged_ICs_TD(EEG,['ICs for ' subj], 'SavePath' ,fullfile(OUTPATH, [subj '_ic_topo.png']), 'PlotOn', false)
    % Remove previously flagged ICs (see tid_psam_ica_preprocessing.m)
    EEG = pop_subcomp( EEG, [], 0);

    % Interpolate bad channels
    if ~isempty(EEG.badchans)
        EEG = pop_interp(EEG, EEG.urchanlocs , 'spherical'); % and interpolate them using urchanlocs
    end

    % Sanity Check: Plot RMS in 10s bins for each electode
    tid_psam_plot_rms_bins_TD(EEG, [subj ' RMS bins'], 'SavePath', fullfile(OUTPATH, [subj '_channel_rms.png']), 'PlotOn', false)

    % Epoching
    EEG = pop_epoch( EEG, EVENTS, [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');

    % Baseline-Removal
    EEG = pop_rmbase( EEG, [BL_FROM 0] ,[]);
    
    % Probability-based removal
    EEG = pop_jointprob(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
    EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);

    % Save dataset
    EEG.setname = [subj '_preprocessed'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_saveset(EEG, 'filename',[subj '_erp_preprocessed.set'],'filepath', OUTPATH);

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
writetable(protocol,fullfile(OUTPATH, 'erp_preprocessing_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_erp_preprocessing_marked_subj.xlsx'))
end

check_done = 'tid_psam_erp_preprocessing_DONE'

close all; delete(wb)












