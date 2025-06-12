% tid_psam_hilbert_preparation.m
%
% Performs a hilbert transformation later used for the SVM analysis and saves data.
%
% Preparation includes the following steps
%
    % Removes EOG electrodes
    % Load data and data set containing ICA weights (see
    %   tid_psam_ica_preprocessing.m)
    % Rename events
    % Remove bad channels as identified in tid_psam_ica_preprocessing.m
    % Attach ICA weights and remove bad components using the ICLabel Plugin
    %   (Pion-Tonachini et al., 2019)
    % Interpolate bad (and removed) channels
    % Loop across frequency bands and apply Hilbert transformation
    % Epoch and remove bad epoch based on the ones identified for the SVM analysis in tid_psam_exclude_trials.m
    %
% Saves data
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
INPATH_RAW = fullfile(MAINPATH, 'data\processed_data\markers_included\');
INPATH_ICA = fullfile(MAINPATH, 'data\processed_data\ica_preprocessed\');
INPATH_EXCLUDED = fullfile(MAINPATH, 'data\processed_data\exclude_trials\');
OUTPATH = fullfile(MAINPATH, 'data\processed_data\hilbert_prepared_clean');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, IDPATH, OUTPATH, INPATH_RAW, INPATH_ICA)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EOG_CHAN = {'E29','E30'}; % Labels of EOG electrodes
EPO_FROM = -1;
EPO_TILL = 0.1;
EVENTS = {'go_act', 'go_pas'};
APLHA_FROM = 8;
APLHA_TILL = 13;
BETA_FROM = 14;
BETA_TILL = 30;
GAMMA_FROM = 31;
GAMMA_TILL = 37;

% Get directory content
dircont_subj = dir(fullfile(IDPATH, 'sub-*.set'));
dircont_subj_exclude_trials = dir(fullfile(INPATH_EXCLUDED, 'sub-*.mat'));

% Sanity Check: Same length of directory contents
if length(dircont_subj) == length(dircont_subj_exclude_trials) 
else
    error('Different number of files')
end

% Initialize sanity check variables
marked_subj = {};
protocol = {};

% Aggregate frequency bands
freq_bands = {[APLHA_FROM APLHA_TILL], [BETA_FROM BETA_TILL], [GAMMA_FROM GAMMA_TILL]};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_hilbert_preparation.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Sanity Check: Same IDs
    if ~strncmp(dircont_subj(subj_idx).name, dircont_subj_exclude_trials(subj_idx).name, 6)
        error('Files from two different subjects')
    end

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Get file of to-be excluded trials based on tid_psam_exclude_trials 
    exclusion_filename = fullfile(INPATH_EXCLUDED,[subj '_excluded_no_probe_trials.mat']);
    load(exclusion_filename) % loads variable 'excluded_trials_all_no_probe'

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_hilbert_preparation.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load ICA data
    EEG = pop_loadset('filename',[subj '_ica_weights.set'],'filepath',INPATH_ICA);

    % Get flagged components as identified in tid_psam_ica_preprocessing.m
    flagged_comps = find(EEG.reject.gcompreject);

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
        if strcmp(EEG.event(event).type, 'S  5') && strncmp(EEG.event(event-1).type, 'con', 3) % Get all go-signals (S  5) during no-probe trials
            if strcmp(EEG.event(event-2).type, 'S 21')
                EEG.event(event).type = 'go_act';
            elseif strcmp(EEG.event(event-2).type, 'S 22')
                EEG.event(event).type = 'go_pas';
            else
                error('Unkown Marker!')
            end
        end
    end

    EEG.setname = [subj '_ready_for_preprocessing'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

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
    EEG.reject.gcompreject = flagged_comps;
    EEG = pop_subcomp( EEG, [], 0);

    % Interpolate bad channels
    if ~isempty(EEG.badchans)
        EEG = pop_interp(EEG, EEG.urchanlocs , 'spherical'); % and interpolate them using urchanlocs
    end

    % Remove EOG channels as they are not used for classification
    EEG = pop_select( EEG, 'rmchannel',EOG_CHAN);
    
    EEG.setname = [subj '_ready_for_Hilbert_transformation'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Loop across frequency bands and apply Hilbert transformation
    for freq_band_num = 1:length(freq_bands)
        EEG = ALLEEG(4);
        CURRENTSET = 4;

        % Apply Filter
        %%EEG = pop_eegfiltnew(EEG, 'locutoff',freq_bands{freq_band_num}(1), 'plotfreqz', 1);
        %%EEG = pop_eegfiltnew(EEG, 'hicutoff',freq_bands{freq_band_num}(2));
        EEG = pop_firws(EEG, 'fcutoff', freq_bands{freq_band_num}(1), 'ftype', 'highpass', 'wtype', 'hamming', 'forder', (2 * ceil((3*(EEG.srate/freq_bands{freq_band_num}(1))) / 2)), 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
        EEG = pop_firws(EEG, 'fcutoff', freq_bands{freq_band_num}(2), 'ftype', 'lowpass', 'wtype', 'hamming', 'forder',(2 * ceil((3*(EEG.srate/freq_bands{freq_band_num}(2))) / 2)), 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

        % Apply Hilbert transformation
        EEG.data = abs(hilbert(EEG.data));

        % Epoch
        EEG = pop_epoch( EEG, EVENTS, [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');
    
        % Reject bad epochs based on the ones identified for the SVM analysis in tid_psam_exclude_trials.m
        EEG.reject.rejglobal = excluded_trials_all_no_probe;
        EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);

        % Store dataset
        EEG.setname = [subj '_after_Hilbert_transformation_clean_epochs_' num2str(freq_bands{freq_band_num}(1)) '_' num2str(freq_bands{freq_band_num}(2))];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        % Store data
        data_hilbert(freq_band_num).freqband = freq_bands{freq_band_num};
        data_hilbert(freq_band_num).hilbert_transformed_epochs_clean = EEG.data;

        % Save data
        save(fullfile(OUTPATH, [subj '_hilbert_preparation_clean.mat']),'data_hilbert', 'freq_bands')
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

check_done = 'tid_psam_hilbert_preparation_DONE'

delete(wb)