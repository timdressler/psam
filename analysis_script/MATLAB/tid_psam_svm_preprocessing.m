% tid_psam_SVM_preprocessing.m
%
% Performs preprocessing for SVM analysis and saves datasets.
%
% Preprocessing includes the following steps
%
    % Load data and data set containing ICA weights (see
    %   tid_psam_ica_preprocessing.m)
    % Rename events
    % Apply a 1 Hz HP-Filter
    % Apply a 30 Hz LP-Filter
    % Remove bad channels as identified in tid_psam_ica_preprocessing.m
    % Attach ICA weights and remove bad components using the ICLabel Plugin
    %   (Pion-Tonachini et al., 2019)
    % Interpolate bad (and removed) channels
    % Epoch data around 'go-signal'
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
% Tim Dressler, 01.05.2025

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
OUTPATH = fullfile(MAINPATH, 'data\processed_data\svm_preprocessed');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_RAW, INPATH_ICA, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EOG_CHAN = {'E29','E30'}; % Labels of EOG electrodes
EPO_FROM = -1;
EPO_TILL = 0.1;
LCF = 1;
HCF = 40;
BL_FROM = -1000;
BL_TILL = -800;
SD_PROB = 3;
SD_PROB_ICA = 3;
EVENTS = {'go_act', 'go_pas'};

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
wb = waitbar(0,'starting tid_psam_svm_preprocessing.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_svm_preprocessing.m'])

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

    % Filter
    EEG = pop_eegfiltnew(EEG, 'locutoff',LCF,'hicutoff',0);
    EEG = pop_eegfiltnew(EEG, 'hicutoff',HCF,'plotfreqz',0);

    % Remove bad channels as identified in tid_psam_ica_preprocessing.m (see above)
    EEG.badchans = chans_to_interp;
    EEG = pop_select(EEG,'nochannel', EEG.badchans);

    EEG.setname = [subj '_ready_for_ICA_weights'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Attach ICA weight to main data
    EEG = pop_editset(EEG,'run', [], 'icaweights','ALLEEG(1).icaweights', 'icasphere','ALLEEG(1).icasphere');
    % Label ICA components with IC Label Plugin (Pion-Tonachini et al., 2019)
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, [0 0;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1]);
    % Sanity Check: Plot flagged ICs
    tid_psam_plot_flagged_ICs_TD(EEG,['ICs for ' subj], 'SavePath' ,fullfile(OUTPATH, [subj '_ic_topo.png']), 'PlotOn', false)
    % Remove flagged ICs
    EEG = pop_subcomp( EEG, [], 0);

    % Interpolate bad channels
    if ~isempty(EEG.badchans)
        EEG = pop_interp(EEG, EEG.urchanlocs , 'spherical'); % and interpolate them using urchanlocs
    end

    % Epoching
    EEG = pop_epoch( EEG, EVENTS, [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');

    % Baseline-Removal
    EEG = pop_rmbase( EEG, [BL_FROM BL_TILL] ,[]);

    % Probability-based removal
    EEG = pop_jointprob(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
    EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);

    % Save dataset
    EEG.setname = [subj '_preprocessed'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_saveset(EEG, 'filename',[subj '_svm_preprocessed.set'],'filepath', OUTPATH);

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
writetable(protocol,fullfile(OUTPATH, 'svm_preprocessing_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_svm_preprocessing_marked_subj.xlsx'))
end

check_done = 'tid_psam_svm_preprocessing_DONE'

close all; delete(wb)












