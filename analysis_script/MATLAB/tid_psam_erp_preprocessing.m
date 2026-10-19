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
%   Load data and data set containing ICA weights (see
%       tid_psam_ica_preprocessing.m)
%   Rename events
%   Apply a 30 Hz LP-Filter
%   Remove bad channels as identified in tid_psam_ica_preprocessing.m
%   Attach ICA weights and remove bad components previously indentified using the ICLabel Plugin
%       (Pion-Tonachini et al., 2019)
%   Interpolate bad (and removed) channels
%   Epoch data around probe events (see above)
%   Apply baseline correction
%   Tag bad epochs based on probability
%
% Saves data
%
% Note. Trials are only marked for exclusion but not excluded yet!
%
% Literature
% Pion-Tonachini, L., Kreutz-Delgado, K., & Makeig, S. (2019).
%   ICLabel: An automated electroencephalographic independent classifier, dataset, and website.
%   NeuroImage, 198, 181-197.
%
% Tim Dressler, 04.04.2025
clear
close all
clc
set(0,'DefaultTextInterpreter','none')

% Set up paths
SCRIPTPATH = cd;
normalizedPath = strrep(SCRIPTPATH, filesep, '/');
expectedSubpath = 'psam/analysis_script/MATLAB';
if contains(normalizedPath, expectedSubpath)
    disp('Path OK')
else
    error('Path not OK')
end
MAINPATH = strrep(SCRIPTPATH, fullfile('analysis_script', 'MATLAB'), '');
INPATH_RAW = fullfile(MAINPATH, 'data');
INPATH_ICA = fullfile(MAINPATH, 'data', 'derivatives', 'preprocessing', 'ica_preprocessing');
OUTPATH = fullfile(MAINPATH, 'data', 'derivatives', 'preprocessing', 'erp_preprocessing');
FUNPATH = fullfile(MAINPATH, 'functions');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_RAW, INPATH_ICA, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
pipeline_params.RANDOM_SEED = 123;
pipeline_params.EOG_CHAN = {'E29','E30'};
pipeline_params.EPO_FROM = -0.2;
pipeline_params.EPO_TILL = 0.400;
pipeline_params.LCF = 0.3;
pipeline_params.HCF = 30;
pipeline_params.BL_FROM = -200;
pipeline_params.SD_PROB = 3;
pipeline_params.EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};

% Unpack and save parameters
rng(pipeline_params.RANDOM_SEED)
tid_psam_unpack_params(pipeline_params);
tid_psam_save_params(pipeline_params, 'tid_psam_erp_preprocessing', OUTPATH);

% Get directory content (BIDS format)
dircont_subj = dir(fullfile(INPATH_RAW, 'sub-*'));
dircont_subj = dircont_subj([dircont_subj.isdir]); % Keep only directories

% Sanity Check: Ensure subjects exist
if isempty(dircont_subj)
    error('No BIDS subject folders found in INPATH_RAW')
end

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_erp_preprocessing.m');

clear subj_idx
for subj_idx = 1:length(dircont_subj)
    % Get current ID and paths
    subj = extractAfter(dircont_subj(subj_idx).name, 'sub-');
    subj_path = fullfile(INPATH_RAW, ['sub-' subj]);
    subj_eeg_path = fullfile(subj_path, 'eeg');
    
    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, ['sub-' subj ' tid_psam_erp_preprocessing.m'])
    
    % Get list of task files for the subject
    dircont_files = dir(fullfile(subj_eeg_path, '*.set'));
    
    for file_idx = 1:length(dircont_files)
        file_name = dircont_files(file_idx).name;
        file_path = dircont_files(file_idx).folder;
        task = extractBefore(extractAfter(file_name,'task-'),'_eeg');
        
        tic;
        % Start eeglab
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        % Load ICA data using dynamic naming matching the reference script
        ica_file_name = ['sub-' subj '_task-' task '_ica_preprocessing.set'];
        EEG = pop_loadset('filename', ica_file_name, 'filepath', INPATH_ICA);
        
        % Get flagged components as identified in tid_psam_ica_preprocessing.m
        flagged_comps = EEG.reject.gcompreject;
        EEG.setname = ['sub-' subj '_task-' task '_ICA_weights'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        
        % Get bad channels as identified in tid_psam_ica_preprocessing.m
        chans_to_interp = EEG.badchans;
        
        % Load raw data
        EEG = pop_loadset('filename', file_name, 'filepath', file_path);
           
        % Add channel locations from BIDS TSV
        chanlocs_file = ['sub-' subj '_electrodes.tsv'];
        EEG = pop_chanedit(EEG, 'load',{fullfile(file_path, chanlocs_file),'filetype','tsv'});
        EEG = eeg_checkset( EEG );
        
        % Add channel type from BIDS TSV
        chantype_file = ['sub-' subj '_task-' task '_channels.tsv'];
        chantype = readtable(fullfile(file_path, chantype_file), 'FileType','text', 'Delimiter','\t');
        
        for i = 1:length(EEG.chanlocs)
            current_label = EEG.chanlocs(i).labels;
            row_idx = find(strcmp(chantype.name, current_label));
            if ~isempty(row_idx)
                EEG.chanlocs(i).type = chantype.type{row_idx};
            else
                EEG.chanlocs(i).type = 'UNKNOWN';
            end
        end
        EEG = eeg_checkset( EEG );
        
        % Remove EOG-Channels
        EEG = pop_select(EEG, 'rmchannel', EOG_CHAN);
        
        % Store channel locations in another field
        EEG.urchanlocs = EEG.chanlocs;
        
        EEG.setname = ['sub-' subj '_task-' task '_ready_for_preprocessing'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        
        % Apply Filter
        % Highpass-Filter % DEV1
        %%LCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth_TD(LCF)); % Get filter order (also see pop_firwsord.m, tid_psam_get_transition_bandwidth_TD.m)
        %%EEG = pop_firws(EEG, 'fcutoff', LCF, 'ftype', 'highpass', 'wtype', 'hamming', 'forder',LCF_ord, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
        
        % Lowpass-Filter
        HCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth_TD(HCF)); % Get filter order (also see pop_firwsord.m, tid_psam_get_transition_bandwidth_TD.m)
        EEG = pop_firws(EEG, 'fcutoff', HCF, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder',HCF_ord, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
        
        % Remove bad channels as identified in tid_psam_ica_preprocessing.m (see above)
        EEG.badchans = chans_to_interp;
        EEG = pop_select(EEG,'nochannel', EEG.badchans);
        EEG.setname = ['sub-' subj '_task-' task '_ready_for_ICA_weights'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        
        % Attach ICA weight to main data
        EEG = pop_editset(EEG,'run', [], 'icaweights','ALLEEG(1).icaweights', 'icasphere','ALLEEG(1).icasphere');
        
        % Remove previously flagged ICs (see tid_psam_ica_preprocessing.m)
        EEG.reject.gcompreject = flagged_comps;
        EEG = pop_subcomp( EEG, [], 0);
        
        % Interpolate bad channels
        if ~isempty(EEG.badchans)
            EEG = pop_interp(EEG, EEG.urchanlocs , 'spherical');
        end
        
        % Sanity Check: Plot RMS in 10s bins for each electode
        tid_psam_plot_rms_bins_TD(EEG, ['sub-' subj '_task-' task ' RMS bins'], 'SavePath', fullfile(OUTPATH, ['sub-' subj '_task-' task '_channel_rms.png']), 'PlotOn', false)
        
        % Epoching
        EEG = pop_epoch( EEG, EVENTS, [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');
        
        % Baseline-Removal
        EEG = pop_rmbase( EEG, [BL_FROM 0] ,[]);
        
        % Probability-based removal
        EEG = pop_jointprob(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
        EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        
        % Save dataset
        EEG.setname = ['sub-' subj '_task-' task '_erp_preprocessed'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        EEG = pop_saveset(EEG, 'filename',['sub-' subj '_task-' task '_erp_preprocessing.set'],'filepath', OUTPATH);
        
        % Update Protocol dynamically for BIDS layout
        file_time = toc;
        if any(strcmp(marked_subj, subj), 'all')
            status_str = 'MARKED';
        else
            status_str = 'OK';
        end
        protocol(end+1, 1:4) = {subj, task, file_time, status_str};
    end
end
% End of processing
if ~isempty(protocol)
    protocol = cell2table(protocol, 'VariableNames',{'subj', 'task', 'time', 'status'});
    writetable(protocol,fullfile(OUTPATH, 'tid_psam_erp_preprocessing_protocol.xlsx'))
end
if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj', 'task', 'issue'});
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_erp_preprocessing_marked_subj.xlsx'))
end
check_done = 'tid_psam_erp_preprocessing_DONE'
close all; delete(wb)