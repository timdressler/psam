% tid_psam_ica_preprocessing_parallel.m
%
% Performs ICA preprocessing for subtasks and saves datasets (in parallel).
% Adapted to match newer architectural standards.
%
% Preprocessing includes the following steps:
%   Load raw data
%   Apply a 1 Hz HP-Filter
%   Identify and remove bad channels
%   Create regular 1s epochs
%   Remove bad epochs based on probability and kurtosis
%   Calculate ICA weights
%   Label bad components using the ICLabel Plugin (Pion-Tonachini et al., 2019)
%
% Literature
% Pion-Tonachini, L., Kreutz-Delgado, K., & Makeig, S. (2019).
%   ICLabel: An automated electroencephalographic independent component classifier.
%   NeuroImage, 198, 181-197.
%
% Tim Dressler, [Current Date]
clear
close all
clc
set(0,'defaulttextInterpreter','none') 

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
INPATH = fullfile(MAINPATH, 'data');
OUTPATH = fullfile(MAINPATH, 'data', 'derivatives', 'preprocessing', 'ica_preprocessing');
TOPO_OUTPATH = fullfile(OUTPATH, 'sanity_plots_ic_topos');
REJECT_OUTPATH = fullfile(OUTPATH, 'sanity_plots_rejected_chans'); 

FUNPATH = fullfile(MAINPATH, 'functions');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH, TOPO_OUTPATH, REJECT_OUTPATH)

% Variables to edit
pipeline_params.RANDOM_SEED = 123;
pipeline_params.TASKS = {'delayedArticulation'}; 
pipeline_params.CHANNELS_TO_REMOVE = {'E29', 'E30'};
pipeline_params.EPOCH_LENGTH_SEC = 1;
pipeline_params.LCF_ICA = 1;
pipeline_params.SD_PROB_ICA = 3;
pipeline_params.THRESH_NUM_BAD_CHANS = 8;
pipeline_params.MAX_BADCHAN_LIMIT = 10;
pipeline_params.ITER_BADCHANS = 5; 
pipeline_params.CHANCORR_CRIT = 0.8; 
pipeline_params.MAX_BROKEN_TIME = 0.4; 
pipeline_params.BADCHAN_FRAC_THRESH = 0.3; 
pipeline_params.FLATLINE_CRIT_BADCHAN = 3;
pipeline_params.LINENOISE_CRIT_BADCHAN = 'off';
pipeline_params.ICLABEL_THRESHOLDS = [0 0; 0.7 1; 0.7 1; 0.7 1; 0.7 1; 0.7 1; 0.7 1];
pipeline_params.NUM_WORKERS = 1; % Set to 1 for laptop/serial mode, increase for HPC

% Check resources
if ismac
    [~, meminfo] = system('sysctl -n hw.memsize');
    total_ram_gb = str2double(strtrim(meminfo)) / 1024^3;
elseif isunix
    [~, meminfo] = system('grep MemTotal /proc/meminfo');
    total_ram_gb = sscanf(meminfo, 'MemTotal: %f') / (1024^2);
else
    [~, sys] = memory;
    total_ram_gb = sys.PhysicalMemory.Total / 1024^3;
end
fprintf('RAM detected: %.1f GB\n', total_ram_gb);
if total_ram_gb < 256 && pipeline_params.NUM_WORKERS > 4
    warning('Less than 256 GB RAM detected (%.1f GB). Running %d workers may be unsafe.', total_ram_gb, pipeline_params.NUM_WORKERS);
    user_input = input('Are you sure you want to continue? (y/n): ', 's');
    if ~strcmpi(user_input, 'y')
        error('Aborted by user. Change number of workers manually in the "Variables to edit" section!');
    end
end

% Setting for parallel processing
delete(gcp('nocreate')) % Close old one if one exists
parpool(pipeline_params.NUM_WORKERS)

tid_psam_unpack_params(pipeline_params);
tid_psam_save_params(pipeline_params, 'tid_psam_ica_preprocessing_parallel', OUTPATH);

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*'));
num_subj = length(dircont_subj);

% Load subject overview
participants = readtable(fullfile(INPATH, 'participants.tsv'), "FileType","text", "Delimiter","\t");

% Load existing logs if they exist
prot_file_path = fullfile(OUTPATH, 'tid_psam_ica_preprocessing_protocol.xlsx');
if exist(prot_file_path, 'file')
    old_protocol = readtable(prot_file_path, 'PreserveVariableNames', true);
else
    old_protocol = table();
end

marked_file_path = fullfile(OUTPATH, 'tid_psam_ica_preprocessing_marked_subj.xlsx');
if exist(marked_file_path, 'file')
    old_marked_subj = readtable(marked_file_path, 'PreserveVariableNames', true);
else
    old_marked_subj = table();
end

flagged_file_path = fullfile(OUTPATH, 'tid_psam_ica_preprocessing_flagged_comps.xlsx');
if exist(flagged_file_path, 'file')
    old_flagged_comps = readtable(flagged_file_path, 'PreserveVariableNames', true);
else
    old_flagged_comps = table();
end

removed_file_path = fullfile(OUTPATH, 'tid_psam_ica_processing_removed_channels.xlsx');
if exist(removed_file_path, 'file')
    old_removed_channels = readtable(removed_file_path, 'PreserveVariableNames', true);
else
    old_removed_channels = table();
end

% Start eeglab
eeglab nogui;

% Preallocate arrays for parallel loop
protocol_all = cell(num_subj, 1);
marked_subj_all = cell(num_subj, 1);
flagged_comps_log_all = cell(num_subj, 1);
removed_channels_log_all = cell(num_subj, 1);

disp('Starting tid_psam_ica_preprocessing_parallel.m...');

% Loop across subjects
parfor subj_idx = 1:num_subj
    set(0, 'DefaultFigureVisible', 'off');
    
    % Get current ID and subject path
    subj = extractAfter(dircont_subj(subj_idx).name, 'sub-');
    subj_path = fullfile(INPATH, ['sub-' subj]);
    subj_eeg_path = fullfile(subj_path, 'eeg');
    
    fprintf('Processing %s (%d / %d)\n', subj, subj_idx, num_subj);
    
    % Get list of to-be-processed .set files
    dircont_files = dir(fullfile(subj_eeg_path, '*.set'));
    dircont_files = dircont_files(contains({dircont_files.name}, TASKS));
    
    % Local variables for this specific worker/subject
    local_protocol = {};
    local_marked = {};
    local_flagged = {};
    local_removed = {};

    % Loop across files
    for file_idx = 1:length(dircont_files)
        file_name = dircont_files(file_idx).name;
        file_path = dircont_files(file_idx).folder;
        task = extractBefore(extractAfter(file_name,'task-'),'_eeg');
        task_has_issue = false;

        subj_file = fullfile(OUTPATH, ['sub-' subj '_task-' task '_ica_preprocessing.set']);
        if exist(subj_file, 'file')
            fprintf('Skipping sub-%s task-%s (already run)\n', subj, task);
            
            % Retrieve original protocol
            old_time = NaN;
            if ~isempty(old_protocol) && ismember('subj', old_protocol.Properties.VariableNames)
                old_idx = find(strcmp(string(old_protocol.subj), string(subj)) & strcmp(string(old_protocol.task), string(task)));
                if ~isempty(old_idx)
                    old_time = old_protocol.time(old_idx(1));
                end
            end
            local_protocol(end+1, 1:4) = {subj, task, old_time, 'SKIPPED '};
            
            % Retrieve original marked subjects
            if ~isempty(old_marked_subj) && ismember('subj', old_marked_subj.Properties.VariableNames)
                old_idx_m = find(strcmp(string(old_marked_subj.subj), string(subj)) & strcmp(string(old_marked_subj.task), string(task)));
                if ~isempty(old_idx_m)
                    for m_idx = 1:length(old_idx_m)
                        local_marked(end+1, 1:3) = {subj, task, char(string(old_marked_subj.issue(old_idx_m(m_idx))))};
                    end
                end
            end
            
            % Retrieve original flagged_comps count and list
            old_n_flagged = NaN;
            old_flagged_str = 'UNKNOWN';
            if ~isempty(old_flagged_comps) && ismember('subj', old_flagged_comps.Properties.VariableNames)
                old_idx_f = find(strcmp(string(old_flagged_comps.subj), string(subj)) & strcmp(string(old_flagged_comps.task), string(task)));
                if ~isempty(old_idx_f)
                    if ismember('n_flagged_comps', old_flagged_comps.Properties.VariableNames)
                        old_n_flagged = old_flagged_comps.n_flagged_comps(old_idx_f(1));
                    end
                    if ismember('flagged_comps', old_flagged_comps.Properties.VariableNames)
                        old_flagged_str = char(string(old_flagged_comps.flagged_comps(old_idx_f(1))));
                    end
                end
            end
            local_flagged(end+1, 1:4) = {subj, task, old_n_flagged, old_flagged_str};
            
            % Retrieve original removed channels
            old_n_removed = NaN;
            old_removed_str = 'UNKNOWN';
            if ~isempty(old_removed_channels) && ismember('subj', old_removed_channels.Properties.VariableNames)
                old_idx_r = find(strcmp(string(old_removed_channels.subj), string(subj)) & strcmp(string(old_removed_channels.task), string(task)));
                if ~isempty(old_idx_r)
                    old_n_removed = old_removed_channels.n_removed_channel(old_idx_r(1));
                    old_removed_str = char(string(old_removed_channels.removed_channels(old_idx_r(1))));
                end
            end
            local_removed(end+1, 1:4) = {subj, task, old_n_removed, old_removed_str};
            continue;
        end

        file_tic = tic;

        % Load data
        EEG = pop_loadset('filename',file_name,'filepath',file_path);

        % Add channel locations
        chanlocs_file = ['sub-' subj '_electrodes.tsv'];
        EEG = pop_chanedit(EEG, 'load',{fullfile(file_path, chanlocs_file),'filetype','tsv'});

        % Add channel type
        chantype_file = ['sub-' subj '_task-' task '_channels.tsv'];
        chantype = readtable(fullfile(file_path, chantype_file), 'FileType','text', 'Delimiter','\t');

        for i = 1:length(EEG.chanlocs)
            current_label = EEG.chanlocs(i).labels;
            row_idx = find(strcmp(chantype.name, current_label));
            if ~isempty(row_idx)
                EEG.chanlocs(i).type = chantype.type{row_idx};
            else
                warning('Channel %s not found. Labeling as UNKNOWN.', current_label);
                EEG.chanlocs(i).type = 'UNKNOWN';
            end
        end

        % Remove Marker-Channel
        EEG = pop_select(EEG, 'rmchannel', CHANNELS_TO_REMOVE);

        % Highpass-Filter
        LCF_ord = pop_firwsord('hamming', EEG.srate, tid_psam_get_transition_bandwidth_TD(LCF_ICA));
        EEG = pop_firws(EEG, 'fcutoff', LCF_ICA, 'ftype', 'highpass', 'wtype', 'hamming', 'forder',LCF_ord, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

        % Find bad channels 
        bemobil_params = struct();
        bemobil_params.chancorr_crit = CHANCORR_CRIT;
        bemobil_params.chan_max_broken_time = MAX_BROKEN_TIME;
        bemobil_params.chan_detect_num_iter = ITER_BADCHANS;
        bemobil_params.chan_detected_fraction_threshold = BADCHAN_FRAC_THRESH;
        bemobil_params.num_chan_rej_max_target = MAX_BADCHAN_LIMIT;
        bemobil_params.flatline_crit = FLATLINE_CRIT_BADCHAN;
        bemobil_params.line_noise_crit = LINENOISE_CRIT_BADCHAN;
        
        [chans_to_interp, rejected_chan_plot_handle] = tid_psam_detect_bad_channels(EEG, ...
            bemobil_params, 'PlotOn', true);
            
        % Save the sanity check plot
        if ~isempty(rejected_chan_plot_handle) && isgraphics(rejected_chan_plot_handle)
            out_file_rej = fullfile(REJECT_OUTPATH, ['sub-' subj '_task-' task '_rejected_chans.png']);           
            exportgraphics(rejected_chan_plot_handle, out_file_rej, 'Resolution', 150);            
            close(rejected_chan_plot_handle);
        end
        
        % Log removed channels
        if ~isempty(chans_to_interp)
            removed_labels = {EEG.chanlocs(chans_to_interp).labels};
            removed_str = strjoin(removed_labels, ', ');
        else
            removed_str = 'none';
        end
        local_removed(end+1, 1:4) = {subj, task, length(chans_to_interp), removed_str};

        % Actually remove the channels
        EEG.badchans = chans_to_interp;
        EEG = pop_select(EEG, 'nochannel', EEG.badchans);

        % Mark subjects if too many bad channels
        if length(chans_to_interp) > THRESH_NUM_BAD_CHANS
            local_marked(end+1, 1:3) = {subj, task, ['large_num_bad_chan_' num2str(length(chans_to_interp))]};
            task_has_issue = true;
        end

        % Epoch, Remove Bad Epochs
        EEG = eeg_regepochs(EEG);
        EEG = pop_jointprob(EEG,1,[1:EEG.nbchan] ,SD_PROB_ICA,0,0,0,[],0);
        EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,SD_PROB_ICA,0,0,0,[],0);
        EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch(EEG, EEG.reject.rejglobal ,0);

        % Run ICA
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
        
        % Label ICA components with IC Label Plugin (Pion-Tonachini et al., 2019)
        EEG = pop_iclabel(EEG, 'default');
        EEG = pop_icflag(EEG, ICLABEL_THRESHOLDS);
        
        flagged_idx = find(EEG.reject.gcompreject);
        n_flagged_comps = length(flagged_idx);
        if ~isempty(flagged_idx)
            flagged_str = strjoin(string(flagged_idx), ', ');
        else
            flagged_str = 'none';
        end
        local_flagged(end+1, 1:4) = {subj, task, n_flagged_comps, flagged_str};

        % Sanity Check: Plot flagged ICs
        tid_psam_plot_flagged_ICs_TD(EEG,['sub-' subj '_task-' task '_ic_topos'], 'SavePath' ,fullfile(TOPO_OUTPATH, ['sub-' subj '_task-' task '_ic_topos.png']), 'PlotOn', false)

        % Save dataset
        EEG.setname = ['sub-' subj '_task-' task '_ica_weights'];
        EEG = pop_saveset(EEG, 'filename',['sub-' subj '_task-' task '_ica_weights.set'],'filepath', OUTPATH);

        % Update Protocol
        file_time = toc(file_tic);
        status_str = 'OK';
        if task_has_issue; status_str = 'MARKED'; end
        local_protocol(end+1, 1:4) = {subj, task, file_time, status_str};

    end % file_idx

    % Store the local worker arrays into the preallocated sliced arrays
    protocol_all{subj_idx} = local_protocol;
    marked_subj_all{subj_idx} = local_marked;
    flagged_comps_log_all{subj_idx} = local_flagged;
    removed_channels_log_all{subj_idx} = local_removed;

end % subj_idx

% Flatten pre-allocated arrays into single cell matrices
protocol = vertcat(protocol_all{:});
marked_subj = vertcat(marked_subj_all{:});
flagged_comps_log = vertcat(flagged_comps_log_all{:});
removed_channels_log = vertcat(removed_channels_log_all{:});

% End of processing
if ~isempty(protocol)
    protocol_tab = cell2table(protocol, 'VariableNames',{'subj', 'task', 'time', 'status'});
    writetable(protocol_tab,fullfile(OUTPATH, 'tid_psam_ica_preprocessing_protocol.xlsx'))
end

if ~isempty(marked_subj)
    marked_subj_tab = cell2table(marked_subj, 'VariableNames',{'subj', 'task', 'issue'});
    writetable(marked_subj_tab,fullfile(OUTPATH, 'tid_psam_ica_preprocessing_marked_subj.xlsx'))
end

if ~isempty(flagged_comps_log)
    flagged_comps_tab = cell2table(flagged_comps_log, 'VariableNames', {'subj', 'task', 'n_flagged_comps', 'flagged_comps'});
    writetable(flagged_comps_tab, fullfile(OUTPATH, 'tid_psam_ica_preprocessing_flagged_comps.xlsx'));
end

if ~isempty(removed_channels_log)
    removed_channels_tab = cell2table(removed_channels_log, 'VariableNames', {'subj', 'task', 'n_removed_channel', 'removed_channels'});
    writetable(removed_channels_tab, fullfile(OUTPATH, 'tid_psam_ica_processing_removed_channels.xlsx'));
end

set(0, 'DefaultFigureVisible', 'on');
delete(gcp('nocreate'))
check_done = 'tid_psam_ica_preprocessing_DONE'