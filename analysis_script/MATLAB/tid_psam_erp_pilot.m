% tid_psam_erp_pilot.m
%
% Performs a 'quick and dirty' ERP analysis (incl. preprocessing).
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
% Epochs data
% Performs baseline correction
% Reject bad epochs using threshold and probability 
%
% Analysis includes the following steps
%
% Load the preprocessed dataset
% Extractes ERP over all conditions
% Extracts ERP for each condition
% Stores ERP amplitudes and latencies
% Creates multiple plots
%
% Tim Dressler, 03.04.2025

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
INPATH_RAW = fullfile(MAINPATH, 'data', 'processed_data', 'markers_included');
INPATH_PROC_CLEAN = fullfile(MAINPATH, 'data', 'processed_data', 'erp_preprocessed_clean');
OUTPATH = fullfile(MAINPATH, 'data', 'analysis_data', 'erp_pilot');
FUNPATH = fullfile(MAINPATH, 'functions');

addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH_RAW, INPATH_PROC_CLEAN,OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
USE_RAW_DATA = 1;

EPO_FROM = -0.2;
EPO_TILL = 0.400;
LCF = 1;
HCF = 40;
BL_FROM = -200;
THRESH = 75;
SD_PROB = 3;
RESAM_ICA = 250;
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};

CHANI = 1;
ERP_FROM = 75;
ERP_TILL = 125;
CB_LIM_LOWER = -4;
CB_LIM_UPPER = 2;

% Get directory content
subj = 96;
dircont_subj = dir(fullfile(INPATH_RAW, ['sub-' num2str(subj) '*']));

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_erp_pilot.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)
    % Start preprocessing

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_erp_preprocessing.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    if USE_RAW_DATA == 0
        % When using preprocessed data

        % Load data
        EEG = pop_loadset('filename',[subj '_erp_preprocessed_clean.set'],'filepath',INPATH_PROC_CLEAN);

        EEG.setname = [subj '_preprocessed_clean'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    elseif USE_RAW_DATA == 1
        % When using rawdata

        % Load data
        EEG = pop_loadset('filename',[subj '_markers_inlcuded.set'],'filepath',INPATH_RAW);

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

        % Add channel locations
        EEG.chanlocs = readlocs( fullfile(MAINPATH,'\config\elec_96ch_adapted.elp')); 

        % Filter
        EEG = pop_firws(EEG, 'fcutoff', LCF, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', (2 * ceil((3*(EEG.srate/LCF)) / 2)), 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
        EEG = pop_firws(EEG, 'fcutoff', HCF, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder',(2 * ceil((3*(EEG.srate/HCF)) / 2)), 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

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

        EEG.setname = [subj '_qad_preprocessed'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        % End of preprocessing
    end

    %% Sart analysis

    % Get all 'Active' condition
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'act_early_alt','act_early_unalt','act_late_alt','act_late_unalt'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [subj '_act_all'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Get ERP
    erp_act_all = mean(EEG.data(:,:,:),3);
    % Setup ERP analysis
    [~,win_start] = min(abs(EEG.times-ERP_FROM));
    [~,win_end] = min(abs(EEG.times-ERP_TILL));
    [~,t_zero] = min(abs(EEG.times));
    % Get ERP amplitude
    maxERPamp = min(erp_act_all(CHANI,win_start:win_end));
    % Get ERP latency
    maxERPsam = find(erp_act_all(CHANI,:) == maxERPamp);
    maxERPlat = EEG.times(maxERPsam);
    % Store ERP
    all_ERP_act_all(:,:, subj_idx) = erp_act_all;

    clear maxERPamp maxERPsam maxERPlat

    % Get all 'Passive' condition
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'pas_early_alt','pas_early_unalt','pas_late_alt','pas_late_unalt'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [subj '_pas_all'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Get ERP
    erp_pas_all = mean(EEG.data(:,:,:),3);
    % Setup ERP analysis
    [~,win_start] = min(abs(EEG.times-ERP_FROM));
    [~,win_end] = min(abs(EEG.times-ERP_TILL));
    [~,t_zero] = min(abs(EEG.times));
    % Get ERP amplitude
    maxERPamp = min(erp_pas_all(CHANI,win_start:win_end));
    % Get ERP latency
    maxERPsam = find(erp_pas_all(CHANI,:) == maxERPamp);
    maxERPlat = EEG.times(maxERPsam);
    % Store ERP
    all_ERP_pas_all(:,:, subj_idx) = erp_pas_all;

    clear maxERPamp maxERPsam maxERPlat

    % Get all 'Active' condition (Control)
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'con_act_early', 'con_act_late'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [subj '_con_act_all'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Get ERP
    erp_con_act_all = mean(EEG.data,3);
    % Store ERP
    all_ERP_con_act_all(:,:, subj_idx) = erp_con_act_all;

    clear maxERPamp maxERPsam maxERPlat

    % Get all 'Passive' condition (Control)
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'con_pas_early', 'con_pas_late'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [subj '_con_pas_all'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Get ERP
    erp_con_pas_all = mean(EEG.data,3);
    % Store ERP
    all_ERP_con_pas_all(:,:, subj_idx) = erp_con_pas_all;

    clear maxERPamp maxERPsam maxERPlat

    % Get all conditions individually
    for cond = 1:length(EVENTS) % Loop over conditions
        EEG = pop_selectevent(ALLEEG(1), 'latency', '-2<=2', 'type', EVENTS{cond}, ...
            'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        EEG.setname = [subj '_' EVENTS{cond}];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        % Get ERP
        erp = mean(EEG.data, 3);

        % Setup ERP analysis
        [~, win_start] = min(abs(EEG.times - ERP_FROM));
        [~, win_end] = min(abs(EEG.times - ERP_TILL));
        [~, t_zero] = min(abs(EEG.times));

        % Get N100 amplitude
        maxERPamp = min(erp(CHANI, win_start:win_end));

        % Get N100 latency
        maxERPsam = find(erp(CHANI, :) == maxERPamp);
        maxERPlat = EEG.times(maxERPsam);

        % Store ERPs dynamically
        all_ERP_cond.(EVENTS{cond})(:, :, subj_idx) = erp;
    end

    % Get control ERPs
    % 'Active - Early' condition
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'con_act_early'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [subj '_con_act_early'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Get ERP
    erp_con_act_early = mean(EEG.data(CHANI,:,:),3);
    % Store ERP
    all_ERP_con_act_early(:,:, subj_idx) = erp_con_act_early;

    % 'Active - Late' condition
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'con_act_late'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [subj '_con_act_late'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Get ERP
    erp_con_act_late = mean(EEG.data(CHANI,:,:),3);
    % Store ERP
    all_ERP_con_act_late(:,:, subj_idx) = erp_con_act_late;

    % 'Passive - Early' condition
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'con_pas_early'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [subj '_con_pas_early'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Get ERP
    erp_con_pas_early = mean(EEG.data(CHANI,:,:),3);
    % Store ERP
    all_ERP_con_pas_early(:,:, subj_idx) = erp_con_pas_early;

    % 'Passive - Late' condition
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'con_pas_late'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [subj '_con_pas_late'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Get ERP
    erp_con_pas_late = mean(EEG.data(CHANI,:,:),3);
    % Store ERP
    all_ERP_con_pas_late(:,:, subj_idx) = erp_con_pas_late;

    % Collapsed 'Early' and 'Late' conditions
    erp_con_act_all = mean([erp_con_act_early; erp_con_act_late],1);
    erp_con_pas_all = mean([erp_con_pas_early; erp_con_pas_late],1);

    % Get corrected ERPs
    % Get all 'Active' condition
    erp_act_all_corrected = erp_act_all - erp_con_act_all;
    all_ERP_act_all_corrected(:,:, subj_idx) = erp_act_all_corrected;

    % Get all 'Passive' condition
    erp_pas_all_corrected = erp_pas_all - erp_con_pas_all;
    all_ERP_pas_all_corrected(:,:, subj_idx) = erp_pas_all_corrected;

    % Get all conditions individually
    % 'Active - Early - Unaltered' condition
    erp_act_early_unalt_corrected = all_ERP_cond.act_early_unalt(:,:,:) - erp_con_act_early;
    all_ERP_act_early_unalt_corrected(:,:, subj_idx) = erp_act_early_unalt_corrected;

    % 'Active - Late - Unaltered' condition
    erp_act_late_unalt_corrected = all_ERP_cond.act_late_unalt(:,:,:) - erp_con_act_late;
    all_ERP_act_late_unalt_corrected(:,:, subj_idx) = erp_act_late_unalt_corrected;

    % 'Active - Early - Altered' condition
    erp_act_early_alt_corrected = all_ERP_cond.act_early_alt(:,:,:) - erp_con_act_early;
    all_ERP_act_early_alt_corrected(:,:, subj_idx) = erp_act_early_alt_corrected;

    % 'Active - Late - Altered' condition
    erp_act_late_alt_corrected = all_ERP_cond.act_late_alt(:,:,:) - erp_con_act_late;
    all_ERP_act_late_alt_corrected(:,:, subj_idx) = erp_act_late_alt_corrected;

    % 'Passive - Early - Unaltered' condition
    erp_pas_early_unalt_corrected = all_ERP_cond.pas_early_unalt(:,:,:) - erp_con_pas_early;
    all_ERP_pas_early_unalt_corrected(:,:, subj_idx) = erp_pas_early_unalt_corrected;

    % 'Passive - Late - Unaltered' condition
    erp_pas_late_unalt_corrected = all_ERP_cond.pas_late_unalt(:,:,:) - erp_con_pas_late;
    all_ERP_pas_late_unalt_corrected(:,:, subj_idx) = erp_pas_late_unalt_corrected;

    % 'Passive - Early - Altered' condition
    erp_pas_early_alt_corrected = all_ERP_cond.pas_early_alt(:,:,:) - erp_con_pas_early;
    all_ERP_pas_early_alt_corrected(:,:, subj_idx) = erp_pas_early_alt_corrected;

    % 'Passive - Late - Altered' condition
    erp_pas_late_alt_corrected = all_ERP_cond.pas_late_alt(:,:,:) - erp_con_pas_late;
    all_ERP_pas_late_alt_corrected(:,:, subj_idx) = erp_pas_late_alt_corrected;

    % Update Protocol
    subj_time = toc;
    protocol{subj_idx,1} = subj;
    protocol{subj_idx,2} = subj_time;
    protocol{subj_idx,3} = 'OK';

end

% Plot ERPs
% Plot 1: ERPs collapsed over 'Active' and 'Passive' Condition
figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
subplot(1,2,1)
plot(EEG.times, mean(all_ERP_act_all(CHANI,:,:),3))
xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Active (collapsed)')

subplot(1,2,2)
plot(EEG.times, mean(all_ERP_pas_all(CHANI,:,:),3))
xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Passive (collapsed)')

sgtitle('Active and Passive Condition ERPs (collapsed)')

saveas(gcf,fullfile(OUTPATH, 'erp_collapsed_pilot.png'))

% Plot 2: ERPs for each condition
figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
colors = lines(4);

subplot(1,2,1)
hold on;
plot(EEG.times, mean(all_ERP_cond.act_early_unalt(CHANI,:,:),3), 'Color', colors(1,:));
plot(EEG.times, mean(all_ERP_cond.act_early_alt(CHANI,:,:),3), 'Color', colors(2,:));
plot(EEG.times, mean(all_ERP_cond.act_late_unalt(CHANI,:,:),3), 'Color', colors(3,:));
plot(EEG.times, mean(all_ERP_cond.act_late_alt(CHANI,:,:),3), 'Color', colors(4,:));
hold off;

xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Active Conditions')
legend({'Early Unaltered', 'Early Altered', 'Late Unaltered', 'Late Altered'}, 'Location', 'northeast')

subplot(1,2,2)
hold on;
plot(EEG.times, mean(all_ERP_cond.pas_early_unalt(CHANI,:,:),3), 'Color', colors(1,:));
plot(EEG.times, mean(all_ERP_cond.pas_early_alt(CHANI,:,:),3), 'Color', colors(2,:));
plot(EEG.times, mean(all_ERP_cond.pas_late_unalt(CHANI,:,:),3), 'Color', colors(3,:));
plot(EEG.times, mean(all_ERP_cond.pas_late_alt(CHANI,:,:),3), 'Color', colors(4,:));
hold off;

xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Passive Conditions')
legend({'Early Unaltered', 'Early Altered', 'Late Unaltered', 'Late Altered'}, 'Location', 'northeast')

sgtitle('Active and Passive Condition ERPs')

saveas(gcf,fullfile(OUTPATH, 'erp_pilot.png'))

% Plot 3: ERPs collapsed over 'Active' and 'Passive' Condition (corrected)
figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
subplot(1,2,1)
plot(EEG.times, mean(all_ERP_act_all_corrected(CHANI,:,:),3))
xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Active (collapsed, corrected)')

subplot(1,2,2)
plot(EEG.times, mean(all_ERP_pas_all_corrected(CHANI,:,:),3))
xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Passive (collapsed, corrected)')

sgtitle('Active and Passive Condition ERPs (collapsed, corrected)')

saveas(gcf,fullfile(OUTPATH, 'erp_collapsed_corrected_pilot.png'))

% Plot 4: ERPs for each condition (corrected)
figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
colors = lines(4);

subplot(1,2,1)
hold on;
plot(EEG.times, mean(all_ERP_act_early_unalt_corrected(CHANI,:,:),3), 'Color', colors(1,:));
plot(EEG.times, mean(all_ERP_act_early_alt_corrected(CHANI,:,:),3), 'Color', colors(2,:));
plot(EEG.times, mean(all_ERP_act_late_unalt_corrected(CHANI,:,:),3), 'Color', colors(3,:));
plot(EEG.times, mean(all_ERP_act_late_alt_corrected(CHANI,:,:),3), 'Color', colors(4,:));
hold off;

xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Active Conditions (corrected)')
legend({'Early Unaltered', 'Early Altered', 'Late Unaltered', 'Late Altered'}, 'Location', 'northeast')

subplot(1,2,2)
hold on;
plot(EEG.times, mean(all_ERP_pas_early_unalt_corrected(CHANI,:,:),3), 'Color', colors(1,:));
plot(EEG.times, mean(all_ERP_pas_early_alt_corrected(CHANI,:,:),3), 'Color', colors(2,:));
plot(EEG.times, mean(all_ERP_pas_late_unalt_corrected(CHANI,:,:),3), 'Color', colors(3,:));
plot(EEG.times, mean(all_ERP_pas_late_alt_corrected(CHANI,:,:),3), 'Color', colors(4,:));
hold off;

xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Passive Conditions (corrected)')
legend({'Early Unaltered', 'Early Altered', 'Late Unaltered', 'Late Altered'}, 'Location', 'northeast')

sgtitle('Active and Passive Condition ERPs (corrected)')

saveas(gcf,fullfile(OUTPATH, 'erp_corrected_pilot.png'))

% Plot 5: ERPs collapsed over 'Active' and 'Passive' Condition (Control) (incl. Topoplots)
grandaverage_ERP_con_act = mean(all_ERP_con_act_all,3);
grandaverage_ERP_con_pas = mean(all_ERP_con_pas_all,3);

[~, win_start] = min(abs(EEG.times - (-200)));
[~, win_end] = min(abs(EEG.times - 300));


figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
subplot(2,2,1)
plot(EEG.times, mean(all_ERP_con_act_all(CHANI,:,:),3))
xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Active (collapsed)')

subplot(2,2,2)
plot(EEG.times, mean(all_ERP_con_pas_all(CHANI,:,:),3))
xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Passive (collapsed)')

subplot(2,2,3)
topoplot(mean(grandaverage_ERP_con_act(:,win_start:win_end),2), EEG.chanlocs, ...
    'emarker2', {CHANI,'o','r',5,1});
colormap("parula")
cb = colorbar;
title(cb, 'Amplitude [µV]')
clim([CB_LIM_LOWER CB_LIM_UPPER])
title('Act')

subplot(2,2,4)
topoplot(mean(grandaverage_ERP_con_pas(:,win_start:win_end),2), EEG.chanlocs, ...
    'emarker2', {CHANI,'o','r',5,1});
colormap("parula")
cb = colorbar;
title(cb, 'Amplitude [µV]')
clim([CB_LIM_LOWER CB_LIM_UPPER])
title('Pas')

sgtitle('Active and Passive Condition ERPs (collapsed, control)')

saveas(gcf,fullfile(OUTPATH, 'erp_collapsed_control_pilot.png'))

% Plot 6: ERPs for each condition (corrected)
figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
colors = lines(4);

subplot(1,2,1)
hold on;
plot(EEG.times, mean(all_ERP_act_early_unalt_corrected(CHANI,:,:),3), 'Color', colors(1,:));
plot(EEG.times, mean(all_ERP_act_early_alt_corrected(CHANI,:,:),3), 'Color', colors(2,:));
plot(EEG.times, mean(all_ERP_pas_early_unalt_corrected(CHANI,:,:),3), 'Color', colors(3,:));
plot(EEG.times, mean(all_ERP_pas_early_alt_corrected(CHANI,:,:),3), 'Color', colors(4,:));
hold off;

xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Early (corrected)')
legend({'Active Early Unaltered', 'Active Early Altered', ' Passive Early Unaltered', 'Passive Early Altered'}, 'Location', 'northeast')

subplot(1,2,2)
hold on;
plot(EEG.times, mean(all_ERP_act_late_unalt_corrected(CHANI,:,:),3), 'Color', colors(1,:));
plot(EEG.times, mean(all_ERP_act_late_alt_corrected(CHANI,:,:),3), 'Color', colors(2,:));
plot(EEG.times, mean(all_ERP_pas_late_unalt_corrected(CHANI,:,:),3), 'Color', colors(3,:));
plot(EEG.times, mean(all_ERP_pas_late_alt_corrected(CHANI,:,:),3), 'Color', colors(4,:));
hold off;

xlim([-200 500])
ylabel('Amplitude [μV]')
xlabel('Time [ms]')
title('Late (corrected)')
legend({'Active Late Unaltered', 'Active Late Altered', ' Passive Late Unaltered', 'Passive Late Altered'}, 'Location', 'northeast')

sgtitle('Early and Late Condition ERPs (corrected)')

saveas(gcf,fullfile(OUTPATH, 'erp_corrected_2_pilot.png'))

% % % Plot 7: Topoplots collapsed over 'Active' and 'Passive' Condition
% % CB_LIM_LOWER = min([mean(all_ERP_act_all(:,ERP_FROM:ERP_TILL),2), mean(all_ERP_pas_all(:,ERP_FROM:ERP_TILL),2)],[],'all');
% % CB_LIM_UPPER = max([mean(all_ERP_act_all(:,ERP_FROM:ERP_TILL),2), mean(all_ERP_pas_all(:,ERP_FROM:ERP_TILL),2)],[],'all');
% % 
% % figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
% % subplot(1,2,1)
% % topoplot(mean(all_ERP_act_all(:,win_start:win_end),2), EEG.chanlocs, 'emarker2', {CHANI,'o','r',5,1})
% % title('Active (collapsed)')
% % colormap("parula")
% % cb = colorbar;
% % title(cb, 'Amplitude [µV]')
% % clim([CB_LIM_LOWER CB_LIM_UPPER])
% % subplot(1,2,2)
% % topoplot(mean(all_ERP_pas_all(:,win_start:win_end),2), EEG.chanlocs, 'emarker2', {CHANI,'o','r',5,1})
% % title('Passive (collapsed)')
% % colormap("parula")
% % cb = colorbar;
% % title(cb, 'Amplitude [µV]')
% % clim([CB_LIM_LOWER CB_LIM_UPPER])
% % 
% % sgtitle('Active and Passive Condition Topoplots (collapsed and averaged ober 70ms-130ms)')
% % 
% % saveas(gcf,fullfile(OUTPATH, 'topo_collapsed_pilot.png'))

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'erp_pilot_protocol.xlsx'))

delete(wb); %%close all;










