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
% Reject bad epochs using threshold and probability (commented out, due to severe artifacts, leaving no epochs unaffected)
%
% Analysis includes the following steps
%
% Load the preprocessed dataset
% Extractes ERP over all conditions
% Extracts ERP for each condition
% Stores ERP amplitudes and latencies
% Plots ERP over all conditions
% Plots ERP for each condition
%
% Tim Dressler, 03.04.2025

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
OUTPATH = fullfile(MAINPATH, 'data\processed_data\erp_pilot');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EPO_FROM = -0.25;
EPO_TILL = 0.750;
LCF = 1;
HCF = 60;
BL_FROM = -250;
THRESH = 75;
SD_PROB = 3;
RESAM_ICA = 250;
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ...
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt', 'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};

CHANI = 1;
ERP_FROM = 75;
ERP_TILL = 125;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*'));

%initialize sanity check variables
marked_subj = {};
protocol = {};
all_epo_lats_OK_ALL = [];

% Setup progress bar
wb = waitbar(0,'starting tid_psam_erp_pilot.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)
    % Start preprocessing

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' pp_erp_pre_proc.m'])

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

    % Epoching
    EEG = pop_epoch( EEG, EVENTS, [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');

    % Baseline-Removal
    EEG = pop_rmbase( EEG, [BL_FROM 0] ,[]);

    % Threshold removal
    % % EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan] ,-THRESH,THRESH,EPO_FROM,EPO_TILL,0,0);

    % Probability-based removal
    % % EEG = pop_jointprob(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
    % % EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
    % % EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    % % EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);

    EEG.setname = [subj '_preprocessed'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % End of preprocessing

    % Sart analysis

    % Get all 'Active' condition
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'act_early_alt','act_early_unalt','act_late_alt','act_late_unalt'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [subj '_act_all'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    % Get ERP
    erp_act_all = mean(EEG.data(CHANI,:,:),3);
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
    erp_pas_all = mean(EEG.data(CHANI,:,:),3);
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
figure;
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

exportgraphics(gcf,fullfile(OUTPATH, 'erp_collapsed_pilot.png'),'Resolution',1000)

% Plot 2: ERPs for each condition
figure;
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
legend({'Early Unaltered', 'Early Altered', 'Late Unaltered', 'Late Altered'}, 'Location', 'best')

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
legend({'Early Unaltered', 'Early Altered', 'Late Unaltered', 'Late Altered'}, 'Location', 'best')

sgtitle('Active and Passive Condition ERPs')

exportgraphics(gcf,fullfile(OUTPATH, 'erp_pilot.png'),'Resolution',1000)

% Plot 3: ERPs collapsed over 'Active' and 'Passive' Condition (corrected)
figure;
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

exportgraphics(gcf,fullfile(OUTPATH, 'erp_collapsed_pilot.png'),'Resolution',1000)

% Plot 4: ERPs for each condition (corrected)
figure;
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
legend({'Early Unaltered', 'Early Altered', 'Late Unaltered', 'Late Altered'}, 'Location', 'best')

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
legend({'Early Unaltered', 'Early Altered', 'Late Unaltered', 'Late Altered'}, 'Location', 'best')

sgtitle('Active and Passive Condition ERPs (corrected)')

exportgraphics(gcf,fullfile(OUTPATH, 'erp_pilot.png'),'Resolution',1000)

% End of processing

protocol = cell2table(protocol, 'VariableNames',{'subj','time', 'status'})
writetable(protocol,fullfile(OUTPATH, 'erp_pilot_protocol.xlsx'))

close(wb)














