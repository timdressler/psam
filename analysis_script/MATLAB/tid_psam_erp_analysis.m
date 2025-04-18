% tid_psam_erp_analysis.m
%
% Performs ERP analysis and exports data for further analysis.

%
% Tim Dressler, 17.04.2025

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
INPATH = fullfile(MAINPATH, 'data\processed_data\erp_preprocessed_clean\');
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\erp_analysis');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EVENTS = {'act_early_unalt', 'act_early_alt', 'act_late_unalt', 'act_late_alt', ... % Real events
    'pas_early_unalt', 'pas_early_alt', 'pas_late_unalt', 'pas_late_alt'};
CON_EVENTS = {'con_act_early', 'con_act_late', ...
    'con_pas_early', 'con_pas_late'};
CHAN = 1;
ERP_FROM = 70;
ERP_TILL = 130;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*_clean.set'));

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
wb = waitbar(0,'starting tid_psam_erp_analysis.m');

clear subj_idx
counter = 1;
cor_counter = 1;
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_erp_analysis.m'])

    tic;
    % Start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Load data
    EEG = pop_loadset('filename',[subj '_erp_preprocessed_clean.set'],'filepath',INPATH);
    EEG.setname = [subj '_all_conds'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    % Get ERPs for each condition
    for cond = 1:length(EVENTS)
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off');

        % Get ERP
        erp = mean(EEG.data, 3);

        % Correct ERP by subtracting the ERP from the equivalent no-probe condition (see Daliri & Max, 2016)
        % Get matching control condition
        con_cond = find(strcmp(erase(EVENTS{cond}, {'_alt', '_unalt'}), erase(string(CON_EVENTS), 'con_')));
        % Get control ERP
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',CON_EVENTS(con_cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off');
        con_erp = mean(EEG.data, 3);
        % Get corrected ERP
        cor_erp = erp - con_erp;

        % Setup ERP analysis
        [~,win_start] = min(abs(EEG.times-ERP_FROM));
        [~,win_end] = min(abs(EEG.times-ERP_TILL));
        [~,t_zero] = min(abs(EEG.times));

        % ERP analysis (uncorrected)
        % Get N100 amplitude (uncorrected)
        erp_amp = min(erp(CHAN,win_start:win_end));
        % Get N100 latency (uncorrected)
        erp_sam = find(erp(CHAN,:) == erp_amp);
        erp_lat = EEG.times(erp_sam);
        % Store ERP
        all_erp{counter,1} = erp;
        % Store ERP data
        all_erp_data{counter,1} = ERP_FROM;
        all_erp_data{counter,2} = ERP_TILL;
        all_erp_data{counter,3} = erp_amp;
        all_erp_data{counter,4} = erp_lat;
        all_erp_data{counter,5} = EVENTS{cond};
        all_erp_data{counter,6} = subj;
        counter = counter+1;

        % ERP analysis (corrected)
        % Get N100 amplitude (corrected)
        cor_erp_amp = min(cor_erp(CHAN,win_start:win_end));
        % Get N100 latency (corrected)
        cor_erp_sam = find(cor_erp(CHAN,:) == cor_erp_amp);
        cor_erp_lat = EEG.times(cor_erp_sam);
        % Store ERP
        all_cor_erp{cor_counter,1} = cor_erp;
        % Store ERP data
        all_cor_erp_data{cor_counter,1} = ERP_FROM;
        all_cor_erp_data{cor_counter,2} = ERP_TILL;
        all_cor_erp_data{cor_counter,3} = cor_erp_amp;
        all_cor_erp_data{cor_counter,4} = cor_erp_lat;
        all_cor_erp_data{cor_counter,5} = EVENTS{cond};
        all_cor_erp_data{cor_counter,6} = subj;
        cor_counter = cor_counter+1;

        % Plot: ERP, control ERP and corrected ERP
        ylim_max = max([erp; con_erp; cor_erp],[], 'all');
        ylim_min = min([erp; con_erp; cor_erp],[], 'all');
        figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
        subplot(131)
        plot(EEG.times, erp(CHAN,:))
        hold on
        scatter(erp_lat, erp_amp)
        hold off
        ylabel('Amplitude [μV]')
        xlabel('Time [ms]')
        ylim([ylim_min ylim_max])
        title('Uncorrected')
        subplot(132)
        plot(EEG.times, con_erp(CHAN,:))
        ylabel('Amplitude [μV]')
        xlabel('Time [ms]')
        ylim([ylim_min ylim_max])
        title('Control')
        subplot(133)
        plot(EEG.times, cor_erp(CHAN,:))
        hold on
        scatter(cor_erp_lat, cor_erp_amp)
        hold off
        ylabel('Amplitude [μV]')
        xlabel('Time [ms]')
        ylim([ylim_min ylim_max])
        title('Corrected')
        sgtitle(['ERPs for ' subj ' and Condition ' EVENTS{cond}])

    end


    % Save dataset ERP data

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
writetable(protocol,fullfile(OUTPATH, 'erp_analysis_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_erp_analysis_marked_subj.xlsx'))
end

check_done = 'tid_psam_erp_preprocessing_DONE'

delete(wb); %%close all;













