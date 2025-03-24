% tid_psam_timingtest.m
%
% Performs timing test.
%
% Tim Dressler, 07.03.2025

clear 
close all
clc

% Start eeglab
eeglab;

% Load data
EEG = pop_loadbv('C:\Users\timdr\OneDrive\Uni_Oldenburg\4_Semester\Master_Thesis\Analysis_Experiment\psam\testing\timingtest\timingtest_240325\sub-99\eeg\', 'timingtest_24032025.vhdr', [], []);
EEG = eeg_checkset( EEG );

% Epoch and BL-Correction
EEG2 = pop_epoch( EEG, {'S 31', 'S 32','S 33','S 34','S 41','S 42','S 43','S 44'},[0  0.2], 'epochinfo', 'yes');
EEG2 = pop_rmbase( EEG2, [ ]);

% Plot: Sound 
fig=figure;
set(fig,'defaultTextFontSize',14);
set(fig,'defaultAxesFontSize',14);
subplot(1,2,1)
for k=1:size(EEG2.data,3)
plot(EEG2.times,(EEG2.data(1,:,k)))
hold on
end
xlabel('time (ms)')
title('Audio marker')

% Find sound onsets
for k=1:size(EEG2.data,3)
   a=find(abs(squeeze(EEG2.data(1,1:end,k))')>1000);  
   sound_onsets(k)=EEG2.times(a(1));
   
end

% Plot: Latency
subplot(1,2,2)
plot(sound_onsets,'bo')
xlabel('trial')
ylabel('sound onset relative to trigger (ms)')

% Display results
disp(['Mean Latency: ' num2str(mean(sound_onsets))])
disp(['SD Latency: ' num2str(std(sound_onsets))])
