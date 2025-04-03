% tid_psam_erp_pilot.m
%
% Performs a 'quick and dirty' ERP analysis (incl. preprocessing).
%
% Preprocessing includes the following steps
%
    % Loads raw data
    % Renames the events
    % Applies band-pass filter
    % Epochs data 
    % Performs baseline correction
    % Reject bad epochs using threshold and probability
    % Stores dataset
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

% Set up paths
SCRIPTPATH = cd;
if regexp(SCRIPTPATH, regexptranslate('wildcard','*psam\experiment_script')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, '\analysis_script\MATLAB');
FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH)















