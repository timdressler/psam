% tid_psam_helper.m
%
% Performs needed helper tasks.
%
% Tim Dressler, 28.07.2025

clear
close all
clc

rng(123)
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
FUNPATH = fullfile(MAINPATH, 'functions');

addpath(FUNPATH);

% Copy individual accuracy heatmaps from the SVM analysis to one folder (for convinience)
tid_psam_copy_files_from_subfolders_TD(...
      'C:\Users\timdr\OneDrive\Uni_Oldenburg\4_Semester\Master_Thesis\Analysis_Experiment\psam\data\analysis_data\svm_analysis', ...
      'sub-*', ...
      'sub-*_accuracy_heatmaps.png', ...
      'individual_accuracy_heatmaps');