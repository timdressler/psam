% Run tid_psam_set_markers.m
%
% Runs MATLAB-based analysis pipeline.
%
% Note. tid_psam_set_markers.m needs to be run manually before running the
%   Praat script as the latter is not able to create the needed folders.
% 
% Tim Dressler, 17.04.2025

% To run the pipeline, run the scripts in the specified order

% Run tid_psam_set_markers.m

% Run tid_psam_beh_preprocessing_1.praat

tid_psam_ica_preprocessing
tid_psam_erp_preprocessing
tid_psam_svm_preprocessing
tid_psam_beh_preprocessing_2
tid_psam_exclude_trials
tid_psam_erp_analysis
tid_psam_hilbert_preparation
tid_psam_svm_preparation
tid_psam_beh_analysis

% Run tid_psam_beh_analysis.R
% Run tid_psam_erp_analysis.R
% Run tid_psam_svm_analysis.py
% Run tid_psam_feature_analysis.py

% Run tid_psam_questionnaire_analysis.R