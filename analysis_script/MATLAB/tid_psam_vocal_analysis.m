% tid_psam_vocal_analysis.m
%
% Performs analysis of vocal data.
% Has to be executed AFTER tid_psam_vocal_analysis.praat.
%
% Tim Dressler, 04.04.2025

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
INPATH = fullfile(MAINPATH, 'data\BIDS\');
INPATH_PRAAT = fullfile(MAINPATH, 'data\analysis_data\vocal_analysis_preparation\');
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\vocal_analysis');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
n_trials = 960;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*'));

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
%%wb = waitbar(0,'starting tid_psam_vocal_analysis.m');

clear subj_idx
for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Update progress bar
    %%waitbar(subj_idx/length(dircont_subj),wb, [subj ' tid_psam_vocal_analysis.m'])

    tic;

    % Get log file & Sanity Check: One log file per subject
    subj_log_filename = dir(fullfile(INPATH, [subj '\beh\*.csv']));
    if numel(subj_log_filename) == 1
        subj_log = readtable(fullfile(subj_log_filename.folder, subj_log_filename.name),'VariableNamingRule', 'preserve');
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'number_of_log_files';
    end

    % Remove all rows from log file which don not resembel experimental trials (e.g. Instruction trials)
    subj_log_cleaned = subj_log(~isnan(subj_log.("mic.started")), :);
    % Sanity Check: Equal db for unaltered and altered probes
    if height(subj_log_cleaned) == n_trials
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'log_dimensions'; 
    end
    % Extract filename from full path for each vocal recording in order to match it with the table produced in Praat
    subj_log_cleaned.filename_tab = cellfun(@(x) strrep(regexp(x, 'recording_mic_.*(?=\.wav)', 'match', 'once'), '.', '_'), subj_log_cleaned.("mic.clip"), 'UniformOutput', false);

    % Get probe properties file & Sanity Check: Correct dimensions
    subj_prope_properties = readtable(fullfile(INPATH, ['stimuli\' subj '\' subj '_probe_properties.xlsx']),'VariableNamingRule', 'preserve');
    if height(subj_prope_properties) == 1 && width(subj_prope_properties) == 7
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'probe_properties_dimensions';
    end
    % Sanity Check: Equal db for unaltered and altered probes
    if round(subj_prope_properties.db_tab_normal) == round(subj_prope_properties.db_tab_pitched)
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'probe_properties_db';
    end

    % Get F0 & RT Table & Sanity Check: Correct dimensions
    subj_f0_rt = readtable(fullfile(INPATH_PRAAT, [subj '_f0_rt_table.csv']),'VariableNamingRule', 'preserve', Delimiter=',');

    if height(subj_f0_rt) == n_trials && width(subj_f0_rt) == 14
    else
        marked_subj{end+1,1} = subj;
        marked_subj{end,2} = 'f0_rt_dimensions';
    end

    % Concatinate log file and F0 & RT Table
    subj_full = join(subj_f0_rt, subj_log_cleaned, 'Keys', 'filename_tab');






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
writetable(protocol,fullfile(OUTPATH, 'tid_psam_vocal_analysis_protocol.xlsx'))

if ~isempty(marked_subj)
    marked_subj = cell2table(marked_subj, 'VariableNames',{'subj','issue'})
    writetable(marked_subj,fullfile(OUTPATH, 'tid_psam_vocal_analysis_marked_subj.xlsx'))
end

%%close(wb)

check_done = 'tid_psam_vocal_analysis_DONE'
