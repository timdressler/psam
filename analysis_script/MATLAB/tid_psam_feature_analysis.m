% tid_psam_feature_analysis.m
%
% Performs .
%
% Processing includes the following steps
%
% Removes EOG electrodes

%
% Saves data
%
% Also see: tid_psam_hjorth_activity_TD, tid_psam_hjorth_mobility_TD, tid_psam_hjorth_complexity_TD.
%
% Tim Dressler, 05.05.2025

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
INPATH = fullfile(MAINPATH, 'data\processed_data\svm_prepared_clean\');
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\feature_analysis');

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

% Variables to edit
EOG_CHANI = [29,30]; % Numbers of EOG electrodes
WINDOWS = {'early', 'late'};
CB_LIM_LOWER = -0.3;
CB_LIM_UPPER = 0.3;

% Get directory content
dircont_subj = dir(fullfile(INPATH, 'sub-*early.csv'));

% Prepare channel locations
% Start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% Add channel locations
EEG.chanlocs = readlocs( fullfile(MAINPATH,'\config\elec_96ch_adapted.elp'));
% Remove EOG channels as they are not used for classification
EEG.chanlocs(EOG_CHANI) = [];

%initialize sanity check variables
marked_subj = {};
protocol = {};

% Setup progress bar
%%wb = waitbar(0,'starting tid_psam_feature_analysis.m');

for subj_idx= 1:length(dircont_subj)

    % Get current ID
    subj = dircont_subj(subj_idx).name;
    subj = regexp(subj, 'sub-\d+', 'match', 'once');

    % Get subject outpath
    subj_outpath = fullfile(OUTPATH,subj);
    mkdir(subj_outpath);

    for win = WINDOWS

        % Load data
        filename = fullfile(INPATH, [subj '_features_' win{1} '.csv']);
        features_labels = readtable(filename);

        % Extract data
        features = features_labels(:,1:end-1); % Get features
        feature_names = features.Properties.VariableNames; % Get feature names
        features = table2array(features);
        labels = table2array(features_labels(:,end)); % Get labels

        % Recode labels (0 = go_pas, 1 = go_act)
        labels_num = zeros(length(labels),1);
        labels_num(strcmp(labels, 'go_act')) = 1;

        % Correlate each feature with the outcome
        for i = 1:size(features,2)
            [r, p] = corrcoef(features(:,i), labels_num);  % Pearson correlation
            corr_results(i) = r(1,2);
            p_values(i) = p(1,2);
        end

        corr_table = table(feature_names', corr_results', p_values', ...
            'VariableNames', {'feature', 'correlation', 'p_value'});

        % Group features by type prefix
        prefixes = cellfun(@(x) regexp(x, '^(.*)_E\d+', 'tokens', 'once'), corr_table.feature, 'UniformOutput', false);
        prefixes = cellfun(@(x) x{1}, prefixes, 'UniformOutput', false);  % Extract nested cell
        unique_prefixes = unique(prefixes);

        for feature = 1:length(unique_prefixes)
            prefix = unique_prefixes{feature};
            mask = startsWith(corr_table.feature, [prefix '_']);
            grouped_table = corr_table(mask, :);

            % Extract electrode names
            electrodes = cellfun(@(x) regexp(x, '(?<=_)(E\d+)', 'match', 'once'), grouped_table.feature, 'UniformOutput', false);
            grouped_table.electrode = electrodes;

            % Store in window-specific struct
            features_grouped_win(feature).feature = prefix;
            features_grouped_win(feature).table = grouped_table;
        end

        % Assign to main struct under current window
        features_grouped.(win{1}) = features_grouped_win;
    end
    clear feature

    % Sanity Check: Same size of features_grouped for early and late windows


    for feature = 1:size(features_grouped.early,2)
        feature_early_name = features_grouped.early(feature).feature;
        feature_late_name = features_grouped.late(feature).feature;

        % Sanity Check: Same feature for early and late windows
        if ~strcmp(feature_early_name, feature_late_name)
            error('Non-matching feature!')
        end

        feature_early_data = table2cell(features_grouped.early(feature).table(:, {'correlation', 'electrode'}));
        feature_late_data = table2cell(features_grouped.late(feature).table(:, {'correlation', 'electrode'}));

        % Sanity Check: Order of features matches EEG.chanlocs
        if ~(strcmp(cell2mat({EEG.chanlocs.labels}'), cell2mat(feature_early_data(:,2))) && strcmp(cell2mat({EEG.chanlocs.labels}'), cell2mat(feature_late_data(:,2))))
            error('Non-matching electrodes!')
        end

        % Topoplots: Correlation with outcome
        figure;
        subplot(121)
        topoplot(cell2mat(feature_early_data(:,1)), EEG.chanlocs);
        title('Early Window')
        colormap("parula")
        cb = colorbar;
        clim([CB_LIM_LOWER CB_LIM_UPPER]);
        cb.Label.String = 'Correlation (r)';

        subplot(122)
        topoplot(cell2mat(feature_late_data(:,1)), EEG.chanlocs);
        title('Late Window')
        colormap("parula")
        cb = colorbar;
        clim([CB_LIM_LOWER CB_LIM_UPPER]);
        cb.Label.String = 'Correlation (r)';

        sgtitle([subj ' - Correlation of ' feature_early_name ' and Outcome'])

        % Save plot
        saveas(gcf,fullfile(subj_outpath, [subj '_topo_correlation_' feature_early_name '.png']))


    end



end

























