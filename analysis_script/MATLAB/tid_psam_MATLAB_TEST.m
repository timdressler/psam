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
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\svm_analysis');
if ~exist(OUTPATH, 'dir'), mkdir(OUTPATH); end

% File selection
files = dir(fullfile(INPATH, 'sub-*early.csv'));

% Grid settings
lower_C = -20;
upper_C = 20;
C_range = logspace(lower_C, upper_C, 32);
gamma_range = logspace(lower_C, upper_C, 32);
nSplits = 5;

protocol = {};

for f = 1:length(files)
    file = files(f).name;
    subj = regexp(file, 'sub-\d+', 'match', 'once');
    fprintf('Processing %s\n', subj);
    tic;

    T = readtable(fullfile(INPATH, file));
    X = T{:, 1:end-1};  % Features
    y = T{:, end};      % Labels

    % Compute means and stds by class
    class_labels = unique(y);
    means = zeros(length(class_labels), size(X, 2));
    sds = zeros(length(class_labels), size(X, 2));

    for c = 1:length(class_labels)
        class = class_labels(c);
        if iscell(y)
            class_data = X(strcmp(y, class), :);
        else
            class_data = X(y == class, :);
        end
        means(c, :) = mean(class_data, 1);
        sds(c, :) = std(class_data, 0, 1);
    end
    disp('Means by class:');
    disp(means);
    disp('Standard deviations by class:');
    disp(sds);

    acc_matrix = zeros(length(gamma_range), length(C_range));
    cv = cvpartition(y, 'KFold', nSplits);

    for i = 1:length(gamma_range)
        gamma_val = gamma_range(i);
        for j = 1:length(C_range)
            C_val = C_range(j);
            fold_accuracies = zeros(nSplits, 1);

            for fold = 1:nSplits
                trainIdx = training(cv, fold);
                testIdx = test(cv, fold);

                X_train = X(trainIdx, :);
                y_train = y(trainIdx);
                X_test = X(testIdx, :);
                y_test = y(testIdx);

                % PCA
                [coeff, score_train, ~, ~, explained] = pca(X_train);
                X_train_pca = score_train(:, 1:120);
                X_test_pca = (X_test - mean(X_train)) * coeff(:, 1:120);

                fprintf("Variance retained: %.2f\n", sum(explained(1:120)));

                % SVM
                template = templateSVM('KernelFunction', 'rbf', ...
                    'BoxConstraint', C_val, 'KernelScale', 1/sqrt(2*gamma_val));
                model = fitcecoc(X_train_pca, y_train, 'Learners', template);
                y_pred = predict(model, X_test_pca);

                if iscell(y_test)
                    acc = mean(strcmp(y_pred, y_test));
                else
                    acc = mean(y_pred == y_test);
                end

                fold_accuracies(fold) = acc;
            end

            acc_matrix(i, j) = mean(fold_accuracies);
            fprintf('C: %.2e, Gamma: %.2e, Acc: %.4f\n', C_val, gamma_val, acc_matrix(i,j));
        end
    end

    elapsed = toc;
    protocol{end+1,1} = subj;
    protocol{end,2} = elapsed;
    protocol{end,3} = 'OK';
    protocol{end,4} = mean(acc_matrix(:));

    % Save accuracy matrix
    csvwrite(fullfile(OUTPATH, [subj '_acc_grid.csv']), acc_matrix);

    % Plot heatmap
    figure('Visible', 'off');
    imagesc(log10(C_range), log10(gamma_range), acc_matrix);
    colorbar;
    xlabel('log_{10}(C)');
    ylabel('log_{10}(Gamma)');
    title(['Validation Accuracy Grid: ' subj]);
    colormap(parula);  % magma isn't available in MATLAB by default
    set(gca, 'YDir', 'normal');
    saveas(gcf, fullfile(OUTPATH, [subj '_heatmap.png']));
    close;
end

% Save protocol
protocol_table = cell2table(protocol, ...
    'VariableNames', {'Subject', 'Time', 'Status', 'MeanAccuracy'});
writetable(protocol_table, fullfile(OUTPATH, 'svm_analysis_protocol.xlsx'));

disp('tid_psam_svm_analysis_DONE');
