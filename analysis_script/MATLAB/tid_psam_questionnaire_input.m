% tid_psam_beh_analysis.m
%
% Performs analysis of vocal data.
%
% Tim Dressler, 07.04.2025

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
INPATH = fullfile(MAINPATH, 'data\questionnaire_data\');
OUTPATH = INPATH;

FUNPATH = fullfile(MAINPATH, '\functions\');
addpath(FUNPATH);

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)

% Get filename
fal_filename = fullfile(INPATH, 'fal_data.xlsx');

% Load or create FAL data
if exist(fal_filename, 'file') ~= 2
    fal_varNames = {
        'subj', 'var1_age', 'var2_handedness', 'var3_sex', 'var4_education', ...
        'var5_occupation', 'var5_hearing_problems', 'var5_hearing_problems_detail', ...
        'var7_ringing_ears', 'var8_sleep', 'var9_sleep_assessment', ...
        'var10_alcohol_yesterday', 'var10_alcohol_yesterday_detail', ...
        'var11_alcohol_today', 'var11_alcohol_today_detail', ...
        'var12_smoking', 'var13_coffee_and_other', 'var14_last_eating', ...
        'var15_currently_neuro_treatment', 'var15_currently_neuro_treatment_detail', ...
        'var16_earlier_neuro_treatment', 'var16_earlier_neuro_treatment_detail', ...
        'var17_other_treatment', 'var17_other_treatment_detail', ...
        'var18_medication', 'var18_medication_detail', ...
        'var19_drugs', 'var19_drugs_detail'
        };
    fal = cell2table(cell(0, numel(fal_varNames)), 'VariableNames', fal_varNames);
else
    fal = readtable(fal_filename);
end

% Get prompt for each variable
fal_promptTexts = containers.Map(...
    fal.Properties.VariableNames, ...
    { ...
    'subj: (Angeben: Nur Zahl)', ...
    'var1_age: (Angeben: Nur Zahl, in Jahren)', ...
    'var2_handedness: (Angeben: left, right or both)', ...
    'var3_sex: (Angeben: m,w,d)', ...
    'var4_education: (Angeben: 0 = kein Abschluss, 1 = Hauptschule, 2 = Mittlere Reife, 3 = Abitur)', ...
    'var5_occupation: (Angeben: 1 = Student, 2 = Berufstätig, 3 = Arbeitslos)', ...
    'var5_hearing_problems: (Angeben: yes or no)', ...
    'var5_hearing_problems_detail: (Spezifizieren)', ...
    'var7_ringing_ears: (Angeben: yes or no)', ...
    'var8_sleep: (Angeben: Nur Zahl, in Stunden)', ...
    'var9_sleep_assessment: (Angeben: 1 = Normal, 2 = Eher lang, 3 = Viel zu kurz)', ...
    'var10_alcohol_yesterday: (Angeben: yes or no)', ...
    'var10_alcohol_yesterday_detail: (Spezifizieren)', ...
    'var11_alcohol_today: (Angeben: yes or no)', ...
    'var11_alcohol_today_detail: (Spezifizieren)', ...
    'var12_smoking: (Angeben: 0 = kein Raucher, 1 = wenig / gar nicht, 2 = normal, 3 = sehr viel)', ...
    'var13_coffee_and_other: (Angeben: 0 = wenig / keine, 2 = normal, 3 = sehr viel)', ...
    'var14_last_eating: (Angeben: Nur Zahl, in Stunden)', ...
    'var15_currently_neuro_treatment: (Angeben: yes or no)', ...
    'var15_currently_neuro_treatment_detail: (Spezifizieren)', ...
    'var16_earlier_neuro_treatment: (Angeben: yes or no)', ...
    'var16_earlier_neuro_treatment_detail: (Spezifizieren)', ...
    'var17_other_treatment: (Angeben: yes or no)', ...
    'var17_other_treatment_detail: (Spezifizieren)', ...
    'var18_medication: (Angeben: yes or no)', ...
    'var18_medication_detail: (Spezifizieren)', ...
    'var19_drugs: (Angeben: yes or no)', ...
    'var19_drugs_detail: (Spezifizieren)' ...
    });

% Get subject ID
subj = inputdlg(fal_promptTexts('subj'), 'Subject Entry', [1 50]);
if isempty(subj)
    disp('Cancelled.');
    return;
end

% Format subject ID
num = str2double(subj{1}); % Convert input to a number
formattedNum = sprintf('%02d', num); % Ensure two-digit format
subj = ['sub-' formattedNum]; % Construct final subject ID

% Check if the subject already exists in the table
if any(strcmp(fal.subj, subj))
    choice = questdlg([subj, ' already exists in FAL. Do you want to overwrite the existing data?'], ...
        'Duplicate Subject', 'Overwrite', 'Cancel', 'Cancel');
    if strcmp(choice, 'Cancel')
        disp('Operation cancelled.');
        fal_skipped = 1;
    elseif strcmp(choice, 'Overwrite')
        disp(['Existing data for ', subj, ' will be overwritten.']);
        % Find the row and overwrite it
        rowIdx = find(strcmp(fal.subj, subj));
        fal(rowIdx, :) = [];
        fal_skipped = 0;
    end
else
    fal_skipped = 0;
end

if fal_skipped ~= 1
    % Create new row
    newRow = cell(1, width(fal));
    newRow{1} = subj;

    % Input values for each variables
    for i = 2:width(fal)
        varName = fal.Properties.VariableNames{i};
        if ~endsWith(varName, '_detail') % Only a
            while true
                answer = inputdlg(fal_promptTexts(varName), 'FAL', [1 70]);
                if isempty(answer)
                    disp('Cancelled.');
                    return;
                end
                % If the answer is empty, keep prompting
                if isempty(answer{1})
                    disp('This field cannot be left empty. Please enter a valid response.');
                else
                    % Otherwise, store the answer and break out of the loop
                    newRow{i} = answer{1};
                    break;
                end
            end
        elseif strcmp(newRow{i-1},'yes') % Only ask for detail if needed
            while true
                answer = inputdlg(fal_promptTexts(varName), 'FAL', [1 70]);
                if isempty(answer)
                    disp('Cancelled.');
                    return;
                end
                % If the answer is empty, keep prompting
                if isempty(answer{1})
                    disp('This field cannot be left empty. Please enter a valid response.');
                else
                    % Otherwise, store the answer and break out of the loop
                    newRow{i} = answer{1};
                    break;
                end
            end
        else
            newRow{i} = 'NaN';
        end
    end

    % Add new row to the table
    newRowTable = cell2table(newRow, 'VariableNames', fal.Properties.VariableNames);
    for i = 1:width(fal)
        if iscell(fal{:, i}) && ~iscell(newRowTable{:, i})
            newRowTable{:, i} = {newRowTable{:, i}};
        end
    end
    fal = [fal; newRowTable];
end

% NASA-TLX
% Get filename
nasatlx_filename = fullfile(INPATH, 'nasatlx_data.xlsx');

% Load or create FAL data
if exist(nasatlx_filename, 'file') ~= 2
    nasatlx_varNames = {
        'subj', 'var1_mental_demand', 'var2_physical_demand', 'var3_performance', 'var4_effort', ...
        'var5_frustration'
        };
    nasatlx = cell2table(cell(0, numel(nasatlx_varNames)), 'VariableNames', nasatlx_varNames);
else
    nasatlx = readtable(nasatlx_filename);
end

% Check if the subject already exists in the table
if any(strcmp(nasatlx.subj, subj))
    choice = questdlg([subj, ' already exists in NASA-TLX. Do you want to overwrite the existing data?'], ...
        'Duplicate Subject', 'Overwrite', 'Cancel', 'Cancel');
    if strcmp(choice, 'Cancel')
        disp('Operation cancelled.');
        nasatlx_skipped = 1;
    elseif strcmp(choice, 'Overwrite')
        disp(['Existing data for ', subj, ' will be overwritten.']);
        % Find the row and overwrite it
        rowIdx = find(strcmp(nasatlx.subj, subj));
        nasatlx(rowIdx, :) = [];
        nasatlx_skipped = 0;
    end
else
    nasatlx_skipped = 0;
end

if nasatlx_skipped ~= 1
    % Create new row
    newRow = cell(1, width(nasatlx));
    newRow{1} = subj;

    % Input values for each variables
    for i = 2:width(nasatlx)
        varName = nasatlx.Properties.VariableNames{i};

        while true
            answer = inputdlg([varName ' (Angeben: Ganze Zahl oder Dezimal mit Punkt (nicht Komma))'], 'NASA-TLX', [1 70]);
            if isempty(answer)
                disp('Cancelled.');
                return;
            end
            % If the answer is empty, keep prompting
            if isempty(answer{1})
                disp('This field cannot be left empty. Please enter a valid response.');
            else
                % Otherwise, store the answer and break out of the loop
                newRow{i} = answer{1};
                break;
            end
        end
    end


    % Add new row to the table
    newRowTable = cell2table(newRow, 'VariableNames', nasatlx.Properties.VariableNames);
    for i = 1:width(nasatlx)
        if iscell(nasatlx{:, i}) && ~iscell(newRowTable{:, i})
            newRowTable{:, i} = {newRowTable{:, i}};
        end
    end
    nasatlx = [nasatlx; newRowTable];
end

% SAM
% Get filename
sam_filename = fullfile(INPATH, 'sam_data.xlsx');

% Load or create FAL data
if exist(sam_filename, 'file') ~= 2
    sam_varNames = {
        'subj', 'var1_mood_break1', 'var2_tiredness_break1', 'var3_mood_break2', 'var4_tiredness_break2', ...
        'var5_mood_break3', 'var6_tiredness_break3','var7_mood_break4', 'var8_tiredness_break4','var9_mood_break5', 'var10_tiredness_break5',...
        'var11_mood_break6', 'var12_tiredness_break6', 'var13_mood_break7', 'var14_tiredness_break7', 'var15_mood_break8', 'var16_tiredness_break8'
        };
    sam = cell2table(cell(0, numel(sam_varNames)), 'VariableNames', sam_varNames);
else
    sam = readtable(sam_filename);
end

% Check if the subject already exists in the table
if any(strcmp(sam.subj, subj))
    choice = questdlg([subj, ' already exists in SAM. Do you want to overwrite the existing data?'], ...
        'Duplicate Subject', 'Overwrite', 'Cancel', 'Cancel');
    if strcmp(choice, 'Cancel')
        disp('Operation cancelled.');
        sam_skipped = 1;
    elseif strcmp(choice, 'Overwrite')
        disp(['Existing data for ', subj, ' will be overwritten.']);
        % Find the row and overwrite it
        rowIdx = find(strcmp(sam.subj, subj));
        sam(rowIdx, :) = [];
        sam_skipped = 0;
    end
else
    sam_skipped = 0;
end

if sam_skipped ~= 1
    % Create new row
    newRow = cell(1, width(sam));
    newRow{1} = subj;

    % Input values for each variables
    for i = 2:width(sam)
        varName = sam.Properties.VariableNames{i};

        while true
            answer = inputdlg([varName ' (Angeben: Ganze Zahl)'], 'SAM', [1 70]);
            if isempty(answer)
                disp('Cancelled.');
                return;
            end
            % If the answer is empty, keep prompting
            if isempty(answer{1})
                disp('This field cannot be left empty. Please enter a valid response.');
            else
                % Otherwise, store the answer and break out of the loop
                newRow{i} = answer{1};
                break;
            end
        end
    end

    % Add new row to the table
    newRowTable = cell2table(newRow, 'VariableNames', sam.Properties.VariableNames);
    for i = 1:width(sam)
        if iscell(sam{:, i}) && ~iscell(newRowTable{:, i})
            newRowTable{:, i} = {newRowTable{:, i}};
        end
    end
    sam = [sam; newRowTable];
end

% Put the old data into an archive
tid_psam_clean_up_folder_TD(OUTPATH)

% Save new FAL data
writetable(fal, fal_filename);
disp([subj, ' FAL has been added and saved.']);

% Save new NASA-TLX data
writetable(nasatlx, nasatlx_filename);
disp([subj, ' NASA-TLX has been added and saved.']);

% Save new SAM data
writetable(sam, sam_filename);
disp([subj, ' SAM has been added and saved.']);


