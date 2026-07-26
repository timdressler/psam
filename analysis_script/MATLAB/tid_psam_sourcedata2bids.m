% participants.tsv DONE

% README DONE

% dataset_description DONE

% task-delayedArticulation DONE

% subj DONE, NEED UPDATE
% eeg DONE
% coordsystem.json DONE
% electrodes.tsv DONE
% channels.tsv DONE
% eeg.vhdr DONE
% eeg.json DONE
% events.json
% events.tsv /includes vocal data + stimuli data


% tid_psam_sourcedata2bids.m
%
% Master-Script that onverts sourcedata (EEG, questionnaire data) to BIDS-conform structure.
%
% Tim Dressler, 01.06.2026

%% Setup
clear
close all
clc
rng(123)
set(0,'DefaultTextInterpreter','none')
eeglab nogui
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
INPATH_SOURCEDATA = fullfile(MAINPATH, 'data', 'sourcedata');
INPATH_QUESTIONNAIRE_SRC = fullfile(INPATH_SOURCEDATA, 'questionnaire_data');
INPATH_TASK_SRC = fullfile(INPATH_SOURCEDATA, 'task_data');
INPATH_VOCAL_SRC = fullfile(MAINPATH, 'data', 'sourcedata', 'processed_data', 'beh_preprocessed_1');
INPATH_EEG_SRC = fullfile(MAINPATH, 'data', 'sourcedata', 'processed_data', 'markers_included');

OUTPATH = fullfile(MAINPATH, 'data');

FUNPATH = fullfile(MAINPATH, 'functions');
addpath(FUNPATH);

% Variables to edit
TASKNAME = 'delayedArticulation';
EXPECTED_SRATE = 1000;
POWERLINE_FREQ = 50;
EEG_NCHANS = 28;
EOG_CHANS = 2;
TRIGGER_NCHANS = 1;
HARDWARE_HP = 0.0159;
HARDWARE_LP = 250;

% Load MANUALLY_EXCLUDED_SUBJ
load(fullfile(INPATH_EEG_SRC, 'manually_excluded_subj.mat')); % loads MANUALLY_EXCLUDED_SUBJ

%% Create root-level files (e.g., participants.tsv, dataset_description_json.json, ...)
h_main_waitbar = waitbar(0, 'Initializing Main BIDS pipeline...', 'Name', 'PSAM BIDS Conversion');
waitbar(1/7, h_main_waitbar, 'Step 1/7: Creating root-level files...');

%% Create root-level files (e.g., participants.tsv, dataset_description_json.json, ...)
% --- Create dataset_description_json.json ---
dataset_description_json = struct();

% Create fields
dataset_description_json.Name = 'PSAM';
dataset_description_json.BIDSVersion = '1.11.1';
dataset_description_json.DatasetType = 'raw';
dataset_description_json.License = 'CC0';
dataset_description_json.Authors = {'Tim Dreßler', 'Andrea Hildebrandt', 'Stefan Debener'};
dataset_description_json.EthicsApprovals = {'Approved under Drs.EK/2025/027 by Kommission für Forschungsfolgenabschätzung und Ethik University of Oldenburg'};
dataset_description_json.Funding = {'tba'};

% Convert to json and write file
dataset_description_json = jsonencode(dataset_description_json, 'PrettyPrint', true);

fid = fopen(fullfile(OUTPATH, 'dataset_description.json'), 'w');
fprintf(fid, '%s', dataset_description_json);
fclose(fid);


disp('--- [OK] dataset_description.json created successfully! ---');

% --- Create participants.tsv ---
% Load and merge questionnaires
fal_data = readtable(fullfile(INPATH_QUESTIONNAIRE_SRC, 'fal_data.xlsx'));
nasatlx_data = readtable(fullfile(INPATH_QUESTIONNAIRE_SRC, 'nasatlx_data.xlsx'));
sam_data = readtable(fullfile(INPATH_QUESTIONNAIRE_SRC, 'sam_data.xlsx'));

participants = innerjoin(fal_data, nasatlx_data, 'Keys', 'subj');
participants = innerjoin(participants, sam_data, 'Keys', 'subj');

% Drop manually excluded subjects
participants = participants(~ismember(participants.subj, MANUALLY_EXCLUDED_SUBJ), :);

% Drop not needed variables
varsToDrop = {'var5_hearing_problems_detail', 'var10_alcohol_yesterday_detail', 'var11_alcohol_today_detail', 'var15_currently_neuro_treatment_detail', ...
    'var16_earlier_neuro_treatment_detail', 'var17_other_treatment_detail', 'var18_medication_detail', 'var19_drugs_detail', ...
    'var10_alcohol_yesterday', 'var11_alcohol_today', 'var12_smoking', 'var13_coffee_and_other', 'var14_last_eating', ...
    'var15_currently_neuro_treatment', 'var16_earlier_neuro_treatment', 'var17_other_treatment', 'var18_medication', ...
    'var19_drugs'};
participants = removevars(participants, varsToDrop);


% Rename existing columns
renameMap = [
    "subj",                            "participant_id" ;
    "var1_age",                        "age" ;
    "var2_handedness",                 "handedness" ;
    "var3_sex",                        "sex" ;
    "var4_education",                  "education" ;
    "var5_occupation",                 "occupation" ;
    "var5_hearing_problems",           "hearing_problems" ;
    "var7_ringing_ears",               "ringing_ears" ;
    "var8_sleep",                      "sleep_duration" ;
    "var9_sleep_assessment",           "sleep_assessment" ;
    "var1_mental_demand",              "mental_demand" ;
    "var2_physical_demand",            "physical_demand" ;
    "var3_performance",                "performance" ;
    "var4_effort",                     "effort" ;
    "var5_frustration",                "frustration" ;
    "var1_mood_break1",                "mood_break1" ;
    "var2_tiredness_break1",           "tiredness_break1" ;
    "var3_mood_break2",                "mood_break2" ;
    "var4_tiredness_break2",           "tiredness_break2" ;
    "var5_mood_break3",                "mood_break3" ;
    "var6_tiredness_break3",           "tiredness_break3" ;
    "var7_mood_break4",                "mood_break4" ;
    "var8_tiredness_break4",           "tiredness_break4" ;
    "var9_mood_break5",                "mood_break5" ;
    "var10_tiredness_break5",          "tiredness_break5" ;
    "var11_mood_break6",               "mood_break6" ;
    "var12_tiredness_break6",          "tiredness_break6" ;
    "var13_mood_break7",               "mood_break7" ;
    "var14_tiredness_break7",          "tiredness_break7" ;
    "var15_mood_break8",               "mood_break8" ;
    "var16_tiredness_break8",          "tiredness_break8"
    ];

% Extract old and new names
oldNames = renameMap(:, 1);
newNames = renameMap(:, 2);

% Apply to your table (assuming your table is named 'participants')
participants = renamevars(participants, oldNames, newNames);

% Add additional columns
% Species
participants.species = repmat("homo sapiens", height(participants), 1); % Everyone is a human here :)

% Reorder table
order = { ...
    'participant_id', ...
    'age', ...
    'sex', ...
    'handedness', ...
    'species', ...
    'education', ...
    'occupation', ...
    'hearing_problems', ...
    'ringing_ears', ...
    'sleep_duration', ...
    'sleep_assessment', ...
    'mental_demand', ...
    'physical_demand', ...
    'performance', ...
    'effort', ...
    'frustration', ...
    'mood_break1', 'tiredness_break1', ...
    'mood_break2', 'tiredness_break2', ...
    'mood_break3', 'tiredness_break3', ...
    'mood_break4', 'tiredness_break4', ...
    'mood_break5', 'tiredness_break5', ...
    'mood_break6', 'tiredness_break6', ...
    'mood_break7', 'tiredness_break7', ...
    'mood_break8', 'tiredness_break8' ...
    };
participants = participants(:, order);

% Write file
writetable(participants, fullfile(OUTPATH,'participants.tsv'), 'FileType', 'text', ...
    'Delimiter', '\t', ...
    'QuoteStrings', false);
disp('--- [OK] participants.tsv created successfully! ---');


% --- Create participants.json ---
participants_json = struct();

% Create fields
participants_json.participant_id.Description = 'participant identifier';

participants_json.species.Description = 'species of the participant';

participants_json.age.Description = 'age of the participant';
participants_json.age.Units = 'year';

participants_json.sex.Description = 'gender of the participant';
participants_json.sex.Levels.val_1 = 'male';
participants_json.sex.Levels.val_2 = 'female';
participants_json.sex.Levels.val_3 = 'diverse';

participants_json.handedness.Description = 'self-reported handedness';
participants_json.handedness.Levels.val_1 = 'right-handed';
participants_json.handedness.Levels.val_2 = 'left-handed';
participants_json.handedness.Levels.val_3 = 'two-handed';

participants_json.education.Description = 'Highest education achieved (German education system)';
participants_json.education.Levels.val_0 = 'no degree';
participants_json.education.Levels.val_1 = 'Hauptschule';
participants_json.education.Levels.val_2 = 'Mittlere-Reife';
participants_json.education.Levels.val_3 = 'Abitur';

participants_json.occupation.Description = 'Occupation of the participant';
participants_json.occupation.Levels.val_1 = 'student';
participants_json.occupation.Levels.val_2 = 'employed';
participants_json.occupation.Levels.val_3 = 'unemployed';

participants_json.hearing_problems.Description = 'Hearing problems';
participants_json.hearing_problems.Levels.val_1 = 'yes';
participants_json.hearing_problems.Levels.val_2 = 'no';

participants_json.ringing_ears.Description = 'Ringing in the ears';
participants_json.ringing_ears.Levels.val_1 = 'yes';
participants_json.ringing_ears.Levels.val_2 = 'no';

participants_json.sleep_duration.Description = 'Sleep of last night';
participants_json.sleep_duration.Units = 'hours';

participants_json.sleep_assessment.Description = 'Self-assessed sleep duration';
participants_json.sleep_assessment.Levels.val_1 = 'normal';
participants_json.sleep_assessment.Levels.val_2 = 'rather long';
participants_json.sleep_assessment.Levels.val_3 = 'way too short';

participants_json.mental_demand.Description = 'Cognitive demand required for processing information and decision-making (0 = low, 10 = high)';
participants_json.physical_demand.Description = 'Physical activity required (e.g., pulling, pushing, steering) (0 = low, 10 = high)';
participants_json.performance.Description = 'Perceived success and satisfaction with task performance (0 = good, 10 = bad)';
participants_json.effort.Description = 'Effort needed to meet task demands (0 = low, 10 = high)';
participants_json.frustration.Description = 'Feelings of stress, irritation, or frustration during the task (0 = low, 10 = high)';

for i = 1:8
    mood_var = sprintf('mood_break%d', i);
    tired_var = sprintf('tiredness_break%d', i);
    participants_json.(mood_var).Description = sprintf('Mood assessed at break %d (1 = very negative, 9 = very positive)', i);
    participants_json.(tired_var).Description = sprintf('Tiredness assessed at break %d (1 = very awake, 9 = very tired)', i);
end

% Convert to json and write file
participants_json = jsonencode(participants_json, 'PrettyPrint', true);

% Fix the MATLAB numeric key issue so the JSON outputs "0" and "1" instead of "val_0"
participants_json = strrep(participants_json, '"val_0"', '"0"');
participants_json = strrep(participants_json, '"val_1"', '"1"');
participants_json = strrep(participants_json, '"val_2"', '"2"');
participants_json = strrep(participants_json, '"val_3"', '"3"');
participants_json = strrep(participants_json, '"val_9"', '"9"');
participants_json = strrep(participants_json, '"val_10"', '"10"');

fid = fopen(fullfile(OUTPATH, 'participants.json'), 'w');
fprintf(fid, '%s', participants_json);
fclose(fid);


disp('--- [OK] participants.json created successfully! ---');

% --- Create README.md ---
readme_lines = {
    '# PSAM: Investigating the Specificity of Pre-Speech Auditory Modulation - From Global Gating to Selective Silence?'
    ''
    '## Overview'
    'The framework of internal forward models proposes that self-generated motor actions are accompanied by an efference copy (EC), which predicts sensory consequences and results in a corollary discharge (CD) that suppresses neural responses in corresponding sensory areas. In the speech-auditory domain, this is known as speaking-induced suppression (SIS). Recent evidence suggests a similar modulatory effect occurs prior to speech onset, termed Pre-Speech Auditory Modulation (PSAM).'
    ''
    'The PSAM dataset investigates whether this pre-speech modulation is specific to the expected vocal outcome (predictive mechanism) or reflects a general attenuation account (gating mechanism). It also explores whether this specificity depends on the phase of speech preparation and how these preparatory mechanisms change over time. The dataset contains comprehensive behavioral metrics, vocal responses, continuous 30-channel EEG data, and longitudinal psychometric questionnaire data collected during a tightly controlled active-passive within-subject speech preparation paradigm.'
    ''
    '## Participants & Screening'
    'Participants undergo an intake screening and initial state evaluation prior to the in-lab experimental testing phase:'
    '1. **Informed Consent:** Written informed consent is obtained from all participants prior to any experimental interaction.'
    '2. **Assessment of Initial State (FAL):** Participants complete the "Fragebogen zur Ausgangslage" (FAL) to record demographic, biographical, educational, and medical history data, as well as baseline sleep details.'
    '3. **Sample Characteristics:** Data were collected from 35 healthy adults at the University of Oldenburg. After excluding participants due to technical issues or insufficient trial quality, the final analyzed sample consists of n=28 participants (20 female, 8 male; mean age = 24.68 years, SD = 3.26).'
    ''
    '## Inclusion / Exclusion Criteria'
    '### Inclusion Criteria'
    '- Age between 18 and 40 years to minimize the risk of including individuals with undiagnosed age-related hearing impairments.'
    '- A minimum German language proficiency of level C1.'
    '- A baseline sleep duration of at least 5 hours the previous night.'
    ''
    '### Exclusion Criteria'
    '- Uncorrected visual impairments or hearing loss.'
    '- Current or historical psychiatric or neurological disorders (including both current and past stuttering).'
    '- Use of medications affecting the central nervous system or illegal drug use.'
    '- Consumption of alcohol on the day of testing.'
    ''
    '## Experimental Design'
    'To systematically track pre-speech modulation, the study utilizes a 2 (Task Condition) x 2 (Probe Type) x 2 (Probe Onset) within-subject factorial design:'
    '- **Task Condition:** Active (participant prepares to vocalize the syllable "/ga/") vs. Passive (participant observes the trial sequence without preparing or responding).'
    '- **Probe Type:** Unaltered (matching the participant''s individual median F0 vocalization) vs. Altered (pitch-shifted down by -4 semitones to introduce an EC mismatch).'
    '- **Probe Onset:** Early (-400 ms relative to the go-signal) vs. Late (-200 ms relative to the go-signal).'
    '- **Trial Configuration:** Participants complete a total of 960 trials across 8 blocks: 240 active no-probe trials, 240 passive no-probe trials, and 60 trials for each of the 8 distinct probe condition combinations.'
    ''
    '## Experimental Procedure'
    'The experiment takes place in a single laboratory session lasting approximately 210 minutes, structured into three consecutive phases:'
    ''
    '### Phase 1: Individual Stimulus Recording'
    '- Participants vocalize a short syllable ("/ga/") 21 times inside a sound-proof chamber following a structured trial sequence.'
    '- Audio recordings are trimmed and processed using Praat.'
    '- The unaltered probe closest to the participant''s median fundamental frequency (F0) is selected to represent a typical vocalization (assumed to match the EC). An altered duplicate is pitch-shifted down by -4 semitones using the Praat Vocal Toolkit to serve as a mismatch.'
    '- Both probes are trimmed to a length of 80 ms (with 10 ms fade-in/fade-out) and normalized to 70 dB SPL.'
    ''
    '### Phase 2: Main Task Block'
    '- Participants complete 8 blocks of 120 trials each (960 total trials) arranged in a pseudo-randomized sequence using a miniblock constraint to prevent pattern learning.'
    '- **Fixation Cross:** Presented for 0.5 to 1.5 seconds (jittered) at the start of each trial.'
    '- **Instruction Cue:** A circle appears containing either "/ga/" (Active trial) or "/xx/" (Passive trial).'
    '- **Delay Period:** Lasts 3 seconds, during which participants actively prepare their speech (in active trials) or remain calm (in passive trials).'
    '- **Auditory Probes:** Played via loudspeakers in 50% of the trials (either early at -400 ms or late at -200 ms). The remaining 50% serve as no-probe baseline control trials.'
    '- **Go-Signal:** The fill color of the circle turns green for 1.5 seconds, prompting the immediate vocalization of "/ga/" in active trials.'
    ''
    '### Phase 3: Follow-Up'
    '- **Self-Assessment Manikin (SAM):** Administered during the 3-minute breaks between every block to longitudinally track valence (mood) and arousal (tiredness).'
    '- **NASA Task Load Index (NASA-TLX):** Completed at the very end of the session to capture subjective metrics of workload across 5 dimensions (mental demand, physical demand, performance, effort, and frustration) on a 0 to 10 scale.'
    ''
    '## Questionnaire & Behavioral Tracking Timeline'
    'Psychometric and state measurements are taken repeatedly across the session, strictly aligned with task events:'
    '- **Pre-Task Intake:** Assessment of Initial State (FAL) capturing biographical profiles and baseline eligibility metrics.'
    '- **Block Breaks (1 to 8):** Custom SAM administered during every 3-minute intermission between blocks to monitor emotional fluctuations and fatigue over time.'
    '- **Post-Task Completion:** Modified NASA-TLX evaluating subjective workload dimensions and perceived strain after all 960 trials are completed.'
    ''
    '## Stimulus Material'
    '- **Target Utterance:** Restricted exclusively to the basic speech sound sound syllable "/ga/".'
    '- **Unaltered Probes:** Individually tailored natural tokens selected at the participant''s median F0 to embody the sensory prediction target.'
    '- **Altered Probes:** Acoustic variants pitch-shifted down by -4 semitones to introduce a controlled sensorimotor prediction error.'
    '- All probes are standard-trimmed to an 80 ms duration, fade-smoothed, and normalized to 70 dB SPL.'
    ''
    '## Data Acquisition & Hardware'
    '- **EEG Data:** Continuous 30-channel cap setup (Easycap) with Ag/AgCl electrodes arranged in a custom equidistant layout with a nose-tip reference. Electrodes E01 to E28 are used for scalp EEG (where E01 corresponds to Cz), and E29-E30 capture EOG. Data are acquired via BrainAmp amplifiers, sampled at 1000 Hz, and band-pass filtered between 0.0159 and 250 Hz with impedances below 20 kΩ.'
    '- **Vocal Responses:** Captured using a Sennheiser MD 43 dynamic microphone positioned ~5 cm from the mouth, routed to a Windows PC via a Focusrite Scarlett 2i2 (3rd Gen) audio interface sampled at 44.1 kHz.'
    '- **Stimulus Presentation:** Managed synchronously via PsychoPy (Version 2024.2.4). Visual assets displayed on a 24-inch screen (1920x1080, 60 Hz) at 110 cm distance. Audio outputs driven by an internal RME HDSP 9632 card, external RME ADI 8 DS MK III D/A converter, Tucker-Davis PA5 attenuators, and HB7 headphone buffers, distributed via Sirocco S30 loudspeakers at ear level (115 cm distance, 45° angle).'
    ''
    '## Dataset Structure'
    'This dataset is formatted according to the Brain Imaging Data Structure (BIDS) standard.'
    ''
    '### Data Protection Notice'
    'Due to strict legal regulations regarding biometric human data protection, the raw vocal audio data (both the initial stimuli recordings and the vocal responses captured during the execution task) **cannot be included** in this public repository.'
    ''
    '### Structure Overview'
    '- `/phenotype`: Contains questionnaire data (`.tsv`) and structural schemas (`.json`) for block-by-block SAM entries, and post-task NASA-TLX indices.'
    '- `/sub-<ID>`: Subject-specific folders containing:'
    '    - `/eeg`: Continuous EEG records (`.set` and `.fdt`), channel layout files, and `_events.tsv` logs.'
    '    - **Note on Speech Metrics:** While raw audio is omitted for data privacy, the fully extracted fundamental frequency (F0) values and vocal onset times for every active trial are directly embedded within the `_events.tsv` table for immediate integration with the neural time-series.'
    };

% Join the lines with newline characters
readme_text = strjoin(readme_lines, '\n');

% Write to file
fid = fopen(fullfile(OUTPATH, 'README.md'), 'w');
fprintf(fid, '%s', readme_text);
fclose(fid);


disp('--- [OK] README.md created successfully! ---');


%% Setup subject-wise folder structure
waitbar(2/7, h_main_waitbar, 'Step 2/7: Creating subject folder structures...');
% Get subjects in sourcedata/task
dircont_subj = dir(fullfile(INPATH_EEG_SRC, 'sub-*.set'));

for subj = 1:length(dircont_subj)
    dircont_subj(subj).name = extractBefore(dircont_subj(subj).name, '_');
end

% Loop over subjects and create folders
for subj = 1:length(dircont_subj)
    % Get subject ID
    subjID = dircont_subj(subj).name;

    % Create subject folder
    path_subject = fullfile(OUTPATH, subjID);
    tid_psam_check_folder_TD(path_subject);
    % Create subfolder for EEG
    psylink_check_folder_TD(fullfile(path_subject, 'eeg'));
end

%% Get and rename EEG files
waitbar(3/7, h_main_waitbar, 'Step 3/7: Copying and renaming EEG files...');

for subj = 1:length(dircont_subj)
    % Get subject ID (e.g., 'sub-95')
    subjID = dircont_subj(subj).name;
    % Define source and destination paths
    src_eeg_dir = INPATH_EEG_SRC;
    dest_eeg_dir = fullfile(OUTPATH, subjID, 'eeg');
    % Base name for BIDS conformity (e.g., 'sub-95_task-delayedArticulation_eeg')
    bids_basename = sprintf('%s_task-%s_eeg', subjID, TASKNAME);
    % File extensions
    extensions = {'.set', '.fdt'};
    for ext_idx = 1:length(extensions)
        ext = extensions{ext_idx};
        % Find the source file with this extension in the subject's sourcedata
        src_file_struct = dir(fullfile(src_eeg_dir, ['*' ext]));
        if ~isempty(src_file_struct)
            % Define full paths
            src_file = fullfile(src_file_struct(1).folder, src_file_struct(1).name);
            dest_file = fullfile(dest_eeg_dir, [bids_basename ext]);
            % Copy the file to the BIDS directory
            copyfile(src_file, dest_file);
            % Update internal references for .vhdr and .vmrk
            if strcmp(ext, '.vhdr') || strcmp(ext, '.vmrk')
                % Read the file contents
                fid = fopen(dest_file, 'r');
                file_content = fread(fid, '*char')';
                fclose(fid);
                % Extract the old base name (removing the extension)
                old_basename = src_file_struct(1).name;
                [~, old_basename_noext, ~] = fileparts(old_basename);
                % Replace the old filename references with the new BIDS base name
                file_content = strrep(file_content, old_basename_noext, bids_basename);
                % Write the updated content back to the destination file
                fid = fopen(dest_file, 'w');
                fwrite(fid, file_content, '*char');
                fclose(fid);
            end
        else
            warning('PSAM:MissingData', 'Could not find %s file for %s in %s', ext, subjID, src_eeg_dir);
        end
    end
    fprintf('--- [OK] Copied and renamed EEG files for %s ---\n', subjID);
end

% Load EEG data and rename event markers and remove not needed ones
for subj = 1:length(dircont_subj)
    % Clean up memory for the next subject
    ALLEEG = [];
    EEG = [];
    CURRENTSET = 0;    
    % Get subject ID (e.g., 'sub-95')
    subjID = dircont_subj(subj).name;    
    % Load EEG file
    subj_eeg_dir = fullfile(OUTPATH, subjID, 'eeg');
    bids_basename = sprintf('%s_task-%s_eeg', subjID, TASKNAME);
    EEG = pop_loadset('filename', [bids_basename '.set'], 'filepath', subj_eeg_dir);    
    % 1. Setup tracking variables
    current_task = ''; % Will hold 'act' or 'pas' based on instruction cues
    events_to_keep = false(1, length(EEG.event)); % Logical index to filter events    
    % 2. Loop through events to rename and flag the ones to keep
    for e = 1:length(EEG.event)        
        % Track the task instruction so we can identify the upcoming Go-Signal
        if strcmp(EEG.event(e).type, 'S 21')
            current_task = 'act';
        elseif strcmp(EEG.event(e).type, 'S 22')
            current_task = 'pas';
        end        
        % Evaluate and rename events
        switch EEG.event(e).type
            % --- AUDIO PROBES (Keep and Rename) ---
            case 'S 931'
                EEG.event(e).type = 'act_early_unalt';
                events_to_keep(e) = true;
            case 'S 932'
                EEG.event(e).type = 'act_early_alt';
                events_to_keep(e) = true;
            case 'S 933'
                EEG.event(e).type = 'act_late_unalt';
                events_to_keep(e) = true;
            case 'S 934'
                EEG.event(e).type = 'act_late_alt';
                events_to_keep(e) = true;
            case 'S 941'
                EEG.event(e).type = 'pas_early_unalt';
                events_to_keep(e) = true;
            case 'S 942'
                EEG.event(e).type = 'pas_early_alt';
                events_to_keep(e) = true;
            case 'S 943'
                EEG.event(e).type = 'pas_late_unalt';
                events_to_keep(e) = true;
            case 'S 944'
                EEG.event(e).type = 'pas_late_alt';
                events_to_keep(e) = true;                
            % --- CONTROL PROBES ---
            case {'con_act_early', 'con_act_late', 'con_pas_early', 'con_pas_late'}
                events_to_keep(e) = true;                
            % --- GO SIGNALS (Keep and Rename based on look-back) ---
            case 'S  5'
                if strcmp(current_task, 'act')
                    EEG.event(e).type = 'go_signal_act';
                elseif strcmp(current_task, 'pas')
                    EEG.event(e).type = 'go_signal_pas';
                end
                events_to_keep(e) = true;                
            % Any other marker (including S 21, S 22, boundary events, etc.) falls through the switch statement and remains events_to_keep = false
        end
    end    
    % 3. Remove all events that were not flagged to be kept
    EEG.event = EEG.event(events_to_keep);    
    % 4. Check dataset 
    EEG = eeg_checkset(EEG);
    pop_saveset(EEG, 'filename', [bids_basename '.set'], 'filepath', subj_eeg_dir);
    
    fprintf('--- [OK] Cleaned and renamed EEG triggers for %s ---\n', subjID);
end




%% EEG Meta data
waitbar(4/7, h_main_waitbar, 'Step 4/7: Creating EEG metadata...');

for subj = 1:length(dircont_subj)
    % Get subject ID (e.g., 'sub-95')
    subjID = dircont_subj(subj).name;
    % --- Create _electrodes.tsv ---
    % Load channel locations
    elp_file = fullfile(MAINPATH, 'config', 'elec_96ch_adapted.elp');
    chanlocs = readlocs(elp_file);
    % Extract labels and coordinates into arrays
    names = {chanlocs.labels}';
    % Check if readlocs successfully generated Cartesian coordinates
    if isfield(chanlocs, 'X') && ~isempty(chanlocs(1).X)
        x = [chanlocs.X]';
        y = [chanlocs.Y]';
        z = [chanlocs.Z]';
    else
        % Fallback to NaNs if coordinate conversion fails
        x = nan(length(names), 1);
        y = nan(length(names), 1);
        z = nan(length(names), 1);
    end
    % Create the initial table
    electrodes_tsv = table(string(names), x, y, z, 'VariableNames', {'name', 'x', 'y', 'z'});
    % Convert NaN to n/a to conform with BIDS standard
    colsToConvert = {'x', 'y', 'z'};
    for i = 1:length(colsToConvert)
        colName = colsToConvert{i};
        data = electrodes_tsv.(colName);
        % Initialize a cell array of strings
        strCol = cell(size(data));
        for j = 1:length(data)
            if isnan(data(j))
                strCol{j} = 'n/a';
            else
                strCol{j} = num2str(data(j), '%.10g');
            end
        end
        % Reassign back to the table as a string array
        electrodes_tsv.(colName) = string(strCol);
    end
    % Write file to the BIDS subject's EEG folder
    dest_eeg_dir = fullfile(OUTPATH, subjID, 'eeg');
    file_name = sprintf('%s_electrodes.tsv', subjID);
    writetable(electrodes_tsv, fullfile(dest_eeg_dir, file_name), 'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);
    fprintf('--- [OK] %s created successfully! ---\n', file_name);
end

for subj = 1:length(dircont_subj)
    % Get subject ID (e.g., 'sub-95')
    subjID = dircont_subj(subj).name;
    % --- Create _coordsystem.json ---
    coordsystem_json = struct();
    % Create fields
    coordsystem_json.EEGCoordinateSystem = 'Other';
    coordsystem_json.EEGCoordinateUnits = 'mm';
    coordsystem_json.EEGCoordinateSystemDescription = 'Lab-specific custom template coordinates applied to all subjects via a custom .elp file. No subject-specific 3D digitization was performed.';
    % Convert to json
    coordsystem_json_text = jsonencode(coordsystem_json, 'PrettyPrint', true);
    % Define destination path and filename (session-less)
    dest_eeg_dir = fullfile(OUTPATH, subjID, 'eeg');
    file_name = sprintf('%s_coordsystem.json', subjID);
    % Write file
    fid = fopen(fullfile(dest_eeg_dir, file_name), 'w');
    fprintf(fid, '%s', coordsystem_json_text);
    fclose(fid);
    fprintf('--- [OK] %s created successfully! ---\n', file_name);
end

% --- Create _channels.tsv ---
for subj = 1:length(dircont_subj)
    % Get subject ID (e.g., 'sub-95')
    subjID = dircont_subj(subj).name;
    % Define BIDS EEG directory and basename
    dest_eeg_dir = fullfile(OUTPATH, subjID, 'eeg');
    bids_basename = sprintf('%s_task-%s_eeg', subjID, TASKNAME);
    set_file = [bids_basename '.set'];
    % Load the EEG data from the BIDS location
    fprintf('Loading %s for channels.tsv extraction...\n', set_file);
    try
        EEG = pop_loadset('filepath', dest_eeg_dir, 'filename', set_file);
    catch ME
        warning('PSAM:LoadError', 'Could not load EEG data for %s: %s', subjID, ME.message);
        continue; % Skip to the next subject if loading fails
    end
    % Convert EEG.chanlocs to a table
    channels_tsv = struct2table(EEG.chanlocs, 'AsArray', true);
    % Ensure 'type' column exists
    if ~ismember('type', channels_tsv.Properties.VariableNames)
        channels_tsv.type = repmat({''}, height(channels_tsv), 1);
    end
    % Filter to only keep 'labels' and 'type'
    channels_tsv = channels_tsv(:, {'labels', 'type'});
    % Rename 'labels' to 'name'
    channels_tsv.Properties.VariableNames{'labels'} = 'name';
    % Convert cell arrays to string arrays for easier text manipulation
    channels_tsv.name = string(channels_tsv.name);
    channels_tsv.type = string(channels_tsv.type);
    % Default everything to EEG first
    channels_tsv.type(:) = "EEG";
    channels_tsv.type(channels_tsv.name == "E29" | channels_tsv.name == "E30") = "EOG";
    channels_tsv.type(channels_tsv.name == "M") = "TRIG";
    % Add the 'units' column
    channels_tsv.units = repmat("uV", height(channels_tsv), 1);
    % Add the 'description' column
    channels_tsv.description = repmat("n/a", height(channels_tsv), 1);
    % Add specific descriptions if needed
    channels_tsv.description(channels_tsv.name == "E01") = "Cz equivalent";
    channels_tsv.description(channels_tsv.name == "E29") = "VEOG under left eye";
    channels_tsv.description(channels_tsv.name == "E30") = "VEOG under right eye";
    channels_tsv.description(channels_tsv.name == "M") = "Trigger channel for audio stimuli";
    % Define file name and path
    file_name = sprintf('%s_task-%s_channels.tsv', subjID, TASKNAME);
    % Write file
    writetable(channels_tsv, fullfile(dest_eeg_dir, file_name), 'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);
    fprintf('--- [OK] %s created successfully! ---\n', file_name);
end

% --- Create _eeg.json ---
for subj = 1:length(dircont_subj)
    % Get subject ID (e.g., 'sub-95')
    subjID = dircont_subj(subj).name;
    dest_eeg_dir = fullfile(OUTPATH, subjID, 'eeg');
    eeg_json = struct();
    % Core fields
    eeg_json.EEGReference = 'nose-tip';
    eeg_json.SamplingFrequency = EXPECTED_SRATE;
    eeg_json.PowerLineFrequency = POWERLINE_FREQ;
    % Hardware and Software filters
    eeg_json.HardwareFilters = struct();
    eeg_json.HardwareFilters.HighpassRC = struct('HalfAmplitudeCutoffHz', HARDWARE_HP);
    eeg_json.HardwareFilters.Lowpass = struct('HalfAmplitudeCutoffHz', HARDWARE_LP);
    eeg_json.SoftwareFilters = 'n/a';    
    eeg_json.TaskName = TASKNAME;
    eeg_json.TaskDescription = 'Delayed articulation paradigm (Active vs. Passive) utilizing unaltered and altered (-4 semitones) auditory probes presented either early or late during the preparatory delay period to investigate Pre-Speech Auditory Modulation (PSAM).';
    % Hardware information
    eeg_json.Manufacturer = 'Brain Products';
    eeg_json.ManufacturersModelName = 'BrainAmp';
    eeg_json.CapManufacturer = 'Easycap';
    eeg_json.EEGPlacementScheme = 'custom equidistant';
    % Institutional information
    eeg_json.InstitutionName = 'University of Oldenburg';
    eeg_json.InstitutionalDepartmentName = 'Department of Psychology, Neuropsychology Lab';
    % Channel counts (using variables defined in your Setup section)
    eeg_json.EEGChannelCount = EEG_NCHANS;
    eeg_json.EOGChannelCount = EOG_CHANS;
    eeg_json.TriggerChannelCount = TRIGGER_NCHANS;
    eeg_json.RecordingType = 'continuous';
    % Further information
    eeg_json.ECGChannelCount = 0;
    eeg_json.EMGChannelCount = 0;
    eeg_json.MISCChannelCount = 0;    
    eeg_json.CapManufacturersModelName = 'EasyCap';
    eeg_json.InstitutionAddress = 'n/a';
    eeg_json.DeviceSerialNumber = 'n/a';
    eeg_json.SoftwareVersions = 'n/a';
    eeg_json.Instructions = 'n/a';
    eeg_json.EEGGround = 'n/a';
    eeg_json.SubjectArtefactDescription = 'n/a';
    % Convert to JSON format
    eeg_json_text = jsonencode(eeg_json, 'PrettyPrint', true);
    % Define file name and path
    file_name = sprintf('%s_task-%s_eeg.json', subjID, TASKNAME);
    % Write file
    fid = fopen(fullfile(dest_eeg_dir, file_name), 'w');
    fprintf(fid, '%s', eeg_json_text);
    fclose(fid);
    fprintf('--- [OK] %s created successfully! ---\n', file_name);
end


%% Create _events.tsv and _events.json files (Merge EEG + Beh + Vocal)
waitbar(5/7, h_main_waitbar, 'Step 5/7: Merging behavioral, vocal, and EEG events...');

for subj = 1:length(dircont_subj)
    % Get subject ID (e.g., 'sub-95')
    subjID = dircont_subj(subj).name;
    
    % Define BIDS EEG directory and basename
    dest_eeg_dir = fullfile(OUTPATH, subjID, 'eeg');
    bids_basename = sprintf('%s_task-%s_eeg', subjID, TASKNAME);
    set_file = [bids_basename '.set'];
    
    % Load the EEG data from the BIDS location
    fprintf('Loading %s for events.tsv extraction...\n', set_file);
    try
        EEG = pop_loadset('filepath', dest_eeg_dir, 'filename', set_file);
    catch ME
        warning('PSAM:LoadError', 'Could not load EEG data for %s: %s', subjID, ME.message);
        continue; % Skip to the next subject if loading fails
    end
    
    % 1. Clean up events_table keeping only latency and type
    events_table = struct2table(EEG.event);
    events_table = events_table(:, {'latency', 'type'});
    
    % Recalculate latency to onset (in seconds) with high-precision formatting
    bids_onset = cellstr(compose('%.11f', (events_table.latency - 1) / EEG.srate));
    % Duration is 0 for all instantaneous EEG events
    bids_duration = cellstr(compose('%.11f', zeros(height(events_table), 1)));
    
    % 2. Load vocal data 
    vocal_table = readtable(fullfile(INPATH_VOCAL_SRC, [subjID '_f0_rt_table.csv']), 'FileType','text', 'Delimiter', ',');
    vocal_table.Properties.VariableNames{'filename_tab'} = 'recording_file'; % Rename for joining
    vocal_table = standardizeMissing(vocal_table, 9999); % Replace 9999 with NaN
    
    % 3. Load log data 
    subj_log_filename = dir(fullfile(INPATH_TASK_SRC, subjID, 'beh', 's*.csv'));
    if numel(subj_log_filename) == 1
        subj_log = readtable(fullfile(subj_log_filename.folder, subj_log_filename.name));
    else
        error('Incorrect number log files for %s!', subjID);
    end
    
    % Clean and format log data (remove instructions, rename conditions)
    subj_log_clean = subj_log(~isnan(subj_log.mic_started), :);
    subj_log_clean.probe_type(strcmp(subj_log_clean.probe_type, 'Normal')) = {'Unaltered'};
    subj_log_clean.probe_type(strcmp(subj_log_clean.probe_type, 'Pitch')) = {'Altered'};
    subj_log_clean.task(strcmp(subj_log_clean.task, '/ga/')) = {'act'};
    subj_log_clean.task(strcmp(subj_log_clean.task, '/xx/')) = {'pas'};
    
    % Extract recording file name from the mic_clip path to match Praat
    subj_log_clean.recording_file = cellfun(@(x) strrep(regexp(x, 'recording_mic_.*(?=\.wav)', 'match', 'once'), '.', '_'), subj_log_clean.mic_clip, 'UniformOutput', false);
    
    % 4. Merge behavioural logfile and vocal data
    merged_beh = innerjoin(subj_log_clean, vocal_table, 'Keys', 'recording_file');
    
    % Recalculate true vocal onset relative to go_signal
    merged_beh.go_stim_started_trial_start = merged_beh.go_stim_started - merged_beh.trial_started;
    merged_beh.recording_vot = merged_beh.rt_tab - (merged_beh.go_stim_started_trial_start - merged_beh.mic_started);
    
    % Add correct_resp column (1 if correct, 0 if incorrect)
    pas_correct = strcmp(merged_beh.task, 'pas') & merged_beh.vocal_response_tab == 0;
    act_correct = strcmp(merged_beh.task, 'act') & merged_beh.vocal_response_tab == 1;
    merged_beh.correct_resp = pas_correct | act_correct;
    
    % 5. Add subject-wise probe properties
    subj_probe_file = fullfile(MAINPATH, 'data', 'sourcedata', 'task_data', 'stimuli', subjID, [subjID '_probe_properties.xlsx']);
    subj_probe_properties = readtable(subj_probe_file);
    sub_unaltered_f0 = num2str(subj_probe_properties.f0_tab_normal(1), '%.2f');
    sub_altered_f0 = num2str(subj_probe_properties.f0_tab_pitched(1), '%.2f');
    
    % 6. Build final BIDS events table
    num_events = height(events_table);
    
    % Initialize columns
    bids_trial = cell(num_events, 1);
    bids_event_type = cell(num_events, 1);
    bids_marker_label = events_table.type;
    
    bids_instruction = cell(num_events, 1); 
    
    bids_vocal_f0 = cell(num_events, 1);
    bids_vocal_rt = cell(num_events, 1);
    bids_vocal_resp = cell(num_events, 1);
    bids_correct_resp = cell(num_events, 1);
    
    bids_probe_onset = cell(num_events, 1);
    bids_probe_type = cell(num_events, 1);
    bids_probe = cell(num_events, 1);
    
    bids_probe_unalt_f0 = repmat({sub_unaltered_f0}, num_events, 1);
    bids_probe_alt_f0 = repmat({sub_altered_f0}, num_events, 1);

    % Loop through chronologically to map trials to triggers
    trial_idx = 1; 
    
    for e = 1:num_events
        marker = events_table.type{e};
        bids_trial{e} = trial_idx;
        
        % Identify Event Type
        if startsWith(marker, 'go_signal')
            bids_event_type{e} = 'go_signal';
        elseif startsWith(marker, 'con_')
            bids_event_type{e} = 'control';
        else
            bids_event_type{e} = 'audio';
        end
        
        % Extract trial-level context (shared for both lines within the trial)
        cur_trial = merged_beh(trial_idx, :);
        
        % Set instruction for both rows in the trial
        if strcmp(cur_trial.task, 'act')
            bids_instruction{e} = 'active';
        else
            bids_instruction{e} = 'passive';
        end
        
        bids_probe_onset{e} = cur_trial.probe_onset_cat{1};
        bids_probe_type{e} = cur_trial.probe_type{1};
        bids_probe{e} = cur_trial.probe{1};
        
        % Clean up values for BIDS standard "n/a"
        if strcmpi(bids_probe_onset{e}, 'None'), bids_probe_onset{e} = 'n/a'; end
        if strcmpi(bids_probe_type{e}, 'None'), bids_probe_type{e} = 'n/a'; end
        
        % Force probe to strict lowercase to match the events.json levels exactly
        if strcmpi(bids_probe{e}, 'None') || strcmpi(bids_probe{e}, 'no')
            bids_probe{e} = 'no'; 
        elseif strcmpi(bids_probe{e}, 'yes')
            bids_probe{e} = 'yes'; 
        end
        
        % Boolean logic for correctness
        if cur_trial.correct_resp
            bids_correct_resp{e} = 'yes';
        else
            bids_correct_resp{e} = 'no';
        end
        
        % Vocal metrics logic
        if strcmp(bids_event_type{e}, 'go_signal')
            if cur_trial.vocal_response_tab == 1
                bids_vocal_resp{e} = 'yes';
                bids_vocal_f0{e} = num2str(cur_trial.f0_tab, '%.2f');
                bids_vocal_rt{e} = num2str(cur_trial.recording_vot, '%.4f');
            else
                bids_vocal_resp{e} = 'no';
                bids_vocal_f0{e} = 'n/a';
                bids_vocal_rt{e} = 'n/a';
            end
            
            % Go-signal is the final event of the trial; move to next trial
            trial_idx = min(trial_idx + 1, height(merged_beh));
        else
            % For audio or control triggers, vocal metrics are purely n/a
            bids_vocal_resp{e} = 'n/a';
            bids_vocal_f0{e} = 'n/a';
            bids_vocal_rt{e} = 'n/a';
        end
    end
    
    % Assemble Final Table
    final_events_tsv = table(bids_onset, bids_duration, bids_trial, bids_event_type, ...
        bids_marker_label, bids_instruction, bids_vocal_f0, bids_vocal_rt, bids_vocal_resp, bids_correct_resp, ...
        bids_probe_onset, bids_probe_type, bids_probe, bids_probe_unalt_f0, bids_probe_alt_f0, ...
        'VariableNames', {'onset', 'duration', 'trial', 'event_type', 'marker_label', ...
        'task_condition', 'vocal_f0', 'vocal_rt', 'vocal_resp', 'correct_resp', 'probe_onset', 'probe_type', ...
        'probe', 'subj_probe_unaltered_f0', 'subj_probe_altered_f0'});
    
    % Save _events.tsv
    tsv_filename = sprintf('%s_task-%s_events.tsv', subjID, TASKNAME);
    writetable(final_events_tsv, fullfile(dest_eeg_dir, tsv_filename), ...
        'FileType', 'text', 'Delimiter', '\t', 'QuoteStrings', false);
    
    fprintf('--- [OK] %s created successfully! ---\n', tsv_filename);



    % --- Create _events.json ---
    events_json = struct();
    
    % Create fields
    events_json.onset.Description = 'Onset time of the event relative to the start of the EEG recording';
    events_json.onset.Units = 's';
    
    events_json.duration.Description = 'Duration of the event';
    events_json.duration.Units = 's';
    
    events_json.trial.Description = 'Chronological trial number within the experiment';
    
    events_json.event_type.Description = 'General category of the experimental event';
    events_json.event_type.Levels.audio = 'Auditory probe presentation';
    events_json.event_type.Levels.go_signal = 'Visual go-signal prompting the vocal response (or silence)';
    events_json.event_type.Levels.control = 'Control/sham auditory marker (no sound presented, timepoint when a marker would have been presented)';
    
    events_json.marker_label.Description = 'Specific condition marker derived from EEG hardware triggers';
    
    events_json.task_condition.Description = 'Task condition (instruction) for the current trial';
    events_json.task_condition.Levels.active = 'Active trial (participant prepares to vocalize /ga/)';
    events_json.task_condition.Levels.passive = 'Passive trial (participant observes without vocalizing)';
    
    events_json.vocal_f0.Description = 'Fundamental frequency (F0) of the vocal response extracted via Praat. Only applicalble if a vocal response was made.';
    events_json.vocal_f0.Units = 'Hz';
    
    events_json.vocal_rt.Description = 'Vocal reaction time calculated relative to the visual go-signal. Only applicalble if a vocal response was made.';
    events_json.vocal_rt.Units = 's';
    
    events_json.vocal_resp.Description = 'Indicates if a vocal response was physically detected during the trial';
    events_json.vocal_resp.Levels.yes = 'Vocal response detected';
    events_json.vocal_resp.Levels.no = 'No vocal response detected';
    
    events_json.correct_resp.Description = 'Whether the participant correctly followed the task instructions';
    events_json.correct_resp.Levels.yes = 'Correct (vocalized during active, or stayed silent during passive)';
    events_json.correct_resp.Levels.no = 'Incorrect (vocalized during passive, or stayed silent during active)';
    
    events_json.probe_onset.Description = 'Timing of the auditory probe relative to the go-signal';
    events_json.probe_onset.Levels.Early = '-400 ms relative to go-signal';
    events_json.probe_onset.Levels.Late = '-200 ms relative to go-signal';
    
    events_json.probe_type.Description = 'Acoustic manipulation of the auditory probe';
    events_json.probe_type.Levels.Unaltered = 'Probe unaltered';
    events_json.probe_type.Levels.Altered = 'Probe pitch-shifted down by -4 semitones';
    
    events_json.probe.Description = 'Indicates whether an auditory probe was actually presented during the trial';
    events_json.probe.Levels.yes = 'Probe presented';
    events_json.probe.Levels.no = 'No probe presented (control/baseline trial)';
    
    events_json.subj_probe_unaltered_f0.Description = 'Baseline F0 of the subject-specific unaltered auditory probe';
    events_json.subj_probe_unaltered_f0.Units = 'Hz';
    
    events_json.subj_probe_altered_f0.Description = 'Pitch-shifted F0 of the subject-specific altered auditory probe (-4 semitones)';
    events_json.subj_probe_altered_f0.Units = 'Hz';

    events_json.StimulusPresentation.OperatingSystem = 'Windows';
    events_json.StimulusPresentation.SoftwareName = 'PsychoPy';
    events_json.StimulusPresentation.SoftwareVersion = '2024.2.4';
    
    % Convert to json
    events_json_text = jsonencode(events_json, 'PrettyPrint', true);
    
    % Write file
    json_filename = sprintf('%s_task-%s_events.json', subjID, TASKNAME);
    fid = fopen(fullfile(dest_eeg_dir, json_filename), 'w');
    fprintf(fid, '%s', events_json_text);
    fclose(fid);
    
    fprintf('--- [OK] %s created successfully! ---\n', json_filename);

end



close(h_main_waitbar)
disp('--- FULL BIDS CONVERSION COMPLETE! ---');