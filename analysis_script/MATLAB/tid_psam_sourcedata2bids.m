% participants.tsv 

% README

% dataset_description

% task-delayedArticulation

% subj
    % eeg
        % coordsystem.json
        % electrodes.tsv
        % channels.tsv
        % eeg.set
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
OUTPATH = fullfile(MAINPATH, 'data');

FUNPATH = fullfile(MAINPATH, 'functions');
addpath(FUNPATH);

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

participants_json.mental_demand.Description = 'Cognitive demand required for processing information and decision-making';
participants_json.mental_demand.Levels.val_0 = 'low';
participants_json.mental_demand.Levels.val_10 = 'high';

participants_json.physical_demand.Description = 'Physical activity required (e.g., pulling, pushing, steering)';
participants_json.physical_demand.Levels.val_0 = 'low';
participants_json.physical_demand.Levels.val_10 = 'high';

participants_json.performance.Description = 'Perceived success and satisfaction with task performance';
participants_json.performance.Levels.val_0 = 'good';
participants_json.performance.Levels.val_10 = 'bad';

participants_json.effort.Description = 'Effort needed to meet task demands';
participants_json.effort.Levels.val_0 = 'low';
participants_json.effort.Levels.val_10 = 'high';

participants_json.frustration.Description = 'Feelings of stress, irritation, or frustration during the task';
participants_json.frustration.Levels.val_0 = 'low';
participants_json.frustration.Levels.val_10 = 'high';

for i = 1:8
    mood_var = sprintf('mood_break%d', i);
    tired_var = sprintf('tiredness_break%d', i);
    
    participants_json.(mood_var).Description = sprintf('Mood assessed at break %d', i);
    participants_json.(mood_var).Levels.val_1 = 'very negative';
    participants_json.(mood_var).Levels.val_9 = 'very positive';
    
    participants_json.(tired_var).Description = sprintf('Tiredness assessed at break %d', i);
    participants_json.(tired_var).Levels.val_1 = 'very awake';
    participants_json.(tired_var).Levels.val_9 = 'very tired';
end

% Convert to json and write file
participants_json = jsonencode(participants_json, 'PrettyPrint', true);

% Fix the MATLAB numeric key issue so the JSON outputs "0" and "1" instead of "val_0"
participants_json = strrep(participants_json, '"val_0"', '"0"');
participants_json = strrep(participants_json, '"val_1"', '"1"');
participants_json = strrep(participants_json, '"val_2"', '"2"');

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
dircont_subj = dir(fullfile(PATH_TASK_SRC, 'sub-*'));

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


disp('--- FULL BIDS CONVERSION COMPLETE! ---');