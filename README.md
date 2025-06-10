# WORK IN PROGRESS!

# Investigating the Specificity of Pre-Speech Auditory Modulation - From Global Gating to Selective Silence?

## In Brief

## This Repository

This repository includes all code related to the 'Investigating the Specificity of Pre-Speech Auditory Modulation - From Global Gating to Selective Silence?' project conducted by Tim Dressler at the University of Oldenburg. The project was supervised by Prof. Dr. Stefan Debener and Prof. Dr. Andrea Hildebrandt.

The repository includes analysis code as well as as code used for running the experiment. No data is included, but can by requested. Next to code, the repository also includes documents such as the participant information. Furthermore, a short theoretical description of the project and its result can be found in this README.

### Structure

```text
psam/
│
├── analysis_script/
│   │
│   ├── MATLAB/
│   │   ├── tid_psam_MATLAB_TEST.m
│   │   ├── tid_psam_beh_analysis.m
│   │   ├── tid_psam_beh_preprocessing_2.m
│   │   ├── tid_psam_erp_analysis.m
│   │   ├── tid_psam_erp_pilot.m
│   │   ├── tid_psam_erp_preprocessing.m
│   │   ├── tid_psam_exclude_trials.m
│   │   ├── tid_psam_feature_analysis.m
│   │   ├── tid_psam_hilbert_preparation.m
│   │   ├── tid_psam_ica_preprocessing.m
│   │   ├── tid_psam_questionnaire_input.m
│   │   ├── tid_psam_set_markers.m
│   │   ├── tid_psam_svm_preparation.m
│   │   └── tid_psam_svm_preprocessing.m
│   │
│   ├── R/
│   │   ├── tid_psam_beh_analysis.R
│   │   ├── tid_psam_erp_analysis.R
│   │   ├── tid_psam_questionnaire_analysis.R
│   │   └── tid_psam_svm_analysis.R
│   │
│   ├── praat/
│   │   │
│   │   ├── plugin_VocalToolkit/
│   │   │
│   │   ├── Praat.exe
│   │   ├── tid_psam_beh_preprocessing_1
│   │   ├── tid_psam_prepare_stimuli
│   │   └── tid_psam_rectangle.wav
│   │
│   ├── python/
│   │   ├── tid_psam_PYTHON_TEST.py
│   │   └── tid_psam_svm_analysis.py
│   │
│   └── tid_psam_run_pipeline.py
│
├── config/
│   ├── elec_96ch.csd
│   ├── elec_96ch.elp
│   ├── elec_96ch_adapted.csd
│   ├── elec_96ch_adapted.csv
│   └── elec_96ch_adapted.elp
│
├── documents/
│
├── experiment_script/
│   ├── tid_psam_create_conditions_file.m
│   ├── tid_psam_determine_loudness.m
│   ├── tid_psam_main_experiment_ALTERNATIVE_NO_CIRCLE.psyexp
│   ├── tid_psam_main_experiment_ALTERNATIVE_NO_CIRCLE.py
│   ├── tid_psam_select_stimuli.m
│   └── tid_psam_stimuli_recording_ALTERNATIVE_NO_CIRCLE_adapted.py
│
├── functions/
│   ├── bemobil_avref.m
│   ├── bemobil_detect_bad_channels.m
│   ├── shadedErrorBar.m
│   ├── tid_psam_archive_subj_TD.m
│   ├── tid_psam_check_folder_TD.R
│   ├── tid_psam_check_folder_TD.m
│   ├── tid_psam_check_folder_clean_up_folder_TD.py
│   ├── tid_psam_check_id_TD.m
│   ├── tid_psam_clean_up_folder_TD.R
│   ├── tid_psam_clean_up_folder_TD.m
│   ├── tid_psam_hjorth_activity_TD.m
│   ├── tid_psam_hjorth_complexity_TD.m
│   ├── tid_psam_hjorth_mobility_TD.m
│   ├── tid_psam_plot_flagged_ICs_TD.m
│   └── tid_psam_plot_rms_bins_TD.m
│
├── testing/
│
├── .gitattributes
├── .gitignore
├── README.md
└── requirements.txt
```
- ```anylsis_code``` contains all code relating to the analysis. Individual scripts can be found in the respective subfolders. To run the entire pipeline (which the raw data is needed for), run ```tid_psam_run_pipeline.py```. 
- ```config``` contains files related to the EEG system. Needed for the analysis.
- ```data``` contains the data. It is not included here be can be requested. If the data is available, it has to be placed in the path. Note that merely the raw data (```BIDS/``` folder) and the questionnaire data (```questionnaire_data/``` folder) is needed for running the pipeline, all other folders are create automatically.
- ```documents``` contains documents used during the data collection such as the participant information. Not needed to run any of the code.
- ```experiment_script``` contains all code related to running the experiment (mainly PsychoPy-based).
- ```functions``` contains all code related to any self-created and fetched functions. Some can be used independently of this project. Check them out.
- ```testing``` contains mutiple files related to testing. Not needed to run any of the code.
- ```requirements.txt``` contains all the libraries used for the python-based analysis scripts. Use ```pip install -r requirements.txt``` to install all needed dependencies.

Important Note. While the ```requirements.txt``` contains all needed dependencies for running the python-based analysis other dependencies have to be installed manually. This includes:
- All used R packages
- All used MATLAB toolboxes

```Praat.exe``` as well as the ```Praat Vocal Processing Toolbox``` are included in the repository. 

### Dependecies

### Usage

0. Make sure to have R, RStudio, MATLAB and VSCode installed.
1. Fetch this repository and clone it.
2. Run ```pip install -r requirements.txt``` while setting ```psam/``` as your current directory.
3. Open all R scripts in RStudio and install the needed packages.
4. Adapt the Path to your ```RScript.exe``` file in ```tid_psam_run_pipeline.py```.
5. Request the data and copy it into ```psam/``` without changing the folder structure.
6. Run ```tid_psam_run_pipeline.py``` using VSCode.


---

## Details

### Theoretical Background


### Methods


### Results


### Discussion




---



