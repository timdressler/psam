# WORK IN PROGRESS!

# Investigating the Specificity of Pre-Speech Auditory Modulation - From Global Gating to Selective Silence?

## In Brief

## This Repository

This repository includes all code related to the 'Investigating the Specificity of Pre-Speech Auditory Modulation - From Global Gating to Selective Silence?' project conducted by Tim Dressler at the University of Oldenburg. The project was supervised by Prof. Dr. Stefan Debener and Prof. Dr. Andrea Hildebrandt.

The repository includes analysis code as well as as code used for running the experiment. No data is included, but can by requested. Next to code, the repository also includes documents such as the participant information. Furthermore, a short theoretical description of the project and its result can be found in this README.

### Structure

```text
psam\
│
├── analysis_script\
│   │
│   ├── MATLAB\
│   │   ├── tid_psam_MATLAB_TEST.m
│   │   ├── tid_psam_beh_analysis.m
│   │   ├── tid_psam_beh_preprocessing_2.m
│   │   ├── tid_psam_check_filters_erp.m
│   │   ├── tid_psam_check_filters_svm.m
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
│   ├── R\
│   │   ├── tid_psam_beh_analysis.R
│   │   ├── tid_psam_erp_analysis.R
│   │   ├── tid_psam_get_package_versions.R
│   │   ├── tid_psam_install_requirements_R.R
│   │   ├── tid_psam_questionnaire_analysis.R
│   │   └── tid_psam_svm_analysis.R
│   │
│   ├── praat\
│   │   │
│   │   ├── plugin_VocalToolkit\
│   │   │
│   │   ├── Praat.exe
│   │   ├── tid_psam_beh_preprocessing_1.praat
│   │   └── tid_psam_prepare_stimuli.praat
│   │
│   ├── python\
│   │   └── tid_psam_svm_analysis.py
│   │
│   └── tid_psam_run_pipeline.py
│
├── config\
│
├── data\
│   │
│   ├── BIDS\
│   │
│   └── questionnaire_data\
│
│
├── documents\
│
├── experiment_script\
│   ├── tid_psam_create_conditions_file.m
│   ├── tid_psam_determine_loudness.m
│   ├── tid_psam_main_experiment_ALTERNATIVE_NO_CIRCLE.psyexp
│   ├── tid_psam_main_experiment_ALTERNATIVE_NO_CIRCLE.py
│   ├── tid_psam_select_stimuli.m
│   └── tid_psam_stimuli_recording_ALTERNATIVE_NO_CIRCLE_adapted.py
│
├── functions\
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
│   ├── tid_psam_get_transition_bandwidth_TD.m
│   ├── tid_psam_hjorth_activity_TD.m
│   ├── tid_psam_hjorth_complexity_TD.m
│   ├── tid_psam_hjorth_mobility_TD.m
│   ├── tid_psam_plot_flagged_ICs_TD.m
│   └── tid_psam_plot_rms_bins_TD.m
│
├── testing\
│
├── .gitattributes
├── .gitignore
├── README.md
└── requirements.txt
```
- ```anylsis_code``` contains all code relating to the analysis. Individual scripts can be found in the respective subfolders. To run the entire pipeline (which the raw data is needed for), run ```tid_psam_run_pipeline.py```. 
- ```config``` contains files related to the EEG system. Needed for the analysis.
- ```data``` contains the data. It is not included here be can be requested. If the data is available, it has to be placed in the path. Note that the raw data (```BIDS/```) and the questionnaire data (```questionnaire_data/```) is needed for running the pipeline. Processed data will be saved automatically in further subfolders under ```data/```.
- ```documents``` contains documents used during the data collection such as the participant information. Not needed to run any of the code.
- ```experiment_script``` contains all code related to running the experiment (mainly PsychoPy-based).
- ```functions``` contains all code related to any self-created and fetched functions. Some can be used independently of this project. Check them out.
- ```testing``` contains mutiple files related to testing. Not needed to run any of the code.
- ```requirements.txt``` contains all the libraries used for the python-based analysis scripts. Use ```pip install -r requirements.txt``` to install all needed python dependencies.

> [!NOTE]  
> While the ```requirements.txt``` contains all needed dependencies for running the python-based analysis and ```tid_psam_install_requirements_R.R``` installs all the requirements needed for the R-based
analysis, MATLAB dependencies have to be installed manually (see below). ```Praat.exe``` as well as the ```Praat Vocal Processing Toolbox``` are included in the repository. 

### Dependecies

This table includes all needed software to run the analysis pipeline. 

### Usage

0. Make sure to have R, RStudio, MATLAB and VSCode installed.
1. Fetch this repository and clone it.
2. Run ```pip install -r requirements.txt``` while setting ```psam/``` as your current directory.
3. Adapt the Path to your ```RScript.exe``` file in ```tid_psam_run_pipeline.py```.
4. Request the data and copy it into ```psam/``` without changing the folder structure.
5. Run ```tid_psam_run_pipeline.py``` using VSCode.






