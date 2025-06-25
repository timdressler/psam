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
├── data\ (not included in this repository, can be requested)
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

| Software         | Version           | Status                     |
|------------------|-------------------|----------------------------|
| MATLAB           | R2024a            | Manual Installation needed |
| R                | 4.5.0             | Manual Installation needed |
| Python           | 3.13.4            | Manual Installation needed |
| Praat            | 6.4.27            | Included in the repository, no installation needed |

This table includes all needed R packages to run the analysis pipeline. 

| R Package        | Version   | Status                |
|------------------|-----------|-----------------------|
| tidyr            | 1.3.1     | Automatically installed |
| afex             | 1.4.1     | Automatically installed |
| emmeans          | 1.11.1    | Automatically installed |
| readxl           | 1.4.5     | Automatically installed |
| car              | 3.1.3     | Automatically installed |
| corrplot         | 0.95      | Automatically installed |
| dplyr            | 1.1.4     | Automatically installed |
| ez               | 4.4.0     | Automatically installed |
| ggplot2          | 3.5.2     | Automatically installed |
| ggstatsplot      | 0.13.1    | Automatically installed |
| DescTools        | 0.99.60   | Automatically installed |
| ggpubr           | 0.6.0     | Automatically installed |
| ggeffects        | 2.2.1     | Automatically installed |
| cowplot          | 1.1.3     | Automatically installed |
| tidyverse        | 2.0.0     | Automatically installed |
| psych            | 2.5.3     | Automatically installed |
| rstatix          | 0.7.2     | Automatically installed |
| rstudioapi       | 0.17.1    | Automatically installed |
| stringr          | 1.5.1     | Automatically installed |
| lme4             | 1.1.37    | Automatically installed |
| performance      | 0.14.0    | Automatically installed |

This table includes all needed Python packages to run the analysis pipeline. 

| Python Package      | Version        | Status                |
|---------------------|----------------|-----------------------|
| colorama            | 0.4.6          | Automatically installed |
| contourpy           | 1.3.2          | Automatically installed |
| cycler              | 0.12.1         | Automatically installed |
| et_xmlfile          | 2.0.0          | Automatically installed |
| fonttools           | 4.58.2         | Automatically installed |
| joblib              | 1.5.1          | Automatically installed |
| kiwisolver          | 1.4.8          | Automatically installed |
| matplotlib          | 3.10.3         | Automatically installed |
| numpy               | 2.3.0          | Automatically installed |
| openpyxl            | 3.1.5          | Automatically installed |
| packaging           | 25.0           | Automatically installed |
| pandas              | 2.3.0          | Automatically installed |
| pillow              | 11.2.1         | Automatically installed |
| pip                 | 25.1           | Automatically installed |
| pyparsing           | 3.2.3          | Automatically installed |
| python-dateutil     | 2.9.0.post0    | Automatically installed |
| pytz                | 2025.2         | Automatically installed |
| scikit-learn        | 1.7.0          | Automatically installed |
| scipy               | 1.15.3         | Automatically installed |
| seaborn             | 0.13.2         | Automatically installed |
| setuptools          | 78.1.1         | Automatically installed |
| six                 | 1.17.0         | Automatically installed |
| threadpoolctl       | 3.6.0          | Automatically installed |
| tqdm                | 4.67.1         | Automatically installed |
| tzdata              | 2025.2         | Automatically installed |
| wheel               | 0.45.1         | Automatically installed |

This table includes all needed Praat packages to run the analysis pipeline. 

| Praat Tool          | Version        | Status                                      |
|---------------------|----------------|---------------------------------------------|
| Praat Vocal Toolkit | 2012–2024      | Included in the repository, no installation needed |



### Usage

> [!IMPORTANT]
> It is assumed that Conda, R and MATLAB (and the MATLAB dependecies) are installed! It is recommended that the Pipeline is run using VSCode.

It is recommended to use a dedicated Python environment (e.g. through conda) to mitigate the risk of potential version conflicts.

```
conda create -n psam python==3.13.4
conda activate psam
pip install -r requirements.txt
```

> [!WARNING]
> As all scripts are called from ```tid_psam_run_pipeline.py```, the path the ```RScript.exe``` needs to be adapted manually in the script!

After requesting the data copy it into the folder as specified above and run ```tid_psam_run_pipeline.py```.





