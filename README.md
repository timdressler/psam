# Investigating the Specificity of Pre-Speech Auditory Modulation - From Global Gating to Selective Silence?

[![DOI](https://zenodo.org/badge/936619418.svg)](https://doi.org/10.5281/zenodo.16632078) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![Lines of Code](https://tokei.rs/b1/github/timdressler/psam)

> [!NOTE]
>  While this repository might be subject to further changes and optimizations, the repository was further published as a time-stamped, immutable, and permanent version under the DOI above which can be used for exact replication attemps.

## This Repository

This repository contains all code related to the research project *Investigating the Specificity of Pre-Speech Auditory Modulation – From Global Gating to Selective Silence?*, conducted by **Tim Dressler** at the **University of Oldenburg**, under the supervision of **Prof. Dr. Stefan Debener** and **Prof. Dr. Andrea Hildebrandt**.

The study was preregistered on [OSF](https://doi.org/10.17605/OSF.IO/23EYF). Any deviations from the preregistration are documented [here](https://osf.io/9quhn).

The repository includes:

- **Analysis code** used for data evaluation  
- **Experimental code** used for stimulus presentation and data collection  
- **Supporting documents**: participant instructions, questionnaires, etc.  
- **Publications** (e.g., posters, papers) – added as they become available  

> **Note**: No data is included in this repository but may be made available upon request.

### Structure

```text
psam\
│
├── analysis_script\
│   │
│   ├── MATLAB\
│   │   ├── tid_psam_beh_analysis.m
│   │   ├── tid_psam_beh_preprocessing_2.m
│   │   ├── tid_psam_check_filters_erp.m
│   │   ├── tid_psam_check_filters_svm.m
│   │   ├── tid_psam_erp_analysis.m
│   │   ├── tid_psam_erp_pilot.m
│   │   ├── tid_psam_erp_preprocessing.m
│   │   ├── tid_psam_exclude_trials.m
│   │   ├── tid_psam_feature_analysis.m
│   │   ├── tid_psam_helper.m
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
├── data\ (NOT INCLUDED IN THIS REPOSITORY, CAN BE REQUESTED)
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
│   ├── tid_psam_copy_files_from_subfolders_TD.m
│   ├── tid_psam_get_transition_bandwidth_TD.m
│   ├── tid_psam_hjorth_activity_TD.m
│   ├── tid_psam_hjorth_complexity_TD.m
│   ├── tid_psam_hjorth_mobility_TD.m
│   ├── tid_psam_plot_flagged_ICs_TD.m
│   └── tid_psam_plot_rms_bins_TD.m
│
├── publications\
│
├── testing\
│
├── .gitattributes
├── .gitignore
├── LICENSE
├── README.md
└── requirements.txt
```
- ```analysis_script``` contains all code relating to the analysis. Individual scripts can be found in the respective subfolders. To run the entire pipeline (which (only) the raw data is needed for), run ```tid_psam_run_pipeline.py```. 
- ```config``` contains files related to the EEG system. Needed for the analysis.
- ```data``` contains the data. It is not included here be can be requested. If the data is available, it has to be placed in the path. Note that the raw data (```BIDS/```) and the questionnaire data (```questionnaire_data/```) is needed for running the pipeline. Processed data will be saved automatically in further subfolders under ```data/```.
- ```documents``` contains documents used during the data collection such as the participant information, instructions and the used questionnaires. Not needed to run any of the code.
- ```experiment_script``` contains all code related to running the experiment (mainly PsychoPy-based).
- ```functions``` contains all code related to any self-created and fetched functions. Some can be used independently of this project. Check them out.
- ```publications``` contains publications (poster, thesis, ...) associated with this research project (once they exist).
- ```testing``` contains mutiple files related to testing. Not needed to run any of the code.
- ```requirements.txt``` contains all the libraries used for the python-based analysis scripts. Use ```pip install -r requirements.txt``` to install all needed python dependencies.

> [!NOTE]  
> While the ```requirements.txt``` contains all needed dependencies for running the python-based analysis and ```tid_psam_install_requirements_R.R``` installs all the requirements needed for the R-based
analysis, MATLAB dependencies have to be installed manually (see below). ```Praat.exe``` as well as the ```Praat Vocal Processing Toolbox``` are included in the repository. 

### Usage

It is recommended to use a dedicated Python environment (e.g., through Conda) to mitigate the risk of potential version conflicts.

```
conda create -n psam python==3.13.4
conda activate psam
```

To install Python dependencies run the following line.

```
pip install -r requirements.txt
```

To install R dependencies, run ```tid_psam_install_requirements_R.R``` in RStudio.

It is assumed that Conda, R and MATLAB (and the MATLAB dependecies, see below) are installed! It is recommended that the pipeline is run using VSCode.

> [!WARNING]
> As all scripts are called from ```tid_psam_run_pipeline.py```, the path the ```RScript.exe``` needs to be adapted manually in the script!

After requesting the data copy it into the folder as specified above and run ```tid_psam_run_pipeline.py``` in VSCode.

### Dependecies

This table includes all needed software to run the analysis pipeline. 

| Type            | Software/Package             | Version         | Status                                           |
|-----------------|------------------------------|------------------|--------------------------------------------------|
| Software        | MATLAB                       | R2024a           | Manual Installation needed                       |
| Software        | R                            | 4.5.0            | Manual Installation needed                       |
| Software        | Python                       | 3.13.4           | Manual Installation needed                       |
| Software        | Praat                        | 6.4.27           | Included in the repository, no installation needed |
| MATLAB Toolbox  | EEGLAB                       | v2023.0          | Manual Installation needed                       |
| MATLAB Toolbox  | ICLabel                      | v1.4             | Manual Installation needed                       |
| MATLAB Toolbox  | Audio Toolbox                | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Control System Toolbox       | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Curve Fitting Toolbox        | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Data Acquisition Toolbox     | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | DSP System Toolbox           | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Image Acquisition Toolbox    | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Image Processing Toolbox     | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Instrument Control Toolbox   | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Mapping Toolbox              | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Parallel Computing Toolbox   | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Signal Processing Toolbox    | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | Statistics and ML Toolbox    | 24.1             | Manual Installation needed                       |
| MATLAB Toolbox  | System Identification Toolbox| 24.1             | Manual Installation needed                       |
| R Package       | tidyr                        | 1.3.1            | Automatically installed                          |
| R Package       | afex                         | 1.4.1            | Automatically installed                          |
| R Package       | emmeans                      | 1.11.1           | Automatically installed                          |
| R Package       | readxl                       | 1.4.5            | Automatically installed                          |
| R Package       | car                          | 3.1.3            | Automatically installed                          |
| R Package       | corrplot                     | 0.95             | Automatically installed                          |
| R Package       | dplyr                        | 1.1.4            | Automatically installed                          |
| R Package       | ez                           | 4.4.0            | Automatically installed                          |
| R Package       | ggplot2                      | 3.5.2            | Automatically installed                          |
| R Package       | ggstatsplot                  | 0.13.1           | Automatically installed                          |
| R Package       | DescTools                    | 0.99.60          | Automatically installed                          |
| R Package       | ggpubr                       | 0.6.0            | Automatically installed                          |
| R Package       | ggeffects                    | 2.2.1            | Automatically installed                          |
| R Package       | cowplot                      | 1.1.3            | Automatically installed                          |
| R Package       | tidyverse                    | 2.0.0            | Automatically installed                          |
| R Package       | psych                        | 2.5.3            | Automatically installed                          |
| R Package       | rstatix                      | 0.7.2            | Automatically installed                          |
| R Package       | rstudioapi                   | 0.17.1           | Automatically installed                          |
| R Package       | stringr                      | 1.5.1            | Automatically installed                          |
| R Package       | lme4                         | 1.1.37           | Automatically installed                          |
| R Package       | performance                  | 0.14.0           | Automatically installed                          |
| R Package       | Rmisc                        | 1.5.1            | Automatically installed                          |
| R Package       | nparLD                       | 2.2              | Automatically installed                          |
| R Package       | RVAideMemoire                | 0.9.83.11        | Automatically installed                          |
| R Package       | effectsize                   | 1.0.1            | Automatically installed                          |
| R Package       | lmerTest                     | 3.1.3            | Automatically installed                          |
| R Package       | smplot2                      | 0.2.5            | Automatically installed                          |
| R Package       | devtools                     | 2.4.5            | Automatically installed                          |
| R Package       | renv                         | 1.1.4            | Automatically installed                          |
| Python Package  | colorama                     | 0.4.6            | Automatically installed                          |
| Python Package  | contourpy                    | 1.3.2            | Automatically installed                          |
| Python Package  | cycler                       | 0.12.1           | Automatically installed                          |
| Python Package  | et_xmlfile                   | 2.0.0            | Automatically installed                          |
| Python Package  | fonttools                    | 4.58.2           | Automatically installed                          |
| Python Package  | joblib                       | 1.5.1            | Automatically installed                          |
| Python Package  | kiwisolver                   | 1.4.8            | Automatically installed                          |
| Python Package  | matplotlib                   | 3.10.3           | Automatically installed                          |
| Python Package  | numpy                        | 2.3.0            | Automatically installed                          |
| Python Package  | openpyxl                     | 3.1.5            | Automatically installed                          |
| Python Package  | packaging                    | 25.0             | Automatically installed                          |
| Python Package  | pandas                       | 2.3.0            | Automatically installed                          |
| Python Package  | pillow                       | 11.2.1           | Automatically installed                          |
| Python Package  | pip                          | 25.1             | Automatically installed                          |
| Python Package  | pyparsing                    | 3.2.3            | Automatically installed                          |
| Python Package  | python-dateutil              | 2.9.0.post0      | Automatically installed                          |
| Python Package  | pytz                         | 2025.2           | Automatically installed                          |
| Python Package  | scikit-learn                 | 1.7.0            | Automatically installed                          |
| Python Package  | scipy                        | 1.15.3           | Automatically installed                          |
| Python Package  | seaborn                      | 0.13.2           | Automatically installed                          |
| Python Package  | setuptools                   | 78.1.1           | Automatically installed                          |
| Python Package  | six                          | 1.17.0           | Automatically installed                          |
| Python Package  | threadpoolctl                | 3.6.0            | Automatically installed                          |
| Python Package  | tqdm                         | 4.67.1           | Automatically installed                          |
| Python Package  | tzdata                       | 2025.2           | Automatically installed                          |
| Python Package  | wheel                        | 0.45.1           | Automatically installed                          |
| Praat Tool      | Praat Vocal Toolkit          | 2012–2024        | Included in the repository, no installation needed |




