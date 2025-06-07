import subprocess
import os
import sys
import re

"""
tid_psam_run_pipeline.py

Master script for running the full analysis pipeline for the PSAM project across all languages.

Pipeline overview:
- Runs initial MATLAB script for setting markers and creating folder structure
- Executes Praat preprocessing script 
- Runs a sequence of MATLAB scripts for preprocessing and analysis
- Runs Python script for SVM classification
- Executes R scripts for behavioral, ERP and SVM analysis

All paths are dynamically resolved relative to the repository root.

Note:
- Praat is called with subprocess.call() due to its requirements
- The Praat script requires hardcoding the location of it

- Execution will stop if any subprocess fails

Tim Dressler, 07.04.2025
"""


# Set up paths
SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
expected_subpath = os.path.join('psam', 'analysis_script')
if not re.search(re.escape(expected_subpath) + r'$', SCRIPTPATH):
    raise EnvironmentError('Path not OK')

MAINPATH = os.path.abspath(os.path.join(SCRIPTPATH, '..'))

# Function: Call scripts
def run_script(command, description):
    print(f"\n ====================== {description} ======================")
    print(f"[DEBUG] Command: {command}")
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error during: {description}")
        print(e)
        sys.exit(1)

# Run pipeline

# MATLAB: tid_psam_set_markers.m
matlab_script = os.path.join(MAINPATH, "analysis_script", "MATLAB", "tid_psam_set_markers.m")
run_script(["matlab", "-batch", f"run('{matlab_script}')"], "Running tid_psam_set_markers.m (MATLAB)") # Has to be run in the beginning as the Praat can't create needed folder structures

# Praat: tid_psam_beh_preprocessing_1.praat
print("\n====================== Running tid_psam_beh_preprocessing_1.praat (Praat) ======================")
praat_exe = os.path.join(MAINPATH, "analysis_script", "praat", "Praat.exe")
praat_script = os.path.join(MAINPATH, "analysis_script", "praat", "tid_psam_beh_preprocessing_1")
subprocess.call([praat_exe, "--run", praat_script]) # Has to be run using subprocess.call, and thus not using run_script

# MATLAB: Multiple scripts (see below)
matlab_scripts = [
    # "tid_psam_ica_preprocessing",
    "tid_psam_erp_preprocessing.m",
    "tid_psam_svm_preprocessing.m",
    "tid_psam_beh_preprocessing_2.m",
    "tid_psam_exclude_trials.m",
    "tid_psam_erp_analysis.m",
    "tid_psam_hilbert_preparation.m",
    "tid_psam_svm_preparation.m",
    "tid_psam_beh_analysis.m",
    "tid_psam_feature_analysis.m"
]

for script in matlab_scripts:
    matlab_path = os.path.join(MAINPATH, "analysis_script", "MATLAB", script)
    run_script(["matlab", "-batch", f"run('{matlab_path}')"], f"Running {script}.m (MATLAB)")

# Python: tid_psam_svm_analysis.py
py_script = os.path.join(MAINPATH, "analysis_script", "python", "tid_psam_svm_analysis.py")
run_script(["python", py_script], "Running tid_psam_svm_analysis.py (Python)")

# R: Multiple scripts (see below)
r_scripts = [
    "tid_psam_beh_analysis.R",
    "tid_psam_erp_analysis.R",
    "tid_psam_questionnaire_analysis.R",
    "tid_psam_scm_analysis.R"
]

for r_script in r_scripts:
    r_path = os.path.join(MAINPATH, "analysis_script", "R", r_script)
    run_script(["Rscript", r_path], f"Running {r_script} (R)")

# End of processing

check_done = "tid_psam_run_pipeline_DONE"
print(check_done)
