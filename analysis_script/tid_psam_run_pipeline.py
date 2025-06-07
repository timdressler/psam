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

# Variables to edit
SKIP_SCRIPTS = [
    # Example: "tid_psam_erp_analysis.m", "tid_psam_questionnaire_analysis.R"
    "tid_psam_erp_preprocessing.m"
]
SKIP_ALREADY_RUN = True # If True, previously run scripts (based on the log file) are not executed again

# Initialize variables
log_file = os.path.join(SCRIPTPATH, "tid_psam_run_pipeline_log.txt")

# Load log file of previously run scripts
if SKIP_ALREADY_RUN and os.path.exists(log_file):
    with open(log_file, "r") as f:
        already_run_scripts = set(line.strip() for line in f)
else:
    already_run_scripts = set()

# Function: Update log file
def log_script(script_name):
    with open(log_file, "a") as f:
        f.write(f"{script_name}\n")

# Function: Call scripts
def run_script(command, description, script_name):
    if script_name in SKIP_SCRIPTS:
        print(f"\n ====================== Skipping {script_name} (manually excluded) ======================")
        return
    if SKIP_ALREADY_RUN and script_name in already_run_scripts:
        print(f"\n ====================== Skipping {script_name} (already run) ======================")
        return

    print(f"\n ====================== {description} ======================")
    print(f"[DEBUG] Command: {command}")
    try:
        subprocess.run(command, check=True)
        log_script(script_name)
    except subprocess.CalledProcessError as e:
        print(f"Error during: {description}")
        print(e)
        sys.exit(1)


# Run pipeline

# MATLAB: tid_psam_set_markers.m
script_name = "tid_psam_set_markers.m"
matlab_script = os.path.join(MAINPATH, "analysis_script", "MATLAB", script_name)
run_script(["matlab", "-batch", f"run('{matlab_script}')"], f"Running {script_name} (MATLAB)", script_name)

# Praat: tid_psam_beh_preprocessing_1.praat
script_name = "tid_psam_beh_preprocessing_1.praat"
if script_name in SKIP_SCRIPTS:
    print(f"\n ====================== Skipping {script_name} (manually excluded) ======================")
elif SKIP_ALREADY_RUN and script_name in already_run_scripts:
    print(f"\n ====================== Skipping {script_name} (already run) ======================")
else:
    print(f"\n====================== Running {script_name} (Praat) ======================")
    praat_exe = os.path.join(MAINPATH, "analysis_script", "praat", "Praat.exe")
    praat_script = os.path.join(MAINPATH, "analysis_script", "praat", "tid_psam_beh_preprocessing_1")
    subprocess.call([praat_exe, "--run", praat_script])
    log_script(script_name)

# MATLAB: Multiple scripts (see below)
matlab_scripts = [
    # "tid_psam_ica_preprocessing.m",
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
    run_script(["matlab", "-batch", f"run('{matlab_path}')"], f"Running {script} (MATLAB)", script)

# Python: tid_psam_svm_analysis.py
script_name = "tid_psam_svm_analysis.py"
py_script = os.path.join(MAINPATH, "analysis_script", "python", script_name)
run_script(["python", py_script], f"Running {script_name} (Python)", script_name)

# R: Multiple scripts (see below)
r_scripts = [
    "tid_psam_beh_analysis.R",
    "tid_psam_erp_analysis.R",
    "tid_psam_questionnaire_analysis.R",
    "tid_psam_scm_analysis.R"
]

for r_script in r_scripts:
    r_path = os.path.join(MAINPATH, "analysis_script", "R", r_script)
    run_script(["Rscript", r_path], f"Running {r_script} (R)", r_script)

# End of processing
check_done = "tid_psam_run_pipeline_DONE"
print(check_done)
