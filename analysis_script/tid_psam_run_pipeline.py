import subprocess
import os
import sys
import re
import time
import pandas as pd
import numpy as np

"""
tid_psam_run_pipeline.py

Master script for running the full analysis pipeline for the PSAM project across all languages.

Pipeline overview:
- tid_psam_set_markers.m (MATLAB)
- tid_psam_beh_preprocessing_1.praat (Praat)
- tid_psam_erp_preprocessing.m (MATLAB)
- tid_psam_svm_preprocessing.m (MATLAB)
- tid_psam_beh_preprocessing_2.m (MATLAB)
- tid_psam_exclude_trials.m (MATLAB)
- tid_psam_erp_analysis.m (MATLAB)
- tid_psam_hilbert_preparation.m (MATLAB)
- tid_psam_svm_preparation.m (MATLAB)
- tid_psam_beh_analysis.m (MATLAB)
- tid_psam_feature_analysis.m (MATLAB)
- tid_psam_svm_analysis.py (Python)
- tid_psam_beh_analysis.R (R)
- tid_psam_erp_analysis.R (R)
- tid_psam_questionnaire_analysis.R (R)
- tid_psam_svm_analysis.R (R)

Logging:
- A plain .txt log ("tid_psam_run_pipeline_log.txt") tracks successfully completed scripts (used for skip logic)
- A .xlsx protocol file ("tid_psam_run_pipeline_protocol.xlsx") logs:
    - script: the filename of each script
    - status: one of "OK", "ERROR", or "SKIPPED"
    - runtime_sec: execution time in seconds (NaN for skipped)

Skipping options:
- SKIP_SCRIPTS: manually list script filenames (including extension) to skip regardless of log status
- SKIP_ALREADY_RUN: if True, skips scripts that appear in the plain-text log file

All paths are dynamically resolved relative to the repository root.

Note:
- Praat is called with subprocess.call() due to its requirements
- The Praat script requires hardcoding the location of it

- The path to Rscript.exe has to be adapted (hardcoded)

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
######################################## For Testing ########################################
# RSCRIPT_EXE = r"C:\Program Files\R\R-4.5.0\bin\Rscript.exe" # Path to Rscript.exe
######################################## For Testing ########################################

RSCRIPT_EXE = r"C:\Program Files\R\R-4.5.0\bin\Rscript.exe" # Path to Rscript.exe

SKIP_SCRIPTS = [ # Manually excluded scripts. Example: "tid_psam_erp_analysis.m", "tid_psam_questionnaire_analysis.R", ...
   "tid_psam_ica_preprocessing.m", "tid_psam_set_markers.m", "tid_psam_beh_preprocessing_2", "tid_psam_beh_analysis.m", "tid_psam_beh_preprocessing_1.praat", "tid_psam_svm_analysis.py", "tid_psam_beh_preprocessing_2.m", "tid_psam_svm_preprocessing", 
   ""
]
SKIP_ALREADY_RUN = True # If True, previously run scripts (based on the log file) are not executed again

# Prepare variables
log_txt_path = os.path.join(SCRIPTPATH, "tid_psam_run_pipeline_log.txt")
log_xlsx_path = os.path.join(SCRIPTPATH, "tid_psam_run_pipeline_last_run_protocol.xlsx")
protocol_df = pd.DataFrame(columns=["script", "status", "runtime_sec"])

# Load log file of previously run scripts
if SKIP_ALREADY_RUN and os.path.exists(log_txt_path):
    with open(log_txt_path, "r") as f:
        already_run_scripts = set(line.strip() for line in f)
else:
    already_run_scripts = set()

# Function: Update log file (.txt)
def log_script_txt(script_name):
    with open(log_txt_path, "a") as f:
        f.write(f"{script_name}\n")

# Function: Update protocol file (.xlsx)
def update_protocol(script_name, status, runtime=None):
    global protocol_df
    row = {
        "script": script_name,
        "status": status,
        "runtime_sec": np.nan if runtime is None else round(runtime, 2)
    }
    protocol_df = pd.concat([protocol_df, pd.DataFrame([row])], ignore_index=True)
    protocol_df.to_excel(log_xlsx_path, index=False)

# Function: Call scripts
def run_script(command, description, script_name, use_call=False):
    if script_name in SKIP_SCRIPTS:
        print(f"\n ====================== Skipping {script_name} (manually excluded) ======================")
        update_protocol(script_name, "SKIPPED")
        return
    if SKIP_ALREADY_RUN and script_name in already_run_scripts:
        print(f"\n ====================== Skipping {script_name} (already run) ======================")
        update_protocol(script_name, "SKIPPED")
        return

    print(f"\n ====================== {description} ======================")
    print(f"[DEBUG] Command: {command}")
    start = time.time()
    try:
        if use_call:
            subprocess.call(command)
        else:
            subprocess.run(command, check=True)
        runtime = time.time() - start
        log_script_txt(script_name)
        update_protocol(script_name, "OK", runtime)
    except subprocess.CalledProcessError as e:
        runtime = time.time() - start
        update_protocol(script_name, "ERROR", runtime)
        print(f"Error during: {description}")
        print(e)
        sys.exit(1)


# Run pipeline (in order)

# MATLAB: tid_psam_set_markers.m
script_name = "tid_psam_set_markers.m"
matlab_script = os.path.join(MAINPATH, "analysis_script", "MATLAB", script_name)
run_script(["matlab", "-batch", f"run('{matlab_script}')"], f"Running {script_name} (MATLAB)", script_name)

# Praat: tid_psam_beh_preprocessing_1.praat
script_name = "tid_psam_beh_preprocessing_1.praat"
praat_exe = os.path.join(MAINPATH, "analysis_script", "praat", "Praat.exe") # Path to Praat.exe
praat_script = os.path.join(MAINPATH, "analysis_script", "praat", "tid_psam_beh_preprocessing_1")
run_script([praat_exe, "--run", praat_script], f"Running {script_name} (Praat)", script_name, use_call=True)

# MATLAB: Multiple scripts (in order, see below)
matlab_scripts = [
    "tid_psam_ica_preprocessing.m",
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

# R: Multiple scripts (in order, see below)
r_scripts = [
    "tid_psam_install_requirements_R.R",
    "tid_psam_beh_analysis.R",
    "tid_psam_erp_analysis.R",
    "tid_psam_questionnaire_analysis.R",
    "tid_psam_svm_analysis.R"
]

for r_script in r_scripts:
    r_path = os.path.join(MAINPATH, "analysis_script", "R", r_script)
    run_script([RSCRIPT_EXE, r_path], f"Running {r_script} (R)", r_script)

# End of processing

check_done = "tid_psam_run_pipeline_DONE"
print(check_done)
