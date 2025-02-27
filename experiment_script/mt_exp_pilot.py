from psychopy import locale_setup
from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, hardware
from psychopy.tools import environmenttools
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER, priority)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

import psychopy.iohub as io
from psychopy.hardware import keyboard
import os
import pandas as pd
import random

# Create a device manager to handle hardware (keyboards, mice, mirophones, speakers, etc.)
deviceManager = hardware.DeviceManager()

# Store info about the experiment session
psychopyVersion = '2024.2.4'
expName = 'PSAM'  
expInfo = {'Subject': ''}
dlg = gui.DlgFromDict(expInfo, title='Experiment Info')
if not dlg.OK:
    core.quit()  # If the user cancels, end the experiment
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion


# Get IDs
subject_id = f"sub-{int(expInfo['Subject']):03d}"

# Set up paths
script_path = os.path.abspath(__file__)
main_path = script_path.replace('\experiment_script\mt_exp_pilot.py', '')
stimuli_path = os.path.join(main_path,'data', subject_id, 'stimuli')
out_path = os.path.join(main_path,'data', subject_id, 'data_experiment')
print(stimuli_path)

# Check paths
paths_to_check = [stimuli_path, out_path]
missing_paths = [path for path in paths_to_check if not os.path.exists(path)]

# Display warning if any paths are missing
if missing_paths:
    warning_dialog = gui.Dlg(title="Path Warning")
    warning_dialog.addText("The following paths do not exist:")
    for path in missing_paths:
        warning_dialog.addText(path)
    warning_dialog.show()

# Exit Experiment if any paths are missing
if missing_paths:
    core.quit()  # or use core.wait(seconds) to pause

# Set up condition file  
conditions_filename = os.path.join(stimuli_path, f"{subject_id}_conditions.xlsx")

import os
import random
import pandas as pd

import os
import random
import pandas as pd

# Set up parameters
n_trials = 16  # Adjust this to 960 if needed

# Define trial types and their counts  
trial_types = {
    'passive_early_normal_sham': n_trials // 16,
    'active_early_normal_sham': n_trials // 16,
    'passive_early_pitch_sham': n_trials // 16,
    'active_early_pitch_sham': n_trials // 16,
    'passive_late_normal_sham': n_trials // 16,
    'active_late_normal_sham': n_trials // 16,
    'passive_late_pitch_sham': n_trials // 16,
    'active_late_pitch_sham': n_trials // 16,
    'passive_early_normal_probe': n_trials // 16,
    'active_early_normal_probe': n_trials // 16,
    'passive_early_pitch_probe': n_trials // 16,
    'active_early_pitch_probe': n_trials // 16,
    'passive_late_normal_probe': n_trials // 16,
    'active_late_normal_probe': n_trials // 16,
    'passive_late_pitch_probe': n_trials // 16,
    'active_late_pitch_probe': n_trials // 16
}

# Generate trials and shuffle
trials = [trial for trial, count in trial_types.items() for _ in range(count)]
random.shuffle(trials)

# Function to get trial properties
def get_trial_params(trial_type):
    """Returns task, probe, probe_type, probe_onset_cat, probe_onset, probe_duration, probe_intensity, and probe_file."""
    
    # Task Type (Active/Passive)
    task = "Active" if "active" in trial_type else "Passive"
    
    # Probe Presence (Yes/No)
    probe = "No" if "sham" in trial_type else "Yes"
    
    # Probe Type (Normal/Pitch) - Sham trials inherit the same probe properties
    probe_type = "Pitch" if "pitch" in trial_type else "Normal"
    
    # Probe Onset Category (Early/Late)
    probe_onset_cat = "Early" if "early" in trial_type else "Late"
    
    # Probe Onset (2.8 for early, 2.9 for late, empty for sham)
    probe_onset = 2.8 if "early" in trial_type and probe == "Yes" else (2.9 if "late" in trial_type and probe == "Yes" else 0)

    # Probe Duration (0.08 for probes, 0 for sham)
    probe_duration = 0.08 if probe == "Yes" else 0

    # Probe Intensity (1 for probes, 0 for sham)
    probe_intensity = 1 if probe == "Yes" else 0

    # Probe File (None for sham, otherwise filename based on probe type)
    probe_file = None if probe == "No" else os.path.join(stimuli_path, f"{subject_id}_{probe_type.lower()}_probe.wav")

    return task, probe, probe_type, probe_onset_cat, probe_onset, probe_duration, probe_intensity, probe_file

# Create DataFrame with required columns
df = pd.DataFrame({
    "subj": subject_id,
    "trial_type": trials,
    "trial_counter": range(len(trials))  # Sequential trial counter
})

# Apply the mapping function
df[[
    "task", "probe", "probe_type", "probe_onset_cat",
    "probe_onset", "probe_duration", "probe_intensity", "probe_file"
]] = df["trial_type"].apply(lambda t: pd.Series(get_trial_params(t)))

# Save DataFrame
conditions_filename = os.path.join(stimuli_path, f"{subject_id}_conditions.xlsx")
df.to_excel(conditions_filename, index=False)




core.quit()


