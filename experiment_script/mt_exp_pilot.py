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

# Set up number of trials  
n_trials = 16 # Set to 960
n_practice_trials = 10 # Set to 20

# Set up practice trials  
practice_trials = ['passive_practice'] * (n_practice_trials // 2) + ['active_practice'] * (n_practice_trials // 2)  
random.shuffle(practice_trials)

# Define trial types and their counts  
trial_types = {
    'passive_no_probe': n_trials // 4,
    'active_no_probe': n_trials // 4,
    'passive_early_normal_probe': n_trials // 16,
    'active_early_normal_probe': n_trials // 16,
    'passive_early_pitch_probe': n_trials // 16,
    'active_early_pitch_probe': n_trials // 16,
    'passive_late_normal_probe': n_trials // 16,
    'active_late_normal_probe': n_trials // 16,
    'passive_late_pitch_probe': n_trials // 16,
    'active_late_pitch_probe': n_trials // 16
}

# Generate trials  
trials = [trial for trial, count in trial_types.items() for _ in range(count)]
random.shuffle(trials)

# Combine practice trials with main trials  
trials = practice_trials + trials  

# Create DataFrame  
df = pd.DataFrame({'trial_type': trials})

# Mapping functions  
def get_probe_params(trial_type):
    """Returns probe_onset, probe_duration, probe_type, and probe_intensity based on trial type."""
    
    if 'no_probe' in trial_type or 'practice' in trial_type:
        return 0, 0, '', 0  # No probe trials → onset = 0, duration = 0, no file, intensity = 0

    # Determine onset time explicitly
    if 'early' in trial_type:
        onset = 2.8
    elif 'late' in trial_type:
        onset = 2.9
    else:
        onset = quit()

    # Determine probe type
    if 'pitch' in trial_type:
        probe_type = f"{subject_id}_pitch_probe.wav"
    else:
        probe_type = f"{subject_id}_normal_probe.wav"

    return onset, 0.08, os.path.join(stimuli_path, probe_type), 1  # Return full parameters


# Apply mapping functions  
df[['probe_onset', 'probe_duration', 'probe_type', 'probe_intensity']] = df['trial_type'].apply(
    lambda t: pd.Series(get_probe_params(t))
)
df['instruction'] = df['trial_type'].apply(lambda t: 'Active' if 'active' in t else 'Passive')

# Save DataFrame  
df.to_excel(conditions_filename, index=False)


##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################

# Start experiment set up
# Set up the window
win = visual.Window([1920,1080], fullscr=True, units='pix')

# Set up stimuli parameters
# Initialize parameters for the cue
initial_radius = 0.7  # Starting size of the cue
final_radius = 0.2    # Final size of the cue
duration = 2.0        # Time in seconds to shrink
hold_time = 1.5       # Time in seconds to hold after reaching final size
    
cue_started = False    # To track when the cue starts shrinking
cue_start_time = None  # To store when the cue started shrinking

# Create Objects
fixation_cross = visual.TextStim(win=win, name='fixation_cross',
    text='+',
    font='Arial',
    pos=(0, 0), draggable=False, height=0.5, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0)

cue = visual.ShapeStim(
    win=win, name='cue',
    size=(0.7, 0.7), vertices='circle',
    ori=0.0, pos=(0, 0), draggable=False, anchor='center',
    lineWidth=1.0,
    colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
    opacity=None, depth=-2.0, interpolate=True)

target = visual.ShapeStim(
    win=win, name='target',
    size=(0.2, 0.2), vertices='circle',
    ori=0.0, pos=(0, 0), draggable=False, anchor='center',
    lineWidth=2.0,
    colorSpace='rgb', lineColor='white', fillColor=[0.0000, 0.0000, 0.0000],
    opacity=1.0, depth=-3.0, interpolate=True)

normal_probe_file = os.path.join(stimuli_path, f"{subject_id}_normal_probe.wav")  
normal_probe = sound.Sound(normal_probe_file,
    stereo=True, 
    hamming=True, 
    speaker='normal_probe',    name='normal_probe'
    )
normal_probe.setSound(normal_probe_file, hamming=True)
normal_probe.setVolume(1.0, log=False)



# End the experiment
win.close()
core.quit()


