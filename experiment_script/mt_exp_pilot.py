from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import os

# Store info about the experiment session
psychopyVersion = 'XXX'
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

# Start experiment set up
# Set up the window
win = visual.Window([1024,768], fullscr=False, units='pix')

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










# End the experiment 
win.close()
core.quit()

