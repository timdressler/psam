from psychopy import prefs
prefs.hardware['audioLib'] = ['PTB']  # Use PTB for precise audio timing
prefs.hardware['audioLatencyMode'] = 3  # High-precision mode

from psychopy import sound, core, visual, event, parallel
import psychtoolbox as ptb
import random
import os

# Create directory for recordings (still used for saving files in future if needed)
recording_dir = "recordings"
if not os.path.exists(recording_dir):
    os.makedirs(recording_dir)

# Initialize PsychoPy window
win = visual.Window(fullscr=True, color=[-1, -1, -1], units='height')  # 'height' ensures true circles

# Fixation cross
fixation = visual.TextStim(win, text='+', color='white', height=0.05)

# Circles (inner and outer)
inner_circle = visual.Circle(win, radius=0.05, edges=128, lineColor='white', fillColor=None)
outer_circle = visual.Circle(win, radius=0.1, edges=128, lineColor='white', fillColor=None)

# Initialize PTB sound system
probe_stim = sound.Sound('C:/Users/timdr/OneDrive/Uni_Oldenburg/4_Semester/Master_Thesis/Analysis_Experiment/psam/testing/stereo_sine_wave.wav', secs=0.2, stereo=True, sampleRate=44100)  # Load a custom sound file

# Initialize Parallel Port
parallel_port = parallel.ParallelPort(address=0x0378)  # Common address for LPT1; change if necessary for your hardware

# Trial Loop
for trial in range(5):  # Run 5 trials as an example
    # Jittered fixation duration (0.5 to 1.5s)
    fixation_time = random.uniform(0.5, 1.5)

    # Send marker for trial start (optional)
    parallel_port.setData(1)  # Send marker value (1), which is a signal for trial start
    core.wait(0.01)  # Ensure the marker is sent properly
    parallel_port.setData(0)  # Reset to 0 after sending the marker

    # Show fixation cross
    fixation.draw()
    win.flip()
    core.wait(fixation_time)

    # Check for ESC key
    if 'escape' in event.getKeys():
        win.close()
        core.quit()

    # Get PTB clock time at trial onset
    trial_start_time = ptb.GetSecs()

    # Schedule sound playback exactly 0.5s after trial onset
    sound_onset_time = trial_start_time + 0.5
    probe_stim.play(when=sound_onset_time)

    # Send marker for sound onset (0.5s after trial onset)
    parallel_port.setData(2)  # Send marker value (2), which signifies sound onset
    core.wait(0.01)  # Wait for a short time to ensure the marker is sent
    parallel_port.setData(0)  # Reset to 0 after sending the marker

    # Show the initial circles
    inner_circle.draw()
    outer_circle.draw()
    win.flip()
    core.wait(1.0)  # Keep them displayed for 1s before shrinking starts

    # Animate shrinking of the outer circle over 2 seconds
    shrink_start_time = ptb.GetSecs()
    shrink_duration = 2.0
    while ptb.GetSecs() - shrink_start_time < shrink_duration:
        elapsed_time = ptb.GetSecs() - shrink_start_time
        new_radius = 0.1 - (elapsed_time / shrink_duration) * (0.1 - 0.05)  # Linearly interpolate radius
        
        # Update and draw circles
        outer_circle.radius = new_radius
        inner_circle.draw()
        outer_circle.draw()
        win.flip()

        # Check for ESC key during animation
        if 'escape' in event.getKeys():
            win.close()
            core.quit()

    # Fill the outer circle with green once shrinking completes
    outer_circle.fillColor = 'green'
    inner_circle.draw()
    outer_circle.draw()
    win.flip()
    core.wait(2.0)  # Keep the final green circle visible for 2s

    # Send marker for circle fully shrunk (filled with green)
    parallel_port.setData(3)  # Send marker value (3), signifying green circle filled
    core.wait(0.01)
    parallel_port.setData(0)  # Reset to 0 after sending the marker

    # Reset outer circle for next trial
    outer_circle.radius = 0.1
    outer_circle.fillColor = None

    # Check for ESC key before next trial
    if 'escape' in event.getKeys():
        win.close()
        core.quit()

# Cleanup
win.close()
core.quit()
