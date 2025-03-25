from psychopy import prefs
prefs.hardware['audioLib'] = ['PTB']  # Use PTB for precise audio timing
prefs.hardware['audioLatencyMode'] = 3  # High-precision mode

from psychopy import sound, core, parallel
import psychtoolbox as ptb

# Initialize PTB sound system
probe_stim = sound.Sound('C:/Experiments/tid_psam/psam/testing/test_rect_44100fs_80ms_2ch.wav', stereo=True, sampleRate=44100)  # Load a custom sound file

# Initialize Parallel Port
parallel_port = parallel.ParallelPort(address=0x3ff8)

# Trial Loop
for trial in range(20):  # Run 20 trials as an example
    # Get PTB clock time at trial onset
    trial_start_time = ptb.GetSecs()

    # Schedule sound playback exactly 0.5s after trial onset
    sound_onset_time = trial_start_time + 0.5

    # Wait for 0.5s relative to trial start, then send the trigger
    core.wait(0.5)  # Wait for 0.5s relative to trial start
    parallel_port.setData(2)  # Send trigger value (2), signifying sound onset
    core.wait(0.01)  # Wait for a short time to ensure the marker is sent
    parallel_port.setData(0)  # Reset the parallel port after sending the marker

    # Play sound at precise timing
    probe_stim.play(when=sound_onset_time)

    # Wait to ensure sound plays before ending trial
    core.wait(1.0)  # Adjust if necessary for the duration of your sound

# Cleanup
core.quit()
