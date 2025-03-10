#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2024.2.4),
    on März 10, 2025, at 14:52
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

# --- Import packages ---
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

# Run 'Before Experiment' code from code_trial
import random
# --- Setup global variables (available in all functions) ---
# create a device manager to handle hardware (keyboards, mice, mirophones, speakers, etc.)
deviceManager = hardware.DeviceManager()
# ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
# store info about the experiment session
psychopyVersion = '2024.2.4'
expName = 'tid_psam_main_experiment'  # from the Builder filename that created this script
# information about this experiment
expInfo = {
    'participant': '',
    'date|hid': data.getDateStr(),
    'expName|hid': expName,
    'psychopyVersion|hid': psychopyVersion,
}

# --- Define some variables which will change depending on pilot mode ---
'''
To run in pilot mode, either use the run/pilot toggle in Builder, Coder and Runner, 
or run the experiment with `--pilot` as an argument. To change what pilot 
#mode does, check out the 'Pilot mode' tab in preferences.
'''
# work out from system args whether we are running in pilot mode
PILOTING = core.setPilotModeFromArgs()
# start off with values from experiment settings
_fullScr = True
_winSize = [2560, 1440]
# if in pilot mode, apply overrides according to preferences
if PILOTING:
    # force windowed mode
    if prefs.piloting['forceWindowed']:
        _fullScr = False
        # set window size
        _winSize = prefs.piloting['forcedWindowSize']

def showExpInfoDlg(expInfo):
    """
    Show participant info dialog.
    Parameters
    ==========
    expInfo : dict
        Information about this experiment.
    
    Returns
    ==========
    dict
        Information about this experiment.
    """
    # show participant info dialog
    dlg = gui.DlgFromDict(
        dictionary=expInfo, sortKeys=False, title=expName, alwaysOnTop=True
    )
    if dlg.OK == False:
        core.quit()  # user pressed cancel
    # return expInfo
    return expInfo


def setupData(expInfo, dataDir=None):
    """
    Make an ExperimentHandler to handle trials and saving.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    dataDir : Path, str or None
        Folder to save the data to, leave as None to create a folder in the current directory.    
    Returns
    ==========
    psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    # remove dialog-specific syntax from expInfo
    for key, val in expInfo.copy().items():
        newKey, _ = data.utils.parsePipeSyntax(key)
        expInfo[newKey] = expInfo.pop(key)
    
    # data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
    if dataDir is None:
        dataDir = _thisDir
    filename = u'../data/BIDS/%s/beh/%s_%s_%s' % (f"sub-{int(expInfo['participant']):02d}", f"sub-{int(expInfo['participant']):02d}", "tid_psam_main_experiment", expInfo['date'])
    # make sure filename is relative to dataDir
    if os.path.isabs(filename):
        dataDir = os.path.commonprefix([dataDir, filename])
        filename = os.path.relpath(filename, dataDir)
    
    # an ExperimentHandler isn't essential but helps with data saving
    thisExp = data.ExperimentHandler(
        name=expName, version='',
        extraInfo=expInfo, runtimeInfo=None,
        originPath='C:\\Users\\timdr\\OneDrive\\Uni_Oldenburg\\4_Semester\\Master_Thesis\\Analysis_Experiment\\psam\\experiment_script\\tid_psam_main_experiment.py',
        savePickle=True, saveWideText=True,
        dataFileName=dataDir + os.sep + filename, sortColumns='alphabetical'
    )
    thisExp.setPriority('thisRow.t', priority.CRITICAL)
    thisExp.setPriority('expName', priority.LOW)
    # return experiment handler
    return thisExp


def setupLogging(filename):
    """
    Setup a log file and tell it what level to log at.
    
    Parameters
    ==========
    filename : str or pathlib.Path
        Filename to save log file and data files as, doesn't need an extension.
    
    Returns
    ==========
    psychopy.logging.LogFile
        Text stream to receive inputs from the logging system.
    """
    # set how much information should be printed to the console / app
    if PILOTING:
        logging.console.setLevel(
            prefs.piloting['pilotConsoleLoggingLevel']
        )
    else:
        logging.console.setLevel('warning')
    # save a log file for detail verbose info
    logFile = logging.LogFile(filename+'.log')
    if PILOTING:
        logFile.setLevel(
            prefs.piloting['pilotLoggingLevel']
        )
    else:
        logFile.setLevel(
            logging.getLevel('info')
        )
    
    return logFile


def setupWindow(expInfo=None, win=None):
    """
    Setup the Window
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    win : psychopy.visual.Window
        Window to setup - leave as None to create a new window.
    
    Returns
    ==========
    psychopy.visual.Window
        Window in which to run this experiment.
    """
    if PILOTING:
        logging.debug('Fullscreen settings ignored as running in pilot mode.')
    
    if win is None:
        # if not given a window to setup, make one
        win = visual.Window(
            size=_winSize, fullscr=_fullScr, screen=0,
            winType='pyglet', allowGUI=False, allowStencil=False,
            monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
            backgroundImage='', backgroundFit='none',
            blendMode='avg', useFBO=True,
            units='height',
            checkTiming=False  # we're going to do this ourselves in a moment
        )
    else:
        # if we have a window, just set the attributes which are safe to set
        win.color = [0,0,0]
        win.colorSpace = 'rgb'
        win.backgroundImage = ''
        win.backgroundFit = 'none'
        win.units = 'height'
    if expInfo is not None:
        # get/measure frame rate if not already in expInfo
        if win._monitorFrameRate is None:
            win._monitorFrameRate = win.getActualFrameRate(infoMsg='Attempting to measure frame rate of screen, please wait...')
        expInfo['frameRate'] = win._monitorFrameRate
    win.hideMessage()
    # show a visual indicator if we're in piloting mode
    if PILOTING and prefs.piloting['showPilotingIndicator']:
        win.showPilotingIndicator()
    
    return win


def setupDevices(expInfo, thisExp, win):
    """
    Setup whatever devices are available (mouse, keyboard, speaker, eyetracker, etc.) and add them to 
    the device manager (deviceManager)
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window in which to run this experiment.
    Returns
    ==========
    bool
        True if completed successfully.
    """
    # --- Setup input devices ---
    ioConfig = {}
    
    # Setup iohub keyboard
    ioConfig['Keyboard'] = dict(use_keymap='psychopy')
    
    # Setup iohub experiment
    ioConfig['Experiment'] = dict(filename=thisExp.dataFileName)
    
    # Start ioHub server
    ioServer = io.launchHubServer(window=win, **ioConfig)
    
    # store ioServer object in the device manager
    deviceManager.ioServer = ioServer
    
    # create a default keyboard (e.g. to check for escape)
    if deviceManager.getDevice('defaultKeyboard') is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='iohub'
        )
    if deviceManager.getDevice('instruction_keys') is None:
        # initialise instruction_keys
        instruction_keys = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='instruction_keys',
        )
    # create speaker 'probe_stim'
    deviceManager.addDevice(
        deviceName='probe_stim',
        deviceClass='psychopy.hardware.speaker.SpeakerDevice',
        index=-1
    )
    # initialise microphone
    deviceManager.addDevice(
        deviceClass='psychopy.hardware.microphone.MicrophoneDevice',
        deviceName='audiointerface',
        index=12,
        maxRecordingSize=240000.0,
        channels=1, 
        sampleRateHz=48000, 
    )
    # return True if completed successfully
    return True

def pauseExperiment(thisExp, win=None, timers=[], playbackComponents=[]):
    """
    Pause this experiment, preventing the flow from advancing to the next routine until resumed.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window for this experiment.
    timers : list, tuple
        List of timers to reset once pausing is finished.
    playbackComponents : list, tuple
        List of any components with a `pause` method which need to be paused.
    """
    # if we are not paused, do nothing
    if thisExp.status != PAUSED:
        return
    
    # start a timer to figure out how long we're paused for
    pauseTimer = core.Clock()
    # pause any playback components
    for comp in playbackComponents:
        comp.pause()
    # make sure we have a keyboard
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        defaultKeyboard = deviceManager.addKeyboard(
            deviceClass='keyboard',
            deviceName='defaultKeyboard',
            backend='ioHub',
        )
    # run a while loop while we wait to unpause
    while thisExp.status == PAUSED:
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=['escape']):
            endExperiment(thisExp, win=win)
        # sleep 1ms so other threads can execute
        clock.time.sleep(0.001)
    # if stop was requested while paused, quit
    if thisExp.status == FINISHED:
        endExperiment(thisExp, win=win)
    # resume any playback components
    for comp in playbackComponents:
        comp.play()
    # reset any timers
    for timer in timers:
        timer.addTime(-pauseTimer.getTime())


def run(expInfo, thisExp, win, globalClock=None, thisSession=None):
    """
    Run the experiment flow.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    psychopy.visual.Window
        Window in which to run this experiment.
    globalClock : psychopy.core.clock.Clock or None
        Clock to get global time from - supply None to make a new one.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    # mark experiment as started
    thisExp.status = STARTED
    # make sure window is set to foreground to prevent losing focus
    win.winHandle.activate()
    # make sure variables created by exec are available globally
    exec = environmenttools.setExecEnvironment(globals())
    # get device handles from dict of input devices
    ioServer = deviceManager.ioServer
    # get/create a default keyboard (e.g. to check for escape)
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='ioHub'
        )
    eyetracker = deviceManager.getDevice('eyetracker')
    # make sure we're running in the directory for this experiment
    os.chdir(_thisDir)
    # get filename from ExperimentHandler for convenience
    filename = thisExp.dataFileName
    frameTolerance = 0.001  # how close to onset before 'same' frame
    endExpNow = False  # flag for 'escape' or other condition => quit the exp
    # get frame duration from frame rate in expInfo
    if 'frameRate' in expInfo and expInfo['frameRate'] is not None:
        frameDur = 1.0 / round(expInfo['frameRate'])
    else:
        frameDur = 1.0 / 60.0  # could not measure, so guess
    
    # Start Code - component code to be run after the window creation
    # Make folder to store recordings from mic
    micRecFolder = filename + '_mic_recorded'
    if not os.path.isdir(micRecFolder):
        os.mkdir(micRecFolder)
    
    # --- Initialize components for Routine "welcome" ---
    welcome_message = visual.TextStim(win=win, name='welcome_message',
        text='Welcome to this experiment',
        font='Arial',
        pos=(0, 0), draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    # Run 'Begin Experiment' code from code_setup
    subject_id = expInfo["participant"] 
    subject_id = f"sub-{subject_id.zfill(2)}"
    thisExp.addData("subject_id", subject_id)
    
    _thisDir2 = os.path.dirname(_thisDir)
    
    filename_cond = u'data\\BIDS\\stimuli\\%s\\%s_conditions.xlsx' % (subject_id, subject_id) # CHANGED
    conditionsFileName = _thisDir2 + os.sep + filename_cond
    thisExp.addData("conditions", conditionsFileName)
    
    
    # --- Initialize components for Routine "instructions" ---
    instructions_message = visual.TextStim(win=win, name='instructions_message',
        text='In this experiment you have two different tasks.\n\nIf the instruction says "Passive" you only have to fixate the moving circle in the middle of the sceen, without doing or saying anything.\n\nIf the instruction says "Active" you are required to produce a short sound (/da/) AFTER the moving cricle reached the size of the inner circle AND the  circle turned green.\n\nIn some trials a sound will be played BEFORE the circle turned green. In this case, still wait until the circle turns green and vocalize /da/ then (in case of Active trials). Try your best to ignore the sound by solely focusing on the circles.',
        font='Arial',
        pos=(0, 0), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    instruction_keys = keyboard.Keyboard(deviceName='instruction_keys')
    
    # --- Initialize components for Routine "ISI" ---
    fixation_cross_ISI = visual.ShapeStim(
        win=win, name='fixation_cross_ISI', vertices='cross',
        size=(0.1, 0.1),
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor='white', fillColor='white',
        opacity=None, depth=0.0, interpolate=True)
    
    # --- Initialize components for Routine "trial" ---
    cue_stim = visual.ShapeStim(
        win=win, name='cue_stim',
        size=(0.7, 0.7), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=4.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    target_stim = visual.ShapeStim(
        win=win, name='target_stim',
        size=(0.2, 0.2), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=4.0,
        colorSpace='rgb', lineColor=[-1.0000, 0.0039, -1.0000], fillColor=[0.0000, 0.0000, 0.0000],
        opacity=1.0, depth=-1.0, interpolate=True)
    probe_stim = sound.Sound(
        'A', 
        secs=-1, 
        stereo=True, 
        hamming=True, 
        speaker='probe_stim',    name='probe_stim'
    )
    probe_stim.setVolume(1.0)
    task_msg = visual.TextStim(win=win, name='task_msg',
        text='',
        font='Arial',
        pos=(0, 0), draggable=False, height=0.04, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-3.0);
    # Run 'Begin Experiment' code from code_trial
    # Import necessary module
    from psychopy import core
    
    # Initialize parameters for the cue
    initial_radius = 0.7  # Starting size of the cue
    final_radius = 0.2    # Final size of the cue
    duration = 2.0        # Time in seconds to shrink
    hold_time = 1.5       # Time in seconds to hold after reaching final size
    
    cue_started = False    # To track when the cue starts shrinking
    cue_start_time = None  # To store when the cue started shrinking
    
    # make microphone object for mic
    mic = sound.microphone.Microphone(
        device='audiointerface',
        name='mic',
        recordingFolder=micRecFolder,
        recordingExt='wav'
    )
    
    # create some handy timers
    
    # global clock to track the time since experiment started
    if globalClock is None:
        # create a clock if not given one
        globalClock = core.Clock()
    if isinstance(globalClock, str):
        # if given a string, make a clock accoridng to it
        if globalClock == 'float':
            # get timestamps as a simple value
            globalClock = core.Clock(format='float')
        elif globalClock == 'iso':
            # get timestamps in ISO format
            globalClock = core.Clock(format='%Y-%m-%d_%H:%M:%S.%f%z')
        else:
            # get timestamps in a custom format
            globalClock = core.Clock(format=globalClock)
    if ioServer is not None:
        ioServer.syncClock(globalClock)
    logging.setDefaultClock(globalClock)
    # routine timer to track time remaining of each (possibly non-slip) routine
    routineTimer = core.Clock()
    win.flip()  # flip window to reset last flip timer
    # store the exact time the global clock started
    expInfo['expStart'] = data.getDateStr(
        format='%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6
    )
    
    # --- Prepare to start Routine "welcome" ---
    # create an object to store info about Routine welcome
    welcome = data.Routine(
        name='welcome',
        components=[welcome_message],
    )
    welcome.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # store start times for welcome
    welcome.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    welcome.tStart = globalClock.getTime(format='float')
    welcome.status = STARTED
    thisExp.addData('welcome.started', welcome.tStart)
    welcome.maxDuration = None
    # keep track of which components have finished
    welcomeComponents = welcome.components
    for thisComponent in welcome.components:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "welcome" ---
    welcome.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 1.5:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *welcome_message* updates
        
        # if welcome_message is starting this frame...
        if welcome_message.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            welcome_message.frameNStart = frameN  # exact frame index
            welcome_message.tStart = t  # local t and not account for scr refresh
            welcome_message.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(welcome_message, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'welcome_message.started')
            # update status
            welcome_message.status = STARTED
            welcome_message.setAutoDraw(True)
        
        # if welcome_message is active this frame...
        if welcome_message.status == STARTED:
            # update params
            pass
        
        # if welcome_message is stopping this frame...
        if welcome_message.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > welcome_message.tStartRefresh + 1.5-frameTolerance:
                # keep track of stop time/frame for later
                welcome_message.tStop = t  # not accounting for scr refresh
                welcome_message.tStopRefresh = tThisFlipGlobal  # on global time
                welcome_message.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'welcome_message.stopped')
                # update status
                welcome_message.status = FINISHED
                welcome_message.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
            )
            # skip the frame we paused on
            continue
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            welcome.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in welcome.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "welcome" ---
    for thisComponent in welcome.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for welcome
    welcome.tStop = globalClock.getTime(format='float')
    welcome.tStopRefresh = tThisFlipGlobal
    thisExp.addData('welcome.stopped', welcome.tStop)
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if welcome.maxDurationReached:
        routineTimer.addTime(-welcome.maxDuration)
    elif welcome.forceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-1.500000)
    thisExp.nextEntry()
    
    # --- Prepare to start Routine "instructions" ---
    # create an object to store info about Routine instructions
    instructions = data.Routine(
        name='instructions',
        components=[instructions_message, instruction_keys],
    )
    instructions.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # create starting attributes for instruction_keys
    instruction_keys.keys = []
    instruction_keys.rt = []
    _instruction_keys_allKeys = []
    # store start times for instructions
    instructions.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    instructions.tStart = globalClock.getTime(format='float')
    instructions.status = STARTED
    thisExp.addData('instructions.started', instructions.tStart)
    instructions.maxDuration = None
    # keep track of which components have finished
    instructionsComponents = instructions.components
    for thisComponent in instructions.components:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "instructions" ---
    instructions.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *instructions_message* updates
        
        # if instructions_message is starting this frame...
        if instructions_message.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            instructions_message.frameNStart = frameN  # exact frame index
            instructions_message.tStart = t  # local t and not account for scr refresh
            instructions_message.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(instructions_message, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'instructions_message.started')
            # update status
            instructions_message.status = STARTED
            instructions_message.setAutoDraw(True)
        
        # if instructions_message is active this frame...
        if instructions_message.status == STARTED:
            # update params
            pass
        
        # *instruction_keys* updates
        waitOnFlip = False
        
        # if instruction_keys is starting this frame...
        if instruction_keys.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            instruction_keys.frameNStart = frameN  # exact frame index
            instruction_keys.tStart = t  # local t and not account for scr refresh
            instruction_keys.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(instruction_keys, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'instruction_keys.started')
            # update status
            instruction_keys.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(instruction_keys.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(instruction_keys.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if instruction_keys.status == STARTED and not waitOnFlip:
            theseKeys = instruction_keys.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
            _instruction_keys_allKeys.extend(theseKeys)
            if len(_instruction_keys_allKeys):
                instruction_keys.keys = _instruction_keys_allKeys[-1].name  # just the last key pressed
                instruction_keys.rt = _instruction_keys_allKeys[-1].rt
                instruction_keys.duration = _instruction_keys_allKeys[-1].duration
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
            )
            # skip the frame we paused on
            continue
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            instructions.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in instructions.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "instructions" ---
    for thisComponent in instructions.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for instructions
    instructions.tStop = globalClock.getTime(format='float')
    instructions.tStopRefresh = tThisFlipGlobal
    thisExp.addData('instructions.stopped', instructions.tStop)
    # check responses
    if instruction_keys.keys in ['', [], None]:  # No response was made
        instruction_keys.keys = None
    thisExp.addData('instruction_keys.keys',instruction_keys.keys)
    if instruction_keys.keys != None:  # we had a response
        thisExp.addData('instruction_keys.rt', instruction_keys.rt)
        thisExp.addData('instruction_keys.duration', instruction_keys.duration)
    thisExp.nextEntry()
    # the Routine "instructions" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler2(
        name='trials',
        nReps=1.0, 
        method='sequential', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=data.importConditions(conditionsFileName), 
        seed=None, 
    )
    thisExp.addLoop(trials)  # add the loop to the experiment
    thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial:
            globals()[paramName] = thisTrial[paramName]
    if thisSession is not None:
        # if running in a Session with a Liaison client, send data up to now
        thisSession.sendExperimentData()
    
    for thisTrial in trials:
        currentLoop = trials
        thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
        if thisTrial != None:
            for paramName in thisTrial:
                globals()[paramName] = thisTrial[paramName]
        
        # --- Prepare to start Routine "ISI" ---
        # create an object to store info about Routine ISI
        ISI = data.Routine(
            name='ISI',
            components=[fixation_cross_ISI],
        )
        ISI.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # store start times for ISI
        ISI.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        ISI.tStart = globalClock.getTime(format='float')
        ISI.status = STARTED
        thisExp.addData('ISI.started', ISI.tStart)
        ISI.maxDuration = random.randint(500, 1500)/1000
        # keep track of which components have finished
        ISIComponents = ISI.components
        for thisComponent in ISI.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "ISI" ---
        # if trial has changed, end Routine now
        if isinstance(trials, data.TrialHandler2) and thisTrial.thisN != trials.thisTrial.thisN:
            continueRoutine = False
        ISI.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # is it time to end the Routine? (based on local clock)
            if tThisFlip > ISI.maxDuration-frameTolerance:
                ISI.maxDurationReached = True
                continueRoutine = False
            
            # *fixation_cross_ISI* updates
            
            # if fixation_cross_ISI is starting this frame...
            if fixation_cross_ISI.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixation_cross_ISI.frameNStart = frameN  # exact frame index
                fixation_cross_ISI.tStart = t  # local t and not account for scr refresh
                fixation_cross_ISI.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixation_cross_ISI, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixation_cross_ISI.started')
                # update status
                fixation_cross_ISI.status = STARTED
                fixation_cross_ISI.setAutoDraw(True)
            
            # if fixation_cross_ISI is active this frame...
            if fixation_cross_ISI.status == STARTED:
                # update params
                pass
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                ISI.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in ISI.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "ISI" ---
        for thisComponent in ISI.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for ISI
        ISI.tStop = globalClock.getTime(format='float')
        ISI.tStopRefresh = tThisFlipGlobal
        thisExp.addData('ISI.stopped', ISI.tStop)
        # the Routine "ISI" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "trial" ---
        # create an object to store info about Routine trial
        trial = data.Routine(
            name='trial',
            components=[cue_stim, target_stim, probe_stim, task_msg, mic],
        )
        trial.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        probe_stim.setSound(stim_file, secs=probe_duration, hamming=True)
        probe_stim.setVolume(probe_intensity, log=False)
        probe_stim.seek(0)
        task_msg.setText(task)
        # Run 'Begin Routine' code from code_trial
        cue_started = False
        aa_green_onset = None  # To store when the cue turns green relative to trial onset
        aa_test_onset = None
        
        # store start times for trial
        trial.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        trial.tStart = globalClock.getTime(format='float')
        trial.status = STARTED
        thisExp.addData('trial.started', trial.tStart)
        trial.maxDuration = 5
        # keep track of which components have finished
        trialComponents = trial.components
        for thisComponent in trial.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "trial" ---
        # if trial has changed, end Routine now
        if isinstance(trials, data.TrialHandler2) and thisTrial.thisN != trials.thisTrial.thisN:
            continueRoutine = False
        trial.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 5.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # is it time to end the Routine? (based on local clock)
            if tThisFlip > trial.maxDuration-frameTolerance:
                trial.maxDurationReached = True
                continueRoutine = False
            
            # *cue_stim* updates
            
            # if cue_stim is starting this frame...
            if cue_stim.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                cue_stim.frameNStart = frameN  # exact frame index
                cue_stim.tStart = t  # local t and not account for scr refresh
                cue_stim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(cue_stim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'cue_stim.started')
                # update status
                cue_stim.status = STARTED
                cue_stim.setAutoDraw(True)
            
            # if cue_stim is active this frame...
            if cue_stim.status == STARTED:
                # update params
                pass
            
            # if cue_stim is stopping this frame...
            if cue_stim.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > cue_stim.tStartRefresh + 4.5-frameTolerance:
                    # keep track of stop time/frame for later
                    cue_stim.tStop = t  # not accounting for scr refresh
                    cue_stim.tStopRefresh = tThisFlipGlobal  # on global time
                    cue_stim.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'cue_stim.stopped')
                    # update status
                    cue_stim.status = FINISHED
                    cue_stim.setAutoDraw(False)
            
            # *target_stim* updates
            
            # if target_stim is starting this frame...
            if target_stim.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                target_stim.frameNStart = frameN  # exact frame index
                target_stim.tStart = t  # local t and not account for scr refresh
                target_stim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(target_stim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'target_stim.started')
                # update status
                target_stim.status = STARTED
                target_stim.setAutoDraw(True)
            
            # if target_stim is active this frame...
            if target_stim.status == STARTED:
                # update params
                pass
            
            # if target_stim is stopping this frame...
            if target_stim.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > target_stim.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    target_stim.tStop = t  # not accounting for scr refresh
                    target_stim.tStopRefresh = tThisFlipGlobal  # on global time
                    target_stim.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'target_stim.stopped')
                    # update status
                    target_stim.status = FINISHED
                    target_stim.setAutoDraw(False)
            
            # *probe_stim* updates
            
            # if probe_stim is starting this frame...
            if probe_stim.status == NOT_STARTED and tThisFlip >= probe_onset-frameTolerance:
                # keep track of start time/frame for later
                probe_stim.frameNStart = frameN  # exact frame index
                probe_stim.tStart = t  # local t and not account for scr refresh
                probe_stim.tStartRefresh = tThisFlipGlobal  # on global time
                # add timestamp to datafile
                thisExp.addData('probe_stim.started', tThisFlipGlobal)
                # update status
                probe_stim.status = STARTED
                probe_stim.play(when=win)  # sync with win flip
            
            # if probe_stim is stopping this frame...
            if probe_stim.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > probe_stim.tStartRefresh + probe_duration-frameTolerance or probe_stim.isFinished:
                    # keep track of stop time/frame for later
                    probe_stim.tStop = t  # not accounting for scr refresh
                    probe_stim.tStopRefresh = tThisFlipGlobal  # on global time
                    probe_stim.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'probe_stim.stopped')
                    # update status
                    probe_stim.status = FINISHED
                    probe_stim.stop()
            
            # *task_msg* updates
            
            # if task_msg is starting this frame...
            if task_msg.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                task_msg.frameNStart = frameN  # exact frame index
                task_msg.tStart = t  # local t and not account for scr refresh
                task_msg.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(task_msg, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'task_msg.started')
                # update status
                task_msg.status = STARTED
                task_msg.setAutoDraw(True)
            
            # if task_msg is active this frame...
            if task_msg.status == STARTED:
                # update params
                pass
            
            # if task_msg is stopping this frame...
            if task_msg.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > task_msg.tStartRefresh + 4.5-frameTolerance:
                    # keep track of stop time/frame for later
                    task_msg.tStop = t  # not accounting for scr refresh
                    task_msg.tStopRefresh = tThisFlipGlobal  # on global time
                    task_msg.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'task_msg.stopped')
                    # update status
                    task_msg.status = FINISHED
                    task_msg.setAutoDraw(False)
            # Run 'Each Frame' code from code_trial
            # Get elapsed time relative to trial onset
            elapsed_time = t  # 't' is automatically provided by PsychoPy
            
            # Ensure cue_stim is always visible from the start
            cue_stim.autoDraw = True  
            
            # Ensure cue_stim stays black and at full size before 1s
            if elapsed_time < 1.0:
                cue_stim.fillColor = [0, 0, 0]  # Keep it black
                cue_stim.size = (0.7, 0.7)  # Ensure initial size
            
            # Start shrinking cue_stim after 1s
            elif elapsed_time >= 1.0:
                if not cue_started:
                    cue_start_time = elapsed_time  # Record start time of shrinking
                    cue_started = True
            
                # Calculate how much time has passed since shrinking started
                time_passed = elapsed_time - cue_start_time
            
                # Compute new size
                new_radius = initial_radius - (time_passed / duration) * (initial_radius - final_radius)
                new_radius = max(new_radius, final_radius)  # Ensure it doesn’t go below final size
            
                # Update cue_stim properties
                cue_stim.size = (new_radius, new_radius)
            
                # Change color when cue reaches final size and record aa_green_onset
                if new_radius <= final_radius:
                    cue_stim.fillColor = [-1.0000, 0.0039, -1.0000]  # Change to green
                    if aa_green_onset is None:  # Record green onset time once
                        aa_green_onset = elapsed_time  
            
            # Remove cue_stim after the hold time
            if aa_green_onset is not None and elapsed_time >= aa_green_onset + hold_time:
                cue_stim.autoDraw = False  # Hide cue_stim
            
            # Track test_stim onset time without modifying its behavior
            if elapsed_time >= 3.0 and aa_test_onset is None:
                aa_test_onset = elapsed_time  # Store exact trial-relative onset time
            
            # Save aa_test_onset and aa_green_onset to output file at the end of the routine
            if aa_test_onset is not None:
                thisExp.addData('aa_test_onset', aa_test_onset)  # Store test_stim onset
            
            if aa_green_onset is not None:
                thisExp.addData('aa_green_onset', aa_green_onset)  # Store cue turning green onset
            
            
            # if mic is starting this frame...
            if mic.status == NOT_STARTED and t >= 2.5-frameTolerance:
                # keep track of start time/frame for later
                mic.frameNStart = frameN  # exact frame index
                mic.tStart = t  # local t and not account for scr refresh
                mic.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(mic, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('mic.started', t)
                # update status
                mic.status = STARTED
                # start recording with mic
                mic.start()
            
            # if mic is active this frame...
            if mic.status == STARTED:
                # update params
                pass
                # update recorded clip for mic
                mic.poll()
            
            # if mic is stopping this frame...
            if mic.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > mic.tStartRefresh + 1.5-frameTolerance:
                    # keep track of stop time/frame for later
                    mic.tStop = t  # not accounting for scr refresh
                    mic.tStopRefresh = tThisFlipGlobal  # on global time
                    mic.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('mic.stopped', t)
                    # update status
                    mic.status = FINISHED
                    # stop recording with mic
                    mic.stop()
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[probe_stim]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                trial.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trial.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "trial" ---
        for thisComponent in trial.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for trial
        trial.tStop = globalClock.getTime(format='float')
        trial.tStopRefresh = tThisFlipGlobal
        thisExp.addData('trial.stopped', trial.tStop)
        probe_stim.pause()  # ensure sound has stopped at end of Routine
        # tell mic to keep hold of current recording in mic.clips and transcript (if applicable) in mic.scripts
        # this will also update mic.lastClip and mic.lastScript
        mic.stop()
        tag = data.utils.getDateStr()
        micClip = mic.bank(
            tag=tag, transcribe='None',
            config=None
        )
        trials.addData(
            'mic.clip', mic.recordingFolder / mic.getClipFilename(tag)
        )
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if trial.maxDurationReached:
            routineTimer.addTime(-trial.maxDuration)
        elif trial.forceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-5.000000)
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'trials'
    
    if thisSession is not None:
        # if running in a Session with a Liaison client, send data up to now
        thisSession.sendExperimentData()
    # save mic recordings
    mic.saveClips()
    
    # mark experiment as finished
    endExperiment(thisExp, win=win)


def saveData(thisExp):
    """
    Save data from this experiment
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    filename = thisExp.dataFileName
    # these shouldn't be strictly necessary (should auto-save)
    thisExp.saveAsWideText(filename + '.csv', delim='auto')
    thisExp.saveAsPickle(filename)


def endExperiment(thisExp, win=None):
    """
    End this experiment, performing final shut down operations.
    
    This function does NOT close the window or end the Python process - use `quit` for this.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window for this experiment.
    """
    if win is not None:
        # remove autodraw from all current components
        win.clearAutoDraw()
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed
        win.flip()
    # return console logger level to WARNING
    logging.console.setLevel(logging.WARNING)
    # mark experiment handler as finished
    thisExp.status = FINISHED
    logging.flush()


def quit(thisExp, win=None, thisSession=None):
    """
    Fully quit, closing the window and ending the Python process.
    
    Parameters
    ==========
    win : psychopy.visual.Window
        Window to close.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    thisExp.abort()  # or data files will save again on exit
    # make sure everything is closed down
    if win is not None:
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed before quitting
        win.flip()
        win.close()
    logging.flush()
    if thisSession is not None:
        thisSession.stop()
    # terminate Python process
    core.quit()


# if running this experiment as a script...
if __name__ == '__main__':
    # call all functions in order
    expInfo = showExpInfoDlg(expInfo=expInfo)
    thisExp = setupData(expInfo=expInfo)
    logFile = setupLogging(filename=thisExp.dataFileName)
    win = setupWindow(expInfo=expInfo)
    setupDevices(expInfo=expInfo, thisExp=thisExp, win=win)
    run(
        expInfo=expInfo, 
        thisExp=thisExp, 
        win=win,
        globalClock='float'
    )
    saveData(thisExp=thisExp)
    quit(thisExp=thisExp, win=win)
