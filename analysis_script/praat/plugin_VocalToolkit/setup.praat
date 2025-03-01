# Praat Vocal Toolkit. Version 2024.12.08
# A Praat plugin with automated scripts for voice processing
# https://www.praatvocaltoolkit.com
# This plugin is open source and can be used for any purpose

if praatVersion < 6420
	beginPause: "Praat Vocal Toolkit - Unsupported Praat version"
		comment: "“Praat Vocal Toolkit” requires Praat version 6.4.20 or higher."
		comment: "Your version of Praat is " + praatVersion$ + ". Please update it at ""https://www.praat.org""."
	endPause: "OK", 1, 1
else
	runScript: "buttons.praat"
	runScript: "updatepresetslists.praat"
endif
