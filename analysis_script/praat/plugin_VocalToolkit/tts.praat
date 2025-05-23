form: "Text to Speech (eSpeak)"
include ttspresetslist.inc
	positive: "Sampling frequency (Hz)", "44100"
	real: "Gap between words (s)", "0.01"
	real: "Pitch multiplier (0.5-2.0)", "1.0"
	real: "Pitch range multiplier (0-2.0)", "1.0"
	real: "Words per minute (80-450)", "175"
	boolean: "Create TextGrid with annotations", 0
	text: "Text", "1 2 3 4 5"
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

s# = selected# ("Sound")
gap_between_words = min(max(gap_between_words, 0), 1)
pitch_multiplier = min(max(pitch_multiplier, 0.5), 2)
pitch_range_multiplier = min(max(pitch_range_multiplier, 0), 2)
words_per_minute = min(max(words_per_minute, 80), 450)

ss = Create SpeechSynthesizer: language$, voice$
Speech output settings: sampling_frequency, gap_between_words, pitch_multiplier, pitch_range_multiplier, words_per_minute, "IPA"

if preview
	Play text: text$
	Remove
	selectObject: s#
else
	ss$ = selected$("SpeechSynthesizer")
	Rename: "tts_" + ss$
	noprogress To Sound: text$, create_TextGrid_with_annotations

	removeObject: ss

	if create_TextGrid_with_annotations
		tg = selected("TextGrid")
		result = selected("Sound")
		selectObject: result
		runScript: "declip.praat"
		plusObject: tg
		View & Edit
		selectObject: result
	else
		runScript: "declip.praat"
	endif
endif
