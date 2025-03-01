form: "Copy duration"
	choice: "Method", 1
		option: "Stretch"
		option: "Cut or add time"
		option: "Change speed"
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

if numberOfSelected("Sound") = 2
	s1 = selected("Sound")
	s1$ = selected$("Sound")
	s2 = selected("Sound", 2)
	s2$ = selected$("Sound", 2)

	selectObject: s1
	dur1 = Get total duration

	selectObject: s2
	dur2 = Get total duration

	wrk = Copy: "wrk"

	if dur1 <> dur2
		if method = 1
			runScript: "changeduration.praat", dur1, "Stretch", 0
		elsif method = 2
			runScript: "fixdc.praat"
			Extract part: 0, dur1, "rectangular", 1, "no"
		elsif method = 3
			runScript: "changespeed.praat", "New duration", 0.5, dur1, 25, 23.976, 0
		endif
		removeObject: wrk
	endif
	result = selected("Sound")

	if preview
include preview.inc
		selectObject: s1, s2
		removeObject: trimmed, pre, result
	else
		Rename: s2$ + "-copyduration-" + s1$
	endif
endif
