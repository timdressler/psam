form: "Hiss filter (low-pass)"
	real: "Frequency (1000-20000 Hz)", "7500"
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

frequency = min(max(frequency, 1000), 20000)

include batch.praat

procedure action
	s = selected("Sound")
	s$ = selected$("Sound")
	int = Get intensity (dB)

	if int <> undefined

include preview1.inc

		wrk = Copy: "wrk"
		runScript: "fixdc.praat"
		Filter (pass Hann band): 0, frequency, 100
		runScript: "fixdc.praat"
		result = selected("Sound")
		removeObject: wrk

include preview2.inc

		if not preview
			Rename: s$ + "-hissfilter_" + string$(frequency)
		endif
	else
		if not preview
			Copy: s$ + "-hissfilter_" + string$(frequency)
		endif
	endif
endproc
