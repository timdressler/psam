# Adapted from the script "Add2_variable" by Chris Darwin, https://groups.io/g/Praat-Users-List/files/Darwin%20scripts

form: "Mix"
	real: "Mix (%)", "50"
	boolean: "Avoid clipping", 0
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

if numberOfSelected("Sound") = 2
	mix = min(max(mix, 0), 100)
	amp1 = min((100 - mix) / 50, 1)
	amp2 = min(mix / 50, 1)

	s1 = selected("Sound")
	s1$ = selected$("Sound")
	s2 = selected("Sound", 2)
	s2$ = selected$("Sound", 2)

	sf1 = 1 / object[s1].dx
	sf2 = 1 / object[s2].dx
	max_sf = max(sf1, sf2)

	selectObject: s1
	dur1 = Get total duration
	ch1 = Get number of channels

	if sf1 = max_sf
		wrk1 = Copy: "wrk1"
	else
		wrk1 = Resample: max_sf, 1
	endif

	runScript: "fixdc.praat"

	selectObject: s2
	dur2 = Get total duration
	ch2 = Get number of channels

	if sf2 = max_sf
		wrk2 = Copy: "wrk2"
	else
		wrk2 = Resample: max_sf, 1
	endif

	runScript: "fixdc.praat"
	max_dur = max(dur1, dur2)
	max_ch = max(ch1, ch2)

	Create Sound from formula: "tmp", max_ch, 0, max_dur, max_sf, "(object[wrk1] * amp1) + (object[wrk2] * amp2)"
	if avoid_clipping
		runScript: "declip.praat"
	endif
	result = selected("Sound")

	removeObject: wrk1, wrk2

	if preview
include preview.inc
		selectObject: s1, s2
		removeObject: trimmed, pre, result
	else
		Rename: s2$ + "-mix-" + s1$
	endif
endif
