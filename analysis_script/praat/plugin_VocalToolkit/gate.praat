form: "Gate"
	real: "Silence threshold (dB)", "-35"
	real: "Minimum silent interval duration (s)", "0.1"
	real: "Minimum sounding interval duration (s)", "0.05"
	boolean: "Mark silences in a TextGrid", 0
	boolean: "Create extra Sound with inverted gate", 0
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

if silence_threshold >= 0
	exitScript: "“Silence threshold” should be a negative number."
endif

include batch.praat

procedure action
	s = selected("Sound")
	s$ = selected$("Sound")
	dur = Get total duration

include preview1.inc

	runScript: "workpre.praat"
	wrk = selected("Sound")

include minmaxf0.praat

	tg_wrk = noprogress To TextGrid (silences): 100, 0, silence_threshold, minimum_silent_interval_duration, minimum_sounding_interval_duration, "silent", "sounding"

	n_intervals = Get number of intervals: 1
	for i to n_intervals
		int_label$[i] = Get label of interval: 1, i
		start_time[i] = Get start time of interval: 1, i
		end_time[i] = Get end time of interval: 1, i
	endfor

	if create_extra_Sound_with_inverted_gate and not preview
		selectObject: wrk
		silent_wrk = Copy: "silent_wrk"
		for i to n_intervals
			if int_label$[i] = "sounding"
				if start_time[i] + 0.005 < end_time[i] - 0.005
					Set part to zero: start_time[i] + 0.005, end_time[i] - 0.005, "at exactly these times"
				endif
			elsif int_label$[i] = "silent"
				Fade in: 0, start_time[i] - 0.005, 0.01, "no"
				Fade out: 0, end_time[i] + 0.005, -0.01, "no"
			endif
		endfor
		runScript: "workpost.praat", dur
		Rename: s$ + "-invertedgate_" + string$(silence_threshold)
		removeObject: silent_wrk
	endif

	selectObject: wrk
	sounding_wrk = Copy: "sounding_wrk"
	for i to n_intervals
		if int_label$[i] = "silent"
			if start_time[i] + 0.005 < end_time[i] - 0.005
				Set part to zero: start_time[i] + 0.005, end_time[i] - 0.005, "at exactly these times"
			endif
		elsif int_label$[i] = "sounding"
			Fade in: 0, start_time[i] - 0.005, 0.01, "no"
			Fade out: 0, end_time[i] + 0.005, -0.01, "no"
		endif
	endfor

	runScript: "workpost.praat", dur
	result = selected("Sound")
	removeObject: wrk, sounding_wrk

include preview2.inc

	if not preview
		Rename: s$ + "-gate_" + string$(silence_threshold)
	endif

	if mark_silences_in_a_TextGrid and not preview
		selectObject: s
		marksilences = Copy: s$ + "-marksilences"
		runScript: "fixdc.praat"

		tg = Create TextGrid: 0, dur, "silences", ""
		Rename: s$ + "-marksilences"
		int_n = 0
		for i to n_intervals - 1
			if end_time[i] - 0.025 > 0 and end_time[i] - 0.025 < dur
				Insert boundary: 1, end_time[i] - 0.025
				int_n += 1
				Set interval text: 1, int_n, int_label$[i]
			endif
		endfor
		int_n += 1
		stt_int = Get start time of interval: 1, int_n

		selectObject: tg_wrk
		int_tim = Get interval at time: 1, stt_int + 0.025
		last_label$ = Get label of interval: 1, int_tim

		selectObject: tg
		Set interval text: 1, int_n, last_label$
		plusObject: marksilences
		View & Edit

		selectObject: result
	endif

	removeObject: tg_wrk
endproc
