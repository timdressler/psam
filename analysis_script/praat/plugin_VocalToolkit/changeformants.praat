form: "Change formants"
	real: "New F1 mean (Hz)", "500.0"
	real: "New F2 mean (Hz)", "1500.0"
	real: "New F3 mean (Hz)", "2500.0"
	real: "New F4 mean (Hz)", "0 (= no change)"
	real: "New F5 mean (Hz)", "0 (= no change)"
	comment: "Formant determination"
	positive: "Maximum formant (Hz)", "5500 (= adult female)"
	comment: "Set 5000 Hz for men, 5500 Hz for women or up to 8000 Hz for children."
	boolean: "Process only voiced parts", 1
	boolean: "Retrieve intensity contour", 1
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

f1 = max(new_F1_mean, 0)
f2 = max(new_F2_mean, 0)
f3 = max(new_F3_mean, 0)
f4 = max(new_F4_mean, 0)
f5 = max(new_F5_mean, 0)

include batch.praat

procedure action
	s = selected("Sound")
	s$ = selected$("Sound")
	original_dur = Get total duration
	int = Get intensity (dB)

	if int <> undefined

include preview1.inc

		if process_only_voiced_parts
			@extractUV
			selectObject: extractUV.s_v
			int_v = Get intensity (dB)
			if int_v <> undefined
				int = int_v
			else
				process_only_voiced_parts = 0
				selectObject: s
				removeObject: extractUV.s_u, extractUV.s_v
			endif
		endif

		runScript: "workpre.praat"
		wrk = selected("Sound")
		sf1 = Get sampling frequency

		runScript: "extractvowels.praat", 0, 0
		vow_tmp = selected("Sound")

		runScript: "workpre.praat"
		vow = selected("Sound")

		formant1 = noprogress nowarn To Formant (robust): 0.005, 5, maximum_formant, 0.025, 50, 1.5, 5, 0.000001
		vf1 = Get mean: 1, 0, 0, "hertz"
		vf2 = Get mean: 2, 0, 0, "hertz"
		vf3 = Get mean: 3, 0, 0, "hertz"
		vf4 = Get mean: 4, 0, 0, "hertz"
		vf5 = Get mean: 5, 0, 0, "hertz"
		df1 = f1 - vf1
		df2 = f2 - vf2
		df3 = f3 - vf3
		df4 = f4 - vf4
		df5 = f5 - vf5

		selectObject: wrk
		hf = Filter (stop Hann band): 0, maximum_formant, 100

		selectObject: wrk
		sf2 = maximum_formant * 2
		rs1 = Resample: sf2, 10

		formant2 = noprogress nowarn To Formant (robust): 0.005, 5, maximum_formant, 0.025, 50, 1.5, 5, 0.000001

		lpc1 = noprogress To LPC: sf2
		plusObject: rs1
		source = Filter (inverse)

		selectObject: formant2
		filtr = Copy: "filtr"

		if f1 <> 0 and abs(df1) < 1000
			Formula (frequencies): "if row = 1 then self + df1 else self fi"
		endif
		if f2 <> 0 and abs(df2) < 2500
			Formula (frequencies): "if row = 2 then self + df2 else self fi"
		endif
		if f3 <> 0 and abs(df3) < 2500
			Formula (frequencies): "if row = 3 then self + df3 else self fi"
		endif
		if f4 <> 0 and abs(df4) < 2500
			Formula (frequencies): "if row = 4 then self + df4 else self fi"
		endif
		if f5 <> 0 and abs(df5) < 2500
			Formula (frequencies): "if row = 5 then self + df5 else self fi"
		endif

		lpc2 = noprogress To LPC: sf2
		plusObject: source
		tmp = Filter: "no"

		rs2 = Resample: sf1, 10
		Formula: "self + object[hf]"

		runScript: "workpost.praat", original_dur
		Scale intensity: int
		runScript: "declip.praat"

		if process_only_voiced_parts
			@mixUV
		endif

		if retrieve_intensity_contour
			tmp3 = selected("Sound")
			if preview
				plusObject: pre1
			else
				plusObject: s
			endif
			runScript: "copyintensitycontour.praat"
			removeObject: tmp3
		endif
		dur = Get total duration
		if dur > 0.5
			Fade in: 0, 0, 0.005, "yes"
			Fade out: 0, dur, -0.005, "yes"
		endif
		result = selected("Sound")

		removeObject: wrk, vow_tmp, vow, formant1, hf, rs1, formant2, lpc1, source, filtr, lpc2, tmp, rs2

include preview2.inc

		if not preview
			Rename: s$ + "-changeformants"
		endif
	else
		if not preview
			Copy: s$ + "-changeformants"
		endif
	endif

endproc

procedure extractUV
	runScript: "voicedunvoiced.praat", 0
	select all
	.s_u = selected("Sound", -2)
	.s_v = selected("Sound", -1)
endproc

procedure mixUV
	.sel_tmp = selected("Sound")
	plusObject: extractUV.s_u
	runScript: "copymix.praat", 50, 0, 0
	removeObject: extractUV.s_u, extractUV.s_v, .sel_tmp
endproc
