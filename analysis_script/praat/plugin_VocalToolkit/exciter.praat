# Adapted from Shekar, Priyanka & Smith, Julius. (2013). Modeling the Harmonic Exciter. https://www.researchgate.net/publication/258333577_Modeling_the_Harmonic_Exciter
# Knoppel, C. (1979). Signal distortion circuit and method of use. Patent No. US4150253 A. https://patents.google.com/patent/US4150253A/en

form: "Exciter"
	positive: "High pass frequency (Hz)", "2000"
	optionmenu: "Phase rotation", 3
		option: "No change"
		option: "+90º Hilbert transform"
		option: "-90º Hilbert transform"
		option: "180º waveform inversion"
	optionmenu: "Type of distortion", 6
		option: "Hard clipping"
		option: "Soft clipping. Quintic"
		option: "Soft clipping. Cubic"
		option: "Soft clipping. Hyperbolic tangent"
		option: "Soft clipping. Algebraic"
		option: "Soft clipping. Arctangent"
	real: "Input gain (dB)", "0.0"
	real: "Positive amplitude limit (Pa)", "1.0"
	real: "Negative amplitude limit (Pa)", "-0.3"
	real: "Mix (dry/wet balance, 0-100 %)", "50"
	boolean: "Scale result to original intensity", 0
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

if positive_amplitude_limit < 0
	exitScript: "“Positive amplitude limit” must be greater than or equal to 0.0"
endif
if negative_amplitude_limit > 0
	exitScript: "“Negative amplitude limit” must be less than or equal to 0.0"
endif

freq = high_pass_frequency
mix = min(max(mix, 0), 100)

include batch.praat

procedure action
	s = selected("Sound")
	s$ = selected$("Sound")
	ch = Get number of channels
	dur = Get total duration
	sf = Get sampling frequency
	int = Get intensity (dB)

	if int <> undefined

include preview1.inc

		nsf = max(sf, 48000)
		if sf <> nsf
			tmp1 = Resample: nsf, 50
		else
			tmp1 = Copy: "tmp1"
		endif
		runScript: "fixdc.praat"

		tmp2 = Copy: "tmp2"
		Scale peak: 0.99
		int2 = Get intensity (dB)
		Scale intensity: int2 + 2

		if phase_rotation = 1
			tmp3 = Copy: "tmp3"
		else
			runScript: "phaserotation.praat", phase_rotation$
			tmp3 = selected("Sound")
		endif

		runScript: "butterworth.praat", "High-pass", freq, 2, "no"
		tmp4 = selected("Sound")
		tmp5 = Filter (stop Hann band): 0, 500, 500

		runScript: "distortion.praat", type_of_distortion$, input_gain, positive_amplitude_limit, negative_amplitude_limit, "No change", "no"
		tmp6 = selected("Sound")
		tmp7 = Filter (stop Hann band): 0, 100, 100

		Scale intensity: int

		x = 0.01 * mix
		amp1 = min(1 + cos(x * pi), 1)
		amp2 = min(1 - cos(x * pi), 1)

		if preview
			dur = min(3, dur)
		endif

		result = Create Sound from formula: "mix" + string$(mix), ch, 0, dur, nsf, "(object[tmp1] * amp1) + (object[tmp7] * amp2)"

		if scale_result_to_original_intensity
			Scale intensity: int
		endif

		runScript: "declip.praat"

		removeObject: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7

include preview2.inc

		if not preview
			Rename: s$ + "-exciter_" + string$(freq) + "_" + phase_rotation$ + "_" + type_of_distortion$ + "_" + string$(input_gain) + "_" + string$(positive_amplitude_limit) + "_" + string$(negative_amplitude_limit) + "_" + string$(mix) + "_" + string$(scale_result_to_original_intensity)
		endif
	else
		if not preview
			Copy: s$ + "-exciter"
		endif

	endif
endproc
