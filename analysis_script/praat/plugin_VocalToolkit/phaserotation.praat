# Parts of code adapted from He and Dellwo (2016) Interspeech. https://www.isca-archive.org/interspeech_2016/he16_interspeech.html

form: "Phase rotation"
	choice: "Rotation", 1
		option: "+90º Hilbert transform"
		option: "-90º Hilbert transform"
		option: "180º waveform inversion"
endform

include batch.praat

procedure action
	s = selected("Sound")
	s$ = selected$("Sound")

	if rotation = 3
		Copy: s$ + "_phaserotation_180º"
		Formula: "-self"
	else
		sf = Get sampling frequency
		ns = Get number of samples
		dur = Get total duration

		if ns > 65535
			if ln(ns - 1) / ln(2) = round(ln(ns - 1) / ln(2))
				dur = (ns + 1) / sf
			elsif ln(ns) / ln(2) = round(ln(ns) / ln(2))
				dur = (ns + 2) / sf
			endif
		endif

		wrk = Create Sound from formula: "wrk", 1, 0, dur, sf, "object[s]"
		sp = To Spectrum (resampled): 50

		if rotation = 1
			ph = Copy: s$ + "_phaserotation_+90º"
			Formula: "if row = 1 then object[sp, 2, col] else -object[sp, 1, col] fi"
		elsif rotation = 2
			ph = Copy: s$ + "_phaserotation_-90º"
			Formula: "if row = 1 then -object[sp, 2, col] else object[sp, 1, col] fi"
		endif

		To Sound (resampled): 50
		removeObject: wrk, sp, ph
	endif
endproc
