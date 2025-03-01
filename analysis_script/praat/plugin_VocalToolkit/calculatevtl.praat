# Adapted from the procedure described at: https://web.archive.org/web/20190616135811/http://www.languagebits.com/phonetics-english/resonant-frequencies-and-vocal-tract-length/
# Original formula from Johnson, Keith. Acoustic and Auditory Phonetics. 2nd ed. Malden, Mass: Blackwell Pub, 2003. p. 96

form: "Calculate vocal tract length"
	comment: "Displays in the Info window the estimated vocal tract length in the"
	comment: "neutral configuration, calculated from a formant frequency value."
	positive: "Formant frequency (Hz)", "3500"
	natural: "Formant number", "4"
	boolean: "Calculate from the selected Sounds", 1
	comment: "Formant determination"
	positive: "Maximum formant (Hz)", "5500 (= adult female)"
	comment: "Set 5000 Hz for men, 5500 Hz for women or up to 8000 Hz for children."
endform

include batch.praat

procedure action
	if i_batch = 1
		if calculate_from_the_selected_Sounds
			appendInfoLine: "Calculate vocal tract length from the selected Sounds..."
		else
			appendInfoLine: "Calculate vocal tract length..."
		endif
	endif

	if calculate_from_the_selected_Sounds
		s = selected("Sound")
		s$ = selected$("Sound")
		int = Get intensity (dB)

		if int <> undefined
			runScript: "workpre.praat"
			tmp1 = selected("Sound")

			runScript: "extractvowels.praat", 0, 0
			tmp2 = selected("Sound")

			runScript: "workpre.praat"
			tmp3 = selected("Sound")

			formant = noprogress nowarn To Formant (robust): 0.005, 5, maximum_formant, 0.025, 50, 1.5, 5, 0.000001
			formant_frequency = Get mean: formant_number, 0, 0, "hertz"

			selectObject: s
			removeObject: tmp1, tmp2, tmp3, formant
		else
			formant_frequency = undefined
		endif

		appendInfoLine: tab$, s, ". Sound ", s$
		t$ = tab$
	else
		t$ = ""
	endif

	prep = 35000 * ((formant_number / 2) - 0.25)
	vtl = prep / formant_frequency

	if calculate_from_the_selected_Sounds or (not calculate_from_the_selected_Sounds and i_batch = 1)
		appendInfoLine: t$, tab$, "Estimated vocal tract length: ", number(fixed$(vtl, 2)), " cm   (mean F", formant_number, " = ", number(fixed$(formant_frequency, 3)), " Hz)", newline$
	endif

	if i_batch = n_batch
		appendInfoLine: "> Adapted from the procedure described at: https://bit.ly/vocal-tract-length-estimation", newline$, newline$, newline$
	endif
endproc
