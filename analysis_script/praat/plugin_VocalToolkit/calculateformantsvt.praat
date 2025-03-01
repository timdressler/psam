# Adapted from the procedure described at: https://web.archive.org/web/20190616135811/http://www.languagebits.com/phonetics-english/resonant-frequencies-and-vocal-tract-length/
# Original formula from Johnson, Keith. Acoustic and Auditory Phonetics. 2nd ed. Malden, Mass: Blackwell Pub, 2003. p. 96

form: "Calculate formants of a vocal tract"
	comment: "Displays in the Info window the average formant values of a vocal tract"
	comment: "in the neutral configuration, calculated from the entered length value."
	positive: "Vocal tract length (cm)", "17.5"
	natural: "Number of formants", "5"
endform

appendInfoLine: "Calculate formants of a vocal tract..."
appendInfoLine: tab$, "Vocal tract length: ", vocal_tract_length, " cm"

for i to number_of_formants
	freq = ((2 * i - 1) * 35000) / (4 * vocal_tract_length)
	appendInfoLine: tab$, "F", i, ": ", number(fixed$(freq, 3)), " Hz"
endfor

appendInfoLine: newline$, "> Adapted from the procedure described at: https://bit.ly/vocal-tract-length-estimation", newline$, newline$, newline$
