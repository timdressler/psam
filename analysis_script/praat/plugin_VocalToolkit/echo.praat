# Adapted from the script "echo.praat" by Ingmar Steiner included in "Automatic Speech Data Processing with Praat"
# https://web.archive.org/web/20221209110042/http://www.coli.uni-saarland.de/~steiner/teaching/2007/winter/praat/praat_lecturenotes.pdf

form: "Echo"
	positive: "Delay (s)", "0.5"
	positive: "Amplitude (Pa)", "0.5"
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

delay = min(max(delay, 0.01), 10)
amplitude = min(max(amplitude, 0.1), 0.9)

include batch.praat

procedure action
	s = selected("Sound")
	s$ = selected$("Sound")
	int = Get intensity (dB)

	if int <> undefined

include preview1.inc

		dur = Get total duration

		wrk = Copy: "wrk"
		runScript: "fixdc.praat"
		tt = dur + (delay * (amplitude * 2)) * 10

		result = Extract part: 0, tt, "rectangular", 1, "no"
		Formula: "self + amplitude * self(x - delay)"
		Fade out: 0, tt, -(tt - dur), "yes"
		runScript: "declip.praat"

		removeObject: wrk

include preview2.inc

		if not preview
			Rename: s$ + "-echo"
		endif
	else
		if not preview
			Copy: s$ + "-echo"
		endif
	endif
endproc
