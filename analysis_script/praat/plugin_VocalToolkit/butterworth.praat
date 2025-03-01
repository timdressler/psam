# Sources: https://www.electronics-tutorials.ws/filter/filter_8.html, https://www.desmos.com/calculator/dxcnvc63av

form: "Butterworth filter"
	choice: "Filter", 1
		option: "High-pass"
		option: "Low-pass"
	positive: "Frequency (Hz)", "1000"
	positive: "Filter order", "2"
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

filter_order = min(max(filter_order, 1), 100)

include batch.praat

procedure action
	s = selected("Sound")
	s$ = selected$("Sound")

	c = frequency
	n = filter_order

include preview1.inc

	if filter = 1
		result = Filter (formula): "self * ((x / c) ^ n) / sqrt(1 + (x / c) ^ (2 * n))"
	elsif filter = 2
		result = Filter (formula): "self * 1 / sqrt(1 + (x / c) ^ (2 * n))"
	endif

include preview2.inc

	if not preview
		Rename: s$ + "_butterworth_" + filter$ + "_" + string$(c) + "_" + string$(n)
	endif
endproc
