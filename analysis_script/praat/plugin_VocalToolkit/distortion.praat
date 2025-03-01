# Adapted from the plugin "Free-Clip" by Jonathan Hyde, available at: https://gitlab.com/JHVenn/Free-Clip

form: "Distortion (clipping)"
	choice: "Type of distortion", 1
		option: "Hard clipping"
		option: "Soft clipping. Quintic"
		option: "Soft clipping. Cubic"
		option: "Soft clipping. Hyperbolic tangent"
		option: "Soft clipping. Algebraic"
		option: "Soft clipping. Arctangent"
	real: "Input gain (dB)", "0.0 (= no change)"
	real: "Positive amplitude limit (Pa)", "0.5"
	real: "Negative amplitude limit (Pa)", "-0.5"
	choice: "Output gain", 1
		option: "No change"
		option: "Scale to original average intensity"
		option: "Normalize (scale peak 0.99)"
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

if positive_amplitude_limit < 0
	exitScript: "“Positive amplitude limit” must be greater than or equal to 0.0"
endif
if negative_amplitude_limit > 0
	exitScript: "“Negative amplitude limit” must be less than or equal to 0.0"
endif

a1 = abs(positive_amplitude_limit)
a2 = abs(negative_amplitude_limit)

include batch.praat

procedure action
	s = selected("Sound")
	s$ = selected$("Sound")
	int = Get intensity (dB)

	if int <> undefined

include preview1.inc

		result = Copy: "tmp"

		if input_gain <> 0
			Scale intensity: int + input_gain
		endif

		Formula: "if self > 0 then if " + string$(a1) + " <> 0 then self * (1 / " + string$(a1) + ") else 0 fi else if " + string$(a2) + " <> 0 then self * (1 / " + string$(a2) + ") else 0 fi fi"

		if type_of_distortion = 1
			Formula: "if self > 0 then min(self, 1) else if self < 0 then max(self, -1) else 0 fi fi"
		elsif type_of_distortion = 2
			Formula: "if abs(self) < 1.25 then self - (256 / 3125) * (self ^ 5) else if self > 0 then 1 else if self < 0 then -1 else 0 fi fi fi"
		elsif type_of_distortion = 3
			Formula: "if abs(self) < 1.5 then self - (4 / 27) * (self ^ 3) else if self > 0 then 1 else if self < 0 then -1 else 0 fi fi fi"
		elsif type_of_distortion = 4
			Formula: "tanh(self)"
		elsif type_of_distortion = 5
			Formula: "self / sqrt(1 + self ^ 2)"
		elsif type_of_distortion = 6
			Formula: "2 / pi * arctan(1.6 * self)"
		endif

		Formula: "if self > 0 then self * " + string$(a1) + " else self * " + string$(a2) + " fi"

		if output_gain = 2
			Scale intensity: int
		elsif output_gain = 3
			Scale peak: 0.99
		endif

include preview2.inc

		if not preview
			Rename: s$ + "-distortion_" + type_of_distortion$ + "__" + string$(input_gain) + "_" + string$(positive_amplitude_limit) + "_" + string$(negative_amplitude_limit) + "_" + output_gain$
		endif
	else
		if not preview
			Copy: "tmp"
			Rename: s$ + "-distortion_" + type_of_distortion$ + "__" + string$(input_gain) + "_" + string$(positive_amplitude_limit) + "_" + string$(negative_amplitude_limit) + "_" + output_gain$
		endif
	endif
endproc
