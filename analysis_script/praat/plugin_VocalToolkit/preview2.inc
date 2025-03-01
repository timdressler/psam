if preview
	pre2 = Extract part: 0, preview_dur, "rectangular", 1, "no"
	nowarn Fade in: 0, 0, 0.025, "yes"
	nowarn Fade out: 0, preview_dur, -0.025, "yes"
	int_pre = Get intensity (dB)
	if int_pre <> undefined
		Play
	endif
	selectObject: s
	removeObject: trimmed, pre1, result, pre2
endif
