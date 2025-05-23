form: "EQ 10 bands"
	real: "Band 1.  31.5 Hz (dB)", "-24"
	real: "Band 2.  63 Hz (dB)", "-24"
	real: "Band 3.  125 Hz (dB)", "-24"
	real: "Band 4.  250 Hz (dB)", "-24"
	real: "Band 5.  500 Hz (dB)", "24"
	real: "Band 6.  1000 Hz (dB)", "24"
	real: "Band 7.  2000 Hz (dB)", "24"
	real: "Band 8.  4000 Hz (dB)", "-24"
	real: "Band 9.  8000 Hz (dB)", "-24"
	real: "Band 10.  16000 Hz (dB)", "-24"
	boolean: "Preview (Apply. Uncheck to publish)", 1
endform

band_1.__31.5_Hz = min(max(round(band_1.__31.5_Hz), -24), 24)
band_2.__63_Hz = min(max(round(band_2.__63_Hz), -24), 24)
band_3.__125_Hz = min(max(round(band_3.__125_Hz), -24), 24)
band_4.__250_Hz = min(max(round(band_4.__250_Hz), -24), 24)
band_5.__500_Hz = min(max(round(band_5.__500_Hz), -24), 24)
band_6.__1000_Hz = min(max(round(band_6.__1000_Hz), -24), 24)
band_7.__2000_Hz = min(max(round(band_7.__2000_Hz), -24), 24)
band_8.__4000_Hz = min(max(round(band_8.__4000_Hz), -24), 24)
band_9.__8000_Hz = min(max(round(band_9.__8000_Hz), -24), 24)
band_10.__16000_Hz = min(max(round(band_10.__16000_Hz), -24), 24)

include batch.praat

procedure action
	s = selected("Sound")
	s$ = selected$("Sound")
	original_dur = Get total duration
	int = Get intensity (dB)

	if int <> undefined

include preview1.inc

		runScript: "workpre.praat"
		wrk = selected("Sound")
		sf = Get sampling frequency
		dur1 = Get total duration

		pointprocess = Create empty PointProcess: "pulse", 0, 0.05
		Add point: 0.025

		pulse = noprogress To Sound (pulse train): sf, 1, 0.05, 2000

		sp_pulse = noprogress To Spectrum: "no"

		buffer = Copy: "buffer"
		Formula: "0"

		sp_eq = Copy: "sp_eq"
		@eqBand: 0, 44.2, band_1.__31.5_Hz, 20
		@eqBand: 44.2, 88.4, band_2.__63_Hz, 20
		@eqBand: 88.4, 177, band_3.__125_Hz, 40
		@eqBand: 177, 354, band_4.__250_Hz, 80
		@eqBand: 354, 707, band_5.__500_Hz, 100
		@eqBand: 707, 1414, band_6.__1000_Hz, 100
		@eqBand: 1414, 2828, band_7.__2000_Hz, 100
		@eqBand: 2828, 5657, band_8.__4000_Hz, 100
		@eqBand: 5657, 11314, band_9.__8000_Hz, 100
		@eqBand: 11314, max(12000, sf / 2), band_10.__16000_Hz, 100
		Filter (stop Hann band): 0, 80, 20
		Filter (pass Hann band): 0, 20000, 100

		pulse_eq_tmp = noprogress To Sound
		dur2 = Get total duration

		pulse_eq = Extract part: (dur2 - 0.05) / 2, dur2 - ((dur2 - 0.05) / 2), "Hanning", 1, "no"
		Scale peak: 0.99
		pulse_sf = Get sampling frequency
		if pulse_sf <> sf
			Override sampling frequency: sf
		endif

		plusObject: wrk
		tmp1 = Convolve: "sum", "zero"

		tmp2 = Extract part: 0.025, dur1 + 0.025, "rectangular", 1, "no"

		runScript: "workpost.praat", original_dur
		Scale intensity: int
		runScript: "declip.praat"
		result = selected("Sound")
		dur3 = Get total duration
		if dur3 > 0.5
			Fade in: 0, 0, 0.005, "yes"
			Fade out: 0, dur3, -0.005, "yes"
		endif

		removeObject: wrk, pointprocess, pulse, sp_pulse, buffer, sp_eq, pulse_eq_tmp, pulse_eq, tmp1, tmp2

include preview2.inc

		if not preview
			Rename: s$ + "-EQ10bands"
		endif
	else
		if not preview
			Copy: s$ + "-EQ10bands"
		endif
	endif
endproc

procedure eqBand: .bnd1, .bnd2, .db, .smoothing
	.amp = 0.00002 * 10 ^ (.db / 20)
	selectObject: buffer
	Formula: "object[sp_pulse]"
	Filter (pass Hann band): .bnd1, .bnd2, .smoothing
	Formula: "self * .amp"
	selectObject: sp_eq
	Formula: "self + object[buffer]"
endproc
