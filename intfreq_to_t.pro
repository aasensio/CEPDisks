function intfreq_to_t, freq, inten, rayleigh=rayleigh
	 PH = 6.62606876d-27
	 PK = 1.3806503d-16
	 PC = 2.99792458d10
	 cte1 = PH * freq / PK
	 cte2 = 2.d0 * PH * freq^3 / (PC^2 * inten)
	 if (keyword_set(rayleigh)) then begin
		return, PC^2/(2.d0*PK*freq^2) * inten
	 endif else begin
		return, cte1 / alog(1.d0+cte2)
	 endelse
end
