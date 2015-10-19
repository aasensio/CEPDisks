; Returns the value of the collisional rate for a given 2-level transition
; so that the value of epsilon is fixed
function collis, nu, aul, T, nh2, epsilon
	epsp = epsilon / (1.d0-epsilon)
	PK = 1.3806503d-16
	PH = 6.62606876d-27
	PC = 2.99792458d10
	
	Eup = nu / PC
; 	print, 'E_up if E_low=0 [cm^-1] : ', Eup	
	return, aul*epsp/(1.d0-exp(-PH*nu/(PK*T))) / nh2
end

pro test
	nu = 3.7240d11
	aul = 1.04d-4
	eps = 10.d0^(-findgen(6)*0.5-0.5)
	for i = 0, 5 do print, eps[i], ' -- Cul = ', collis(nu, aul, 15.d0, 5.d5, eps[i])
	stop
end