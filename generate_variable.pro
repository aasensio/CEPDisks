; Locate the positions of the disk where I will set the values of
; epsilon. These values of epsilon are put at R=R_in and at
; different heights. We will generate molecular models with collisional
; rates so that the value of epsilon at these points are the ones we want
pro locate_positions
	a = ddread('disk_models/standard_400au.dat',offset=3,/countal)
	
	rr = reform(a(0,*))
	hh = reform(a(1,*))

	ind = uniq(rr)
	nr = n_elements(ind)
	
	i = 0 
	equal = where(rr eq rr[ind[i]])
	n = n_elements(equal)

	h = hh(equal)
	T = reform(a[2,equal])
	nh2 = reform(a[3,equal])
	
	hdup = [-reverse(h[1:n-1]),h]*1.5d13
	Tdup = [reverse(T[1:n-1]),T]
	nhdup = [reverse(nh2[1:n-1]),nh2]
	
	nh = n_elements(nhdup)
		
	kul = 2.9016223e-12
	aul = 1.04d-4
	nu = 3.7240e+011
	epsp = kul*nhdup/aul*(1.d0-exp(-6.626d-27*nu/1.381d-16/tdup))	
	eps = epsp / (1.d0+epsp)
	
	col_dens = fltarr(nh)
	for i = 1, nh-1 do begin
		col_dens[i] = int_tabulated(hdup[0:i],nhdup[0:i])
	endfor
	Av = col_dens / 5.d21
	
	!p.multi = [0,1,2]
	plot, hdup/1.5d13, Av, xran=[-20,0],/ylog
	very,2.
	verx,-7.2
	verx,-4.6
	plot, hdup/1.5d13, Tdup, xran=[-20,0],yran=[18,60],ysty=1
	very, 22.
	verx,-4.6
	!p.multi=0
	
	kul = 5.585d-13
	aul = 1.04d-4
	nu = 3.7240e+011
	epsp = kul*nhdup/aul*(1.d0-exp(-6.626d-27*nu/1.381d-16/tdup))	
	eps = epsp / (1.d0+epsp)
	plot, hdup/1.5d13, eps, xran=[-10,0], /ylog
	verx, -7.2
	verx, -4.6
; 	very, 1.d-3, line=1
; 	very, 1.d-2, line=1
; 	very, 1.d-1, line=1
	very, 0.5, line=1
	print, kul
	
	stop
end

pro build_all_variable
; 	generate_variable, 1, 1.d-5
; 	generate_variable, 1, 1.d-6
; 	generate_variable, 1, 1.d-7
; 	generate_variable, 1, 1.d-8
; 	generate_variable, 1, 1.d-9
; 	generate_variable, 1, 1.d-10
	
	generate_variable, 2, 1.d-5
	generate_variable, 2, 1.d-6
	generate_variable, 2, 1.d-7
	generate_variable, 2, 1.d-8
	generate_variable, 2, 1.d-9
	generate_variable, 2, 1.d-10
end

;******************************************
; Modify these quantities to change the total optical depth in the vertical 
; and radial directions
;******************************************
; 	Rout_Rin = 10.d0   ; Increasing it, we increase the optical depth in the radial
; 	                   ; direction without affecting the vertical direction
;  h_Rin = 1.d0       ; Increasing it, we increase the overall optical depth of the disk


pro generate_variable, disk_case, xmol

; These values are fixed so that the 2-level models give the expected
; values for epsilon
	
	a = ddread('disk_models/standard_400au.dat',offset=3,/countal)
	
	rr = reform(a[0,*])
	hh = reform(a[1,*])

	ind = uniq(rr)
	nr = n_elements(ind)
	
	root = 'models_variable/MODEL_standard400_CASE'+$
		strtrim(string(disk_case,FORMAT='(I1)'),2)+$
		'_xmol'+strtrim(string(xmol,FORMAT='(E7.1)'),2)+$
		'/'
		
	test = file_test(root,/directory)
	if (test eq 0) then begin
		file_mkdir, root
	endif
		
	for i = 0, nr-1 do begin
		print, i
		equal = where(rr eq rr[ind[i]])
		n = n_elements(equal)

		h = hh(equal)
		T = reform(a[2,equal])
		nh2 = reform(a[3,equal])
		nhd2p = reform(a[4,equal]*nh2*0.25)

		hdup = [-reverse(h[1:n-1]),h]*1.5d13
		Tdup = [reverse(T[1:n-1]),T]
		nhdup = [reverse(nh2[1:n-1]),nh2]
		
		nh = n_elements(hdup)
		
		dzdup = fltarr(nh)
		for j = 0, nh-2 do begin
			dzdup[j] = hdup[j+1]-hdup[j]
		endfor
		dzdup[nh-1] = dzdup[nh-2]
				
		col_dens = fltarr(nh)
		for j = 1, nh-1 do begin
			col_dens[j] = int_tabulated(hdup[0:j],nhdup[0:j])
		endfor
		Av = col_dens / 5.d21
		Av[nh/2+1:nh-1] = reverse(Av[0:nh/2-1])
		
		case (disk_case) of
			1 : begin
					index = where(Av le 2.d0)
				 end
			2 : begin
					index = where(Av gt 2.d0) ;where(Tdup gt 15.d0 and Av gt 1.d0)
				 end
; 			3 : begin
; 					index = where(Tdup le 22.d0)
; 				 end
		endcase
		nmol = replicate(1.d-8,nh)
		if (index[0] ne -1) then begin
			nmol[index] = xmol * nhdup[index]
		endif
								
		openw,2,root+'radius_'+strtrim(string(i),2)+'.model'
		printf,2,'MODEL at radius = '+strtrim(string(rr[ind[i]]),2)+' AU'
		printf,2,' h [cm]   T [K]   n(H2) [cm^-3]   n(mol) [cm^-3]'
		printf,2,'--------------------------------------------------'
		printf,2, 2*n-1
		printf,2,FORMAT='(4(3X,E10.2))',transpose([[dzdup],[Tdup],[nhdup],[nmol]])
		close,2
		
	endfor
	
; 	stop
	
end
