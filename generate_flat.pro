pro build_all_flat
	generate_flat, 10.d0, 1.d0, 1.d-3
	generate_flat, 10.d0, 1.d0, 1.d-2
end

;******************************************
; Modify these quantities to change the total optical depth in the vertical 
; and radial directions
;******************************************
; 	Rout_Rin = 10.d0   ; Increasing it, we increase the optical depth in the radial
; 	                   ; direction without affecting the vertical direction
;  h_Rin = 1.d0       ; Increasing it, we increase the overall optical depth of the disk


pro generate_flat, Rout_Rin, h_Rin, nh2dp_cte

; These values are fixed so that the 2-level models give the expected
; values for epsilon

	Tcte = 15.d0
	nh2_cte = 5.d5	
	Rin = 40.d0	
				
	Rout = Rout_Rin * Rin
	h_tot = h_Rin * Rin	
	
	nh = 100
	nr = 80
			
	rr = findgen(nr)/(nr-1.d0)*(Rout-Rin) + Rin
	deltah = 2*h_tot / nh
	hh = replicate(deltah,nh)
	
	root = 'models_flat/MODEL_RoutRin_'+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
		'_hRin_'+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
		'_nmol_'+strtrim(string(nh2dp_cte,FORMAT='(E7.1)'))+'/'
		
	test = file_test(root,/directory)
	if (test eq 0) then begin
		file_mkdir, root
	endif

	for i = 0, nr-1 do begin
		print, i
				
		T = replicate(Tcte,nh)
		nh2 = replicate(nh2_cte,nh)		
		nhd2p = replicate(nh2dp_cte,nh)

		
		openw,2,root+'radius_'+strtrim(string(i),2)+'.model'
		printf,2,'CONSTANT MODEL at radius = '+strtrim(string(rr[i]),2)+' AU'
		printf,2,' h [cm]   T [K]   n(H2) [cm^-3]   n(H2D+) [cm^-3]'
		printf,2,'--------------------------------------------------'
		printf,2, nh
		printf,2,FORMAT='(4(3X,E15.5))',transpose([[hh*1.5d13],[T],[nh2],[nhd2p]])
		close,2
					
	endfor
	
end
