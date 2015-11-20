pro dibuja_profiles_paper
		
; Rout/Rin=10     h/Rin=1
	Rout_Rin = 10.d0
	h_Rin = 0.3d0
	
	restore,'profiles.idl'
	
	abre_ps,'profiles_grid_v_inclination.eps',/mitad,/enca
	
	erase
	multiplot2, [4,4], mxtitle='x', mytitle='Normalized flux'
	
	epsilons = reverse(['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1'])
	neps = n_elements(epsilons)
	
	mols = [1.d-3,1.d-2,1.d-1,1.d0,1.d1,1.d2,1.d3]
	nmols = n_elements(mols)
	which_mols = [2,3,4,5,6]
			
	angles = [1] ;90,60,30,1]
	nangles = n_elements(angles)
	
	vels = [0.d0] ;,1.d0,3.d0,5.d0,7.d0,10.d0]
	nvels = n_elements(vels)
	which_vel = [0,2,3,5]
	
	labels_angles = ['Face-on','i=30!M%!7','i=60!M%!7','Edge-on']
	
	max_y = [0.0699,0.0699,0.0699,0.0699]
	
	taus = fltarr(5)
				
; loop_eps,loop_mols,loop_angle,loop_vels,freq
	
	for loop_ang = n_elements(angles)-1, 0, -1 do begin
		for loop_v = 0, n_elements(vels)-1 do begin
			plot, [0,0], /nodata, xran=[-5,5], xsty=1, yran=[0,max_y[3-loop_ang]], ysty=1, charsize=0.8
			if (loop_v eq 0) then begin
				xyouts,-4.5,0.060,labels_angles[3-loop_ang],/data,charsize=0.8
			endif
			if (loop_ang eq 3) then begin
				xyouts,1.5,0.060,'v/!7D!6v='+$
					strtrim(string(vels[which_vel[loop_v]],FORMAT='(I2)'),2),/data,charsize=0.8
			endif
			for loop_tau = 0, 4 do begin
				oplot, v_axis,$
					profs[2,which_mols[loop_tau],loop_ang,which_vel[loop_v],*],line=loop_tau,$
					thick=3
				taus[loop_tau] = optdepth[2,which_mols[loop_tau],3,0]
			endfor
			if (loop_ang eq 0 and loop_v eq 3) then begin
				legend,'!7s!6='+strtrim(string(taus,FORMAT='(I4)'),2),line=findgen(5),$
					/right,/top,charsize=0.6,spacing=0,thick=replicate(3,5)
			endif
			multiplot2
		endfor
	endfor
	
	cierra_ps
	
	multiplot2,/default		
	
	stop
end

pro dibuja2_tau
		
; Rout/Rin=10     h/Rin=1
	Rout_Rin = 10.d0
	h_Rin = 10.d0
	
	!p.multi = [0,2,2]
	
	epsilons = reverse(['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1'])
	neps = n_elements(epsilons)
	
	mols = [1.d-3,1.d-2,1.d-1,1.d0]
	nmols = n_elements(mols)
	
	angles = [90,60,30,1]
	nangles = n_elements(angles)
	
	vels = [0.d0]
	nvels = n_elements(vels)
	
	for loop_angle = 0, nangles-1 do begin
		
		angle = angles[loop_angle]
		
		abre_ps,"tau_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
			"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+"_angle_"+$
			strtrim(string(angle,FORMAT='(I2)'),2)+".ps"
				
		for loop_vels = 0, nvels-1 do begin
		
			for loop_eps = 0, 5 do begin
			
				epsilon = epsilons[loop_eps]
														
				for loop_mols = 0, 3 do begin
					nmol = mols[loop_mols]
						
					root_output = $
						"MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
						"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
						"_nmol_"+strtrim(string(nmol,FORMAT='(E7.1)'))+"_eps_"+epsilon
				
					restore,root_output+'/disk_angle'+strtrim(string(angle,FORMAT='(I2)'),2)+$
						'_vel_'+strtrim(string(vels[loop_vels],FORMAT='(F4.1)'),2)+'.idl'
									
					tvframe,tot_opt_depth[*,*,50], /sample, /aspect, /bar, $
					tit='!7e!6='+epsilon+' - i='+strtrim(string(angle,FORMAT='(I2)'),2)+'!M%!6',$
						btit='Optical depth'
					
				endfor
				
			endfor
						
		endfor
		
		cierra_ps
		
	endfor
	
	!p.multi=0
	
	stop
end

pro dibuja3_prof
		
; Rout/Rin=10     h/Rin=1
	Rout_Rin = 10.d0
	h_Rin = 0.3d0
	
	!p.multi = [0,3,2]
	
	epsilons = reverse(['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1'])
	neps = n_elements(epsilons)
	
	mols = [1.d-3,1.d-2,1.d-1,1.d0,1.d1,1.d2,1.d3]
	nmols = n_elements(mols)
	
	angles = [90,60,30,1]
	nangles = n_elements(angles)
	
	vels = [0.d0,1.d0,3.d0,5.d0,7.d0,10.d0]
	nvels = n_elements(vels)
	
	optdepth = fltarr(nmols)
	
	for loop_angle = 0, nangles-1 do begin
		
		angle = angles[loop_angle]
		
		abre_ps,"disk_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
			"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+"_angle_"+$
			strtrim(string(angle,FORMAT='(I2)'),2)+".ps"
				
		for loop_vels = 0, nvels-1 do begin
		
			for loop_eps = 0, neps-1 do begin
			
				epsilon = epsilons[loop_eps]
			
				plot, [0,0], /nodata, xran=[-5,5], xsty=1, yran=[0,0.12], ysty=1, $
					tit='!7e!6='+epsilon+' - i='+strtrim(string(angle,FORMAT='(I2)'),2)+'!M%!6',$
					xtit='x',ytit='Normalized flux'
				xyouts, -4.5,0.11,'v/!7D!6v='+strtrim(string(vels[loop_vels],FORMAT='(F4.1)'),2),$
					charsize=0.7
			
				for loop_mols = 0, nmols-1 do begin
					nmol = mols[loop_mols]
						
					root_output = $
						"MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
						"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
						"_nmol_"+strtrim(string(nmol,FORMAT='(E7.1)'))+"_eps_"+epsilon
				
					restore,root_output+'/disk_angle'+strtrim(string(angle,FORMAT='(I2)'),2)+$
						'_vel_'+strtrim(string(vels[loop_vels],FORMAT='(F4.1)'),2)+'.idl'
				
					n_disk = n_elements(where(image3[*,*,49] ne 0.d0))
				
					total_flux = total(total(image3,1),1) / n_disk
					area_total_flux = total(total_flux)
				
					oplot, v_axis, total_flux/area_total_flux, line=loop_mols, thick=3
					
					optdepth[loop_mols] = max(tot_opt_depth[*,*,50])
					
				endfor
				
				legend,'!7s!6='+strtrim(string(optdepth,FORMAT='(F8.2)'),2),line=findgen(nmols),$
					/right,/top,charsize=0.6,spacing=0,thick=replicate(3,nmols)
								
;  				tvframe,tot_opt_depth[*,*,50], /sample, /bar
				
			endfor
			
						
		endfor
		
		cierra_ps
		
	endfor
	
	!p.multi=0
	
	stop
end


pro dibuja3_tau
		
; Rout/Rin=10     h/Rin=1
	Rout_Rin = 10.d0
	h_Rin = 0.3d0
	
	!p.multi = [0,2,2]
	
	epsilons = reverse(['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1'])
	neps = n_elements(epsilons)
	
	mols = [1.d-3,1.d-2,1.d-1,1.d0]
	nmols = n_elements(mols)
	
	angles = [90,60,30,1]
	nangles = n_elements(angles)
	
	vels = [0.d0]
	nvels = n_elements(vels)
	
	for loop_angle = 0, nangles-1 do begin
		
		angle = angles[loop_angle]
		
		abre_ps,"tau_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
			"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+"_angle_"+$
			strtrim(string(angle,FORMAT='(I2)'),2)+".ps"
				
		for loop_vels = 0, nvels-1 do begin
		
			for loop_eps = 0, 5 do begin
			
				epsilon = epsilons[loop_eps]
														
				for loop_mols = 0, 3 do begin
					nmol = mols[loop_mols]
						
					root_output = $
						"MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
						"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
						"_nmol_"+strtrim(string(nmol,FORMAT='(E7.1)'))+"_eps_"+epsilon
				
					restore,root_output+'/disk_angle'+strtrim(string(angle,FORMAT='(I2)'),2)+$
						'_vel_'+strtrim(string(vels[loop_vels],FORMAT='(F4.1)'),2)+'.idl'
									
					tvframe,tot_opt_depth[*,*,50], /sample, /aspect, /bar, $
					tit='!7e!6='+epsilon+' - i='+strtrim(string(angle,FORMAT='(I2)'),2)+'!M%!6',$
						btit='Optical depth'
					
				endfor
				
			endfor
						
		endfor
		
		cierra_ps
		
	endfor
	
	!p.multi=0
	
	stop
end