function disk_prof, h, eps, mol, angle, vel, peakpos
	restore,'profiles.idl'
	restore,'peakpos_fixed.idl'
	
	indh = where(abs(h_disks-h) lt 1.d-3)
	indeps = where(abs(float(epsilons)-eps) lt 1.d-3)
	indmol = where(abs(mols-mol) lt 1.d-3)
	indang = where(abs(angles-angle) lt 1.d-3)
	indvel = where(abs(vels-vel) lt 1.d-3)
	
	if (indh[0] ne -1 and indeps[0] ne -1 and indmol[0] ne -1 and $
		indang[0] ne -1 and indvel[0] ne -1) then begin
		res = profs[indh,indeps,indmol,indang,indvel,*]
		
		peakpos = peakpos[indh,indeps,indmol,indang,indvel,0]
		if (peakpos eq -1) then peakpos = 0.d0
	endif else begin
		res = 0.d0		
		print, 'Not found'
	endelse
	return, res
end

pro read_all
			
	epsilons = ['5.d-1','1.d-2']
	nepsilons = n_elements(epsilons)
	
	xmols = [1.d-5,1.d-6,1.d-7,1.d-8,1.d-9,1.d-10]
	nxmols = n_elements(xmols)
	
	angles = [1,30,60,90]
	nangles = n_elements(angles)
	
	vels = [0.d0,3.d0,5.d0,10.d0]
	nvels = n_elements(vels)
		
	disk_cases = [1,2]
	ncases = n_elements(disk_cases)
	
	optdepth = fltarr(ncases,nepsilons,nxmols,nangles,nvels)
	profs = fltarr(ncases,nepsilons,nxmols,nangles,nvels,100)
	
	for loop_cases = 0, ncases-1 do begin
		disk_case = disk_cases[loop_cases]
		
		for loop_xmols = 0, nxmols-1 do begin
			xmol = xmols[loop_xmols]
			print, xmol
			
			for loop_eps = 0, nepsilons-1 do begin
				epsilon = epsilons[loop_eps]
				
				root_output = 'MODEL_standard400_CASE'+$
					strtrim(string(disk_case,FORMAT='(I1)'),2)+$
					'_xmol'+strtrim(string(xmol,FORMAT='(E7.1)'),2)+'/'
								
				for loop_vels = 0, nvels-1 do begin
					for loop_angles = 0, nangles-1 do begin
						
						angle = angles[loop_angles]
						vel = vels[loop_vels]
						
						restore,root_output+'/disk_angle'+strtrim(string(angle,FORMAT='(I2)'),2)+$
							'_vel_'+strtrim(string(vel,FORMAT='(F4.1)'),2)+'_eps'+epsilon+'.idl'
							
						n_disk = n_elements(where(image3[*,*,49] ne 0.d0))
						
						total_flux = total(total(image3,1),1) / n_disk
						area_total_flux = total(total_flux)
									
						optdepth[loop_cases,loop_eps,loop_xmols,loop_angles,loop_vels] = $
							max(tot_opt_depth[*,*,50])
						profs[loop_cases,loop_eps,loop_xmols,loop_angles,loop_vels,*] = $
							total_flux/area_total_flux

					endfor
				endfor
						
			endfor
		endfor
	endfor
		
	save,epsilons,xmols,angles,vels,disk_cases,v_axis,optdepth,profs,filename='profiles.idl'
	stop
end