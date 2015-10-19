; Calculate the main beam temperature using a square beam
pro plot_flux_square_beam, dir
	
; Restore the calculations performed for one of the disks and
; any of the positions. You only have to change the file with 
; the data
	restore,dir+'/disk_angle90.idl'
	
; Generate a circular mask where the disk is placed
	nx = n_elements(image3[*,0,0])
	ny = n_elements(image3[0,*,0])
	mask = fltarr(nx,ny)
	ind = where(image3[*,*,50] ne 0.d0)
	mask[ind] = 1.d0
	mask[nx/2-4:nx/2+4,ny/2-4:ny/2+4] = 1.d0
										
; flux_center contains the flux, while v_center contains the 
; velocity axis (transformed from the wavelength axis)
	flux_center = fltarr(4,nfreq)
	v_center = fltarr(4,nfreq)
	 
; For each wavelength, we average the image of the disk, thus
; obtaining the wavelength variation of the flux
 	for i = 0, nfreq-1 do begin
	 	
 		ind = where(mask eq 1)
		im = reform(image3[*,*,i])
		flux_center[i] = mean(im[ind])
			
 	endfor
 	
; This is the frequency of the H2D+ transition that we are dealing with
 	nu = 372.4d9
 	
; Transform to velocity units
 	v_center = (freq_axis-nu)/nu*3.d5
	 		 	 
	plot,v_center,flux_center
		
 	stop
end