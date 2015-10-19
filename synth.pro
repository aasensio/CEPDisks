;***********************************
; Set some nice colors for plotting
;***********************************
pro cwpal
	r = [0,255,234,0  ,0  ,255,0  ,255,255,0  ,175,255]
	g = [0,255,0  ,140,216,235,235,0  ,148,126,0  ,67 ]
	b = [0,255,0  ,234,0  ,0  ,228,200,0  ,0  ,201,67 ]
	tvlct, r, g, b
end
;***********************************
; Returns psim and psi0 of linear short-characteristics scheme
; This routine is special for the vectorized formal solution
;***********************************
function linear_polycoef, dtm, n_freq    	
	u = dblarr(n_freq,2)
	salida = dblarr(n_freq,2)
	exu = dblarr(n_freq)
	inds = where( (dtm ne 0.d0) and (dtm gt 0.001))
	if (inds(0) ne -1) then begin
	    exu(inds) = exp(-dtm(inds))
	    u(inds,0) = 1.d0 - exu(inds)
	    u(inds,1) = dtm(inds) - 1.d0 + exu(inds)
       	    salida(inds,0) = u(inds,1)/dtm(inds)
	    salida(inds,1) = u(inds,0) - salida(inds,0)
	endif
	
	inds = where( (dtm ne 0.d0) and (dtm lt 0.001))	
	if (inds(0) ne -1) then begin	    
	    salida(inds,0) = dtm(inds)/2.d0 - dtm(inds)^2/6.d0
	    salida(inds,1) = dtm(inds)/2.d0 - dtm(inds)^2/3.d0
	endif
			
	return, salida
end

;***********************************
; Returns psim and psi0 of linear short-characteristics scheme
;***********************************
function linearcoef, dtm
	u = dblarr(2)
	salida = dblarr(2)
	if (dtm ne 0.d0) then begin

		if (dtm > 0.001) then begin
			exu = exp(-dtm)
			u(0) = 1.d0 - exu
			u(1) = dtm - 1.d0 + exu
			salida(0) = u(1)/dtm
			salida(1) = u(0) - salida(0)
		endif else begin
			salida(0) = dtm/2.d0 - dtm^2/6.d0
			salida(1) = dtm/2.d0 - dtm^2/3.d0
		endelse
		
	endif
		
	return, salida
end

;***********************************
; This routine solves the radiative transfer equation using the
; linear short-characteristics method. Vectorized form for simultaneously
; solving the RT equation for many wavelengths.
; Let the source function and the opacity be known at two consecutive 
; points M and O. Assuming a linear behavior of the source function between
; both points, the intensity at point O can be calculated using an 
; analytical integration of the formal solution of the RT equation.
; I_O = I_M*exp(-dtm) + psi_M*S_m + psi_0*S_0
; 
;  INPUT
;  x : geometrical axis along the ray
;  chi : absorption coefficient for each point along the ray and each wavelength 
;          chi(n_freq,npoints)
;  S : source function for each point along the ray and each wavelength S(n_freq,npoints)
;  I0 : boundary condition for each wavelength I0(n_freq)
;  mu : inclination angle (not used in this program and set always to 1)
;  OUTPUT
;  This function returns the intensity at all points along the ray for all wavelengths
;  tot_optical_depth : total optical depth along the ray
;***********************************
function shortcar_source_linear_poly, x, chi, S, mu, I0, tot_optical_depth
	
; Calculate the number of points along the ray and the number
; of wavelength points
	tamano = n_elements(x)
	nf = n_elements(chi(*,0))
	
; Set the boundary condition	
	Inten = dblarr(nf,tamano)
	Inten(*,0) = I0
	
	tot_optical_depth = fltarr(nf)

; Loop over the points along the ray
	for k=1, tamano-1 do begin

		chitot = dblarr(2,nf)
		Stot = dblarr(2,nf)
		exu = dblarr(nf)
		
		km = k-1
		
		dm = abs((abs(x(km))-abs(x(k))) / mu)
		
		chitot(0,*) = chi(*,km)
		chitot(1,*) = chi(*,k)
				
		Stot(0,*) = S(*,km)
		Stot(1,*) = S(*,k)
				
; Calculate optical depth
		dtm = 0.5d0 * (chitot(0,*) + chitot(1,*)) * dm
	   inds = where(dtm gt 0.001)
		if (inds(0) ne -1) then exu(inds) = 1.d0*exp(-dtm(inds))
		inds = where(dtm lt 0.001)
		if (inds(0) ne -1) then exu(inds) = 1.d0 - dtm(inds) + dtm(inds)^2/2.d0
		
; Apply the linear short-characteristics formula
		Inten(*,k) = Inten(*,k-1) * exu + total(linear_polycoef(dtm, nf) * transpose(Stot),2)
		
; Add the contribution to the total optical depth of this step
		tot_optical_depth = tot_optical_depth + dtm
		
	endfor
	
	return, Inten
end


;***********************************
; This is the main routine that solves the RT equation along many rays crossing the
; disk.
; It reads a result file placed on a directory (for using different calculations)
; INPUT
;  dir : directory where the file "results.idl" is placed. This file contains the
;        the height, radius, population and kinetic temperature for all points of
;        the disk. (See routine gen_results for details how to build such a file)
;  v_Rin : velocity in units of the Doppler width (which is in fact 1 km/s), so that
;          it is indeed the velocity in the inner radius in units of km/s
;  angle : angle between the plane of the disk and the line-of-sight. angle=90 is 
;          for face-on observation
;***********************************
pro synth, dir, v_Rin, angle, clean=clean, fileplus=fileplus, npix=npix
    
; Set the color palette
    loadct,0
	 cwpal
	 	 
; Restore the results of the NLTE calculation for each position on the disk
	 if (not keyword_set(fileplus)) then begin
	 	fileplus = ''
	 endif
	 restore,DIR+'/results'+fileplus+'.idl'
	 	 
	 nhmax = n_elements(height[0,*])
	 

; Plot the disk structure
	 hmax = max(height)
	 rmax = max(radius)
	 size_max = max([hmax,rmax])
	 	 
	 if (not windowavailable(1)) then window,1
	 wset,1
	 plot, radius[*,nhmax-1], height[*,nhmax-1], yran=[-size_max*1.01,size_max*1.01], $
	 	xran=[-size_max*1.01,size_max*1.01], xsty=1, ysty=1, xtit='X',ytit='Z'
	 oplot, radius[*,nhmax-1], height[*,0]
	 oplot, -radius[*,nhmax-1], height[*,nhmax-1]
	 oplot, -radius[*,nhmax-1], height[*,0]
	 
	 if (not windowavailable(2)) then window,2
	 wset,2
	 plot, radius[*,nhmax-1], height[*,nhmax-1], yran=[-size_max*1.01,size_max*1.01], $
	 	xran=[-size_max*1.01,size_max*1.01], xsty=1, ysty=1, xtit='Y',ytit='Z'
	 oplot, radius[*,nhmax-1], height[*,0]
	 oplot, -radius[*,nhmax-1], height[*,nhmax-1]
	 oplot, -radius[*,nhmax-1], height[*,0]
	 
	 theta = findgen(100)/99.*2.d0*!DPI
	 if (not windowavailable(3)) then window,3
	 wset,3
 	 plot, rmax*cos(theta), rmax*sin(theta), yran=[-rmax*1.01,rmax*1.01], $
	 	xran=[-rmax*1.01,rmax*1.01], xsty=1, xtit='X',ytit='Y', ysty=1
 	 	 
; Number of radius points and height points
	 nr = n_elements(radius[*,0])
	 nh = n_elements(radius[0,*])
	 
; Viewing angles
	 theta = angle * !DPI / 180.d0
	 z0 = max(height) + 5.105d0
	 	 
	 rmax = max(radius)
	 rmin = min(radius)
	 hmax = max(height)
	 hmin = min(height)
	 	 	 
	 print, 'Keplerian velocity at R_in = ', v_Rin
	 print, 'Keplerian velocity at R_out = ', v_Rin / sqrt(rmax/rmin)	 	 
	 
; Desired size of the steps along the ray
; We choose to have at least 40 along each ray

; 	 ray_step = 2*hmax / sin(theta) / 40.d0
; 	 ray_step = min([hmax,rmax]) / 20.d0
	 
; Measure the maximum length of the ray. Normally, it will go through
; the disk because h << R. However, test for all situations
	 length = 2*hmax / sin(theta)
	 xcut = length * cos(theta)	 
	 if (xcut gt 2*rmax) then begin
	 	length = 2*rmax
	 endif
	 ray_step = length / 30.d0
	 print, 'Ray step = ', ray_step

	 	 
; Number of pixels where to perform the ray-tracing (npx)
; Since the pixel size is always constant, increasing this number reduces the
; size of the disk on the final image
	 if (keyword_set(npix)) then begin
	 	npx = npix
	 endif else begin
	 	npx = 40
	 endelse
	 step = 2.1d0 * rmax / npx
	 image = dblarr(npx,npx)
	 image2 = dblarr(npx,npx)	 

; Upper and lower boundaries of the disk
	 upper_boundary = dblarr(nr)
	 lower_boundary = dblarr(nr)

; Set the maximum (and minimum) height for each radial position
	 for i = 0, nr-1 do begin
	 	  upper_boundary[i] = max(height[i,*])
		  lower_boundary[i] = min(height[i,*])
	 endfor

; Calculate a parabolic fit to the boundaries of the disk
	 coef_upper = poly_fit(radius[*,0],upper_boundary,1)
	 coef_lower = poly_fit(radius[*,0],lower_boundary,1)
	 
	 up = poly(radius[*,0],coef_upper)
	 low = poly(radius[*,0],coef_lower)
	 
; Verify that the fits work properly
	 wset,1
	 oplot, radius[*,0], up, col=3
	 oplot, radius[*,0], low, col=3
	 
	 wset,2
	 oplot, radius[*,0], up, col=3
	 oplot, radius[*,0], low, col=3
	 	 	 
; When the disk is inclined, we have to make sure that we see the
; complete disk in order not to loose any part of it
; Positions of the extreme rays of the integration
; The following formulae calculate some simple geometrical factors in order
; to make sure that the complete disk is seen
; 	 D1 = -(-(rmax+50.d0) + (z0 - hmax) * tan(!DPI/2.d0-theta)) 
; 	 D2 = -((rmax+50.d0) + 2.d0 * hmax * tan(!DPI/2.d0-theta))
; 	 print, D1, D2
	 
	 dd = z0 - hmax
	 D1 = rmax - dd * tan(!DPI/2.d0-theta) + 50.d0
	 D2 = -rmax - (2*hmax+dd) * tan(!DPI/2.d0-theta) - 50.d0
	 print, D1, D2	 
	 
; We generate the position of the pixels of the "camera" that we use
; to observe the disk. This camera is built such that, when sending inclined
; rays passing through its pixels, we cover the whole disk
	 x0_arr = dindgen(npx)/(npx-1.d0) * (D2-D1) + D1
	 y0_arr = dindgen(npx)/(npx-1.d0) * 2*rmax - rmax ;(D2-D1) + D1 ;2*hmax - hmax
	 	 	 
	 print, 'R(max) : ', rmax
	 print, 'R(min) : ', rmin

; Some fundamental data to describe the atomic properties of the transition
; THEY HAVE TO BE CONSISTENT WITH THE VALUES USED FOR THE NLTE CALCULATION
	 PC = 2.99792458d10
	 PH = 6.62606876d-27
	 PK = 1.3806503d-16
	 UMA = 1.66053873d-24
	 H2DP_MASS = 4.d0
	 AU = 1.49598d13
	 G = 6.673d-8
	 M = 1.98892d33 * 0.5d0
	 
	 Aul = 1.04d-4
	 nu = 372.4d9          ; 
	 vtherm = 1.d5        ; Microturbulent velocity
	 gu = 3.d0
	 gl = 1.d0
	 
	 Bul = PC^2 / (2.d0*PH*nu^3) * Aul
	 Blu = gu/gl * Bul
	 
; Einstein coefficients
	 print, 'Aul : ', Aul
	 print, 'Bul : ', Bul
	 print, 'Blu : ', Blu
	 
; We generate a frequency axis to cover the spectral line. The range [-7,7]
; should be enough
	 nfreq = 100
	 v_axis = dindgen(nfreq) / (nfreq-1.d0) * 14.d0 - 7.d0
	 freq_axis = v_axis * 1.d5 * nu / PC + nu
	 
; This array will contain the final disk image for each wavelength
	 image3 = dblarr(npx,npx,nfreq)
; and this will contain the total optical depth for each wavelength
	 tot_opt_depth = dblarr(npx,npx,nfreq)
	 
	 xpos = dblarr(npx,npx)
	 ypos = dblarr(npx,npx)
	 
; The population of the upper and lower level of the transition is known
; at an irregular grid because the disk model was given in (radius,height) coordinates.
; For interpolation purposes, it is better to have everything on a regular grid
	 r_interpol = radius
	 h_interpol = height
	 pop_lower = reform(popul[*,*,0])
	 pop_upper = reform(popul[*,*,1])
	 
; This generates a regular grid
	 triangulate,r_interpol,h_interpol,tr,b
	 
; And we interpolate the population of the lower and upper level, together with
; the kinetic temperature at each point
	 ninterp = 200
	 nl_surf = trigrid(r_interpol,h_interpol,pop_lower,tr,$
	 	xgrid=rgrid,ygrid=hgrid,nx=ninterp,ny=ninterp)
	 nu_surf = trigrid(r_interpol,h_interpol,pop_upper,tr,$
	 	xgrid=rgrid,ygrid=hgrid,nx=ninterp,ny=ninterp)
	 t_surf = trigrid(r_interpol,h_interpol,tkin,tr,$
	 	xgrid=rgrid,ygrid=hgrid,nx=ninterp,ny=ninterp)
	 	
	 
	 velo = dblarr(npx,npx)
	 velo_noproject = dblarr(npx,npx)
	 	 	 	 
; This is the core of the calculation
; Loop over all the points in the final image (pixels)
	 for loopx = 0, npx-1 do begin	 
	 	  
	 	  for loopy = 0, npx-1 do begin
	 	   	
				x0 = x0_arr[loopx]
				y0 = y0_arr[loopy] ;-1.05*rmax + loopy * step
				
				xpos[loopx,loopy] = x0
				ypos[loopx,loopy] = y0
				
				n = 0
				col = 0
				rayout = 0
				i = 0
				chi_vector = dblarr(nfreq,max([hmax,rmax]*2)/ray_step + 100)
				source_vector = dblarr(nfreq,max([hmax,rmax]*2)/ray_step + 100)
								
; We need to calculate all the properties along the ray and finally
; carry out the formal solution
				while (rayout eq 0) do begin
				
; Calculate the position along the ray
	 				 t = i*ray_step
					 x = t * cos(theta) + x0
					 y = y0
					 z = -t * sin(theta) + z0
					 
; Calculate the radius position of the ray inside the disk
					 rp = sqrt(x^2+y^2)
					 phi = atan(y,x)
					 
					 h_upper = poly(rp,coef_upper)
					 h_lower = poly(rp,coef_lower)

; Is the ray marching inside the disk?
					 if (rp le rmax and rp ge rmin and z le h_upper and z ge h_lower) then begin

; Plot the position of the present point
						  wset,1
						  plots,x,z,psym=1
						  
						  wset,2
						  plots,y,z,psym=1
						  
						  wset,3
						  plots,x,y,psym=1
						  
						  col = rp
						  
						  ix = (rp-rmin) / (rmax-rmin) * ninterp
						  iy = (z-hmin) / (hmax-hmin) * ninterp
						  
; Calculate the value of the lower and upper levels population and the kinetic temperature
; Since everything is now in a regular grid, we can use a bi-linear interpolation
						  n_l = bilinear(nl_surf,ix,iy)
						  n_u = bilinear(nu_surf,ix,iy)
						  Tk = bilinear(T_surf,ix,iy)
						  						  
						  deltanu = nu / PC * vtherm

; This constant includes the area normalization of the line profile (1/sqrt(pi)*deltanu)
	 					  factor = PH * nu / (4.d0 * !DPI * sqrt(!DPI) * deltanu)
						  
; Calculate absorption and emission coefficients for the line
						  eta = factor * (n_l*Blu - n_u*Bul)
						  eps = factor * n_u * Aul
						  
; Projection of the velocity on the LOS
; 						  v = sqrt(G*M / (rp*AU))
						  
						  v = v_Rin / sqrt(rp/rmin) * 1.d5
						  						  						  
						  phi_vector = [-sin(phi),cos(phi),0.d0]
						  ray_vector = [x,y,z] - [x0,y0,z0]
						  ray_vector = ray_vector / sqrt(total(ray_vector*ray_vector))
						  		  
						  proy_los = total(phi_vector*ray_vector)
						  
						  v_los = v * proy_los
						  
						  velo[loopx,loopy] = v_los
						  velo_noproject[loopx,loopy] = v
						  						  						  
; Normalized line profile. We choose a gaussian profile with a thermal broadening
						  prof = exp(-(freq_axis-nu-nu*v_los/PC)^2 / deltanu^2)
						  
; The opacity multiplied by the line profile
						  chi_vector[*,n] = eta * prof
; and the source function
						  source_vector[*,n] = replicate(eps/eta,nfreq)
						  
						  n = n + 1
					 endif else begin
		   			  wset,1
		   			  plots,x,z,psym=4,col=2
		   			  wset,2
		   			  plots,y,z,psym=4,col=2
		   			  wset,3
		   			  plots,x,y,psym=4,col=2
					 endelse

; Has the ray reached the boundary?
					 if (z le hmin) then begin
						  rayout = 1
					 endif
					 i = i + 1

				endwhile
				
; Build an artificial image of the disk with a color proportional to the
; radius position (it is only because it is nice)
				image(loopx,loopy) = col
				
; And now, the formal solution along the ray that we have built
				n = n_elements(where(chi_vector(nfreq/2,*) ne 0.d0))
							
				if (n gt 1) then begin
					 ind = where(chi_vector(nfreq/2,*) ne 0.d0)
					 nind = n_elements(ind)
					 
; Try to avoid where the opacity is zero because it gives a non-defined source function
; because there is no background continuum. Only solve for those points having
; non-zero opacity
					 chi_vector = chi_vector(*,ind)
					 chi_vector = reverse(chi_vector,2)
					 
					 source_vector = source_vector(*,ind)
					 source_vector = reverse(source_vector,2)
					 
					 step_distance = ray_step * AU
					 distance = dindgen(nind) * step_distance
					 
; Formal solution of the RT equation
					 res = shortcar_source_linear_poly(distance, chi_vector, source_vector, $
					 	1.d0, 0.d0, tot_optical_depth)
					 	
; Save the last point along the rays.
					 image3(loopx,loopy,*) = res(*,nind-1)
					 
					 tot_opt_depth[loopx,loopy,*] = tot_optical_depth
					 					 	 	   		 
				endif
				
		  endfor
		  
	 endfor
	 
; The following lines convert the images to temperature units but may be deleted
	 temper = dblarr(npx,npx,nfreq)
	 Ta = dblarr(npx,npx,nfreq)
	 
	 for i = 0, nfreq-1 do begin
	 	  intensity = reform(image3(*,*,i))
	 	  ind = where(intensity ne 0.d0)
	 	  temp = intensity
; Image in temperature units	 
	 	  if (ind(0) ne -1) then begin
		   	temp(ind) = intfreq_to_t(nu,intensity(ind))
		   	temper(*,*,i) = temp
		  endif
	 endfor
	 
; And finally, save the results	 
	 save,image3,tot_opt_depth,radius,height,xpos,ypos,nfreq,npx,freq_axis,v_axis,velo,$
	 	filename=DIR+'/disk_angle'+strtrim(string(angle,FORMAT='(I2)'),2)+$
	 	'_vel_'+strtrim(string(v_Rin,FORMAT='(F4.1)'),2)+fileplus+'.idl'
	 		 	 
	 if (keyword_set(clean)) then delwins
	 		 	 
end