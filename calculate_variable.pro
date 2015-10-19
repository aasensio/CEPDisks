;************************************
; This script solves the radiative transfer problem for all the
; radial shells
;************************************
pro calculate_variable, disk_case, xmol, epsilon

	root_model = 'models_variable/MODEL_standard400_CASE'+$
		strtrim(string(disk_case,FORMAT='(I1)'),2)+$
		'_xmol'+strtrim(string(xmol,FORMAT='(E7.1)'),2)+$
		'/'
	root_output = 'output_variable/MODEL_standard400_CASE'+$
		strtrim(string(disk_case,FORMAT='(I1)'),2)+$
		'_xmol'+strtrim(string(xmol,FORMAT='(E7.1)'),2)+$
		'/'
		
	nr = 79

	file = strarr(35)	
	test = file_test(root_output,/directory)
	if (test eq 0) then begin
		file_mkdir, root_output
	endif

; Open the original config file
	openr,2,'config/2-level.config.original'
	readf,2,file
	close,2

; Loop over all the shells
	for i = 0, nr-1 do begin
		print, i
		
		file[7] = "'"+root_output+"radius_"+strtrim(string(i),2)+"_eps"+epsilon+"'"
		file[10] = "'"+root_model+"radius_"+strtrim(string(i),2)+".model'"
		file[13] = "'molecules/VARIABLE/2-level_eps"+epsilon+"_case"+$
			strtrim(string(disk_case,FORMAT='(I1)'),2)+".model'"

		openw,2,'config/2-level.config',width=132
		for k = 0, 34 do begin
			printf,2,file[k]
		endfor
		close,2

; Run the cep code
		spawn,'./cep'
				
	endfor
end

pro lb, u, n
	caca = ''
	for i = 0, n-1 do begin
		readf,u,caca		
	endfor
end

;*********************************
; Read the results of the calculations and generate the file results.idl
;*********************************
pro read_variable, disk_case, xmol, epsilon

	root_model = 'models_variable/MODEL_standard400_CASE'+$
		strtrim(string(disk_case,FORMAT='(I1)'),2)+$
		'_xmol'+strtrim(string(xmol,FORMAT='(E7.1)'),2)+$
		'/'
	root_output = 'output_variable/MODEL_standard400_CASE'+$
		strtrim(string(disk_case,FORMAT='(I1)'),2)+$
		'_xmol'+strtrim(string(xmol,FORMAT='(E7.1)'),2)+$
		'/'
				
	a = ddread('disk_models/standard_400au.dat',offset=3,/countal)
	rr = reform(a[0,*])
	hh = reform(a[1,*])
		
	ind = uniq(rr)
   radius = rr[ind]
   
   nr = n_elements(radius)
   nh = 99
   
	pop = dblarr(nr,nh,2)
	depart = dblarr(nr,nh,2)
	tex = dblarr(nr,nh)
	h = dblarr(nr,nh)
	r = dblarr(nr,nh)
	tkin = dblarr(nr,nh)
	temp = dblarr(3,nh)
	temp2 = dblarr(4,nh)

	for i = 0, nr-1 do begin
		file = root_output+'radius_'+strtrim(string(i),2)+'_eps'+epsilon+'.out'
			
		openr,2,file
		lb,2,118
		readf,2,temp
		pop[i,*,*] = transpose(temp[1:2,*])
		lb,2,3
		lb,2,105
		readf,2,temp
		depart[i,*,*] = transpose(temp[1:2,*])
		lb,2,236
		readf,2,temp
		tex[i,*] = reform(temp[2,*])
		close,2
		
		if (i eq 0) then begin
			print, 'tau_vert_max = ', max(temp[0,*])			
		endif
						
		file = root_model+'radius_'+strtrim(string(i),2)+'.model'
		openr,2,file
		lb,2,4
		readf,2,temp2		
		for j = nh/2+1, nh-1 do h[i,j] = h[i,j-1] + temp2[0,j] / 1.5d13		
		h[i,0:nh/2-1] = -reverse(h[i,nh/2+1:nh-1],2)
		close,2		
		r[i,*] = replicate(radius[i],nh)
		tkin[i,*] = reform(temp2[1,*])
				
	endfor
	
	minr = min(r)
	maxr = max(r)
	minh = min(h)
	maxh = max(h)

	nh = n_elements(r[0,*])
	
	height = h
	radius = r
	popul = pop

	save,height,radius,popul,tkin,filename=root_output+'results_eps'+epsilon+'.idl'		
	
end