pro calculate_all_flat
	calculate, 10.d0, 1.d0, 1.d-3, '3.d-3'
	read_flat, 10.d0, 1.d0, 1.d-3, '3.d-3'
	
	calculate, 10.d0, 1.d0, 1.d-2, '3.d-3'
	read_flat, 10.d0, 1.d0, 1.d-2, '3.d-3'
end

;************************************
; This script solves the radiative transfer problem for all the
; radial shells
;************************************
pro calculate, Rout_Rin, h_Rin, nh2dp_cte, epsilon

	root_model = "models_flat/MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
		"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
		"_nmol_"+strtrim(string(nh2dp_cte,FORMAT='(E7.1)'))+"/"
	root_output = "output_flat/MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
		"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
		"_nmol_"+strtrim(string(nh2dp_cte,FORMAT='(E7.1)'))+"_eps_"+epsilon+"/"
		
	nr = 80

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
	for i = 0, 0 do begin
		print, i
		
		file[7] = "'"+root_output+"radius_"+strtrim(string(i),2)+"'"
		file[10] = "'"+root_model+"radius_"+strtrim(string(i),2)+".model'"
		file[13] = "'molecules/2-level_eps"+epsilon+".model'"

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
pro read_flat, Rout_Rin, h_Rin, nh2dp_cte, epsilon

	root_model = "models_flat/MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
		"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
		"_nmol_"+strtrim(string(nh2dp_cte,FORMAT='(E7.1)'))+"/"
	root_output = "output_flat/MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
		"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
		"_nmol_"+strtrim(string(nh2dp_cte,FORMAT='(E7.1)'))+"_eps_"+epsilon+"/"
		
; These values are fixed so that the 2-level models give the expected
; values for epsilon

	Tcte = 15.d0
	nh2_cte = 5.d5	
	Rin = 40.d0	
				
	Rout = Rout_Rin * Rin
	h_tot = h_Rin * Rin
	
	nr = 80
	nh = 100
			
	radius = dindgen(nr)/(nr-1.d0)*(Rout-Rin) + Rin

	pop = dblarr(nr,nh,2)
	depart = dblarr(nr,nh,2)
	tex = dblarr(nr,nh)
	h = dblarr(nr,nh)
	r = dblarr(nr,nh)
	tkin = dblarr(nr,nh)
	temp = dblarr(3,nh)
	temp2 = dblarr(4,nh)

	for i = 0, nr-1 do begin
		file = root_output+'radius_0.out'
		openr,2,file
		lb,2,119
		readf,2,temp
		pop[i,*,*] = transpose(temp[1:2,*])
		lb,2,3
		lb,2,106
		readf,2,temp
		depart[i,*,*] = transpose(temp[1:2,*])
		lb,2,238
		readf,2,temp
		tex[i,*] = reform(temp[2,*])
		close,2
		
		if (i eq 0) then begin
			print, 'tau_vert_max = ', max(temp[0,*])
			print, 'approx. tau_horiz_max (if constant) = ', $
				max(temp[0,*]) * (Rout_Rin-1.d0)/h_Rin
		endif
				

		file = root_model+'radius_'+strtrim(string(i),2)+'.model'
		openr,2,file
		lb,2,4
		readf,2,temp2
		htot = total(temp2[0,*])
		h[i,0] = htot / 1.5d13 / 2.d0
		for j = 1, nh-1 do begin
			h[i,j] = h[i,j-1] - temp2[0,j] / 1.5d13
		endfor
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

	save,height,radius,popul,tkin,filename=root_output+'results.idl'
	
end