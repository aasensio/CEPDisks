@generate_flat
@calculate
@synth
;************************************
; This script runs several disk models
; The paramters are:
;    Rout_Rin, h_Rin, nh2dp_cte, [epsilon]
;************************************
pro batch_calculate1

	epsilons = ['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1']
; 	nmols = [1.d-3,1.d-2,1.d-1,1.d0]
	nmols = [1.d1,1.d2,1.d3]

	Rout_Rin = 10.d0
	h_Rin = 1.d0
	
	for j = 0, n_elements(nmols)-1 do begin
		nmol = nmols[j]
		for i = 0, 5 do begin
			generate_flat, Rout_Rin, h_Rin, nmol
			calculate, Rout_Rin, h_Rin, nmol, epsilons[i]
			read_flat, Rout_Rin, h_Rin, nmol, epsilons[i]
		endfor
	endfor
	
	stop
end

pro batch_calculate2

	epsilons = ['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1']
; 	nmols = [1.d-3,1.d-2,1.d-1,1.d0]
	nmols = [1.d1,1.d2,1.d3]

	Rout_Rin = 10.d0
	h_Rin = 10.d0
	
	for j = 0, n_elements(nmols)-1 do begin
		nmol = nmols[j]
		for i = 0, 5 do begin
			generate_flat, Rout_Rin, h_Rin, nmol
			calculate, Rout_Rin, h_Rin, nmol, epsilons[i]
			read_flat, Rout_Rin, h_Rin, nmol, epsilons[i]
		endfor
	endfor
	
	stop
end

pro batch_calculate3

	epsilons = ['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1']
; 	nmols = [1.d-3,1.d-2,1.d-1,1.d0]
	nmols = [1.d1,1.d2,1.d3]

	Rout_Rin = 10.d0
	h_Rin = 0.3d0
	
	for j = 0, n_elements(nmols)-1 do begin
		nmol = nmols[j]
		for i = 0, 5 do begin
			generate_flat, Rout_Rin, h_Rin, nmol
			calculate, Rout_Rin, h_Rin, nmol, epsilons[i]
			read_flat, Rout_Rin, h_Rin, nmol, epsilons[i]
		endfor
	endfor
	
	stop
end

pro batch_calculate4

	epsilons = ['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1']
 	nmols = [1.d-3,1.d-2,1.d-1,1.d0,1.d1,1.d2,1.d3]*10.d0

	Rout_Rin = 10.d0
	h_Rin = 0.03d0
	
	for j = 0, n_elements(nmols)-1 do begin
		nmol = nmols[j]
		for i = 0, 5 do begin
			generate_flat, Rout_Rin, h_Rin, nmol
			calculate, Rout_Rin, h_Rin, nmol, epsilons[i]
			read_flat, Rout_Rin, h_Rin, nmol, epsilons[i]
		endfor
	endfor
	
	stop
end


pro batch_synth1
	
; Rout/Rin=10     h/Rin=1
	Rout_Rin = 10.d0
	h_Rin = 1.d0
	
	epsilons = ['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1']
	neps = n_elements(epsilons)
	
	mols = [1.d-3,1.d-2,1.d-1,1.d0]
	mols = [1.d1,1.d2,1.d3]
	nmols = n_elements(mols)
	
	angles = [90,60,30,1]
	angles = [90,70,60,30,10]
	nangles = n_elements(angles)
	
	vels = [0.d0,1.d0,3.d0,5.d0,7.d0,10.d0]	
	nvels = n_elements(vels)

	for loop_mols = 0, nmols-1 do begin
		n_mol = mols[loop_mols]
		for loop_eps = 0, 5 do begin
			epsilon = epsilons[loop_eps]
	
			root_output = $
				"OUTPUT_FLAT/MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
				"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
				"_nmol_"+strtrim(string(n_mol,FORMAT='(E7.1)'))+"_eps_"+epsilon
			
			for loop_vel = 0, nvels-1 do begin
				for loop_angles = 0, nangles-1 do begin
					synth, root_output, vels[loop_vel], angles[loop_angles]
				endfor
			endfor
					
		endfor
	endfor
		
	stop
end

pro batch_synth2
	
; Rout/Rin=10     h/Rin=10
	Rout_Rin = 10.d0
	h_Rin = 10.d0
	
	epsilons = ['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1']
	neps = n_elements(epsilons)
	
	mols = [1.d-3,1.d-2,1.d-1,1.d0]
	mols = [1.d1,1.d2,1.d3]
	nmols = n_elements(mols)
	
	angles = [90,60,30,1]
	angles = [90,70,60,30,10]
	nangles = n_elements(angles)
	
	vels = [0.d0,1.d0,3.d0,5.d0,7.d0,10.d0]	
	nvels = n_elements(vels)

	for loop_mols = 0, nmols-1 do begin
		n_mol = mols[loop_mols]
		for loop_eps = 0, 5 do begin
			epsilon = epsilons[loop_eps]
	
			root_output = $
				"OUTPUT_FLAT/MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
				"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
				"_nmol_"+strtrim(string(n_mol,FORMAT='(E7.1)'))+"_eps_"+epsilon
			
			for loop_vel = 0, nvels-1 do begin
				for loop_angles = 0, nangles-1 do begin
					synth, root_output, vels[loop_vel], angles[loop_angles]
				endfor
			endfor
					
		endfor
	endfor
		
	stop
end

pro batch_synth3
	
; Rout/Rin=10     h/Rin=10
	Rout_Rin = 10.d0
	h_Rin = 0.3d0
	
	epsilons = ['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1']
	neps = n_elements(epsilons)
	
 	mols = [1.d-3,1.d-2,1.d-1,1.d0,1.d1,1.d2,1.d3]
	nmols = n_elements(mols)
	
	angles = [90,70,60,50,30,10,1]
; 	angles = [50] ;90,70,50,30,10]
	nangles = n_elements(angles)
	
	vels = [0.d0,1.d0,3.d0,5.d0,7.d0,10.d0]	
	nvels = n_elements(vels)

	for loop_mols = 0, nmols-1 do begin
		n_mol = mols[loop_mols]
		for loop_eps = 0, 5 do begin
			epsilon = epsilons[loop_eps]
	
			root_output = $
				"OUTPUT_FLAT/MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
				"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
				"_nmol_"+strtrim(string(n_mol,FORMAT='(E7.1)'))+"_eps_"+epsilon
			
			for loop_vel = 0, nvels-1 do begin
				for loop_angles = 0, nangles-1 do begin
					synth, root_output, vels[loop_vel], angles[loop_angles]
				endfor
			endfor
					
		endfor
	endfor
		
	stop
end

; Moshe wanted a model with h_Rin=0.03. In order to keep the same populations
; and optical depths, we divide the height by a factor 10 and multiply
; Aul by a factor 10
pro batch_synth4
	
; Rout/Rin=10     h/Rin=10
	Rout_Rin = 10.d0
	h_Rin = 0.03d0
	
	epsilons = ['1.d-3','3.d-3','1.d-2','3.d-2','1.d-1','3.d-1']
	neps = n_elements(epsilons)
	
 	mols = [1.d-3,1.d-2,1.d-1,1.d0,1.d1,1.d2,1.d3]*10.d0
	nmols = n_elements(mols)
	
	angles = [90,70,60,50,30,10,1]
; 	angles = [50] ;90,70,50,30,10]
	nangles = n_elements(angles)
	
	vels = [0.d0,1.d0,3.d0,5.d0,7.d0,10.d0]	
	nvels = n_elements(vels)

	for loop_mols = 0, nmols-1 do begin
		n_mol = mols[loop_mols]
		for loop_eps = 0, 5 do begin
			epsilon = epsilons[loop_eps]
	
			root_output = $
				"OUTPUT_FLAT/MODEL_RoutRin_"+strtrim(string(Rout_Rin,FORMAT='(E7.1)'),2)+$
				"_hRin_"+strtrim(string(h_Rin,FORMAT='(E7.1)'),2)+$
				"_nmol_"+strtrim(string(n_mol,FORMAT='(E7.1)'))+"_eps_"+epsilon
			
			for loop_vel = 0, nvels-1 do begin
				for loop_angles = 0, nangles-1 do begin
					synth, root_output, vels[loop_vel], angles[loop_angles]
				endfor
			endfor
					
		endfor
	endfor
		
	stop
end