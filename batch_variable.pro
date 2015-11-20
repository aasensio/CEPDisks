@generate_variable
@calculate_variable
@synth

;************************************
; This script runs several disk models
; for case 1
;************************************
pro batch_variable_calculate1

	epsilons = ['1.d-1','1.d-2']
	xmols = [1.d-5,1.d-6,1.d-7,1.d-8,1.d-9,1.d-10]
	
	for j = 0, n_elements(xmols)-1 do begin
		xmol = xmols[j]
		for i = 0, n_elements(epsilons)-1 do begin			
			calculate_variable, 1, xmol, epsilons[i]
			read_variable, 1, xmol, epsilons[i]
		endfor
	endfor
	
	stop
end

;************************************
; This script runs several disk models
; for case 2
;************************************
pro batch_variable_calculate2
	
	epsilons = ['1.d-1','1.d-2']
	xmols = [1.d-5,1.d-6,1.d-7,1.d-8,1.d-9,1.d-10]
	
	for j = 0, n_elements(xmols)-1 do begin
		xmol = xmols[j]
		for i = 0, n_elements(epsilons)-1 do begin			
			calculate_variable, 2, xmol, epsilons[i]
			read_variable, 2, xmol, epsilons[i]
		endfor
	endfor
	
	stop
end


pro batch_synth1
		
	epsilons = ['1.d-1'] ;,'1.d-2']
	nepsilons = n_elements(epsilons)
	
	xmols = [1.d-5,1.d-6,1.d-7,1.d-8,1.d-9,1.d-10]
	nxmols = n_elements(xmols)
	
	angles = [1,30,60,90]
	nangles = n_elements(angles)
	
	vels = [0.d0,3.d0,5.d0,10.d0]
	nvels = n_elements(vels)
	
	disk_case = 1
	
	for loop_xmol = 0, nxmols-1 do begin
		xmol = xmols[loop_xmol]
		for loop_eps = 0, nepsilons-1 do begin
			epsilon = epsilons[loop_eps]
			
			root_output = 'output_variable/MODEL_standard400_CASE'+$
				strtrim(string(disk_case,FORMAT='(I1)'),2)+$
				'_xmol'+strtrim(string(xmol,FORMAT='(E7.1)'),2)+'/'
							
			for loop_vel = 0, nvels-1 do begin
				for loop_angles = 0, nangles-1 do begin					
					synth, root_output, vels[loop_vel], angles[loop_angles], $
						fileplus='_eps'+epsilon, npix=40
				endfor
			endfor
					
		endfor
	endfor
		
	stop
end

pro batch_synth2
		
	epsilons = ['1.d-1'] ;,'1.d-2']
	nepsilons = n_elements(epsilons)
	
	xmols = [1.d-5,1.d-6,1.d-7,1.d-8,1.d-9,1.d-10]
	nxmols = n_elements(xmols)
	
	angles = [1,30,60,90]
	nangles = n_elements(angles)
	
	vels = [0.d0,3.d0,5.d0,10.d0]
	nvels = n_elements(vels)
	
	disk_case = 2
	
	for loop_xmol = 0, nxmols-1 do begin
		xmol = xmols[loop_xmol]
		for loop_eps = 0, nepsilons-1 do begin
			epsilon = epsilons[loop_eps]
			
			root_output = 'output_variable/MODEL_standard400_CASE'+$
				strtrim(string(disk_case,FORMAT='(I1)'),2)+$
				'_xmol'+strtrim(string(xmol,FORMAT='(E7.1)'),2)+'/'
							
			for loop_vel = 0, nvels-1 do begin
				for loop_angles = 0, nangles-1 do begin					
					synth, root_output, vels[loop_vel], angles[loop_angles], $
						fileplus='_eps'+epsilon, npix=40
				endfor
			endfor
					
		endfor
	endfor
		
	stop
end