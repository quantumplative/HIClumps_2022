import pandas as pd 
import numpy as np 
import pyfof_HI_clumps
outputs=['003840', '003840', '003968', '003840', '003968', '003840', '003840', '003840', '003840', '003840' ]
#models=['GM2', 'GM2noBHs', "GM2SI1", "GM3", "GM3noBHs", "GM3SI1", "P0", "P0noBHs","GM1", "GM1noBHs"]
#models=['P0', "P0noBHs"]
models = ['GM2']

for i in range(len(models)): 
	print('working on ' + models[i])
	halo, output = pyfof_HI_clumps.get_model_outputs_halos(models[i])
	
	for j in range(len(output)):
		if(len(str(output[j])) == 4):
			try:
				pyfof_HI_clumps.clump_grp(models[i], '00' + str(output[j]), 1, 4e-08, 2, halo[j])
				print('saving clump.grp and clump.grp.mass files')
			except: 
				print('something up with ' + models[i] + ' output ' + str(output[j]))
		else: 
			try:
				pyfof_HI_clumps.clump_grp(models[i], '000' + str(output[j]), 1, 4e-08, 2, halo[j])
				print('saving clump.grp and clump.grp.mass files')
			except: 
				print('something up with ' + models[i] + ' output ' + str(output[j]))
