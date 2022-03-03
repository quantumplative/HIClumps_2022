import os
import numpy as np 
import pandas as pd 
import pyfof_HI_clumps 
import matplotlib.pyplot as plt
from statistics import mode

def adjacent_timestep_timescales(iord_data, model_timetrace_info):
	adjacent_clumpiord_below = []
	delta_t = []

	output = np.array(iord_data['output'])
	iord_data_clumpID = np.array(iord_data['clumpID'])
	iord_data_clumpgr = np.array(iord_data['clump.grp'])
	iord_data_clumpM  = np.array(iord_data['clump.grp.mass'])
	#print('iord_data_clumpM ', iord_data_clumpM)
	grps = []
	clumpmasses_4iord = [] 
	for j in range(len(iord_data)):
		if(iord_data_clumpID[j] == 1.0 and j != len(iord_data)-1):
			grps.append(iord_data_clumpgr[j]) 
			if(iord_data_clumpID[j+1] == 1.0):
				adjacent_clumpiord_below.append(1)
			else:
				adjacent_clumpiord_below.append(0)
			clumpmasses_4iord.append(iord_data_clumpM[j])
			delta_t.append(abs(float(model_timetrace_info[model_timetrace_info['output'] == output[j]]['time']) - float(model_timetrace_info[model_timetrace_info['output'] == output[j+1]]['time'])))
		elif(iord_data_clumpID[j] == 1.0 and j == len(iord_data)-1):
			grps.append(iord_data_clumpgr[j]) 
			clumpmasses_4iord.append(iord_data_clumpM[j])
			delta_t.append(abs(float(model_timetrace_info[model_timetrace_info['output'] == output[j]]['time']) - float(model_timetrace_info[model_timetrace_info['output'] == output[j-1]]['time'])))
		else:
			adjacent_clumpiord_below.append(0)
			delta_t.append(0)
	clumpmasses_4iord = np.array(clumpmasses_4iord)
	mean_mass = clumpmasses_4iord.mean()
	min_mass  = min(clumpmasses_4iord)
	max_mass  = max(clumpmasses_4iord)
	std_mass  = np.std(clumpmasses_4iord)
	print('clumpmasses_4iord ', clumpmasses_4iord)
	return np.sum(delta_t), len(np.unique(grps)), mean_mass, std_mass, min_mass, max_mass

def model_adjacent_timescales_upperbound(model): 
	model_timetrace_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')	
	df = pd.read_pickle(model + '_alliords_ever_inclumps_ID.pkl')
	iords = np.unique(df['iord'])
	outfn = model + '_alliords_ever_inclumps_timescales.dat'
	
	#check if outfile already exists	
	if os.path.exists(outfn):
		print('outfn ' + outfn + ' already exists')
		outfile = open(outfn, 'r')
		lines = [line for line in outfile]
		outfile.close()
		
		#finds what the last group in the file is 	
		if len(lines) == 1 :
			#file only has column names wriitten
			last_iord_in_file = 0
		else:
			last_iord_in_file = len(lines) - 1
	#creates an outfile if one doesnt already exist
	else:
		print('creating outfile ' + outfn)
		last_iord_in_file = 0
		with open(outfn, 'a') as outfile:
			outfile.write('model iord timeinclumps nclumps mean_mass std_mass min_mass max_mass\n')

	for i in range(last_iord_in_file, len(iords)):
		print('working on iord = ' +  str(iords[i]) + ' ' + str(i+1) + ' out of ' + str(len(iords)))
		iord_i_data = df[df['iord'] == iords[i]]
		t, ngrps, mean_mass, std_mass, min_mass, max_mass = adjacent_timestep_timescales(iord_i_data, model_timetrace_info)
		#iords_timeinclumps.append(t)
		#grps_per_iords.append(ngrps)
		with open(outfn, 'a') as outfile: 
			outfile.write('%s %i %.2f %i %.2f %.2f %.2f %.2f\n'%(model, iords[i], t, ngrps, mean_mass, std_mass,  min_mass, max_mass))
