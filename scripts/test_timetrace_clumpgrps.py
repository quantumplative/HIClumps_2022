import os
import pyfof
import math 
import pyfof_HI_clumps 
import pynbody 
import numpy as np
import pandas as pd 
from scipy.stats import skew
import track_halo_number
import seaborn as sns
import matplotlib.pyplot as plt 
from statistics import mode	
import BullockMuller2004_profiles as BMP 

mp = pynbody.array.SimArray(1.67262e-24)
mp.units = 'g'
kb = pynbody.array.SimArray(1.38 * 1e-16)
kb.units =  'erg K**-1'
G = pynbody.array.SimArray(6.67*1e-8) 
G.units = 'cm**3 g**-1 s**-2'

def grpclump_masstrace(model, si, ei): 
	halo, output = pyfof_HI_clumps.get_model_outputs_halos(model)
	start_output = output[si]
	output = output[si:ei]
	df = open_timetrace_grpclump_dat(model, si, ei)
	print('df = ', df)
	output_df = df[df['output1'] == start_output]
	print('output_df = ', output_df)
	clumpgrp_df = output_df[output_df['grpmatch'] > 0]
	grpmatch = clumpgrp_df['grpmatch'].to_numpy()
	print('grpmatch = ', grpmatch)
	grpmass  = clumpgrp_df['mass[Msol]'].to_numpy()
	grpdist  = clumpgrp_df['r[kpc]'].to_numpy()
	trace_mass     = []
	trace_dist     = []
	trace_grpmatch = []
	trace_output   = [] 
	for i in range(len(grpmatch)):
		mass_ij = [] 
		grp_ij  = [] 
		output_ij = []
		dist_ij = []
		#for output index 0 
		#the initial group you are tracing from the start_output is the first entry in grp_ij 
		grp_ij.append(grpmatch[i]) 
		mass_ij.append(grpmass[i])
		dist_ij.append(grpdist[i])
		output_ij.append(start_output) 
		grpmatch_i = grpmatch[i]
		j = 0 
		print('grpmatch[' + str(i) + '] = ' + str(grpmatch[i]))
		while(grpmatch_i > 0 and j < ei-si-1): 
			traceoutput = output[1+j] 
			traceout_df = df[df['output1'] == traceoutput] 
			#should have len 1. every grp is unique in a given output 
			#grptrace_df = traceout_df[traceout_df['grp'] == grpmatch[i]]
			grptrace_df = traceout_df[traceout_df['grp'] == grpmatch_i]
			print('grptrace_df = ', grptrace_df)
			try:
				grp_ij.append(grptrace_df['grpmatch'].to_numpy()[0])
				mass_ij.append(grptrace_df['mass[Msol]'].to_numpy()[0])
				dist_ij.append(grptrace_df['r[kpc]'].to_numpy()[0])
				output_ij.append(traceoutput)
				grpmatch_i = grptrace_df['grpmatch'].to_numpy()[0]
			except: 
				print('weird grp number')
			j = j + 1
		trace_mass.append(mass_ij)
		#array of arryas with group numbers. each array is a group traced across timesteps to new group numbers  
		trace_grpmatch.append(grp_ij)
		trace_output.append(output_ij)
		trace_dist.append(dist_ij)
	return np.array(trace_mass), np.array(trace_dist), np.array(trace_grpmatch), np.array(trace_output)

def open_timetrace_grpclump_dat(model, si, ei, metals):
	halo, output = pyfof_HI_clumps.get_model_outputs_halos(model)
	output = output[si:ei]
	halo = halo[si:ei]
	fn2 = pyfof_HI_clumps.get_fn(model)
	for j in range(len(output)-1):
		#df_fn = fn2 + 'timetrace_grps_lowtohigh_outputs' + str(output[j]) + '_' + str(output[j+1]) + '_wmass.dat'
		if metals: 
			df_fn = fn2 + 'timetrace_grps_lowtohigh_outputs' + str(output[j]) + '_' + str(output[j+1]) + '_clump_properties_metals_fixvr.dat'
		else:
			df_fn = fn2 + 'timetrace_grps_lowtohigh_outputs' + str(output[j]) + '_' + str(output[j+1]) + '_clump_properties.dat'
		try:
			df = pd.read_table(df_fn, sep = '\s+')
		
			if(j==0):
				dfs = df
			else:
				dfs = pd.concat([df, dfs])
		except: 
			print(df_fn + ' is empty or DNE')
	return dfs

def timetrace_grpclumps(model, link_len, HI_cut, min_m, si): 
	fn2 = pyfof_HI_clumps.get_fn(model)
	
	#get progenitor halos and connected outputs for model from tangos	
	halo, output = pyfof_HI_clumps.get_model_outputs_halos(model)
	halo = halo[si:si+2]
	output = output[si:si+2] 
	print('output = ', output) 
	print('halo = ', halo) 
	mp = pynbody.array.SimArray(1.67262e-27)
	mp.units = 'kg'
	for j in range(len(output)-1):
		outfn = fn2 + 'timetrace_grps_lowtohigh_outputs' + str(output[j]) + '_' + str(output[j+1]) + '_clump_properties.dat'
		with open(outfn, 'a') as outfile:
			outfile.write('model output1 output2 grp grpmatch fraction mass[Msol] r[kpc] log10avg_pressure[Msols**-2kpc**-1] temp[K] log10rho[kgcm**-3] log10nHI size[kpc]\n')
		print('j = ', j)
		print('----------------------working on output ' + str(output[j]) + ' halo ' + str(halo[j]) + '-------------------------')
		#try: 
		if(len(str(output[j+1])) == 4):
			tsoutfn = fn2 + '00' + str(output[j+1])
		else: 
			tsoutfn = fn2 + '000' + str(output[j+1])
		print('tsoutfn = ' + str(tsoutfn))
		if(len(str(output[j])) == 4): 
			print('trying to load grps and iords for ' + '00' + str(output[j]))
			pyfof_grps, grp_iords, s2 = clumps_and_iords(model, '00' + str(output[j]), link_len, HI_cut, min_m, halo[j])
		else:
			print('trying to load grps and iords for ' + '000' + str(output[j]))
			pyfof_grps, grp_iords, s2 = clumps_and_iords(model, '000' + str(output[j]), link_len, HI_cut, min_m, halo[j])
		print('got ' + str(len(pyfof_grps)) + ' grps and iords for output1 = ' + str(output[j]))
		#clump grp ID associated with every gas particle in a given time step 
		clump_grp = np.loadtxt(tsoutfn + '.clump.grp')
		s = pynbody.load(tsoutfn)
		print('loaded s for output2 = ' + str(output[j+1]))
		print('working on tracing particles from output1 = ' + str(output[j]) + ' into output2 = ' + str(output[j+1]))	
		s.physical_units()
		s.g['clump.grp'] = clump_grp
		print('assigned clump.grp array to s')
		#go through all of the groups in a given output j
		print('len of grp_iords = ' + str(len(grp_iords))) 
		print('grp_iords = ' + str(grp_iords)) 
		for i in range(len(grp_iords)):
			print('working on grp with index ' + str(i))
			print('grp_iords[' + str(i) + '] = ' + str(grp_iords[i]))
			grpmap = np.in1d(s.g['iord'], grp_iords[i])
			unigrp = np.unique(s.g[grpmap]['clump.grp'])
			s2grpmap = np.in1d(s2.g['iord'], grp_iords[i])
			print('working on grp with index ' + str(i))
			print("partices in adj ts grps " + str(s.g[grpmap]['clump.grp']))
			print("unique grps in adj ts =  " + str(unigrp))	
			try:
				matching_grp = mode(s.g[grpmap]['clump.grp'])
				grp = i+1
			except: 
				if(len(unigrp)==len(s.g[grpmap])):
					if(len(s.g[grpmap]) > 1): 
						matching_grp = min(unigrp)
					else: 
						matching_grp = -1.0 #iords dont exist in the traced to ts
					grp = i+1
				elif(len(unigrp)==2):
					matching_grp = min(unigrp)
					grp = i+1
				elif(len(unigrp)==1):
					matching_grp = unigrp[0]
					grp = i+1
				elif(len(unigrp) > 1): 
					matching_grp = max(unigrp)  
					grp = i+1 
			print('assigned matching grp')
			print('matching_grp = ', matching_grp)
			if(len(s.g[grpmap]) > 0):
				fraction  = len(s.g[grpmap]['clump.grp'][s.g[grpmap]['clump.grp'] == matching_grp])/len(s.g[grpmap])
			else: 
				fraction = 1.0
			print('fraction  = ' + str(fraction))
			print('mass of original grp = ' + str(s.g[grpmap]['mass']))
			print('particle mass with matching_grp number = ' + str(s.g[grpmap]['mass'][s.g[grpmap]['clump.grp'] == matching_grp]))
			pyfof_grpi = s2.g[s2grpmap]
			avg_temp = pyfof_grpi['temp'].in_units('K').mean()
		
			grp_rho = np.array(pyfof_grpi['rho'].in_units('kg cm**-3'))
			avg_rho = np.log10(grp_rho).mean()
			grp_n = (pyfof_grpi['HI'] * grp_rho / mp)
			avg_n = np.log10(grp_n).mean()
			
			mass = sum(s2.g[s2grpmap]['mass'])
			pressure_i = np.log10(s2.g[s2grpmap]['p'].in_units('Msol s**-2 kpc**-1')).mean()
			print('pressure = ', s2.g[s2grpmap]['p'])
			print('log 10 avg pressure = ', np.log10(s2.g[s2grpmap]['p'].in_units('Msol s**-2 kpc**-1')).mean())
			pyfof_comi = pynbody.analysis.halo.center_of_mass(s2.g[s2grpmap])
			pyfof_posi = s2.g[s2grpmap]['pos']
			pyfof_pos_mcom = pyfof_posi - pyfof_comi
			dist_from_grpcom = max(((pyfof_pos_mcom[:, 0])**2 + (pyfof_pos_mcom[:, 1])**2 + (pyfof_pos_mcom[:, 2])**2)**(0.5))
			print('size = ' + str(dist_from_grpcom) + ' kpc')
			cgm_com_dist   = ((pyfof_comi[0])**2 + (pyfof_comi[1])**2 + (pyfof_comi[2])**2)**(0.5)  
			print('mass = ' + str(mass))
			with open(outfn, 'a') as outfile: 
				outfile.write('%s %s %s %i %i %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n' % (model, output[j], output[j+1], grp, matching_grp, fraction, mass, cgm_com_dist, pressure_i, avg_temp, avg_rho, avg_n, dist_from_grpcom))	
		print('done tracing grps from output1 = ' + str(output[j])) 

def timetrace_grpclumps_metals_vr(model, link_len, HI_cut, min_m, si): 
	'''
	match and trace groups through time. find mass of clumps and radial COM velocity and average metallicity of the clumps. This data should be used in combo with data produced using timetrace_grpclumps
	'''
	fn2 = pyfof_HI_clumps.get_fn(model)
	
	#get progenitor halos and connected outputs for model from tangos	
	halo, output = pyfof_HI_clumps.get_model_outputs_halos(model)
	halo = halo[si:si+2]
	output = output[si:si+2] 
	print('output = ', output) 
	print('halo = ', halo) 
	mp = pynbody.array.SimArray(1.67262e-27)
	mp.units = 'kg'
	for j in range(len(output)-1):
		outfn = fn2 + 'timetrace_grps_lowtohigh_outputs' + str(output[j]) + '_' + str(output[j+1]) + '_clump_properties_metals_fixvr.dat'
		with open(outfn, 'a') as outfile:
			outfile.write('model output1 output2 grp grpmatch fraction mass[Msol] r[kpc] vr[kms^-1] logmetals size[kpc] oxymass[Msol] femass[Msol]\n')
		print('j = ', j)
		print('----------------------working on output ' + str(output[j]) + ' halo ' + str(halo[j]) + '-------------------------')
		#try: 
		if(len(str(output[j+1])) == 4):
			tsoutfn = fn2 + '00' + str(output[j+1])
		else: 
			tsoutfn = fn2 + '000' + str(output[j+1])
		print('tsoutfn = ' + str(tsoutfn))
		if(len(str(output[j])) == 4): 
			print('trying to load grps and iords for ' + '00' + str(output[j]))
			pyfof_grps, grp_iords, s2 = clumps_and_iords(model, '00' + str(output[j]), link_len, HI_cut, min_m, halo[j])
		else:
			print('trying to load grps and iords for ' + '000' + str(output[j]))
			pyfof_grps, grp_iords, s2 = clumps_and_iords(model, '000' + str(output[j]), link_len, HI_cut, min_m, halo[j])
		print('got ' + str(len(pyfof_grps)) + ' grps and iords for output1 = ' + str(output[j]))
		#clump grp ID associated with every gas particle in a given time step 
		clump_grp = np.loadtxt(tsoutfn + '.clump.grp')
		s = pynbody.load(tsoutfn)
		print('loaded s for output2 = ' + str(output[j+1]))
		print('working on tracing particles from output1 = ' + str(output[j]) + ' into output2 = ' + str(output[j+1]))	
		s.physical_units()
		s.g['clump.grp'] = clump_grp
		print('assigned clump.grp array to s')
		#go through all of the groups in a given output j
		print('len of grp_iords = ' + str(len(grp_iords))) 
		print('grp_iords = ' + str(grp_iords)) 
		for i in range(len(grp_iords)):
			print('working on grp with index ' + str(i))
			print('grp_iords[' + str(i) + '] = ' + str(grp_iords[i]))
			grpmap = np.in1d(s.g['iord'], grp_iords[i])
			unigrp = np.unique(s.g[grpmap]['clump.grp'])
			s2grpmap = np.in1d(s2.g['iord'], grp_iords[i])
			print('working on grp with index ' + str(i))
			print("partices in adj ts grps " + str(s.g[grpmap]['clump.grp']))
			print("unique grps in adj ts =  " + str(unigrp))	
			try:
				matching_grp = mode(s.g[grpmap]['clump.grp'])
				grp = i+1
			except: 
				if(len(unigrp)==len(s.g[grpmap])):
					if(len(s.g[grpmap]) > 1): 
						matching_grp = min(unigrp)
					else: 
						matching_grp = -1.0 #iords dont exist in the traced to ts
					grp = i+1
				elif(len(unigrp)==2):
					matching_grp = min(unigrp)
					grp = i+1
				elif(len(unigrp)==1):
					matching_grp = unigrp[0]
					grp = i+1
				elif(len(unigrp) > 1): 
					matching_grp = max(unigrp)  
					grp = i+1 
			print('assigned matching grp')
			print('matching_grp = ', matching_grp)
			if(len(s.g[grpmap]) > 0):
				fraction  = len(s.g[grpmap]['clump.grp'][s.g[grpmap]['clump.grp'] == matching_grp])/len(s.g[grpmap])
			else: 
				fraction = 1.0
			print('fraction  = ' + str(fraction))
			print('mass of original grp = ' + str(s.g[grpmap]['mass']))
			print('particle mass with matching_grp number = ' + str(s.g[grpmap]['mass'][s.g[grpmap]['clump.grp'] == matching_grp]))
			pyfof_grpi = s2.g[s2grpmap]
			
			mass = sum(s2.g[s2grpmap]['mass'])
			oxmass = sum(s2.g[s2grpmap]['OxMassFrac']*s2.g[s2grpmap]['mass'])
			femass = sum(s2.g[s2grpmap]['FeMassFrac']*s2.g[s2grpmap]['mass'])
			
			mean_metals = sum(pyfof_grpi['mass']*pyfof_grpi['metals'])/sum(pyfof_grpi['mass'])

			pyfof_comi = pynbody.analysis.halo.center_of_mass(s2.g[s2grpmap])
			pyfof_posi = s2.g[s2grpmap]['pos']
			
			pyfof_pos_mcom = pyfof_posi - pyfof_comi
			dist_from_grpcom = max(((pyfof_pos_mcom[:, 0])**2 + (pyfof_pos_mcom[:, 1])**2 + (pyfof_pos_mcom[:, 2])**2)**(0.5))
			pyfof_com_veli = pynbody.analysis.halo.center_of_mass_velocity(s2.g[s2grpmap]).in_units('km s^-1')
			
			cgm_com_dist   = ((pyfof_comi[0])**2 + (pyfof_comi[1])**2 + (pyfof_comi[2])**2)**(0.5)  
			
			#rad_comvel     = ((pyfof_com_veli[0])**2 + (pyfof_com_veli[1])**2 + (pyfof_com_veli[2])**2)**(0.5)  
			rad_comvel = (pyfof_comi * pyfof_com_veli).sum(axis=0) / cgm_com_dist
			#print('com times com_vel = ', (pyfof_comi * pyfof_com_veli))
			#print('com dot com_vel = ', (pyfof_comi * pyfof_com_veli).sum(axis=0))
			print('mass = ' + str(mass))
			
			with open(outfn, 'a') as outfile: 
				outfile.write('%s %s %s %i %i %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n' % (model, output[j], output[j+1], grp, matching_grp, fraction, mass, cgm_com_dist, rad_comvel, np.log10(mean_metals), dist_from_grpcom, oxmass, femass))	
		print('done tracing grps from output1 = ' + str(output[j])) 

def clumps_and_iords(model, output, link_len, HI_cut, min_m, halo):
	fn = pyfof_HI_clumps.get_fn(model)
	s2 = pynbody.load(fn + output)
	s2.physical_units()
	pyfof_HI_clumps.add_NHI(s2)
	if (model == 'GM3noBHs' or model =='GM3'):
		h2 = s2.halos(write_fpos = False)
	else:
		h2 = s2.halos()
	h21 = h2[halo]
	pynbody.analysis.halo.center(h21, mode = 'com')
	h21_rcut = h21.g[h21.g['r'].in_units('kpc') > 15]
	
	s2 = pyfof_HI_clumps.cut_gal(s2)
	z = s2.properties['z']

	h21_HIcut = h21_rcut[h21_rcut['HIn']  > HI_cut]
	
	pyn_data =  np.array(h21_HIcut['pos'])
	#print('pyn_data = ', pyn_data)
	pyn_groups = pyfof.friends_of_friends(pyn_data, link_len)
	pyfof_grp_minm = pyfof_HI_clumps.grp_min_m(pyn_groups, min_m)
	pyfof_iord = pyfof_HI_clumps.get_iords(s2.g, h21_HIcut, pyfof_grp_minm)
	return pyfof_grp_minm, pyfof_iord, s2 

def plot_grpclump_massloss(model, si, ei):
	model_timetrace_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
	trace_mass, trace_dist, trace_grpmatch, trace_output = grpclump_masstrace(model, si, ei) 
	for i in range(len(trace_grpmatch)):
		#if(len(trace_grpmatch[i])>=3):
		time = []
		for t in trace_output[i]:
			time.append(model_timetrace_info[model_timetrace_info['output'] == t]['time'].to_numpy()[0])
		print('time = ', time)
		fig = plt.figure(figsize=(12,10))
		plt.plot(time, trace_mass[i]/max(trace_mass[i]), lw=2, c='k')
		plt.xlabel(r'cosmic time[Gyr]', fontsize=16)
		plt.ylabel(r'Cold clump mass / Max cold clump mass', fontsize=16)
		plt.tick_params(labelsize=16)
		plt.hlines(0.5, plt.xlim()[0]-.25, plt.xlim()[1]+.25)
		plt.xlim([plt.xlim()[0]+.25, plt.xlim()[1]-0.25])
		plt.annotate(r'M$_{\rm clump}$(t=t$_{0}$) = ' + "{:.2e}".format(trace_mass[i][0]) + ' M$_{\odot}$', xy=(0.6, 0.9), xycoords='axes fraction', fontsize=16)
		plt.savefig(model + '_massloss_startoutput_' + str(trace_output[si]) + '_tracegrp' + str(trace_grpmatch[i][0]) + '.pdf')
		plt.show()

def plot_grpclump_integratedtime(model, si, ei):
	model_timetrace_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
	trace_mass, trace_dist, trace_grpmatch, trace_output = grpclump_masstrace(model, si, ei) 
	lower_time = []
	lower_mass = []
	full_time  = []
	full_mass  = []
	for i in range(len(trace_grpmatch)):
		#if(len(trace_grpmatch[i])>=3):
		time = []
		z = [] 
		for t in trace_output[i]:
			time.append(model_timetrace_info[model_timetrace_info['output'] == t]['time'].to_numpy()[0])
			z.append(model_timetrace_info[model_timetrace_info['output'] == t]['z'].to_numpy()[0])
		print('time = ', time)
		max_mass  = max(trace_mass[i])
		norm_mass = np.array(trace_mass[i])/max_mass
		if(len(norm_mass)==2):
			#can only get a lower time bound with 2 timesteps! 
			ncross = 0  
			if(norm_mass[1] < 0.5 and norm_mass[0] > 0.5):
				slope = (norm_mass[1] - norm_mass[0]) / (time[1] - time[0])
				b = norm_mass[1]
				xoff = time[1]
				lower_time.append((0.5-b)/slope + xoff - time[0])
				lower_mass.append(max(trace_mass[i]))
				ncross = ncross + 1
				print('intercept time = ', (0.5-b)/slope + xoff)
				print('integrated lower time = ', (0.5-b)/slope + xoff - time[0])
			elif(norm_mass[1] > 0.5 and norm_mass[0] > 0.5):
				lower_time.append(time[1]-time[0])
				lower_mass.append(max_mass)
				ncross = ncross + 1
				print('integrated lower time = ', time[1]-time[0]) 
			else:
				print('there is a combo you have not considered.')
				print('norm_mass[0] = ', norm_mass[0])
				print('norm_mass[1] = ', norm_mass[1])
		else: 
			print('tracing group number through more than 2 timesteps')
			ncross = 0  
			for j in range(len(norm_mass)-1):
				ti = 0 
				if(norm_mass[j] < 0.5 and norm_mass[j+1] > 0.5):
					xoff = time[j]
					b = norm_mass[j]
					slope = (norm_mass[j+1] - norm_mass[j]) / (time[j+1] - time[j])  
					print('b = ', b)
					print('slope = ', slope) 
					print('intercept time = ', (0.5-b)/slope + time[j])
					if(ncross == 0): 
						t1 = (0.5-b)/slope + xoff
					else:
						ti = (0.5-b)/slope + xoff
					ncross = ncross + 1 
				elif(norm_mass[j] > 0.5 and norm_mass[j+1] < 0.5): 
					xoff = time[j+1]
					b = norm_mass[j+1]
					slope = (norm_mass[j+1] - norm_mass[j]) / (time[j+1] - time[j])  
					print('b = ', b)
					print('slope = ', slope) 
					print('intercept time = ', (0.5-b)/slope + time[j+1])
					if(ncross == 0): 
						t1 = (0.5-b)/slope + xoff
					elif(j == len(norm_mass) - 1):
						tf = (0.5-b)/slope + xoff
					elif(j != len(norm_mass)-1 and ncross > 2):
						ti = (0.5-b)/slope + xoff
					ncross = ncross + 1 
				print('ncross = ', ncross)
			if(ncross > 1):
				tf = (0.5-b)/slope + xoff
			elif(ncross == 0): 
				tf = time[-1]
			elif(ncross == 1): 
				tf = t1
				t1 = time[0]
			#if(ncross > 1 and norm_mass[0] > min(norm_mass) and norm_mass[-1] > min(norm_mass)):
				#print('found a clump that is losing mass')
				#lower_time.append((time[-1]-ti) + (t1 - time[0]))
				#lower_mass.append(max_mass) 
			#else: 	
			full_time.append((tf-ti) + (ti-t1))
			full_mass.append(max_mass)
			print('integrated time = ', (tf-ti) + ti-t1)			

		
	fig = plt.figure(figsize=(12,10))
	plt.scatter(full_mass, full_time, s=200, c='k', marker='o')
	plt.scatter(lower_mass, lower_time, s=200, c='indianred', marker='^')
	plt.semilogx()
	plt.tick_params(labelsize=16)
	plt.xlabel(r'log$_{10}$ M$_{\rm clump}$ [M$_{\odot}$]', fontsize=16)
	plt.ylabel(r't$_{\rm clump}$ [Gyr]', fontsize=16)
	plt.annotate(r'$z_{\rm start}$ = ' + str(z[0]), xy=(0.8, 0.9), xycoords='axes fraction', fontsize=16)
	plt.savefig(model + '_z_' + str(z[0]) + '_clumpmass_vsclumplifetime.pdf', fontsize=16)
	plt.show()

def plot_grpclump_massdef_massloss(model, si, ei):
	#output time information
	model_timetrace_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
	trace_mass, trace_dist, trace_grpmatch, trace_output = grpclump_masstrace(model, si, ei) 
	dMsolUnit  = 1.5928853e16
	dMinGasMass= 2.66289e-12
	lower_time = []
	lower_mass = []
	full_time  = []
	full_mass  = []
	for i in range(len(trace_grpmatch)):
		#if(len(trace_grpmatch[i])>=3):
		time = []
		z = [] 
		for t in trace_output[i]:
			time.append(model_timetrace_info[model_timetrace_info['output'] == t]['time'].to_numpy()[0])
			z.append(model_timetrace_info[model_timetrace_info['output'] == t]['z'].to_numpy()[0])
		fig = plt.figure(figsize=(12,10))
		plt.plot(time, trace_mass[i], lw=2, c='k')
		plt.xlabel(r'cosmic time[Gyr]', fontsize=16)
		plt.ylabel(r'Cold clump mass', fontsize=16)
		plt.tick_params(labelsize=16)
		plt.hlines(2*dMinGasMass*dMsolUnit, plt.xlim()[0]-.25, plt.xlim()[1]+.25)
		plt.xlim([plt.xlim()[0]+.25, plt.xlim()[1]-0.25])
		plt.annotate(r'M$_{\rm clump}$(t=t$_{0}$) = ' + "{:.2e}".format(trace_mass[i][0]) + ' M$_{\odot}$', xy=(0.6, 0.9), xycoords='axes fraction', fontsize=16)
		plt.savefig(model + '_massloss_minmassdef_startoutput_' + str(trace_output[si]) + '_tracegrp' + str(trace_grpmatch[i][0]) + '.pdf')
		plt.semilogy()
		plt.show()

def plot_grpclump_massdef_intergratedtime(model, si, ei):
	#for timescale calculations 
	model_timetrace_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
	#for min mass mask 
	dMsolUnit  = 1.5928853e16
	dMinGasMass= 2.66289e-12
	outfn = model + '_allz_massdef_intergratedtimes_wdist_andTS.dat' 
	outfile = open(outfn, 'w')
	outfile.write('model tform[Gyr] lifetime[Gyr] NTS maxmass[Msol] r[kpc] type\n')
	lower_time = []
	lower_mass = []
	lowerform_t= []
	lower_dist = [] 
	full_time  = []
	full_mass  = []
	fullform_t = []
	full_dist  = [] 
	for j in range(ei-si-1): 	
		trace_mass, trace_dist, trace_grpmatch, trace_output = grpclump_masstrace(model, si+j, ei) 
		if(j!=0):
			trace_massjm1, trace_distjm1, trace_grpmatchjm1, trace_outputjm1 = grpclump_masstrace(model, si+j-1, ei)
			previous_grps = [] 
			for k in range(len(trace_grpmatchjm1)):
				#get the group associated with this time step which is the 2nd entry (index 1 of trace_grpmatchjm1[k] for the kth grp
				if(len(trace_grpmatchjm1[k]) > 1):
					previous_grps.append(trace_grpmatchjm1[k][1])
			tracegrps = [] 
			for k in range(len(trace_grpmatch)):
				#new trace grps for current timestep 
				tracegrps.append(trace_grpmatch[k][0]) 
			newgrps_map = np.logical_not(np.in1d(tracegrps, previous_grps))
			#redefine trace arrays to only look at newly formed grps 
			trace_mass = trace_mass[newgrps_map]
			trace_grpmatch = trace_grpmatch[newgrps_map]
			trace_output   = trace_output[newgrps_map]
			trace_dist     = trace_dist[newgrps_map]
		else: 
			print('j should be 0. j = ', j)
		for i in range(len(trace_grpmatch)):
			time = []
			z = [] 
			for t in trace_output[i]:
				time.append(model_timetrace_info[model_timetrace_info['output'] == t]['time'].to_numpy()[0])
				z.append(model_timetrace_info[model_timetrace_info['output'] == t]['z'].to_numpy()[0])
			time = np.array(time)
			z    = np.array(z)
			clumpmap = np.array(trace_mass[i]) > 2*dMsolUnit*dMinGasMass
			clumptimes = time[clumpmap]
			trace_mass_i = np.array(trace_mass[i])
			trace_dist_i = np.array(trace_dist[i])
			if(len(clumptimes) == 2):
				#lower bound-- cant see grow to peak or destruction
				lower_time.append(clumptimes[-1] - clumptimes[0])
				lower_mass.append(max(trace_mass_i[clumpmap]))	
				lower_dist.append(max(trace_dist_i[clumpmap]))
				lowerform_t.append(clumptimes[0])	
				outfile.write('%s %.2f %.2f %i %.2f %.2f %s\n' % (model, clumptimes[0], clumptimes[-1] - clumptimes[0], len(clumptimes), max(trace_mass_i[clumpmap]), max(trace_dist_i[clumpmap]), 'lower'))
			elif(len(clumptimes) > 2):
				full_time.append(clumptimes[-1] - clumptimes[0])
				full_mass.append(max(trace_mass_i[clumpmap])) 
				full_dist.append(max(trace_dist_i[clumpmap]))
				fullform_t.append(clumptimes[0])	
				outfile.write('%s %.2f %.2f %i %.2f %.2f %s\n' % (model, clumptimes[0], clumptimes[-1] - clumptimes[0],len(clumptimes), max(trace_mass_i[clumpmap]),max(trace_dist_i[clumpmap]) ,'full'))
	outfile.close()

def plot_cosmictime_vs_deltat(model):
	model_timetrace_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
	time = model_timetrace_info['time'].to_numpy()
	delta_t = []
	for i in range(len(time)-1):
    		delta_t.append(time[i] - time[i+1])	
	delta_t.append(time[-1])
	fig = plt.figure(figsize=(12,10))
	plt.scatter(time, delta_t, s=200, c='k')
	plt.xlabel(r'cosmic time[Gyr]', fontsize=16)
	plt.ylabel(r'$\Delta$ t [Gyr]', fontsize=16)
	plt.tick_params(labelsize=16)
	plt.savefig(model + '_time_vsdeltatime.pdf')

def plot_grpclump_massfrac_integratedtime(model, si, ei, plotmassloss):
	model_timetrace_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
	lower_time = []
	lower_mass = []
	full_time  = []
	full_mass  = []
	
	#go through all timesteps 
	for j in range(ei-si-1): 	
		trace_mass, trace_dist, trace_grpmatch, trace_output = grpclump_masstrace(model, si+j, ei) 
		if(j!=0):
			trace_massjm1, trace_distjm1, trace_grpmatchjm1, trace_outputjm1 = grpclump_masstrace(model, si+j-1, ei)
			previous_grps = [] 
			for k in range(len(trace_grpmatchjm1)):
				#get the group associated with this time step which is the 2nd entry (index 1 of trace_grpmatchjm1[k] for the kth grp
				if(len(trace_grpmatchjm1[k]) > 1):
					previous_grps.append(trace_grpmatchjm1[k][1])
			tracegrps = [] 
			for k in range(len(trace_grpmatch)):
				#new trace grps for current timestep 
				tracegrps.append(trace_grpmatch[k][0]) 
			newgrps_map = np.logical_not(np.in1d(tracegrps, previous_grps))
			#redefine trace arrays to only look at newly formed grps 
			trace_mass = trace_mass[newgrps_map]
			trace_grpmatch = trace_grpmatch[newgrps_map]
			trace_output   = trace_output[newgrps_map]
		else: 
			print('j should be 0. j = ', j)
		for i in range(len(trace_grpmatch)):
			time = []
			z = [] 
			for t in trace_output[i]:
				time.append(model_timetrace_info[model_timetrace_info['output'] == t]['time'].to_numpy()[0])
				z.append(model_timetrace_info[model_timetrace_info['output'] == t]['z'].to_numpy()[0])
			time = np.array(time)
			z    = np.array(z)
			trace_mass_i = np.array(trace_mass[i])
			max_mass  = max(trace_mass_i)
			norm_mass = np.array(trace_mass_i)/max_mass
			if(len(norm_mass)==2):
				#can only get a lower time bound with 2 timesteps! 
				ncross = 0  
				if(norm_mass[1] < 0.5 and norm_mass[0] > 0.5):
					slope = (norm_mass[1] - norm_mass[0]) / (time[1] - time[0])
					b = norm_mass[1]
					xoff = time[1]
					lower_time.append((0.5-b)/slope + xoff - time[0])
					lower_mass.append(max(trace_mass_i))
					ncross = ncross + 1
					print('intercept time = ', (0.5-b)/slope + xoff)
					print('integrated lower time = ', (0.5-b)/slope + xoff - time[0])
				elif(norm_mass[1] > 0.5 and norm_mass[0] > 0.5):
					lower_time.append(time[1]-time[0])
					lower_mass.append(max_mass)
					ncross = ncross + 1
					print('integrated lower time = ', time[1]-time[0]) 
				else:
					print('there is a combo you have not considered.')
					print('norm_mass[0] = ', norm_mass[0])
					print('norm_mass[1] = ', norm_mass[1])
			else: 
				print('tracing group number through more than 2 timesteps')
				ncross = 0  
				for j in range(len(norm_mass)-1):
					ti = 0 
					if(norm_mass[j] < 0.5 and norm_mass[j+1] > 0.5):
						xoff = time[j]
						b = norm_mass[j]
						slope = (norm_mass[j+1] - norm_mass[j]) / (time[j+1] - time[j])  
						print('b = ', b)
						print('slope = ', slope) 
						print('intercept time = ', (0.5-b)/slope + time[j])
						if(ncross == 0): 
							t1 = (0.5-b)/slope + xoff
						else:
							ti = (0.5-b)/slope + xoff
						ncross = ncross + 1 
					elif(norm_mass[j] > 0.5 and norm_mass[j+1] < 0.5): 
						xoff = time[j+1]
						b = norm_mass[j+1]
						slope = (norm_mass[j+1] - norm_mass[j]) / (time[j+1] - time[j])  
						print('b = ', b)
						print('slope = ', slope) 
						print('intercept time = ', (0.5-b)/slope + time[j+1])
						if(ncross == 0): 
							t1 = (0.5-b)/slope + xoff
						elif(j == len(norm_mass) - 1):
							tf = (0.5-b)/slope + xoff
						elif(j != len(norm_mass)-1 and ncross > 2):
							ti = (0.5-b)/slope + xoff
						ncross = ncross + 1 
					print('ncross = ', ncross)
				if(ncross > 1):
					tf = (0.5-b)/slope + xoff
				elif(ncross == 0): 
					tf = time[-1]
				elif(ncross == 1): 
					tf = t1
					t1 = time[0]
				#if(ncross > 1 and norm_mass[0] > min(norm_mass) and norm_mass[-1] > min(norm_mass)):
					#print('found a clump that is losing mass')
					#lower_time.append((time[-1]-ti) + (t1 - time[0]))
					#lower_mass.append(max_mass) 
				#else: 	
				full_time.append((tf-ti) + (ti-t1))
				full_mass.append(max_mass)
				print('integrated time = ', (tf-ti) + ti-t1)			
			if(plotmassloss):
				fig = plt.figure(figsize=(12,10))
				plt.plot(time, trace_mass_i/max(trace_mass_i), lw=2, c='k')
				plt.xlabel(r'cosmic time[Gyr]', fontsize=16)
				plt.ylabel(r'Cold clump mass / Max cold clump mass', fontsize=16)
				plt.tick_params(labelsize=16)
				plt.hlines(0.5, plt.xlim()[0]-.25, plt.xlim()[1]+.25)
				plt.xlim([plt.xlim()[0]+.25, plt.xlim()[1]-0.25])
				plt.annotate(r'M$_{\rm clump}$(t=t$_{0}$) = ' + "{:.2e}".format(trace_mass[i][0]) + ' M$_{\odot}$', xy=(0.6, 0.9), xycoords='axes fraction', fontsize=16)
				plt.savefig(model + '_massloss_startoutput_' + str(trace_output[i][0]) + '_tracegrp' + str(trace_grpmatch[i][0]) + '.pdf')
				#plt.show()
				plt.close()

		
	fig = plt.figure(figsize=(12,10))
	plt.scatter(full_mass, full_time, s=200, c='k', marker='o')
	plt.scatter(lower_mass, lower_time, s=200, c='indianred', marker='^')
	plt.semilogx()
	#plt.loglog()
	plt.tick_params(labelsize=16)
	plt.xlabel(r'log$_{10}$ M$_{\rm clump}$ [M$_{\odot}$]', fontsize=16)
	plt.ylabel(r't$_{\rm clump}$ [Gyr]', fontsize=16)
	plt.annotate(r'$z_{\rm start}$ = ' + str(z[0]), xy=(0.8, 0.9), xycoords='axes fraction', fontsize=16)
	plt.savefig(model + '_z_' + str(z[0]) + '_clumpmass_vsclumplifetime.pdf', fontsize=16)

def get_pyfof_com(h1_r, pyfof_iord_i):
	iord_map = np.in1d(h1_r['iord'], pyfof_iord_i)
	pyfof_grpi = h1_r[iord_map]
	pyfof_comi = pynbody.analysis.halo.center_of_mass(pyfof_grpi)
	return pyfof_comi

def redshift_rvir_grpclump_dat(model, si, ei, metals): 
	dfs = open_timetrace_grpclump_dat(model, si, ei, metals)
	z_array = []
	t_array = []
	rvir    = []
	Mvir    = [] 
	Nclumps = [] 
	mgas    = []
	mstar   = [] 
	halo_track_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
	halo_properties = pd.read_table('../datfiles/' + model + '_galaxy_properties.dat', sep='\s+')
	output1 = np.array(dfs['output1'])
	print('dfs[output1] = ', output1)
	for ts in output1:
		print('ts = ', ts)
		print('z = ', halo_track_info[halo_track_info['output'] == ts]['z'])
		t_array.append(float(halo_track_info[halo_track_info['output'] == ts]['time']))
		z_array.append(float(halo_track_info[halo_track_info['output'] == ts]['z']))
		rvir.append(float(halo_properties[halo_properties['output'] == ts]['rvir[kpc]']))
		Mvir.append(float(halo_properties[halo_properties['output'] == ts]['dmass[Msol]']))
		mgas.append(float(halo_properties[halo_properties['output'] == ts]['gmass[Msol]']))
		mstar.append(float(halo_properties[halo_properties['output'] == ts]['smass[Msol]']))
		Nclumps.append(len(dfs[dfs['output1'] == ts]))
	G.units = 'cm**3 g**-1 s**-2'
	dfs['Vvir[kms^-1]'] = (G.in_units('kpc km**2 Msol**-1 s**-2')*Mvir/rvir)**(0.5)
	dfs['z'] = z_array
	dfs['t[Gyr]'] = t_array
	dfs['rvir[kpc]'] = rvir 
	dfs['Mvir[Msol]'] = Mvir 
	dfs['Mstar[Msol]'] = mstar
	dfs['Mgas[Msol]']  = mgas
	dfs['rbyrvir'] = dfs['r[kpc]'] / dfs['rvir[kpc]']
	if metals: 
		dfs['vrbyvvir'] = dfs['vr[kms^-1]'] / dfs['Vvir[kms^-1]'] 
	dfs['Nclumps'] = Nclumps
	return dfs 

def COMdist_hist(model,  si, ei, cmap, show):
	if model == 'P0': 
		cmap = 'rocket_r'
	dfs = redshift_rvir_grpclump_dat(model, si, ei, False) 
	fig = plt.figure(figsize=(12,10))
	sns.histplot(data=dfs, x='rbyrvir', hue='z', multiple='stack', log_scale=True, palette=cmap)
	plt.xlabel(r'r$_{\rm clump}$/r$_{\rm vir}$', fontsize=16)
	plt.ylabel(r'log$_{10}$ Count', fontsize=16)
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.annotate(model, xy=(0.1, 0.9),xycoords='axes fraction', fontsize=16)
	plt.savefig(model + '_COMdist_rvirnorm_hist.pdf')
	if show: 
		plt.show() 
	else: 
		plt.close()

def COMdist_kde(model, si, ei, cmap, show): 
	if model == 'P0': 
		cmap = 'rocket_r'
	dfs = redshift_rvir_grpclump_dat(model, si, ei) 
	fig = plt.figure(figsize=(12,10))
	sns.kdeplot(data=dfs, x='rbyrvir', hue='z', log_scale=True, palette=cmap, lw=3)
	plt.xlabel(r'r$_{\rm clump}$/r$_{\rm vir}$', fontsize=16)
	plt.ylabel('Density', fontsize=16)
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.xlim(1e-2, 2e0)
	plt.ylim(0, 0.65)
	plt.annotate(model, xy=(0.1, 0.9),xycoords='axes fraction', fontsize=16)
	plt.savefig(model + '_COMdist_rvirnorm_kde.pdf')
	if show: 
		plt.show() 
	else: 
		plt.close()

def plot_vdisp_r(model, output, z, show): 
	sizelim = 0.25
	infile = '../datfiles/' + model + '_'+ output + '_HI4e-08_1kpc_clump_data.dat'
	df = pd.read_table(infile, sep='\s+')
	df = df[df['max_r[kpc]'] > sizelim]	
	fig = plt.figure(figsize=(12, 10))
	plt.scatter(df['r[kpc]'], df['vdist[km/s]'], s=200, color=pyfof_HI_clumps.model_colors(model), edgecolor='k')
	plt.title(model + r' $z = $' + str(z), fontsize=16)
	plt.loglog()
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.xlabel(r'log$_{10}$ r$_{\rm clump}$ [kpc]', fontsize=16)
	plt.ylabel(r'log$_{10} \sigma_v$ [km/s]', fontsize=16)
	plt.savefig(model + '_' + str(output[0]) + '_vdisp_clump_profile.pdf')
	if show:
		plt.show()
	else:
		plt.close()

def plot_size_r(model, output, z, show): 
	sizelim = 0.25
	infile = '../datfiles/' + model + '_'+ output + '_HI4e-08_1kpc_clump_data.dat'
	df = pd.read_table(infile, sep='\s+')
	df = df[df['max_r[kpc]'] > sizelim]	
	fig = plt.figure(figsize=(12, 10))
	plt.scatter(df['r[kpc]'], df['max_r[kpc]'], s=200, color=pyfof_HI_clumps.model_colors(model), edgecolor='k')
	plt.title(model + r' $z = $' + str(z), fontsize=16)
	plt.loglog()
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.xlabel(r'log$_{10}$ r$_{\rm clump}$ [kpc]', fontsize=16)
	plt.ylabel(r'log$_{10}$ size [kpc]', fontsize=16)
	plt.savefig(model + '_' + str(output[0]) + '_size_clump_profile.pdf')
	if show:
		plt.show()
	else:
		plt.close()

def pressure_confinement_plot(model, si, adiabatic, show):
	halo, output = pyfof_HI_clumps.get_model_outputs_halos(model)
	halo = halo[si:si+2]
	output = output[si:si+2] 
	fn2 = pyfof_HI_clumps.get_fn(model)
	infile = fn2 + 'timetrace_grps_lowtohigh_outputs' + str(output[0]) + '_' + str(output[1]) + '_clump_properties.dat'
	df = pd.read_table(infile, sep='\s+')
	
	halo_track_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
	halo_properties = pd.read_table('../datfiles/' + model + '_galaxy_properties.dat', sep='\s+')
	halo_df = halo_properties[halo_properties['output'] == output[0]]
	z = float(halo_track_info[halo_track_info['output'] == output[0]]['z'])
	rvir  = float(halo_df['rvir[kpc]'])
	Mstar = float(halo_df['smass[Msol]'])
	Mvir  = float(halo_df['dmass[Msol]']) 	
	if(len(str(output[0])) == 4):	
		h1 = pyfof_HI_clumps.get_mainHalo(model, '00' + str(output[0]), halo[0])
	else:
		h1 = pyfof_HI_clumps.get_mainHalo(model, '000' + str(output[0]), halo[0])
	try: 
		h1 = h1[h1['amiga.grp']==halo[0]]
		print('got rid of substructure')
		outfile = model + '_z' + str(z) + '_pressure_profile_coldclump_pressure_profile_nosub_wquant' + '.pdf'
	except: 
		print('amiga.grp files dont exist for model ' + model)
		outfile = model + '_z' + str(z) + '_pressure_profile_coldclump_pressure_profile_wquant' + '.pdf'
	CGM = pyfof_HI_clumps.cut_gal(h1)
	cold_CGM = CGM[CGM['temp'].in_units('K') < 10**5]
	hot_CGM = CGM[CGM['temp'].in_units('K') >= 10**5]
	cold_CGMprofile = pynbody.analysis.profile.Profile(cold_CGM, min='10 kpc', max='250 kpc')
	hot_CGMprofile = pynbody.analysis.profile.Profile(hot_CGM, min='10 kpc', max='250 kpc')
	dfpressure = pynbody.array.SimArray(10**df['log10avg_pressure[Msols**-2kpc**-1]'])
	dfpressure.units = 'Msol s**-2 kpc**-1'

	#make profile error bars 
	cold_CGMquant = pynbody.analysis.profile.QuantileProfile(cold_CGM, rmin='10 kpc', rmax='250 kpc') 
	hot_CGMquant  = pynbody.analysis.profile.QuantileProfile(hot_CGM, rmin='10 kpc', rmax='250 kpc') 
	
	#print('hot_CGMquant["p"] = ', hot_CGMquant['p'])
	
	fig = plt.figure(figsize=(12, 10))
	plt.scatter(df['r[kpc]'], dfpressure.in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), s=200, color=pyfof_HI_clumps.model_colors(model), edgecolor='k', label='fof cold clumps')
	plt.plot(cold_CGMquant['rbins'], cold_CGMquant['p'][:, 1].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), 'k', ls='--', lw=3, label=r'cool CGM (T < 10$^{5}$ K')
	plt.plot(hot_CGMquant['rbins'], hot_CGMquant['p'][:, 1].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), 'k', ls='-', lw=3, label=r'hot CGM (T $\geq$ 10$^5$ K)')
	#plt.plot(cold_CGMquant['rbins'], cold_CGMquant['p'][:, 1].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), 'k', ls='--', lw=3, label=r'cool CGM (T < 10$^{5}$ K')
	#plt.plot(hot_CGMquant['rbins'], hot_CGMquant['p'][:, 1].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), 'k', ls='-', lw=3, label=r'hot CGM (T $\geq$ 10$^5$ K)')
	#plt.fill_between(cold_CGMquant['rbins'], cold_CGMquant['p'][:,0].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), cold_CGMquant['p'][:,2].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), color='k', alpha=0.5)
	#plt.fill_between(hot_CGMquant['rbins'], hot_CGMquant['p'][:,0].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), hot_CGMquant['p'][:,2].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), color='k', alpha=0.5)
	plt.fill_between(cold_CGMquant['rbins'], cold_CGMquant['p'][:,0].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), cold_CGMquant['p'][:,2].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), color='k', alpha=0.5)
	plt.fill_between(hot_CGMquant['rbins'], hot_CGMquant['p'][:,0].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), hot_CGMquant['p'][:,2].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), color='k', alpha=0.5)
	plt.axhspan(1, 100, alpha=0.5, color='gray', label='Stocke et al. 2013')
	if adiabatic: 
		T = 1e5
		plt.plot(cold_CGMprofile['rbins'], BMP.adiabatic_pressure_profile(cold_CGMprofile['rbins'], BMP.Mvir, BMP.z)/kb.in_units('g cm**2 s**-2 K**-1'), 'mediumblue', lw=3, label=r'adiabatic NFW Maller & Bullock 2004')
		#plt.plot(cold_CGMprofile['rbins'], (BMP.Thalo/BMP.Tcool)*BMP.adiabatic_pressure_profile(cold_CGMprofile['rbins'], BMP.Mvir, BMP.z)/kb.in_units('g cm**2 s**-2 K**-1'), 'b', lw=3, label=r'cold adiabatic NFW Maller & Bullock 2004')
	plt.loglog()
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.xlabel(r'log$_{10}$ r$_{\rm clump}$ [kpc]', fontsize=16)
	plt.ylabel(r'P/k$_{\rm b}$ [cm$^{-3}$ K]', fontsize=16)
	plt.ylim(1e-1, 2e5)
	plt.xlim(10, rvir + 10)
	plt.legend(fontsize=14)
	plt.title(model + r' $z = $' + str(z), fontsize=16)
	plt.savefig(outfile)
	if show:
		plt.show()
	else:
		plt.close()


def allGMs_pressure_confinement_plot(models, sinds, adiabatic, show):
	fig = plt.figure(figsize=(12, 10))
	for i in range(len(models)):
		halo, output = pyfof_HI_clumps.get_model_outputs_halos(models[i])
		halo = halo[sinds[i]:sinds[i]+2]
		output = output[sinds[i]:sinds[i]+2] 
		fn2 = pyfof_HI_clumps.get_fn(models[i])
		infile = fn2 + 'timetrace_grps_lowtohigh_outputs' + str(output[0]) + '_' + str(output[1]) + '_clump_properties.dat'
		df = pd.read_table(infile, sep='\s+')
	
		halo_track_info = pd.read_table('../datfiles/' + models[i] + '_halo_track_info.dat', sep='\s+')
		halo_properties = pd.read_table('../datfiles/' + models[i] + '_galaxy_properties.dat', sep='\s+')
		halo_df = halo_properties[halo_properties['output'] == output[0]]
		z = float(halo_track_info[halo_track_info['output'] == output[0]]['z'])
		rvir  = float(halo_df['rvir[kpc]'])
		Mstar = float(halo_df['smass[Msol]'])
		Mvir  = float(halo_df['dmass[Msol]']) 	
		print('rvir = ', rvir)
	
		h1 = pyfof_HI_clumps.get_mainHalo(models[i], '00' + str(output[0]), halo[0])
		try: 
			h1 = h1[h1['amiga.grp']==halo[0]]
			print('got rid of substructure')
			outfile = 'allGMs_z' + str(z) + '_pressure_profile_coldclump_pressure_profile_nosub' + '.pdf'
		except: 
			print('amiga.grp files dont exist for model ' + models[i])
			outfile = 'allGMs_z' + str(z) + '_pressure_profile_coldclump_pressure_profile' + '.pdf'
		CGM = pyfof_HI_clumps.cut_gal(h1)
		cold_CGM = CGM[CGM['temp'].in_units('K') < 10**5]
		hot_CGM = CGM[CGM['temp'].in_units('K') > 10**6]
		cold_CGMprofile = pynbody.analysis.profile.Profile(cold_CGM, min='10 kpc', max=str(rvir) + ' kpc')
		hot_CGMprofile = pynbody.analysis.profile.Profile(hot_CGM, min='10 kpc', max=str(rvir) + ' kpc')
		dfpressure = pynbody.array.SimArray(10**df['log10avg_pressure[Msols**-2kpc**-1]'])
		dfpressure.units = 'Msol s**-2 kpc**-1'

		plt.scatter(df['r[kpc]'], dfpressure.in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), s=200, color=pyfof_HI_clumps.model_colors(models[i]), edgecolor='k', label=models[i] + ' fof cold clumps')
		plt.plot(cold_CGMprofile['rbins'], cold_CGMprofile['p'].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), color=pyfof_HI_clumps.model_colors(models[i]), ls='--', lw=3, label=models[i] + r'cold CGM (T < 10$^{5}$ K')
		plt.plot(hot_CGMprofile['rbins'], hot_CGMprofile['p'].in_units('g cm**-1 s**-2')/kb.in_units('g cm**2 s**-2 K**-1'), color=pyfof_HI_clumps.model_colors(models[i]), ls='-', lw=3, label=models[i] + r'hot CGM (T > 10$^6$ K)')
		plt.axhspan(1, 100, alpha=0.5, color='gray', label='Stocke et al. 2013')
		if adiabatic: 
			T = 1e5
			plt.plot(cold_CGMprofile['rbins'], BMP.adiabatic_pressure_profile(cold_CGMprofile['rbins'], BMP.Mvir, BMP.z)/kb.in_units('g cm**2 s**-2 K**-1'), 'maroon', lw=3, label=r'hot adiabatic NFW Maller & Bullock 2004')
			plt.plot(cold_CGMprofile['rbins'], (BMP.Thalo/BMP.Tcool)*BMP.adiabatic_pressure_profile(cold_CGMprofile['rbins'], BMP.Mvir, BMP.z)/kb.in_units('g cm**2 s**-2 K**-1'), 'b', lw=3, label=r'cold adiabatic NFW Maller & Bullock 2004')
	plt.loglog()
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.xlabel(r'log$_{10}$ r$_{\rm clump}$ [kpc]', fontsize=16)
	plt.ylabel(r'P/k$_{\rm b}$ [cm$^{-3}$ K]', fontsize=16)
	plt.legend(fontsize=14)
	plt.ylim(1e-1, 2e5)
	plt.title(r' $z = $' + str(z), fontsize=16)
	plt.savefig(outfile)
	if show:
		plt.show()
	else:
		plt.close()


def com_pos_ridgeplot(model, si, ei, sample,  sample_i):
	if sample: 
		i = si 
		j = 1
		while i < ei: 
			df = redshift_rvir_grpclump_dat(model, i, i+2, False)
			if i==si: 
				dfs = df 
			else:
				dfs = pd.concat([df, dfs])
			i = sample_i*j
			j = j + 1
	else:
		dfs = redshift_rvir_grpclump_dat(model, si, ei, False)
	sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
	z = np.unique(dfs['z'])
	N = np.unique(dfs['Nclumps'])
	print('number of redshift outputs ', len(z))
	print('number of N ', len(N))
	pal = sns.color_palette(palette='Reds', n_colors=len(z))
	g = sns.FacetGrid(dfs, row='z', hue='z', aspect=15, height=0.7, palette=pal)
	g.map(sns.kdeplot, 'rbyrvir', common_norm=False,
      	bw_adjust=1, clip_on=False,
      	fill=True, alpha=1, linewidth=1.5)
	g.map(sns.kdeplot, 'rbyrvir', common_norm=False,
      	bw_adjust=1, clip_on=False,
      	color="w", lw=2)
	g.map(plt.axhline, y=0,
      	lw=2, clip_on=False)
	for i, ax in enumerate(g.axes.flat):
		ax.text(-0.175, 0.02, r'$z$ = ' + str(z[i]),
		fontweight='bold', fontsize=15,
		color=ax.lines[-1].get_color())
    		
		ax.text(1.2, 0.02, r'N = ' + str(N[i]),
		fontweight='bold', fontsize=15,
		color=ax.lines[-1].get_color())
	g.fig.subplots_adjust(hspace=-0.3)
	g.set_titles("")
	g.set(yticks=[])
	g.despine(bottom=True, left=True)

	plt.setp(ax.get_xticklabels(), fontsize=15, fontweight='bold')
	plt.xlabel(r'r$_{\rm clump}$/r$_{\rm vir}$', fontweight='bold', fontsize=15)
	plt.gcf().subplots_adjust(bottom=0.15)
	g.fig.suptitle(model + ' clump radial position',
               ha='center',
               fontsize=20,
               fontweight=20)
	plt.savefig(model + '_ridgeplot_com_dist.pdf')
	plt.show()

def Jeans_mass(P):
	cw = pynbody.array.SimArray(11.5)
	cw.units = 'km s**-1'
	G = pynbody.array.SimArray(6.6743*10**(-11))
	G.units = 'm**3 kg**-1 s**-2'
	return (9*cw.in_units('kpc s**-1')**4)/(5*(G.in_units('kpc**3 Msol**-1 s**-2'))**(3/2)*P**(1/2))

def Jeans_mass_vs_mass(model, si, ei, show): 
	dfs = redshift_rvir_grpclump_dat(model, si, ei, False)
	size_lim = 0.25
	dfs_sizelim = dfs[dfs['size[kpc]'] > size_lim]
	z = np.unique(dfs_sizelim['z'])
	dfs_lowz = dfs_sizelim[dfs_sizelim['z'] == min(z)]
	dfpressure = pynbody.array.SimArray(10**dfs_lowz['log10avg_pressure[Msols**-2kpc**-1]'])
	dfpressure.units = 'Msol s**-2 kpc**-1'
	jeans_mass = Jeans_mass(dfpressure)
	massmin = 0.7*min(dfs_lowz['mass[Msol]'])
	massmax = 1.3*max(dfs_lowz['mass[Msol]'])
	masslims = np.linspace(massmin, massmax, 100)
	ones = np.ones(100)	
	
	#plot 
	fig = plt.figure(figsize=(12,10))
	plt.scatter(dfs_lowz['mass[Msol]'], jeans_mass/dfs_lowz['mass[Msol]'], s=100, c=dfs_lowz['rbyrvir'], cmap='magma', edgecolor='k')
	plt.xlabel(r'M$_{\rm clump}$ [M$_{\odot}$]', fontsize=16)
	plt.ylabel(r'M$_{\rm J}$/M$_{\rm clump}$', fontsize=16)
	cb = plt.colorbar()
	cb.set_label(r'r$_{\rm clump}$/r$_{\rm vir}$', fontsize=16)
	plt.loglog()
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.plot(masslims, ones, lw=2, c='k')
	plt.annotate(model + r'; $z = $' + str(min(z)) , xy=(0.8, 0.9),xycoords='axes fraction', fontsize=16)
	cb.ax.tick_params(labelsize=16)
	plt.xlim(massmin, massmax)
	plt.savefig(model + '_jeansmassbymassvsmass_rclumpcolor.pdf')
	if show:
		plt.show()
	else:
		plt.close()

def make_comdist_moments(models, si_array, ei_array):
	size_lim = 0.25 #kpc 
	for i in range(len(models)):
		dfs = redshift_rvir_grpclump_dat(models[i], si_array[i], ei_array[i], True) 
		z = np.unique(dfs['z'])
		dfs_sizelim = dfs[dfs['size[kpc]'] > size_lim]
		rbyrvir_mean = [] 
		rbyrvir_std  = [] 
		rbyrvir_skew = []
		for j in range(len(z)): 
			dfs_sizelim_z = dfs_sizelim[dfs_sizelim['z'] == z[j]]
			rbyrvir_z = dfs_sizelim_z['rbyrvir'].to_numpy()
			rbyrvir_mean.append(rbyrvir_z.mean())
			rbyrvir_std.append(rbyrvir_z.std())
			rbyrvir_skew.append(skew(rbyrvir_z))
		np.savetxt(models[i] + '_rbyrvir_skew_metals.txt', rbyrvir_skew)		
		np.savetxt(models[i] + '_rbyrvir_std_metals.txt',  rbyrvir_std)		
		np.savetxt(models[i] + '_rbyrvir_mean_metals.txt', rbyrvir_mean)		
		np.savetxt(models[i] + '_rbyrvir_redshift_metals.txt', z)		
		print('done w ' + models[i])

def make_radcomvel_moments(models, si_array, ei_array):
	size_lim = 0.25 #kpc 
	for i in range(len(models)):
		df  = redshift_rvir_grpclump_dat(models[i], si_array[i], ei_array[i], True)
		z = np.unique(df['z'])
		df_sizelim = df[df['size[kpc]'] > size_lim]
		vr_mean = [] 
		vr_std  = [] 
		vr_skew = []
		for j in range(len(z)): 
			df_sizelim_z = df_sizelim[df_sizelim['z'] == z[j]]
			vr_z = df_sizelim_z['vr[kms^-1]'].to_numpy()
			vr_mean.append(vr_z.mean())
			vr_std.append(vr_z.std())
			vr_skew.append(skew(vr_z))
		np.savetxt(models[i] + '_fixvr_skew.txt', vr_skew)		
		np.savetxt(models[i] + '_fixvr_std.txt',  vr_std)		
		np.savetxt(models[i] + '_fixvr_mean.txt', vr_mean)		
		np.savetxt(models[i] + '_fixvr_redshift.txt', z)		

def mass_prof(halo):
	p = pynbody.analysis.profile.Profile(halo, vmin =.01, ndim=3)	
	mass = p['mass_enc']
	rbin = p['rbins']
	return rbin, mass 

def arrayval_ind(array, val):
	return np.abs(array - val).argmin()

def lin_interp(xarr, yarr, x):
	ind = arrayval_ind(xarr, x)
	if ind == len(xarr) - 1:
		y = yarr[ind]
	elif xarr[ind] > x and ind != 0:
		y =  yarr[ind - 1] + (x - xarr[ind - 1]) * ((yarr[ind] - yarr[ind - 1])/(xarr[ind] - xarr[ind - 1]))
	elif xarr[ind] > x and ind == 0:
		y = yarr[0]
	elif xarr[ind] < x:
		y = yarr[ind] + (x - xarr[ind]) * ((yarr[ind + 1] - yarr[ind])/(xarr[ind + 1] -  xarr[ind]))
	return y
	
def make_migration_free_fall_time(model, si, ei): 
	size_lim = 0.25 #kpc
	G = pynbody.array.SimArray(6.6743 * 1e-11) 
	G.units = 'm**3 kg**-1 s**-2'
	outfn = model + '_migration_fullden_freefall_times_z.dat' 
	outfile = open(outfn, 'w')
	outfile.write('model t[Gyr] z migration_time[Gyr] freefall[Gyr] meanmass[Msol] meansize[kpc] avgr[kpc] avgvr[kms^-1]\n')
	#this free fall time is wrong 
	for i in range(ei-si): 
		marr = []
		try:
			halo_track_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
			df  = redshift_rvir_grpclump_dat(model, si+i, si+i+2, True)
			df_sizelim = df[df['size[kpc]'] > size_lim]
			output = df_sizelim['output1'].to_numpy()[0]
			print('halo output = ', output)
			h = halo_track_info[halo_track_info['output'] == output]['halo'].to_numpy()[0]
			print('halo number = ', h)
			if(len(str(output)) == 3):  
				halo = pyfof_HI_clumps.get_mainHalo(model, '000' + str(output), h)
			else: 
				halo = pyfof_HI_clumps.get_mainHalo(model, '00' + str(output), h)
			rbins, pmass = mass_prof(halo)
			t = df_sizelim['t[Gyr]'][0]
			z = df_sizelim['z'][0]
			vr   = pynbody.array.SimArray(df_sizelim['vr[kms^-1]'].to_numpy())
			vr.units = 'km s**-1'
			r    = pynbody.array.SimArray(df_sizelim['r[kpc]'].to_numpy())
			r.units = 'kpc'
			size = pynbody.array.SimArray(df_sizelim['size[kpc]'].to_numpy())
			size.units = 'kpc'
			mass = pynbody.array.SimArray(df_sizelim['mass[Msol]'].to_numpy())
			mass.units = 'Msol'
			Mvir = pynbody.array.SimArray(df_sizelim['Mvir[Msol]'].to_numpy()[0])
			Mvir.units = 'Msol'
			Mgas = pynbody.array.SimArray(df_sizelim['Mgas[Msol]'].to_numpy()[0])
			Mgas.units = 'Msol'
			Mstar = pynbody.array.SimArray(df_sizelim['Mstar[Msol]'].to_numpy()[0])
			Mstar.units = 'Msol'
			Rvir = pynbody.array.SimArray(df_sizelim['rvir[kpc]'].to_numpy()[0])
			Rvir.units = 'kpc'
			tmig = r/vr 
			tmig_avg = tmig.mean()
			#tmig_avg = r.mean() / vr.mean()
			tmig_gyr = tmig_avg.in_units('Gyr')
			#mean_den =  (Mvir + Mstar + Mgas) / (4*np.pi*(Rvir)**3/3)
			#tdyn2 = 0.767/(G*mean_den) #B+T eqn 2.39
			for j in range(len(r)):
				marr.append(lin_interp(rbins, pmass, r[j]))
			Mgal = pynbody.array.SimArray(marr)
			Mgal.units = 'Msol'
			print('Mgal at r = ', Mgal)
			print('total M = ', Mvir + Mstar + Mgas) 
			tdyn = (np.pi*r**(3/2))/(G.in_units('kpc**3 Msol**-1 Gyr**-2')*Mgal)**(0.5)
			print('tdyn = ', tdyn)
			#tdyn = (tdyn2.in_units('Gyr**2'))**(0.5)
			outfile.write('%s %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n' %(model, t, z, np.abs(tmig_gyr), tdyn.mean(), mass.mean(), size.mean(), r.mean(), vr.mean()))
		except: 
			print('something up with index model ' + model + ' index ' + str(si+i))	
	outfile.close()

def make_clump_count_datfile(model, si, ei):
	size_lim = 0.25 #kpc 
	outfn = '../datfiles/' + model + '_Nclumps_time_z_metals_amiga.dat' 
	outfile = open(outfn, 'w')
	outfile.write('model t[Gyr] z Nclump\n')
	df  = redshift_rvir_grpclump_dat(model, si, ei, True)
	df_sizelim = df[df['size[kpc]'] > size_lim]
	z = np.unique(df_sizelim['z'])
	for i in z: 
		Nclumps = len(df_sizelim[df_sizelim['z'] == i])
		print('z = ', z) 
		print('Nclumps = ', Nclumps)
		t = np.array(df_sizelim[df_sizelim['z'] == i]['t[Gyr]'])[0]
		outfile.write('%s %.2f %.2f %i\n' %(model, t, i, Nclumps))
	outfile.close()

def rbyrvir_moment_plot(models, moment):
	fig = plt.figure(figsize=(12,10))
	for i in range(len(models)):
		moment_array = np.loadtxt(models[i] + '_rbyrvir_' + moment + '_metals.txt')
		z = np.loadtxt(models[i] + '_rbyrvir_redshift_metals.txt')
		#for arrays made with clump_prop.dat files which used AHF...
		#moment_array = np.loadtxt(models[i] + '_rbyrvir_' + moment + '.txt')
		#z = np.loadtxt(models[i] + '_rbyrvir_redshift.txt')
		
		plt.scatter(z, moment_array, s=200, edgecolor='k', facecolor=pyfof_HI_clumps.model_colors(models[i]), label=models[i])
	plt.xlabel(r'$z$', fontsize=16)
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.legend(fontsize=14)
	plt.axvline(x=1, color='k', lw=2)
	if(moment == 'std'):
		plt.ylabel(r'$\sigma$ r/r$_{\rm vir}$', fontsize=16)
	elif(moment == 'mean'):
		plt.ylabel(r' $\bar{\rm{r}/\rm{r}_{\rm vir}}$', fontsize=16)
	elif(moment == 'skew'):
		plt.ylabel(r'skew r/r$_{\rm vir}$', fontsize=16)
	#using metals dat files 
	plt.savefig('amiga_GMs_comdist_moment_' + moment + '.pdf')
	plt.show()

def rbyrvir_moment_subplots(models):
	f, ax = plt.subplots(figsize=(12,10), nrows=2, ncols=1, sharex=True)
	moments = ['mean', 'std']
	moment_ylabel = [r' $\overline{\rm{r}/\rm{r}_{\rm vir}}$', r'$\sigma$ r/r$_{\rm vir}$']
	#f, ax = plt.subplots(figsize=(10,12), nrows=3, ncols=1, sharex=True)
	#moments = ['mean', 'std', 'skew']
	#moment_ylabel = [r' $\overline{\rm{r}/\rm{r}_{\rm vir}}$', r'$\sigma$ r/r$_{\rm vir}$', r'skew r/r$_{\rm vir}$']
	for i in range(len(models)):
		for m in range(len(moments)): 	
			moment_array = np.loadtxt('../txtfiles/' + models[i] + '_rbyrvir_' + moments[m] + '_metals.txt')
			z = np.loadtxt('../txtfiles/' + models[i] + '_rbyrvir_redshift_metals.txt')
			#using clump dat files which used AHF 
			#moment_array = np.loadtxt('../txtfiles/' + models[i] + '_rbyrvir_' + moments[m] + '.txt')
			#z = np.loadtxt('../txtfiles/' + models[i] + '_rbyrvir_redshift.txt')
		
			ax[m].scatter(z, moment_array, s=200, edgecolor='k', facecolor=pyfof_HI_clumps.model_colors(models[i]), label=models[i])
			ax[m].axvline(x=1, color='k', lw=2)
			ax[m].set_ylabel(moment_ylabel[m], fontsize=24)
			ax[m].tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=24, which='both')
			#ax[m].semilogx()
	plt.xlabel(r'$z$', fontsize=24)
	ax[1].legend(fontsize=22, ncol=2)
	plt.tight_layout()	
	
	plt.savefig('../pdfs/Hubble2022_amiga_GMs_comdist_suplot_moments_tight.pdf')
	plt.show()

def rbyrvir_moment_BHvsnoBHs_subplots(compare_models, compare_model_color):
	f, ax = plt.subplots(figsize=(10,30), nrows=3, ncols=1, sharex=True)
	moments = ['mean', 'std', 'skew']
	moment_ylabel = [r' $\overline{\rm{r}/\rm{r}_{\rm vir}}$', r'$\sigma$ r/r$_{\rm vir}$', r'skew r/r$_{\rm vir}$']
	for i in range(len(compare_models)):
		for m in range(len(moments)): 	
			compare_model = compare_models[i]
			moment_array = np.loadtxt(compare_model[0] + '_rbyrvir_' + moments[m] + '.txt')
			z = np.loadtxt(compare_model[0] + '_rbyrvir_redshift.txt')
			moment_array2 = np.loadtxt(compare_model[1] + '_rbyrvir_' + moments[m] + '.txt')
			z2 = np.loadtxt(compare_model[1] + '_rbyrvir_redshift.txt')
			model1_map = np.in1d(z, z2) 
			model2_map = np.in1d(z2, z) 
			moment_array = moment_array[model1_map]
			z = z[model1_map]
			moment_array2 = moment_array2[model2_map]
			z2 = z2[model2_map]
			moment_ratio = moment_array/moment_array2
			ax[m].scatter(z, moment_ratio, s=200, edgecolor='k', facecolor=compare_model_color[i], label=compare_model[0] + '/' + compare_model[1])
			ax[m].axvline(x=1, color='k', lw=2)
			ax[m].axhline(y=1, color='k', lw=2)
			ax[m].set_ylabel(moment_ylabel[m], fontsize=16)
			ax[m].tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
			ax[m].semilogx()
	plt.xlabel(r'log$_{10} z$', fontsize=16)
	ax[0].legend(fontsize=12)
	
	plt.savefig('logz_GMs_BHvsnoBHs_comdist_suplot_moments.pdf')
	plt.show()

def vr_moment_subplots(models):
	f, ax = plt.subplots(figsize=(10,12), nrows=3, ncols=1, sharex=True)
	moments = ['mean', 'std', 'skew']
	moment_ylabel = [r' $\overline{\rm{v}_{r}}$ [km/s]', r'$\sigma \rm{v}_r$ [km/s]', r'skew v$_r$']
	for i in range(len(models)):
		for m in range(len(moments)): 	
			moment_array = np.loadtxt('../txtfiles/' + models[i] + '_fixvr_' + moments[m] + '.txt')
			z = np.loadtxt('../txtfiles/' + models[i] + '_fixvr_redshift.txt')
		
			ax[m].scatter(z, moment_array, s=200, edgecolor='k', facecolor=pyfof_HI_clumps.model_colors(models[i]), label=models[i])
			ax[m].axvline(x=1, color='k', lw=2)
			ax[m].set_ylabel(moment_ylabel[m], fontsize=16)
			ax[m].tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
			if(moments[m] == 'mean'):
				ax[m].axhline(y=0, color='k', lw=2)
				
	plt.xlabel(r'$z$', fontsize=16)
	ax[0].legend(fontsize=12)
	plt.tight_layout()	
	plt.savefig(models[0] + '_fixvr_suplot_moments.pdf')
	plt.show()

def plot_avg_timescales(models, fullden):	
	fig = plt.figure(figsize=(12, 10))
	for i in range(len(models)):
		if fullden:
			df = pd.read_table(models[i] + '_migration_fullden_freefall_times_z.dat', sep='\s+')
		else:
			df = pd.read_table(models[i] + '_migration_freefall_times_z.dat', sep='\s+')
		plt.scatter(df['z'], df['freefall[Gyr]'], s=200, color=pyfof_HI_clumps.model_colors(models[i]), marker='s', label = models[i] + r' t$_{\rm dyn}$')
		plt.scatter(df['z'], df['migration_time[Gyr]'], s=200, color=pyfof_HI_clumps.model_colors(models[i]), label=models[i] + r' t$_{\rm migr}$')
	plt.semilogy()
	plt.axvline(x=1, color='k', lw=2)
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.xlabel(r'$z$', fontsize=16)
	plt.ylabel(r'log$_{10}$ timescale [Gyr]', fontsize=16)
	plt.legend(fontsize=14, ncol=2)
	#plt.ylim(1e-1, 2e5)
	if fullden:
		plt.savefig(models[0] + '_' + models[-1] + '_plot_fullden_freefall_vs_migrtime.pdf')
	else: 
		plt.savefig(models[0] + '_' + models[-1] + '_plot_freefall_vs_migrtime.pdf')
	plt.show()

def plot_avg_timescales_flowcolors(models, fullden):	
	fig = plt.figure(figsize=(12, 10))
	for i in range(len(models)):
		if fullden:
			df = pd.read_table(models[i] + '_migration_fullden_freefall_times_z.dat', sep='\s+')
		else:
			df = pd.read_table(models[i] + '_migration_freefall_times_z.dat', sep='\s+')
		outflow_df = df[df['avgvr[kms^-1]'] > 0]
		inflow_df = df[df['avgvr[kms^-1]'] < 0]
		plt.scatter(df['z'], df['freefall[Gyr]'], s=200, color='k', marker='s')
		plt.scatter(outflow_df['z'], outflow_df['migration_time[Gyr]'], s=200, color='r')
		plt.scatter(inflow_df['z'], inflow_df['migration_time[Gyr]'], s=200, color='b')
		if(i == 0):
			plt.scatter(df['z'], df['freefall[Gyr]'], s=200, color='k', marker='s', label = r't$_{\rm dyn}$')
			plt.scatter(outflow_df['z'], outflow_df['migration_time[Gyr]'], s=200, color='r', label=r'v$_{r} > 0$ t$_{\rm migr}$')
			plt.scatter(inflow_df['z'], inflow_df['migration_time[Gyr]'], s=200, color='b', label=r'v$_{r} < 0$ t$_{\rm migr}$')
	plt.semilogy()
	plt.axvline(x=1, color='k', lw=2)
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.xlabel(r'$z$', fontsize=16)
	plt.ylabel(r'log$_{10}$ timescale [Gyr]', fontsize=16)
	plt.legend(fontsize=14)
	plt.ylim(8e-3, 2e2)
	if fullden:
		plt.savefig(models[0] + '_' + models[-1] + '_plot_fullden_freefall_vs_migrtime_flowcolors.pdf')
	else:
		plt.annotate('spherical', xy=(0.8, 0.1), xycoords='axes fraction', fontsize=16) 
		plt.savefig(models[0] + '_' + models[-1] + '_plot_freefall_vs_migrtime_flowcolors.pdf')
	plt.show()

def plot_r_vr_zcolor(models):	
	fig = plt.figure(figsize=(12, 10))
	for i in range(len(models)):
		df = pd.read_table(models[i] + '_migration_fullden_freefall_times_z.dat', sep='\s+')
		plt.scatter(df['avgvr[kms^-1]'], df['avgr[kpc]'], s=50, c=df['z'].to_numpy(), marker='o', ec='white', cmap='Reds')
		print('done w model ', models[i])
	plt.loglog()
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	cb = plt.colorbar()
	cb.set_label(r'r$z$', fontsize=16)
	plt.ylabel(r'log$_{10}$ $\bar{r}$ [kpc]', fontsize=16)
	plt.xlabel(r'log$_{10}$  $\bar{v_r}$ [km s$^{-1}$]', fontsize=16)
	#plt.legend(fontsize=14)
	#plt.ylim(8e-3, 2e2)
	plt.savefig('../pdfs/' + models[0] + '_' + models[-1] + '_plot_r_vs_vr_zcolor.pdf')
	plt.show()

def GMs_metallicity_plot(models, sinds, show, makeprofs, scatter):
	Zsun = 0.012
	fig = plt.figure(figsize=(12, 10))
	for i in range(len(models)):
		halo, output = pyfof_HI_clumps.get_model_outputs_halos(models[i])
		halo = halo[sinds[i]:sinds[i]+2]
		output = output[sinds[i]:sinds[i]+2] 
		fn2 = pyfof_HI_clumps.get_fn(models[i])
		infile = fn2 + 'timetrace_grps_lowtohigh_outputs' + str(output[0]) + '_' + str(output[1]) + '_clump_properties_metals_fixvr.dat'
		df  = pd.read_table(infile, sep='\s+')
		df = df[df['size[kpc]'] > 0.25]	
		logmetals = np.array(df['logmetals'])
		r = np.array(df['r[kpc]'])
		halo_track_info = pd.read_table('../datfiles/' + models[i] + '_halo_track_info.dat', sep='\s+')
		z = float(halo_track_info[halo_track_info['output'] == output[0]]['z'])
		if makeprofs:
			h1 = pyfof_HI_clumps.get_mainHalo(models[i], '00' + str(output[0]), halo[0])
			try: 
				h1 = h1[h1['amiga.grp']==halo[0]]
				print('got rid of substructure')
				outfile = 'GMs_z' + str(z) + '_metallicity_profile_nosub' + '.pdf'
			except: 
				print('amiga.grp files dont exist for model ' + models[i])
				outfile = 'GMs_z' + str(z) + '_metallicity_profile' + '.pdf'
			CGM = pyfof_HI_clumps.cut_gal(h1)
			CGMprofile = pynbody.analysis.profile.Profile(CGM, min='10 kpc', max='250 kpc')
			cold_CGM = CGM[CGM['temp'].in_units('K') < 10**5]
			cold_CGMprofile = pynbody.analysis.profile.Profile(cold_CGM, min='10 kpc', max='250 kpc')
			plt.plot(cold_CGMprofile['rbins'], cold_CGMprofile['metals']/Zsun, color=pyfof_HI_clumps.model_colors(models[i]), ls='-', lw=3, label=models[i] + r'background CGM')
			np.savetxt('../txtfiles/' + models[i] + '_coldCGM_prof_rbins.txt', cold_CGMprofile['rbins'])
			np.savetxt('../txtfiles/' + models[i] + '_coldCGM_metalsbyZsum.txt', cold_CGMprofile['metals']/Zsun)
		else: 
			rbins = np.loadtxt('../txtfiles/' + models[i] + '_coldCGM_prof_rbins.txt')
			metalsbyzsun = np.loadtxt('../txtfiles/' + models[i] + '_coldCGM_metalsbyZsum.txt')
			plt.plot(rbins, metalsbyzsun, color=pyfof_HI_clumps.model_colors(models[i]), ls='-', lw=3, label=models[i] + r'background CGM')
			outfile = 'GMs_z' + str(z) + '_metallicity_madeprofile' + '.pdf'
		if scatter:
			plt.scatter(r, 10**logmetals/Zsun, s=50, alpha=0.5, color=pyfof_HI_clumps.model_colors(models[i]), edgecolor='k', label=models[i] + ' fof cold clumps')
			outfile = outfile[:-3] + '_scatter.pdf'
	plt.loglog()
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.xlabel(r'log$_{10}$ r$_{\rm clump}$ [kpc]', fontsize=16)
	plt.ylabel(r'$Z/Z_{\odot}$', fontsize=16)
	plt.legend(fontsize=14)
	#plt.ylim(1e-1, 2e5)
	plt.title(r' $z = $' + str(z), fontsize=16)
	plt.savefig(outfile)
	if show:
		plt.show()
	else:
		plt.close()

def make_metal_profile_jointplts(models, sinds): 
	rc={'font.size': 14, 'axes.labelsize': 14, 'legend.fontsize': 14, 'legend.title_fontsize': 14,
    	'axes.titlesize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'xtick.direction': 'in', 'ytick.direction': 'in'}
	sns.set(rc=rc)
	sns.set_style("ticks", rc=rc)
	Zsun = 0.012
	sizelim = 0.25
	pcolor = []
	for i in range(len(models)):
		halo, output = pyfof_HI_clumps.get_model_outputs_halos(models[i])
		halo = halo[sinds[i]:sinds[i]+2]
		output = output[sinds[i]:sinds[i]+2] 
		fn2 = pyfof_HI_clumps.get_fn(models[i])
		infile = fn2 + 'timetrace_grps_lowtohigh_outputs' + str(output[0]) + '_' + str(output[1]) + '_clump_properties_metals_fixvr.dat'
		df  = pd.read_table(infile, sep='\s+')
		pcolor.append(pyfof_HI_clumps.model_colors(models[i]))
		if(i==0):
			dfs = df
		else:
			dfs = pd.concat([df, dfs])
	dfs = dfs[dfs['size[kpc]'] > sizelim]
	dfs['logr[kpc]'] = np.log10(dfs['r[kpc]'])
	dfs['logmetalsbyzsun'] = np.log10(10**dfs['logmetals']/Zsun)
	g = sns.jointplot(data=dfs, x='r[kpc]', y='logmetalsbyzsun', hue='model', height=10, palette=pcolor, legend=True)
	g.set_axis_labels(r'log$_{10}$ r [kpc]', r'log$_{10} Z/Z_{\odot}$')
	plt.show() 
