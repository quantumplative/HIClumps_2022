import pynbody
import matplotlib.pyplot as plt 
import pandas as pd 
import numpy as np 
import pyfof_HI_clumps
import pynbody.plot as pp

def subplot_SFH_Nclumps(models):
	fig = plt.figure(figsize=(12, 8))
	grid = plt.GridSpec(4, 1, wspace=.25, hspace=.25)
	plt.subplots_adjust(wspace=.25, hspace=.25)
	ax1 = fig.add_subplot(grid[:2, 0])
	ax2 = fig.add_subplot(grid[2:, 0], sharex=ax1)
	plt.setp(ax2.get_xticklabels(), visible=False)
	for i in range(len(models)):
		halo, output = pyfof_HI_clumps.get_model_outputs_halos(models[i])
		df = pd.read_table(models[i] + '_Nclumps_time_z.dat' , sep='\s+')
		h1 = pyfof_HI_clumps.get_mainHalo(models[i], '00' + str(output[-1]), halo[-1])
		#plt.subplot(grid[:2, 0]) 
		sf, rbins = pp.sfh(h1,trange=[0, 12], lw=2,color=pyfof_HI_clumps.model_colors(models[i]), label=models[i], subplot=ax1, legend=True)
		plt.subplot(grid[2:, 0])
		plt.plot(df['t[Gyr]'], df['Nclump'], lw=2, ls='-', color=pyfof_HI_clumps.model_colors(models[i]))
	plt.subplot(grid[:2, 0])
	plt.ylabel(r'log$_{10}$ SFR [M$_{\odot}$ yr$^{-1}$]', fontsize=16)
	plt.legend(fontsize=12)
	plt.tick_params(axis='both', labelsize=16, direction='in')
	plt.xlim(0, 12)
	plt.subplot(grid[2:, 0])
	plt.xlim(0, 12)
	plt.tick_params(axis='both', labelsize=16, direction='in')
	plt.xlabel('time [Gyr]', fontsize=16)
	plt.ylabel(r'log$_{10}$ N$_{clumps}$', fontsize=16)
	plt.semilogy()
	plt.savefig('subplots_12Gyr_' + models[0] + '-' + models[-1] + '_SFHandNClumps.pdf')
	plt.show()	

def plot_SFH_Nclumps(models):
	fig, ax = plt.subplots(figsize= (15, 10))
	ax2 = ax.twinx() 
	for i in range(len(models)):
		halo, output = pyfof_HI_clumps.get_model_outputs_halos(models[i])
		df = pd.read_table(models[i] + '_Nclumps_time_z.dat' , sep='\s+')
		h1 = pyfof_HI_clumps.get_mainHalo(models[i], '00' + str(output[-1]), halo[-1])
		sf, rbins = pp.sfh(h1,trange=[0, 14], lw=2,color=pyfof_HI_clumps.model_colors(models[i]), subplot=ax, label=models[i], legend=True,log=True)
		ax2.plot(df['t[Gyr]'], df['Nclump'], lw=2, ls='--', color=pyfof_HI_clumps.model_colors(models[i]))
		print('done with ', models[i])
	ax.set_ylabel(r'SFR [M$_{\odot}$/yr]', fontsize=16) 
	ax.set_xlabel(r't [Gyr]', fontsize=16)
	ax2.set_ylabel(r'N$_{\rm clump}$', fontsize=16)
	ax2.semilogy() 	
	ax.tick_params(axis='both', labelsize=16, direction='in')
	ax2.tick_params(axis='both', labelsize=16, direction='in')
	plt.savefig(models[0] + '-' + models[-1] + '_SFHandNClumps.pdf')
	plt.show()	
