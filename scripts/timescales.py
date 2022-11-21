import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

def subplot_avg_timescales_ratio_flowcolors(models, fullden):	
	fig, ax = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
	for i in range(len(models)):
		if fullden:
			df = pd.read_table('../datfiles/' + models[i] + '_migration_fullden_freefall_times_z.dat', sep='\s+')
		else:
			df = pd.read_table(models[i] + '_migration_freefall_times_z.dat', sep='\s+')
		outflow_df = df[df['avgvr[kms^-1]'] > 0]
		inflow_df = df[df['avgvr[kms^-1]'] < 0]
		ax[0].scatter(df['z'], df['freefall[Gyr]'], s=200, color='k', marker='s')
		ax[0].scatter(outflow_df['z'], outflow_df['migration_time[Gyr]'], s=200, color='r')
		ax[0].scatter(inflow_df['z'], inflow_df['migration_time[Gyr]'], s=200, color='b')
		ax[1].scatter(outflow_df['z'], outflow_df['migration_time[Gyr]']/outflow_df['freefall[Gyr]'], s=200, color='r')
		ax[1].scatter(inflow_df['z'], inflow_df['migration_time[Gyr]']/inflow_df['freefall[Gyr]'], s=200, color='b')
		
		if(i == 0):
			ax[0].scatter(df['z'], df['freefall[Gyr]'], s=200, color='k', marker='s', label = r't$_{\rm dyn}$')
			ax[0].scatter(outflow_df['z'], outflow_df['migration_time[Gyr]'], s=200, color='r', label=r'v$_{r} > 0$ t$_{\rm migr}$')
			ax[0].scatter(inflow_df['z'], inflow_df['migration_time[Gyr]'], s=200, color='b', label=r'v$_{r} < 0$ t$_{\rm migr}$')
	ax[0].semilogy()
	ax[1].semilogy()
	ax[0].axvline(x=1, color='k', lw=2)
	ax[1].axvline(x=1, color='k', lw=2)
	ax[1].axhline(y=1, color='k', lw=2)
	ax[0].tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	ax[1].tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	ax[1].set_xlabel(r'$z$', fontsize=16)
	ax[0].set_ylabel(r'log$_{10}$ timescale [Gyr]', fontsize=16)
	ax[1].set_ylabel(r'log$_{10}$ t$_{\rm migr}$/t$_{\rm dyn}$', fontsize=16)
	ax[0].legend(fontsize=14)
	if fullden:
		plt.savefig(models[0] + '_' + models[-1] + '_subplot_fullden_freefall_vs_migrtime_flowcolors.pdf')
	else:
		plt.annotate('spherical', xy=(0.8, 0.1), xycoords='axes fraction', fontsize=16) 
		plt.savefig(models[0] + '_' + models[-1] + '_subplot_freefall_vs_migrtime_flowcolors.pdf')
	
	plt.tight_layout()
	plt.show()

