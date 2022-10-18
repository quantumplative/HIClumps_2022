import matplotlib.cm as cm
import numpy as np
import seaborn as sns
import clump_mass_function
import sns_clump_phase_plots
import pandas as pd
import matplotlib.pyplot as plt

rc={'font.size': 14, 'axes.labelsize': 14, 'legend.fontsize': 14, 'legend.title_fontsize': 14, 'axes.titlesize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14}	
sns.set(rc=rc)
sns.set_style("ticks", rc=rc)

def combine_dat_models(models, output):
	rows = 0
	for i in range(len(models)):
		dat_tab = clump_mass_function.open_pd_table(models[i], '00' + str(output), 4e-8, 1)
		dat_tab = clump_mass_function.pdtable_put_output(dat_tab, output)
		rows += len(dat_tab)
		print('model ' +  models[i] + ' output ' + str(output) + ' has ' + str(len(dat_tab)) + ' rows')
		if(i == 0):
			dat_tabf = dat_tab
		else:
			dat_tabf = combine_pdtables(dat_tabf, dat_tab)
	print('total rows = ' + str(rows))
	return dat_tabf

def combine_pdtables(t1, t2):
	return pd.concat([t1, t2])

def clump_qtys_jointplt(models, output): 
	rc={'font.size': 20, 'axes.labelsize': 20, 'legend.fontsize': 18, 'legend.title_fontsize': 20, 'axes.titlesize': 20, 'xtick.labelsize': 20, 'ytick.labelsize': 20, 'xtick.direction' : 'in', 'ytick.direction': 'in'}	
	sns.set(rc=rc)
	sns.set_style("ticks", rc=rc)

	z = {3840: 0.06, 3456: 0.17, 1536: 1.18, 384: 4.58}	
	
	colmap = []
	for i in range(len(models)):
		colmap.append(sns_clump_phase_plots.model_colors(models[i]))	
	dat_models = combine_dat_models(models, output)
	#cmaps = ['plasma_r', 'cividis_r', 'magma_r', 'twilight_r', 'inferno_r']
	#qtys = ['clump_mass[Msol]', 'avg_temp[K]', 'avg_pressure[Pa]',
       #'avg_nHden[cm^-3]', 'vdist[km/s]']
	#ylabels = [r'log$_{10}$ Clump Mass [M$_{\odot}$]', r'log$_{10}$ Temperature [K]', r'log$_{10}$ Pressure [Pa]', 'log$_{10}$ n$_{HI}$ [cm$^{-3}$]', 'log$_{10}$ $\sigma_V$ [km/s]']
	#ylimits = [[5, 10], [3.5, 6], [5, 12], [-8, 1], [-1, 2.5]]
	qtys = ['clump_mass[Msol]']
	ylabels = [r'log$_{10}$ Clump Mass [M$_{\odot}$]']
	ylimits = [[5, 10]]
	dat_models['logSize[kpc]'] = np.log10(dat_models['clump_r_avg[kpc]'])
	cut_clumps = dat_models[dat_models['logSize[kpc]'] < np.log10(.25)]
	print('cutting ' + str(len(cut_clumps)) + ' clumps out of ' + str(len(dat_models)))
	for i in range(len(qtys)):
		dat_models['log' + qtys[i]] = np.log10(dat_models[qtys[i]])
		cut_clumps['log' + qtys[i]] = np.log10(cut_clumps[qtys[i]])
		g = sns.jointplot(data=dat_models, x='logSize[kpc]', y='log' + qtys[i], hue='model', s=100, palette=colmap, legend=True, height=10)
		if(qtys[i] == 'avg_temp[K]'): 
			g.ax_joint.annotate(r'z = ' + str(z[output]), xy=(0.15, 0.9), xycoords='axes fraction', fontsize=20)
		else:
			g.ax_joint.annotate(r'z = ' + str(z[output]), xy=(0.85, 0.9), xycoords='axes fraction', fontsize=20)
		g.set_axis_labels(r'log$_{10}$ Size [kpc]', ylabels[i])
		g.ax_marg_y.set_ylim(ylimits[i][0], ylimits[i][1])
		g.ax_marg_x.set_xlim(-2.5, 2.5)
		g.ax_joint.axvline(x=np.log10(0.25), color='gray', lw=2, ls='--') # resolution limit in kpc - grav softening 
		g.ax_joint.scatter(cut_clumps['logSize[kpc]'], cut_clumps['log' + qtys[i]], s=100, fc='gray', ec='white')
		g.ax_joint.annotate(r'$\epsilon$',color='gray', xy=(-0.5, ylimits[i][1]-0.1), xycoords='data')
		#g.ax_marg_x.axvline(x=np.log10(0.25), color='k', lw=2) # for smoothed histograms 
		#g.ax_joint.legend(loc='upper left')
		if(qtys[i] == 'vdist[km/s]'): 
			outfn = '../plots/' + str(output) + '_vdist[kms^-1]' + '_jointplt.png'
		else: 
			outfn = '../plots/' + str(output) + '_' + qtys[i] + '_jointplt.png'
		plt.savefig(outfn)
		plt.show()
