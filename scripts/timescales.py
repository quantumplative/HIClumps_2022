import pynbody
import palettable
import numpy as np 
import pandas as pd
import seaborn as sns
import pyfof_HI_clumps
import matplotlib.pyplot as plt 
from sns_clump_phase_plots import model_colors 
from test_timetrace_clumpgrps import lin_interp, mass_prof, redshift_rvir_grpclump_dat

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

def make_timescales_accrection_df(model, si, ei): 
	size_lim = 0.25 #kpc
	G = pynbody.array.SimArray(6.6743 * 1e-11) 
	G.units = 'm**3 kg**-1 s**-2'
	outfn = '../datfiles/' + model + '_allclumps_migration_fullden_freefall_times_z.csv' 
	
	halo_track_info = pd.read_table('../datfiles/' + model + '_halo_track_info.dat', sep='\s+')
	for i in range(ei - si):
		try:
			df  = redshift_rvir_grpclump_dat(model, si+i, si+i+2, True)
			output = df['output1'].to_numpy()[0]
			print('halo output = ', output)
			h = halo_track_info[halo_track_info['output'] == output]['halo'].to_numpy()[0]
			marr = []
			fn2 = pyfof_HI_clumps.get_fn(model)
			h = halo_track_info[halo_track_info['output'] == output]['halo'].to_numpy()[0]
			print('halo number = ', h)
			if(len(str(output)) == 3):  
				halo = pyfof_HI_clumps.get_mainHalo(model, '000' + str(output), h)
			else: 
				halo = pyfof_HI_clumps.get_mainHalo(model, '00' + str(output), h)
			rbins, pmass = mass_prof(halo)
			vr   = pynbody.array.SimArray(df['vr[kms^-1]'].to_numpy())
			vr.units = 'km s**-1'
			r    = pynbody.array.SimArray(df['r[kpc]'].to_numpy())
			r.units = 'kpc'
			size = pynbody.array.SimArray(df['size[kpc]'].to_numpy())
			size.units = 'kpc'
			mass = pynbody.array.SimArray(df['mass[Msol]'].to_numpy())
			mass.units = 'Msol'
			tmig = r/vr 
			tmig_avg = tmig.mean()
			tmig_gyr = tmig_avg.in_units('Gyr')
			for j in range(len(r)):
				marr.append(lin_interp(rbins, pmass, r[j]))
			Mgal = pynbody.array.SimArray(marr)
			Mgal.units = 'Msol'
			print('Mgal at r = ', Mgal)
			tdyn = (np.pi*r**(3/2))/(G.in_units('kpc**3 Msol**-1 Gyr**-2')*Mgal)**(0.5)
			print('tdyn = ', tdyn)
			df['tmig[Gyr]'] = tmig.in_units('Gyr')
			df['tdyn[Gyr]'] = tdyn
			df['clump_acc[Msolyr^-1]'] = (mass / tmig.in_units('yr'))
			if(i==0):
				dfs = df
			else:
				dfs = pd.concat([df, dfs])
		except: 
			print('something up w index si + i ', si + i)
	dfs = dfs[dfs['size[kpc]'] > size_lim]
	dfs.to_csv(outfn)

def sns_displot_clumpacc(model): 
	infile = '../datfiles/' + model + '_allclumps_migration_fullden_freefall_times_z.csv' 
	df = pd.read_csv(infile)
	rc={'font.size': 14, 'axes.labelsize': 14, 'legend.fontsize': 14, 'legend.title_fontsize': 14, 'axes.titlesize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14}
	sns.set(rc=rc)
	sns.set_style("ticks", rc=rc)
	df['ab_clump_acc'] = abs(df['clump_acc[Msolyr^-1]'])	
	df['log_clump_acc'] = np.log10(df['ab_clump_acc'])
	g = sns.displot(data=df, x='log_clump_acc', hue='z', kind='kde', height=10, palette='magma')
	g.set_xlabels(r'log$_{10}$ $\dot{\rm M}_{\rm clump}$ [M$_{\odot}$ yr$^{-1}$]')
	g.savefig('../plots/' + model + '_kde_clump_acc.pdf')
	plt.show()


def sns_displot_tmigbytdyn(model): 
	infile = '../datfiles/' + model + '_allclumps_migration_fullden_freefall_times_z.csv' 
	df = pd.read_csv(infile)
	rc={'font.size': 14, 'axes.labelsize': 14, 'legend.fontsize': 14, 'legend.title_fontsize': 14, 'axes.titlesize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14}
	sns.set(rc=rc)
	sns.set_style("ticks", rc=rc)
	df['logtmigbytdyn'] = np.log10(abs(df['tmig[Gyr]']/df['tdyn[Gyr]']))
	n = len(np.unique(df['output1']))	
	evenly_spaced_interval = np.linspace(0.1, 1, n)
	colors = [palettable.cartocolors.sequential.Sunset_7.mpl_colormap(x) for x in evenly_spaced_interval]
	g = sns.displot(data=df, x='logtmigbytdyn', hue='z', kind='kde', height=10, palette=colors)
	g.set_xlabels(r'log$_{10}$ t$_{\rm mig}$/t$_{\rm dyn}$')
	g.savefig('../plots/' + model + '_kde_tmigbytdyn_palettable.pdf')
	plt.show()

def sns_catplot_tmigbytdyn(model):
	infile = '../datfiles/' + model + '_allclumps_migration_fullden_freefall_times_z.csv' 
	df = pd.read_csv(infile)
	rc={'font.size': 14, 'axes.labelsize': 14, 'legend.fontsize': 14, 'legend.title_fontsize': 14, 'axes.titlesize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14}
	sns.set(rc=rc)
	sns.set_style("ticks", rc=rc)
	df['logtmigbytdyn'] = np.log10(abs(df['tmig[Gyr]']/df['tdyn[Gyr]']))
	n = len(np.unique(df['output1']))	
	evenly_spaced_interval = np.linspace(0.1, 1, n)
	colors = [palettable.cartocolors.sequential.Sunset_7.mpl_colormap(x) for x in evenly_spaced_interval]
	g = sns.catplot(data=df, y='logtmigbytdyn', x='z', kind='boxen', height=10, palette='Reds')
	g.set_ylabels(r'log$_{10}$ t$_{\rm mig}$/t$_{\rm dyn}$')
	g.set_xlabels(r'z')
	g.savefig('../plots/' + model + '_catplot_tmigbytdyn_Reds.pdf')
	plt.show()

def comparemodel_sns_catplot_tmigbytdyn(models, outputs): 
	rc={'font.size': 14, 'axes.labelsize': 14, 'legend.fontsize': 14, 'legend.title_fontsize': 14, 'axes.titlesize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14}
	sns.set(rc=rc)
	sns.set_style("ticks", rc=rc)
	pcolors = []
	for i in range(len(models)):
		infile = '../datfiles/' + models[i] + '_allclumps_migration_fullden_freefall_times_z.csv'
		df = pd.read_csv(infile)
		pcolors.append(model_colors(models[i]))
		print('model ' + models[i])
		for j in range(len(outputs)):
			print('output ' + str(outputs[j]))
			if(i == 0 and j == 0):
				dfs = df[df['output1'] == outputs[j]]
			else:
				dfs = pd.concat([dfs, df[df['output1'] == outputs[j]]])
	dfs['logtmigbytdyn'] = np.log10(abs(dfs['tmig[Gyr]']/dfs['tdyn[Gyr]']))
	g = sns.catplot(data=dfs, y='logtmigbytdyn', x='z', kind='box', hue='model', height=10, palette=pcolors)
	g.set_ylabels(r'log$_{10}$ t$_{\rm mig}$/t$_{\rm dyn}$')
	g.set_xlabels(r'z')
	g.savefig('../plots/compareGMs' + models[0] + '_' + models[-1] + '_catplot_tmigbytdyn_box.pdf')
	plt.show()

def stat_clump_acc_rate(model):
	infile = '../datfiles/' + model + '_allclumps_migration_fullden_freefall_times_z.csv'
	df = pd.read_csv(infile)
	df_inflow = df[df['tmig[Gyr]'] < 0]
	df_inflow['abtmig[Gyr]'] = abs(df_inflow['tmig[Gyr]'])
	df_inflow = df_inflow[df_inflow['abtmig[Gyr]'] < 1e6]	
	outputs = np.unique(df_inflow['output1'])
	z = np.unique(df_inflow['z'])
	mean_clumpaccr = [] 
	median_clumpaccr = [] 
	print('outputs ', outputs)
	for i in range(len(outputs)):
		mean_tmig = df_inflow[df_inflow['output1'] == outputs[i]]['abtmig[Gyr]'].mean()
		median_tmig = df_inflow[df_inflow['output1'] == outputs[i]]['abtmig[Gyr]'].median()	
		sum_mass = df_inflow[df_inflow['output1'] == outputs[i]]['mass[Msol]'].sum()	
		mean_clumpaccr.append(sum_mass/(mean_tmig * 1e9))
		median_clumpaccr.append(sum_mass/(median_tmig * 1e9))
	np.savetxt('../datfiles/' + model + '_clumpacc_z.txt', z)
	np.savetxt('../datfiles/' + model + '_mean_clumpacc.txt', mean_clumpaccr)
	np.savetxt('../datfiles/' + model + '_median_clumpacc.txt', median_clumpaccr)

def plot_mean_clump_acc(models): 
	fig = plt.figure(figsize = (10, 10))
	for m in models:
		z = np.loadtxt('../datfiles/' + m + '_clumpacc_z.txt')
		mean_clumpaccr = np.loadtxt('../datfiles/' + m + '_mean_clumpacc.txt')
		plt.plot(z, mean_clumpaccr, lw=2, color=model_colors(m))
	plt.xlabel(r'$z$', fontsize=16)
	plt.ylabel(r'$\dot{\rm{M}}_{\rm clump}$ [M$_{\odot}$ yr$^{-1}$]', fontsize=16)
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.savefig('../plots/' + models[0] + models[-1] + '_mean_stat_clumpaccr.pdf')
	plt.show()

def plot_median_clump_acc(models): 
	fig = plt.figure(figsize = (10, 10))
	for m in models:
		z = np.loadtxt('../datfiles/' + m + '_clumpacc_z.txt')
		median_clumpaccr = np.loadtxt('../datfiles/' + m + '_median_clumpacc.txt')
		plt.plot(z, median_clumpaccr, lw=2, color=model_colors(m))
	plt.xlabel(r'$z$', fontsize=16)
	plt.ylabel(r'$\dot{\rm{M}}_{\rm clump}$ [M$_{\odot}$ yr$^{-1}$]', fontsize=16)
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.savefig('../plots/' + models[0] + models[-1] + '_median_stat_clumpaccr.pdf')
	plt.show()
	

