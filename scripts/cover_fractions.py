import pandas as pd
import numpy as np 
import sns_clump_phase_plots
import pyfof_HI_clumps
import matplotlib.pyplot as plt 

def plt_NHI_vs_coverfrac(gal): 
	models = ['P0', 'P0noBHs', 'GM1', 'GM1noBHs', 'GM2', 'GM2noBHs', 'GM3', 'GM3noBHs', 'GM3SI1']
	NHI_data = pd.read_table('GMs_covering_frac_003456_gal' + str(gal) + '.dat', sep='\s+')
	h148_NHI_data = pd.read_table('h148_covering_frac_003456_gal' + str(gal) + '.dat', sep='\s+')
	CGM2_data = pd.read_csv('akaxia_cov_frac_data.csv')
	fig = plt.figure(figsize=(12,10))
	for i in range(len(models)):
    		data = NHI_data[NHI_data['model'] == models[i]]
    		plt.plot(data['NHI_threshold'], data['covering_frac'], color=sns_clump_phase_plots.model_colors(models[i]), lw=2, label=models[i])
	#plt CGM data 
	plt.errorbar(10**CGM2_data['logN_thresh'], CGM2_data['cov_frac'], yerr=[CGM2_data['cf_lolims'], CGM2_data['cf_uplims']], c='k', label=r'CGM$^2$')
	#plt.errorbar(10**CGM2_data['logN_thresh'], CGM2_data['cov_frac'], yerr=CGM2_data['cf_uplims'], c='k')
	plt.plot(h148_NHI_data['NHI_threshold'], data['covering_frac'], color='r', lw=2, label='h148')
	plt.semilogx()
	plt.ylabel('Covering fraction < 150 kpc', fontsize = 20)
	plt.xlabel(r'Log$_{10}$ N$_{\rm HI}$ [cm$^{-2}$]', fontsize=20)
	plt.tick_params(labelsize=18)
	plt.legend(fontsize=18)
	#plt.ylim(0,1)
	#plt.show()
	plt.savefig('GMs_coveringfrac_vs_NHIthreshold_gal' + str(gal) + '_semilogxi_withCGM2_wh148.pdf')
