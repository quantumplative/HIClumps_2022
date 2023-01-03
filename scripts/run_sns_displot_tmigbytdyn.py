import timescales

#models = ['P0', 'P0noBHs', 'GM1', 'GM1noBHs', 'GM2']
models = ['P0']
for m in models: 
	timescales.sns_displot_tmigbytdyn(m)
