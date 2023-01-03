import timescales

models = ['P0', 'P0noBHs', 'GM1', 'GM1noBHs', 'GM2']

for m in models: 
	timescales.sns_displot_clumpacc(m)
